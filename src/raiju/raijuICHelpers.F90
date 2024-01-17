module raijuICHelpers
    use kdefs
    use XML_Input
    use planethelper

    use raijudefs
    use raijutypes
    use raijugrids
    use raijuetautils

    implicit none

    contains

    subroutine initRaijuIC_DIP(Model, Grid, State, iXML)
        type(raijuModel_T) , intent(in)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in)    :: iXML

        real(rp) :: cpcp
            !! Things we get from the xml file
        real(rp) :: L, bVol
            !! Internal quantities
        integer :: i,j
            !! Loops and indices

        write(*,*) "Initializing a dipole field with perscribed cpcp"

        ! Using L-dependent pressure profile from Liemohn 2003
        call iXML%Set_Val(cpcp,"prob/cpcp", 50.0)  ! Cross-polar cap potential [kV]
        ! Start by setting all variables we init to zero
        State%Bmin = 0.0
        State%xyzMin = 0.0
        State%topo = 1  ! Everything closed
        State%active = 1 ! Everything active
        State%espot = 0.0
        State%thcon = 0.0  ! Conjugate points are the same
        State%phcon = 0.0
        State%bvol = 0.0

        ! Ni+1, Nj+1 variables
        do i=Grid%shGrid%isg,Grid%shGrid%ieg+1
            L = DipColat2L(Grid%shGrid%th(i))
            bVol = DipFTV_L(L, Model%planet%magMoment) ! [Rx/nT]
            do j=Grid%shGrid%jsg,Grid%shGrid%jeg+1
                State%xyzMin(i,j,XDIR) = L*cos(Grid%shGrid%ph(j))
                State%xyzMin(i,j,YDIR) = L*sin(Grid%shGrid%ph(j))
            ! Conjugate points
                State%thcon(i,j) = Grid%shGrid%th(i)
                State%phcon(i,j) = Grid%shGrid%ph(j)

            ! Bmin surface vars
                State%Bmin(i,j,ZDIR) = Model%planet%magMoment/L**3.0
                State%bvol(i,j) = bVol

            ! Electrostatic potential [kV]
                ! Taken from Toffoletto's rcm.x
                State%espot(i,j) = -cpcp/2. * sin(Grid%shGrid%ph(j)) &
                                    * Grid%shGrid%th(1)/Grid%shGrid%th(i) !& 
                                    !+ Model%planet%psiCorot/L*sin(Grid%shGrid%th(i))**2
            enddo
        enddo

        ! Ni, Nj variables
        !do i=Grid%shGrid%isg,Grid%shGrid%ieg
        !    L = DipColat2L(Grid%shGrid%thc(i))
        !    do j=Grid%shGrid%jsg,Grid%shGrid%jeg
        !    enddo  
        !enddo

        ! Finally, init etas
        call initEtaPresets(Model, Grid, State, iXML)

        ! Calc moments
        call EvalMoments(Grid, State)


    end subroutine initRaijuIC_DIP


!------
! Hard-coded eta distribution presets
!------

    subroutine initEtaPresets(Model, Grid, State, iXML)
        type(raijuModel_T) , intent(in)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in)    :: iXML

        character(len=strLen) :: epType
            !! Eta preset identifier
        real(rp) :: D, kT_i, kT_e, L_peak, dL, tiote, phiMin, phiMax
            !! Quantities from xml
        real(rp) ::  D0, L, bVol, vm, dist
            !! Internal quantities
        real(rp) :: a0, a1, a2, a3
            !! For bell
        integer :: i,j,sIdx
        
        
        call iXML%Set_Val(epType,"prob/etaPreset" , "RING")

        call iXML%Set_Val(D0,"prob/D" , 10.0 )  ! [#/cc]
        call iXML%Set_Val(kT_i,"prob/T" , 30.0 )  ! [keV]
        call iXML%Set_Val(tiote,"prob/tiote", 4.0)  ! Ratio of ion temp over electron temp
        kT_e = kT_i/tiote
        !call iXML%Set_Val(phiMin,"prob/phiMin", 0.0  )  ! Min phi to populate etas
        !call iXML%Set_Val(phiMax,"prob/phiMax", 360.0)  ! Min phi to populate etas

        State%eta = 0.0

        select case(trim(epType))
            case("RING") ! Fully symmetric ring parameterized by Leimohn 2003
                call iXML%Set_Val(L_peak,"prob/L" , 4.0  )  ! [Re]
                call iXML%Set_Val(dL,"prob/dL", 0.625    )  ! [Re]
        
                do i=Grid%shGrid%isg,Grid%shGrid%ieg
                    ! Fully symmetric, so only need to set certain things per i
                    L = DipColat2L(Grid%shGrid%thc(i))
                    bVol = DipFTV_L(L, Model%planet%magMoment) ! [Rx/nT]
                    vm = bVol**(-2./3.)

                    D = D0*exp(-abs(L-L_peak)/dL)
                    do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                        call applyDkT2Eta(Model, Grid, State, D, kT_i, kT_e, vm)
                    enddo
                enddo
                State%eta_last = State%eta
            
            case("BELL")  ! Bell-like density pulse at L0 along midnight

                a0 = 0.355768
                a1 = 0.487396
                a2 = 0.144232
                a3 = 0.012604
                call iXML%Set_Val(L_peak,"prob/L" , 4.0  )  ! [Re], location of peak
                call iXML%Set_Val(dL    ,"prob/dL", 1.0  )  ! [Re], full width of bell

                ! Using Nuttall window cause it has continuous first derivative (https://en.wikipedia.org/wiki/Window_function)
                ! w[n] = a0 - a1*cos(2x) + a2*cos(4x) - a3*cos(6x), where x = PI*n/N
                ! n is location, N is width (dL for us)
                ! n = distance from L_peak, Done using coordinates in the equatorial plane
                ! L_peak assumed to be along midnight axis, so x_peak = -L_peak and y_peak = 0
                do i=Grid%shGrid%isg,Grid%shGrid%ieg
                    L = DipColat2L(Grid%shGrid%thc(i))
                    bVol = DipFTV_L(L, Model%planet%magMoment) ! [Rx/nT]
                    vm = bVol**(-2./3.)
                    do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                        dist = sqrt( (State%xyzMin(i,j,XDIR) + L_peak)**2 + (State%xyzMin(i,j,YDIR) - 0.0)**2 ) + dL/2
                        
                        if (dist > dL) then
                            !! Leave eta here at zero
                            cycle
                        endif
                        D = D0*(a0 &
                               -a1*cos(2*PI*dist/dL) &
                               +a2*cos(4*PI*dist/dL) &
                               -a3*cos(6*PI*dist/dL))
                        call applyDkT2Eta(Model, Grid, State, D, kT_i, kT_e, vm)
                    enddo
                enddo
                State%eta_last = State%eta
            
            case("LINE")  ! Line of constant eta for all channels

                call iXML%Set_Val(D0,"prob/constEta" , 1.0)
                j = Grid%shGrid%Np/2
                i = Grid%shGrid%Nt/2
                State%eta(i-6:i+6,j-6:j+6,:) = D0
                State%eta_last = State%eta

            case default
                write(*,*) "ERROR: Pick a valid eta preset in raijuICHelpers.F90:initEtaPresets"
                stop
       
        end select

         ! TODO: Implement psphere IC
        if (Model%doPlasmasphere .and. spcExists(Grid, F_PSPH)) then
            write(*,*) "RAIJU Warning: plasmasphere on but no IC is implemented in initRaijuIC_DIP."
        endif

        contains

        subroutine applyDkT2Eta(Model, Grid, State, D, kT_i, kT_e, vm)
            type(raijuModel_T) , intent(in)    :: Model
            type(raijuGrid_T)  , intent(in)    :: Grid
            type(raijuState_T) , intent(inout) :: State
            real(rp), intent(in) :: D, kT_i, kT_e, vm
            ! Protons first
            sIdx = spcIdx(Grid, F_HOTP)
            call DkT2SpcEta(Model, Grid%spc(sIdx), &
                        State%eta(i,j,Grid%spc(sIdx)%kStart:Grid%spc(sIdx)%kEnd), &
                        D, kT_i, vm)

            ! Then electrons
            sIdx = spcIdx(Grid, F_HOTE)
            call DkT2SpcEta(Model, Grid%spc(sIdx), &
                        State%eta(i,j,Grid%spc(sIdx)%kStart:Grid%spc(sIdx)%kEnd), &
                        D, kT_e, vm)
        end subroutine applyDkT2Eta

    end subroutine initEtaPresets
    

end module raijuICHelpers