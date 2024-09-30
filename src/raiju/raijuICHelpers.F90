module raijuICHelpers
    use kdefs
    use XML_Input
    use planethelper
    use math
    use earthhelper

    use raijudefs
    use raijutypes
    use raijugrids
    use raijuetautils

    implicit none

    contains

!------
! Routines responsible for making sure full state is defined
!------

    subroutine initRaijuIC_DIP(Model, Grid, State, iXML)
        !! Most general and customizable configuration
        !! Will set dipole B field but get rest of definition (eta and electrostatic potential) from xml input
        type(raijuModel_T) , intent(in)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in)    :: iXML

        write(*,*) "Initializing a dipole field with XML-determined electrostatic potential and eta distribution"

        call initBDipole(Model, Grid, State)
        call initEspotPresets(Model, Grid, State, iXML)
        call initEtaPresets(Model, Grid, State, iXML)
        call EvalMoments(Grid, State)

    end subroutine initRaijuIC_DIP


    subroutine initRaijuIC_testThetaAdvection(Model, Grid, State, iXML)
        !! Dipole field, purely radial convection
        !! V will vary, up to postprocessing to check and make sure its right everywhere
        type(raijuModel_T) , intent(in)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in)    :: iXML

        integer :: i,j
        integer :: jMid, jWidth, iWidth
        real(rp) :: pDrop, delPhi

        ! B-field
        call initBDipole(Model, Grid, State)

        !Espot
        call iXML%Set_Val(jWidth,"test/jWidth" , 6)
        call iXML%Set_Val(pDrop,"prob/cpcp" , 50.0_rp)  ! keV
        jMid = Grid%shGrid%Np/2
        delPhi = Grid%shGrid%ph(jMid+jWidth+1) - Grid%shGrid%ph(jMid-jWidth)

        State%espot(:,:jMid-jWidth) = 0
        State%espot(:,jMid+jWidth+1:) = pDrop
        do j=jMid-jWidth, jMid+jWidth
            State%espot(:,j) = pDrop * (Grid%shGrid%ph(j) - Grid%shGrid%ph(jMid-jWidth)) / delPhi 
        enddo

        ! Eta
        iWidth = 4
        State%eta(Grid%shGrid%is:Grid%shGrid%is+iWidth,jMid-jWidth:jMid+jWidth,:) = 1.0
        State%eta_last = State%eta

        call EvalMoments(Grid, State)


    end subroutine initRaijuIC_testThetaAdvection

!------
! Hard-coded magnetic field configurations
!------

    subroutine initBDipole(Model, Grid, State)
        type(raijuModel_T) , intent(in)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State

        real(rp) :: L
        real(rp), dimension(3) :: xyzIono
        integer :: i,j
        

        State%Bmin   = 0.0
        State%xyzMin = 0.0
        State%thcon  = 0.0  ! Conjugate points are the same
        State%phcon  = 0.0
        State%bvol   = 0.0
        State%topo        = RAIJUCLOSED  ! Everything closed
        State%active      = RAIJUACTIVE  ! Everything active
        State%active_last = RAIJUACTIVE  ! Everything active

        ! Ni+1, Nj+1 corner variables
        do i=Grid%shGrid%isg,Grid%shGrid%ieg+1
            do j=Grid%shGrid%jsg,Grid%shGrid%jeg+1
                xyzIono(1) = Model%planet%ri_m/Model%planet%rp_m * sin(Grid%shGrid%th(i)) * cos(Grid%shGrid%ph(j))
                xyzIono(2) = Model%planet%ri_m/Model%planet%rp_m * sin(Grid%shGrid%th(i)) * sin(Grid%shGrid%ph(j))
                xyzIono(3) = Model%planet%ri_m/Model%planet%rp_m * cos(Grid%shGrid%th(i))

                L = DipoleL(xyzIono)

                State%xyzMin(i,j,XDIR) = L*cos(Grid%shGrid%ph(j))
                State%xyzMin(i,j,YDIR) = L*sin(Grid%shGrid%ph(j))
            ! Conjugate points
                State%thcon(i,j) = Grid%shGrid%th(i)
                State%phcon(i,j) = Grid%shGrid%ph(j)

            ! Bmin surface vars
                State%Bmin(i,j,ZDIR) = Model%planet%magMoment*G2nT/L**3.0  ! [nT]
                State%bvol(i,j) = DipFTV_L(L, Model%planet%magMoment) ! [Rp/nT]
            enddo
        enddo

        ! Ni, Nj cell-center variables
        do i=Grid%shGrid%isg,Grid%shGrid%ieg
            do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                xyzIono(1) = Model%planet%ri_m/Model%planet%rp_m * sin(Grid%shGrid%thc(i)) * cos(Grid%shGrid%phc(j))
                xyzIono(2) = Model%planet%ri_m/Model%planet%rp_m * sin(Grid%shGrid%thc(i)) * sin(Grid%shGrid%phc(j))
                xyzIono(3) = Model%planet%ri_m/Model%planet%rp_m * cos(Grid%shGrid%thc(i))

                L = DipoleL(xyzIono)
                State%xyzMincc(i,j,XDIR) = L*cos(Grid%shGrid%phc(j))
                State%xyzMincc(i,j,YDIR) = L*sin(Grid%shGrid%phc(j))
                State%bvol_cc(i,j) = DipFTV_L(L, Model%planet%magMoment) ! [Rp/nT]
            enddo
        enddo

    end subroutine initBDipole


!------
! Hard-coded eta distribution presets
!------

    subroutine initEtaPresets(Model, Grid, State, iXML)
        !! Based on xml input, creates desired eta distribution
        type(raijuModel_T) , intent(in)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in)    :: iXML

        character(len=strLen) :: epType
            !! Eta preset identifier
        real(rp) :: D, kT_i, kT_e, L_peak, dL, phi_peak, dPhi, tiote, phiMin, phiMax
            !! Quantities from xml
        real(rp) ::  D0, L, bVol, vm, dist
            !! Internal quantities
        real(rp) :: a0, a1, a2, a3
            !! For bell
        integer :: i,j,sIdx,kPsph
        
        
        call iXML%Set_Val(epType,"prob/etaPreset" , "RING")

        call iXML%Set_Val(D0,"prob/D" , 10.0 )  ! [#/cc]
        call iXML%Set_Val(kT_i,"prob/T" , 30.0 )  ! [keV]
        call iXML%Set_Val(tiote,"prob/tiote", 4.0)  ! Ratio of ion temp over electron temp
        kT_e = kT_i/tiote
        !call iXML%Set_Val(phiMin,"prob/phiMin", 0.0  )  ! Min phi to populate etas
        !call iXML%Set_Val(phiMax,"prob/phiMax", 360.0)  ! Min phi to populate etas

        State%eta = 0.0
        State%eta_last = 0.0

        select case(trim(epType))
            case("RING") ! Fully symmetric ring parameterized by Leimohn 2003
                call iXML%Set_Val(L_peak,"prob/L" , 4.0  )  ! [Re]
                call iXML%Set_Val(dL,"prob/dL", 0.625    )  ! [Re]
        
                ! ! $OMP PARALLEL DO default(shared) &
                ! ! $OMP schedule(dynamic) &
                ! ! $OMP private(i,j,L,bVol,vm,D)
                do i=Grid%shGrid%isg,Grid%shGrid%ieg
                    ! Fully symmetric, so only need to set certain things per i
                    L = DipColat2L(Grid%shGrid%thc(i))
                    bVol = DipFTV_L(L, Model%planet%magMoment) ! [Rx/nT]
                    vm = bVol**(-2./3.)

                    D = D0*exp(-abs(L-L_peak)/dL)
                    do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                        call applyDkT2Eta(Model, Grid, State, i, j, D, kT_i, kT_e, vm)
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
                        call applyDkT2Eta(Model, Grid, State, i, j, D, kT_i, kT_e, vm)
                    enddo
                enddo
                State%eta_last = State%eta
            
            case("WEDGE")  ! Small wedge of constant eta for all channels

                call iXML%Set_Val(D0,"prob/constEta" , 1.0)
                call iXML%Set_Val(L_peak  ,"prob/L" , 4.0 )  ! [Re], center L of wedge
                call iXML%Set_Val(dL      ,"prob/dL", 1.0 )  ! [Re], L extent of wedge
                call iXML%Set_Val(phi_peak,"prob/phi", 180.0)  ! [deg], center phi of wedge
                call iXML%Set_Val(dPhi    ,"prob/dPhi", 10.0)  ! [deg], phi width of wedge
                phi_peak = phi_peak*deg2rad
                dPhi = dPhi*deg2rad
                do i=Grid%shGrid%isg,Grid%shGrid%ieg
                    L = DipColat2L(Grid%shGrid%thc(i))
                    do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                        if (L > (L_peak - dL/2) .and. L < (L_peak + dL/2) .and. &
                            Grid%shGrid%phc(j) > (phi_peak - dPhi/2) .and. &
                            Grid%shGrid%phc(j) < (phi_peak + dPhi/2) ) then
                                State%eta(i,j,:) = D0
                        endif
                    enddo
                enddo
                
                State%eta_last = State%eta

            case("YLINE")  ! Line in the tail along Y-SM direction

                call iXML%Set_Val(D0,"prob/constEta" , 1.0)
                call iXML%Set_Val(L,"prob/Xval" , -6.0)  ! Note, just using L here to not make a whole bunch of variables for each option
                do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                    do i=Grid%shGrid%is,Grid%shGrid%ie
                        if ( any(State%xyzmin(i:i+1,j:j+1,XDIR) .le. L) .and. any(State%xyzmin(i:i+1,j:j+1,XDIR) .ge. L) ) then
                            State%eta(i,j,:) = D0
                        endif
                    enddo
                enddo
                State%eta_last = State%eta

            case default
                write(*,*) "(RAIJU) ERROR: Pick a valid eta preset in raijuICHelpers.F90:initEtaPresets"
                stop
       
        end select

        call setRaijuInitPsphere(Model, Grid, State, Model%psphInitKp)

        contains

        subroutine applyDkT2Eta(Model, Grid, State, i, j, D, kT_i, kT_e, vm)
            type(raijuModel_T) , intent(in)    :: Model
            type(raijuGrid_T)  , intent(in)    :: Grid
            type(raijuState_T) , intent(inout) :: State
            integer :: i, j
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

!------
! Electrostatic potential distributions
!------

    subroutine initEspotPresets(Model, Grid, State, iXML)
        !! Based on xml input, chooses desired electrostatic potential distribution
        !! Actual espot setting is done in different routines so that they are more generally accessible
        type(raijuModel_T) , intent(in)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in)    :: iXML

        character(len=strLen) :: espotType
            !! Eta preset identifier
        integer :: i,j
        
        
        call iXML%Set_Val(espotType,"prob/espotPreset" , "TOFFO")

        select case(trim(espotType))
            case("TOFFO") ! Simple distribution from Toffoletto's rcm.x
                call espot_TOFFO(Grid, State, iXML)
            case("BT") ! Simple distribution from Toffoletto's rcm.x
                call espot_BT(Grid, State, iXML)
            case default
                write(*,*) "(RAIJU) ERROR: Pick a valid electrostatic potential preset in raijuICHelpers.F90:initEspotPresets"
                stop
        end select

    end subroutine initEspotPresets


    subroutine espot_TOFFO(Grid, State, iXML)
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T)  , intent(in)    :: iXML
        
        real(rp) :: cpcp
        integer :: i,j

        State%espot = 0.0
        call iXML%Set_Val(cpcp,"prob/cpcp", 50.0)  ! Cross-polar cap potential [kV]

        ! Ni+1, Nj+1 variables
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg+1
            do i=Grid%shGrid%isg,Grid%shGrid%ieg+1
                ! Taken from Toffoletto's rcm.x
                State%espot(i,j) = -cpcp/2. * sin(Grid%shGrid%ph(j)) &
                                    * Grid%shGrid%th(1)/Grid%shGrid%th(i)
            enddo
        enddo
    end subroutine espot_TOFFO


    subroutine espot_BT(Grid, State, iXML)
        !! Electrostatic potential described in "Basic Space Plasma Physics" - Baumjohann and Treumann
        !! (page 99, Eqs 5.14-5.15)
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T)  , intent(in)    :: iXML
        
        real(rp) :: cpcp, gamma
        real(rp) :: Agamma
        real(rp) :: dy, L
        integer :: i,j

        !call iXML%Set_Val(E0,"prob/E0", 50.0)  ! Electric field in mV/m
        call iXML%Set_Val(cpcp,"prob/cpcp", 50.0)  ! cross polar cap potential in kV
        call iXML%Set_Val(gamma,"prob/fShield", 1.0)  ! Shielding factor. 1 = no shielding, 2-3 is more realistic

        ! Get deltaY
        j = Grid%shGrid%Np/4
        dy = State%xyzmin(Grid%shGrid%is,j,YDIR) - State%xyzmin(Grid%shGrid%ie,j,YDIR)  ! Rp

        Agamma = 0.5*cpcp*dy**(-1.0*gamma)
        
        do j = Grid%shGrid%jsg,Grid%shGrid%jeg+1
            do i = Grid%shGrid%isg,Grid%shGrid%ieg+1
                L = sqrt( State%xyzmin(i,j,XDIR)**2 + State%xyzmin(i,j,YDIR)**2 )  ! Rp
                State%espot(i,j) = -1*Agamma * L**gamma * sin(Grid%shGrid%ph(j))
            enddo
        enddo
        

    end subroutine espot_BT


!------
! Plasmasphere initialization
!------

    subroutine setRaijuInitPsphere(Model, Grid, State, Kp)
        type(raijuModel_T) , intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        real(rp) :: Kp

        integer :: psphIdx

        if (Model%doPlasmasphere .and. spcExists(Grid, F_PSPH)) then
            psphIdx = spcIdx(Grid, F_PSPH)
            State%eta     (:,:,Grid%spc(psphIdx)%kStart) = getInitPsphere(Grid, State, Kp)
            State%eta_last(:,:,Grid%spc(psphIdx)%kStart) = State%eta(:,:,Grid%spc(psphIdx)%kStart)
        endif
    
    end subroutine


    function getInitPsphere(Grid, State, Kp) result(etaPsph)
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp) :: Kp

        integer :: i,j
        real(rp) :: den, vm
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: etaPsph

        write(*,*) "RAIJU initializing plasmasphere with Kp =",Kp

        etaPsph = 0.0

        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,den,vm)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                den = GallagherXY(State%xyzMincc(i,j,XDIR), State%xyzMincc(i,j,YDIR), Kp)  ! [#/cc]
                vm = State%bvol_cc(i,j)**(-2./3.)
                etaPsph(i,j) = den/(vm**1.5)*sclEta
            enddo
        enddo

    end function getInitPsphere

end module raijuICHelpers