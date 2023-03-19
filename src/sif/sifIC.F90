module sific
    use kdefs
    use XML_Input
    use planethelper

    use sifdefs
    use siftypes
    use sifgrids
    use sifetautils

    implicit none

    contains

    subroutine initSifIC_DIP(Model, Grid, State, iXML)
        type(sifModel_T) , intent(in)    :: Model
        type(sifGrid_T)  , intent(in)    :: Grid
        type(sifState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in)    :: iXML

        real(rp) :: D, kT_i, kT_e, L_peak, dL, tiote, cpcp
            !! Things we get from the xml file
        real(rp) :: D0, L, bVol, vm
            !! Calculated quantities
        integer :: i,j,sIdx
            !! Loops and indices

        write(*,*) "Initializing a Maxwellian in a dipole field"

        ! Using L-dependent pressure profile from Liemohn 2003
        call iXML%Set_Val(D0,"prob/D" , 10.0 )  ! [#/cc]
        call iXML%Set_Val(kT_i,"prob/T" , 30.0 )  ! [keV]
        call iXML%Set_Val(L_peak,"prob/L" , 4.0  )  ! [Re]
        call iXML%Set_Val(dL,"prob/dL", 0.625)  ! [Re]
        call iXML%Set_Val(tiote,"prob/tiote", 4.0)  ! Ratio of ion temp over electron temp
        call iXML%Set_Val(cpcp,"prob/cpcp", 50.0)  ! Cross-polar cap potential [kV]
        kT_e = kT_i/tiote
        ! Start by setting all variables we init to zero
        State%eta = 0.0
        State%Bmin = 0.0
        State%xyzMin = 0.0
        State%topo = 1  ! Everything closed
        State%active = 1 ! Everything active
        State%espot = 0.0
        State%thc = 0.0  ! Conjugate points are the same
        State%phc = 0.0
        State%bvol = 0.0

        ! Ni+1, Nj+1 variables
        do i=1,Grid%shGrid%Nt+1
            L = DipColat2L(Grid%shGrid%th(i))
            do j=1,Grid%shGrid%Np+1
                State%xyzMin(i,j,XDIR) = L*cos(Grid%shGrid%ph(j))
                State%xyzMin(i,j,YDIR) = L*sin(Grid%shGrid%ph(j))
            enddo
        enddo

        ! Ni, Nj variables
        do i=1,Grid%shGrid%Nt
            ! Fully symmetric, so only need to set certain things per i
            L = DipColat2L(Grid%shGrid%thc(i))
            D = D0*exp(-abs(L-L_peak)/dL)

            bVol = DipFTV_L(L, Model%planet%magMoment) ! [Rx/nT]
            vm = bVol**(-2./3.)
            write(*,*) "L=",L
            do j=1,Grid%shGrid%Np
            ! Init etas
                ! Protons first
                sIdx = spcIdx(Grid, F_HOTP)
                associate(spc => Grid%spc(sIdx))
                    call DkT2SpcEta(spc, State%eta(i,j,spc%kStart:spc%kEnd), D, kT_i, vm)
                end associate

                ! Then electrons
                sIdx = spcIdx(Grid, F_HOTE)
                associate(spc => Grid%spc(sIdx))
                    call DkT2SpcEta(spc, State%eta(i,j,spc%kStart:spc%kEnd), D, kT_e, vm)
                end associate

                ! TODO: Implement psphere IC
                if (Model%doPlasmasphere .and. spcExists(Grid, F_PSPH)) then
                    write(*,*) "SIF Warning: plasmasphere on but no IC is implemented in initSIFIC_DIP."
                endif

            ! Bmin surface vars
                State%Bmin(i,j,ZDIR) = Model%planet%magMoment/L**3.0
                State%bvol(i,j) = bVol
            
            ! Conjugate points
                State%thc(i,j) = Grid%shGrid%thc(i)
                State%phc(i,j) = Grid%shGrid%phc(j)
    
            ! Electrostatic potential [kV]
                ! Taken from Toffoletto's rcm.x
                State%espot(i,j) = -cpcp/2. * sin(Grid%shGrid%phc(j)) &
                                    * Grid%shGrid%thc(1)/Grid%shGrid%thc(i) & 
                                    + Model%planet%psiCorot/L*sin(Grid%shGrid%thc(i))**2
            enddo  
        enddo


    end subroutine initSifIC_DIP

end module sific