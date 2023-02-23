module sific
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

        integer :: i,j,k
        real(rp) :: D, kT, L_peak, dL, tiote
            !! Things we get from the xml file
        real(rp) :: D0, L, bVol

        write(*,*) "Initializing a Maxwellian in a dipole field"

        ! Using L-dependent pressure profile from Liemohn 2003
        call iXML%Set_Val(D0,"prob/D" , 10.0 )  ! [#/cc]
        call iXML%Set_Val(kT ,"prob/T" , 30.0 )  ! [keV]
        call iXML%Set_Val(L_peak,"prob/L" , 4.0  )  ! [Re]
        call iXML%Set_Val(dL,"prob/dL", 0.625)  ! [Re]
        call iXML%Set_Val(tiote,"prob/tiote", 4.0)  ! Ratio of ion temp over electron temp

        !P0 = DkT2P(D,kT)  ! [nPa]

        ! Start by setting all species etas to zero
        State%eta = 0.0

        ! Fully symmetric, so only need to set certain things per i
        do i=1,Grid%shGrid%Nt
            L = DipColat2L(Grid%shGrid%thc(i))
            D = D0*exp(-(L-L_peak)/dL)

            bVol = DipFTV_L(L, Model%planet%magMoment) ! [Rx/nT]

            ! Protons first
            call DkT2SpcEta(Model, State, Grid%spc(HOTE)%flav, D, kT)
            
        enddo

    end subroutine initSifIC_DIP

end module sific