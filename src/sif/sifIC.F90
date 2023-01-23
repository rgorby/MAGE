module sific
    use XML_Input

    use sifdefs
    use siftypes
    use sifgrids

    implicit none

    contains

    subroutine initSifIC_DIP(Model, Grid, State, iXML)
        type(sifModel_T) , intent(in)    :: Model
        type(sifGrid_T)  , intent(in)    :: Grid
        type(sifState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in) :: iXML

        real(rp) :: D, T, L, dL, tiote

        write(*,*) "Initializing a Maxwellian in a dipole field"

        call iXML%Set_Val(D ,"prob/D" , 10.0 )  ! [#/cc]
        call iXML%Set_Val(T ,"prob/T" , 30.0 )  ! [keV]
        call iXML%Set_Val(L ,"prob/L" , 4.0  )  ! [Re]
        call iXML%Set_Val(dL,"prob/dL", 0.625)  ! [Re]
        call iXML%Set_Val(tiote,"prob/tiote", 4.0)  ! Ratio of ion temp over electron temp


    end subroutine initSifIC_DIP

end module sific