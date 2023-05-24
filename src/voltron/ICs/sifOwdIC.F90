module sifuseric
    !! IC for one way driving
    use XML_Input

    use sifTypes
    use sifCplTypes

    implicit none

    contains


    subroutine SIFinitStateUserIC(Model,Grid,State,inpXML)
        type(sifModel_T) , intent(in)    :: Model
        type(sifGrid_T)  , intent(in)    :: Grid
        type(sifState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        write(*,*)"lol"
        stop
    end subroutine SIFinitStateUserIC


    subroutine SIFinitCplUserIC(Model,Grid,State,cplBase)
        type(sifModel_T) , intent(in)    :: Model
        type(sifGrid_T)  , intent(in)    :: Grid
        type(sifState_T) , intent(in) :: State
        type(sif_cplBase_T), intent(inout) :: cplBase

        write(*,*)"lol"
        stop
    end subroutine SIFinitCplUserIC

end module sifuseric