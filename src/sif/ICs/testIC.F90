module sifuseric
    use XML_Input

    use siftypes

    implicit none

    contains

    subroutine initSifUserIC(Model,Grid,State,inpXML)
        type(sifModel_T) , intent(in)    :: Model
        type(sifGrid_T)  , intent(in)    :: Grid
        type(sifState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        write(*,*)"Hey"
    end subroutine initSifUserIC

end module sifuseric