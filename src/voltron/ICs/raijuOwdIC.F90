module raijuuseric
    !! IC for one way driving
    use XML_Input
    use volttypes

    use raijuTypes

    implicit none

    contains


    subroutine RAIJUinitStateUserIC(Model,Grid,State,inpXML)
        type(raijuModel_T) , intent(in)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        write(*,*)"lol"
        stop
    end subroutine RAIJUinitStateUserIC

end module raijuuseric