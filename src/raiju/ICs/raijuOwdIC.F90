module raijuuseric
    !! IC for one way driving
    use XML_Input
    use volttypes

    use raijuTypes
    use raijuICHelpers

    implicit none

    contains


    subroutine raijuInitState_useric(Model,Grid,State,inpXML)
        type(raijuModel_T) , intent(inout)    :: Model
        type(raijuGrid_T)  , intent(in)    :: Grid
        type(raijuState_T) , intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML

        ! Just do simple dipole step0 for now. Later we will probably do the whole empirical starter here
        write(*,*) "raijuInitState_useric called from raijuOwdIC but I'm just using default DIP settings"
        Model%initState => initRaijuIC_DIP
    end subroutine raijuInitState_useric

end module raijuuseric