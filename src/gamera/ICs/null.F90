module usergamic
    use types
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use background

    implicit none

    contains
    
    subroutine initUser(Model,Grid,State,xmlInp)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        write(*,*) 'Attempting to use null user routine, bailing out!'
        stop
    end subroutine initUser

end module usergamic
