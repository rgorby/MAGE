! Defining the base class for boundary conditions
module basebc
    use types

    implicit none

    ! Base class for Boundar Conditions. All other BCs will inherit from this.
    type :: baseBC_T
        contains

        procedure baseInit
        procedure :: baseFailBC

        ! functions which will be over-written by sub-classes
        procedure :: doInit => baseInit
        procedure :: doBC => baseFailBC

    end type baseBC_T

    contains

    !A null initialization function for BCs that don't require initialization
    subroutine baseInit(bc,Model,Grid,State,inpXML)
        class(baseBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State
        type(XML_Input_T), intent(in) :: inpXML
    end subroutine baseInit

    ! this should never be called
    subroutine baseFailBC(bc,Model,Grid,State)
        class(baseBC_T), intent(in) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        print *, "Base BC function called. This should never happen"
        stop
    end subroutine baseFailBC

end module basebc

