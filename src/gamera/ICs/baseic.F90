! base type for ICs to inherit from and extend
module baseic

    implicit none

    !ID for different ICs
    integer, parameter :: SPRINKLER_IC_ID=1,REMIX_IC_ID=2

    type icHolder
        class(baseIC_T), allocatable :: p
    end type icHolder

    type baseIC_T
        contains

        ! functions which will be over-written by sub-classes
        procedure :: doInit => baseInit
        procedure :: getType => baseGetType
        procedure :: doBC => null
        procedure :: doAdjustE => null
        procedure :: doAdjustStep => null
        procedure :: doAdjustFlux => null

    end type baseIC_T

contains

    function baseInit(baseIC)
        class(baseIC_T), intent(inout) :: baseIC)
    end function baseInit

    function baseGetType(baseIC)
        integer :: baseGetType
        class(baseIC_T), intent(inout) :: baseIC
        print *,"Base version of getType called in baseic, this should not happen"
        baseGetType = -1
    end function baseGetType


end module baseic

