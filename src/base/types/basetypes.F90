! Module containing base types from which other models will inherit

module basetypes

    use gdefs
    use ioclock
    use ioH5
    use xml_input

    implicit none

    type :: BaseModel_T
        !Output info
        type (IOClock_T) :: IO

        !Unit information
        !type (Units_T) :: Units
        ! should this be a type which is extended by other apps?

        !Whether this model is controlled by another model
        logical :: isSub = .false.

        contains
    end type BaseModel_T

    type :: BaseGrid_T
        contains
        ! this base class contains nothing, but its existence is useful for interfaces
    end type BaseGrid_T

    type, abstract :: BaseState_T
        !Time info
        real(rp) :: time = -HUGE
        real(rp) :: MJD = -HUGE

        contains
    end type BaseState_T

    type, abstract :: BaseApp_T
        contains
        ! the base App intentionally does NOT contain instances of the base model/grid/state
        ! this allows the specific sub-apps to include their specific instances of model/grid/state
        ! including the option to not use them, or have multiple copies of them

        ! functions which will be over-written by sub-apps
        procedure(InitModel_interface), deferred :: InitModel
        procedure(InitIO_interface), deferred :: InitIO
        procedure(WriteRestart_interface), deferred :: WriteRestart
        procedure(ReadRestart_interface), deferred :: ReadRestart
        procedure(WriteConsoleOutput_interface), deferred :: WriteConsoleOutput
        procedure(WriteFileOutput_interface), deferred :: WriteFileOutput
        procedure(WriteSlimFileOutput_interface), deferred :: WriteSlimFileOutput
        procedure(AdvanceModel_interface), deferred :: AdvanceModel

    end type BaseApp_T

    type, abstract :: BaseOptions_T
        contains
        ! this base class contains nothing, but will be extended and passed into Init by other classes
    end type BaseOptions_T

    abstract interface
        subroutine InitModel_interface(App, Xml)
            import BaseApp_T ! I don't want to talk about it
            import XML_Input_T
            class(BaseApp_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine InitModel_interface

        subroutine InitIO_interface(App, Xml)
            import BaseApp_T
            import XML_Input_T
            class(BaseApp_T), intent(inout) :: App
            type(XML_Input_T), intent(inout) :: Xml
        end subroutine InitIO_interface

        subroutine WriteRestart_interface(App, resFile)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
            character(len=*), intent(in) :: resFile
        end subroutine WriteRestart_interface

        subroutine ReadRestart_interface(App, resFile)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
            character(len=*), intent(in) :: resFile
        end subroutine ReadRestart_interface

        subroutine WriteConsoleOutput_interface(App)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
        end subroutine WriteConsoleOutput_interface

        subroutine WriteFileOutput_interface(App, outFile, nStep)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
            character(len=*), intent(in) :: outFile
            integer, intent(in) :: nStep
        end subroutine WriteFileOutput_interface

        subroutine WriteSlimFileOutput_interface(App, outFile, nStep)
            import BaseApp_T
            class(BaseApp_T), intent(inout) :: App
            character(len=*), intent(in) :: outFile
            integer, intent(in) :: nStep
        end subroutine WriteSlimFileOutput_interface

        subroutine AdvanceModel_interface(App, dt)
            import BaseApp_T
            import rp
            class(BaseApp_T), intent(inout) :: App
            real(rp), intent(in) :: dt
        end subroutine AdvanceModel_interface
    end interface

end module

