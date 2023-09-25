submodule (gamtypes) gamtypessub
    use gamapp

    contains

    ! procedures for gamApp_T, expanded on from src/base/types/gamtypes.F90
    module subroutine gamInitModel(App, Xml)
        class(gamApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call initGamera(App, Xml)

    end subroutine gamInitModel

    module subroutine gamInitIO(App, Xml)
        class(gamApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call initGamIO(App%Model, Xml)

    end subroutine gamInitIO

    module subroutine gamWriteRestart(App)
        class(gamApp_T), intent(inout) :: App

        call resOutput(App%Model, App%Grid, App%oState, App%State)

    end subroutine gamWriteRestart

    module subroutine gamReadRestart(App, resId, nRes)
        class(gamApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        character(len=strLen) :: inH5

        call getRestart(App%Model,App%Grid,resId,nRes,inH5)
        call readH5Restart(App%Model, App%Grid, App%State, App%oState,inH5)

    end subroutine gamReadRestart

    module subroutine gamWriteConsoleOutput(App)
        class(gamApp_T), intent(inout) :: App

        call consoleOutput(App%Model, App%Grid, App%State)

    end subroutine gamWriteConsoleOutput

    module subroutine gamWriteFileOutput(App)
        class(gamApp_T), intent(inout) :: App

        call fOutput(App%Model, App%Grid, App%State)

    end subroutine gamWriteFileOutput

    module subroutine gamWriteSlimFileOutput(App)
        class(gamApp_T), intent(inout) :: App

        call App%WriteFileOutput()

    end subroutine gamWriteSlimFileOutput

    module subroutine gamAdvanceModel(App, dt)
        class(gamApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt

        call stepGamera(App)

    end subroutine gamAdvanceModel

end submodule

