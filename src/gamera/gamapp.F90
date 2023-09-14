! implementation of subroutines in GameraApp type

module gamapp
    use gamtypes
    use step
    use init
    use mhdgroup
    use output

    implicit none

    type, extends(BaseOptions_T) :: gamOptions_T
        procedure(StateIC_T), pointer, nopass :: userInitFunc

        contains
    end type gamOptions_T

    type, extends(BaseApp_T) :: gamApp_T
        type(Model_T)  :: Model
        type(Grid_T)   :: Grid
        type(State_T)  :: State, oState
        type(Solver_T) :: Solver

        type(gamOptions_T) :: gOptions

        contains

        procedure :: InitModel => gamInitModel
        procedure :: InitIO => gamInitIO
        procedure :: WriteRestart => gamWriteRestart
        procedure :: ReadRestart => gamReadRestart
        procedure :: WriteConsoleOutput => gamWriteConsoleOutput
        procedure :: WriteFileOutput => gamWriteFileOutput
        procedure :: WriteSlimFileOutput => gamWriteSlimFileOutput
        procedure :: AdvanceModel => gamAdvanceModel

    end type gamApp_T


    contains

    ! procedures for gamApp_T
    subroutine gamInitModel(App, Xml)
        class(gamApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call initGamera(App, Xml)

    end subroutine gamInitModel

    subroutine gamInitIO(App, Xml)
        class(gamApp_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call initGamIO(App%Model, Xml)

    end subroutine gamInitIO

    subroutine gamWriteRestart(App)
        class(gamApp_T), intent(inout) :: App

        call resOutput(App%Model, App%Grid, App%oState, App%State)

    end subroutine gamWriteRestart

    subroutine gamReadRestart(App, resId, nRes)
        class(gamApp_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        character(len=strLen) :: inH5

        call getRestart(App%Model,App%Grid,resId,nRes,inH5)
        call readH5Restart(App%Model, App%Grid, App%State, App%oState,inH5)

    end subroutine gamReadRestart

    subroutine gamWriteConsoleOutput(App)
        class(gamApp_T), intent(inout) :: App

        call consoleOutput(App%Model, App%Grid, App%State)

    end subroutine gamWriteConsoleOutput

    subroutine gamWriteFileOutput(App)
        class(gamApp_T), intent(inout) :: App

        call fOutput(App%Model, App%Grid, App%State)

    end subroutine gamWriteFileOutput

    subroutine gamWriteSlimFileOutput(App)
        class(gamApp_T), intent(inout) :: App

        call App%WriteFileOutput()

    end subroutine gamWriteSlimFileOutput

    subroutine gamAdvanceModel(App, dt)
        class(gamApp_T), intent(inout) :: App
        real(rp), intent(in) :: dt

        call stepGamera(App)

    end subroutine gamAdvanceModel

    ! implementation
    subroutine initGamIO(Model, xmlInp)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: xmlInp

        logical :: doFatIO

        !Output/Restart (IOCLOCK)
        call Model%IO%init(xmlInp,Model%t,Model%ts)
        call xmlInp%Set_Val(Model%doDivB ,'output/DivB'    ,.true. )
        call xmlInp%Set_Val(doFatIO      ,'output/doFatIO' ,.false.)
        if (doFatIO) then
            call SetFatIO()
        endif
    end subroutine

    subroutine initGamera(gameraApp, xmlInp)
        class(gamApp_T), intent(inout) :: gameraApp
        type(XML_Input_T), intent(inout) :: xmlInp

        character(len=strLen) :: kaijuRoot

        gameraApp%Model%isLoud = .true.

        call xmlInp%SetRootStr('Kaiju/Gamera')
        call xmlInp%SetVerbose(.true.)

        ! try to verify that the XML file has "Kaiju" as a root element
        kaijuRoot = ""
        call xmlInp%Get_Key_Val("/Gamera/sim/H5Grid",kaijuRoot,.false.)
        if(len(trim(kaijuRoot)) /= 0) then
            write(*,*) "The input XML appears to be of an old style."
            write(*,*) "As of June 12th, 2021 it needs a root element of <Kaiju>."
            write(*,*) "Please modify your XML config file by adding this line at the top:"
            write(*,*) "<Kaiju>"
            write(*,*) "and this line at the bottom:"
            write(*,*) "</Kaiju>"
            write(*,*) "OR (preferred) convert your configuration to an INI file and use"
            write(*,*) " the XMLGenerator.py script to create conforming XML files."
            write(*,*) "Please refer to the python script or"
            write(*,*) " the [Generating XML Files] wiki page for additional info."
            stop
        endif

        ! read debug flags
        call xmlInp%Set_Val(writeGhosts,"debug/writeGhosts",.false.)
        call xmlInp%Set_Val(writeMagFlux,"debug/writeMagFlux",.false.)

        !Initialize Grid/State/Model (Hatch Gamera)
        !Will enforce 1st BCs, caculate 1st timestep, set oldState
        call Hatch(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Solver,xmlInp,gameraApp%gOptions%userInitFunc)
        call cleanClocks()

        if (.not. gameraApp%Model%isSub) then
            if (.not. gameraApp%Model%isRestart) call fOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
            call consoleOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        
    end subroutine initGamera

    subroutine stepGamera(gameraApp)
        class(gamApp_T), intent(inout) :: gameraApp

        !update the state variables to the next timestep
        call UpdateStateData(gameraApp)

        !Calculate new timestep
        call Tic("DT")
        gameraApp%Model%dt = CalcDT(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call Toc("DT")

        !Enforce BCs
        call Tic("BCs")
        call EnforceBCs(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call Toc("BCs")
        
    end subroutine stepGamera

    subroutine UpdateStateData(gameraApp)
        class(gamApp_T), intent(inout) :: gameraApp

        call Tic("Gamera",.true.)
        !Advance system
        call AdvanceMHD(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Solver,gameraApp%Model%dt)
        call Toc("Gamera",.true.)

        !Call user-defined per-step function
        !NOTE: Do this before updating time
        if (associated(gameraApp%Model%HackStep)) then
            call Tic("HackStep")
            call gameraApp%Model%HackStep(gameraApp%Model,gameraApp%Grid,gameraApp%State)
            call Toc("HackStep")
        endif

        !Update info
        gameraApp%Model%ts = gameraApp%Model%ts+ 1
        gameraApp%Model%t  = gameraApp%Model%t + gameraApp%Model%dt

    end subroutine UpdateStateData

end module

