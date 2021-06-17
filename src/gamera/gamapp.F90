! Main data objects and functions to perform a gamera simulation

module gamapp
    use gamtypes
    use step
    use init
    use mhdgroup

    implicit none

    type gamApp_T
        type(Model_T)  :: Model
        type(Grid_T)   :: Grid
        type(State_T)  :: State, oState
        type(Solver_T) :: Solver

    end type gamApp_T

    contains

    subroutine initGamera(gameraApp, userInitFunc, optFilename,doIO)
        type(gamApp_T), intent(inout) :: gameraApp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        character(len=*), optional, intent(in) :: optFilename
        logical, optional, intent(in) :: doIO

        character(len=strLen) :: inpXML,kaijuRoot
        type(XML_Input_T) :: xmlInp
        logical :: doIOX

        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
        else
            !Find input deck
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

        if (present(doIO)) then
            doIOX = doIO
        else
            doIOX = .true.
        endif

        !Create XML reader
        write(*,*) 'Reading input deck from ', trim(inpXML)
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Gamera',.true.)

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
        call Hatch(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Solver,xmlInp,userInitFunc)
        call cleanClocks()

        if (doIOX) then
            if (.not. gameraApp%Model%isRestart) call fOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
            call consoleOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        
    end subroutine initGamera

    subroutine stepGamera(gameraApp)
        type(gamApp_T), intent(inout) :: gameraApp

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

        call Tic("Gamera")
        !Advance system
        call AdvanceMHD(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Solver,gameraApp%Model%dt)
        call Toc("Gamera")

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

end module gamapp

