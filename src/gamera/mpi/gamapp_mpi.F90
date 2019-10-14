! Main data objects and functions to perform a gamera simulation

module gamapp_mpi
    use gamtypes
    use step
    use init
    use mhdgroup
    use gamapp
    use mpi

    implicit none

    type gamAppMpi_T extends gamApp_T
        integer :: NumRi=1,NumRj=1,NumRk=1
        integer :: Ri=1,Rj=1,Rk=1
        type(MPI_Comm) :: gamMpiComm
    end type gamAppMpi_T

    contains

    subroutine initGamera_mpi(gamAppMpi, userInitFunc, optFilename)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        
        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
            call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")
        else
            !Find input deck
            call getIDeckStr(inpXML)

        endif

        write(*,*) 'Reading input deck from ', trim(inpXML)
        inquire(file=inpXML,exist=fExist)
        if (.not. fExist) then
            write(*,*) 'Error opening input deck, exiting ...'
            write(*,*) ''
            stop
        endif
        
        !Create XML reader
        xmlInp = New_XML_Input(trim(inpXML),'Gamera',.true.)

        !Initialize Grid/State/Model (Hatch Gamera)
        !Will enforce 1st BCs, caculate 1st timestep, set oldState
        call Hatch_mpi(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Solver,xmlInp,userInitFunc)
        call cleanClocks()

        if (.not. gameraApp%Model%isRestart) call fOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call consoleOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)

    end subroutine initGamera_mpi

    subroutine stepGamera_mpi(gameraApp)
        type(gamApp_T), intent(inout) :: gameraApp

        call Tic("Gamera")
        !Advance system
        call AdvanceMHD(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Solver,gameraApp%Model%dt)
        call Toc("Gamera")

        !Update info
        gameraApp%Model%ts = gameraApp%Model%ts+1
        gameraApp%Model%t = gameraApp%Model%t+gameraApp%Model%dt

        !Call user-defined per-step function
        if (associated(gameraApp%Model%HackStep)) then
            call Tic("HackStep")
            call gameraApp%Model%HackStep(gameraApp%Model,gameraApp%Grid,gameraApp%State)
            call Toc("HackStep")
        endif

        !Update ghost cells
        call Tic("Halos")
        call HaloUpdate(gameraApp)
        call Toc("Halos")

        !Calculate new timestep
        call Tic("DT")
        gameraApp%Model%dt = CalcDT(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call Toc("DT")

        !Enforce BCs
        call Tic("Halos")
        call EnforceBCs(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call Toc("Halos")


        !Output if necessary
        call Tic("IO")
        if (modulo(gameraApp%Model%ts,gameraApp%Model%tsOut) ==0) then
            call consoleOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        if (gameraApp%Model%t >= gameraApp%Model%tOut) then
            call fOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        if (gameraApp%Model%doResOut .and. (gameraApp%Model%t >= gameraApp%Model%tRes)) then
            !print *,"RESTART :: ", Model%doResOut, Model%t >= Model%tRes
            call resOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        call Toc("IO")

    end subroutine stepGamera_mpi

    subroutine haloUpdate(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
    end subroutine haloUpdate

end module gamapp_mpi

