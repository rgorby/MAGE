!Driver for mpi decomposed Gamera coupled with Voltron (remix only)

program voltron_mpix
    use clocks
    use gamCouple_mpi_G2V
    use voltapp_mpi
    use output
    use voltio
    use uservoltic
    use mpi_f08
    use xml_input

    implicit none

    ! create allocatable objects for the possible app types
    ! only one will be used by each executable
    type(voltAppMpi_T), allocatable :: vApp
    type(gamCouplerMpi_gam_T), allocatable :: gApp

    procedure(StateIC_T), pointer :: userInitFunc => initUser
    integer :: ierror, length, provided, worldSize, worldRank, numHelpers
    type(MPI_Comm) :: voltComm
    integer :: required=MPI_THREAD_MULTIPLE
    character( len = MPI_MAX_ERROR_STRING) :: message
    character(len=strLen) :: inpXML, helpersBuf
    logical :: useHelpers, helperQuit
    integer(KIND=MPI_AN_MYADDR) :: tagMax
    logical :: tagSet
    type(XML_Input_T) :: xmlInp
    real(rp) :: nextDT

    ! initialize MPI
    !Set up MPI with or without thread support
#ifdef _OPENMP
    call MPI_INIT_THREAD(required,provided,ierror)
    if(provided < required) then
        print *,"Not support for MPI_THREAD_MULTIPLE, aborting!"
        call abort()
    end if
    !print *,"MPI + OpenMP !!!"
#else
    !print *," MPI without threading"
    call MPI_INIT(ierror)
#endif

    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    ! initialize mpi data type
    call setMpiReal()

    call initClocks()

    call mpi_comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tagMax, tagSet, ierror)
    !print *, 'Tag Upper-Bound = ', tagMax

    call getIDeckStr(inpXML)
    call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

    !Create XML reader
    write(*,*) 'Reading input deck from ', trim(inpXML)
    xmlInp = New_XML_Input(trim(inpXML),'Kaiju',.true.)

    ! need to know how many voltron helpers there are
    call getIDeckStr(inpXML)
    call ReadXmlImmediate(trim(inpXML),'/Kaiju/Voltron/Helpers/useHelpers',helpersBuf,'F',.false.)
    read(helpersBuf,*) useHelpers
    call ReadXmlImmediate(trim(inpXML),'/Kaiju/Voltron/Helpers/numHelpers',helpersBuf,'0',.false.)
    read(helpersBuf,*) numHelpers
    if(.not. useHelpers) numHelpers = 0

    ! create a new MPI communicator for just Gamera
    !    for now this is always all ranks excep the last one (which is reserved for voltron)
    call MPI_Comm_Size(MPI_COMM_WORLD, worldSize, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    call MPI_Comm_Rank(MPI_COMM_WORLD, worldRank, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    if(worldRank .ge. (worldSize-1-numHelpers)) then
        ! voltron rank, including helpers (FOR NOW)
        allocate(vApp)
        call MPI_Comm_Split(MPI_COMM_WORLD, voltId, worldRank, voltComm, ierror)
        if(ierror /= MPI_Success) then
            call MPI_Error_string( ierror, message, length, ierror)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
        end if

        ! voltron
        vApp%vOptions%gamUserInitFunc => initUser
        vApp%vOptionsMpi%allComm = MPI_COMM_WORLD
        vApp%vOptionsMpi%allVoltComm = voltComm
        call initVoltron_mpi(vApp)

        if(vApp%amHelper) then
            ! do helper loop
            helperQuit = .false.
            do while(.not. helperQuit)
                call Tic("Omega", .true.)
                call helpVoltron(vApp, helperQuit)
                call Toc("Omega", .true.)
            end do
        else
            ! voltron run loop
            do while (vApp%time < vApp%tFin)
                !Start root timer
                call Tic("Omega", .true.)

                !Advance Voltron models one coupling step
                call Tic("StepVoltron")
                call stepVoltron_mpi(vApp)
                call Toc("StepVoltron")

                !IO checks
                call Tic("IO", .true.)
                !Console output
                if (vApp%IO%doConsole(vApp%ts)) then
                    call consoleOutputV(vApp,vApp%gApp)
                    if (vApp%IO%doTimerOut) then
                        call printClocks()
                    endif
                    call cleanClocks()
                elseif (vApp%IO%doTimer(vApp%ts)) then
                    if (vApp%IO%doTimerOut) then
                        call printClocks()
                    endif
                    call cleanClocks()
                endif

                !Restart output
                if (vApp%IO%doRestart(vApp%time)) then
                    call resOutputV(vApp,vApp%gApp)
                endif
                !Data output
                if (vApp%IO%doOutput(vApp%time)) then
                    call fOutputV(vApp,vApp%gApp)
                endif

                call Toc("IO", .true.)
                call Toc("Omega", .true.)
            enddo
            ! if using helpers, tell them to quit
            if(vApp%useHelpers) call vhReqHelperQuit(vApp)
        endif ! end all voltron loops

        call endVoltronWaits(vApp)
    else
        ! gamera rank
        allocate(gApp)
        call MPI_Comm_Split(MPI_COMM_WORLD, gamId, worldRank, gApp%gOptionsMpi%gamComm, ierror)
        if(ierror /= MPI_Success) then
            call MPI_Error_string( ierror, message, length, ierror)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
        end if

        gApp%gOptionsCplMpiG%allComm = MPI_COMM_WORLD
        gApp%gOptionsMpi%doIO = .false.
        gApp%gOptions%userInitFunc => initUser

        call gApp%InitModel(xmlInp)
        call gApp%InitIO(xmlInp)

        do while (gApp%Model%t < gApp%Model%tFin)
            call Tic("Omega") !Start root timer

            !Step model
            nextDT = min(gApp%Model%tFin-gApp%Model%t, gApp%Model%IO%nextIOTime(gApp%Model%t, gApp%Model%ts, gApp%Model%dt)-gApp%Model%t)
            call gApp%AdvanceModel(nextDT)

            !Output if necessary
            call Tic("IO")

            if (gApp%Model%IO%doConsole(gApp%Model%ts)) then
                call gApp%WriteConsoleOutput()

                !Timing info
                if (gApp%Model%IO%doTimerOut) call printClocks()
                call cleanClocks()

            elseif (gApp%Model%IO%doTimer(gApp%Model%ts)) then
                if (gApp%Model%IO%doTimerOut) call printClocks()
                call cleanClocks()
            endif

            if (gApp%Model%IO%doOutput(gApp%Model%t)) then
                call gApp%WriteFileOutput(gApp%Model%IO%nOut)
            endif

            if (gApp%Model%IO%doRestart(gApp%Model%t)) then
                call gApp%WriteRestart(gApp%Model%IO%nRes)
            endif

            call Toc("IO")
            call Toc("Omega")
        end do

    endif

    call MPI_FINALIZE(ierror)
    write(*,*) "Fin Voltron Mpi"
    

end program voltron_mpix

