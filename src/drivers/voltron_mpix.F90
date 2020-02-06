!Driver for mpi decomposed Gamera coupled with Voltron (remix only)

program voltron_mpix
    use clocks
    use gamapp_mpi
    use voltapp_mpi
    use gam2VoltComm_mpi
    use output
    use voltio
    use uservoltic
    use mpi

    implicit none

    ! this driver requires that each rank be either Gamera OR Voltron, never both
    logical :: isGamera = .true.

    ! gamera data
    type(gamAppMpi_T) :: gApp
    type(gam2VoltCommMpi_T) :: g2vComm

    ! voltron data
    type(voltAppMpi_T) :: vApp

    procedure(StateIC_T), pointer :: userInitFunc => initUser

    integer :: ierror, length, provided, worldSize, worldRank, gamComm
    integer :: required=MPI_THREAD_MULTIPLE
    character( len = MPI_MAX_ERROR_STRING) :: message
    logical :: deepPerformed

    ! initialize MPI
    !Set up MPI with or without thread support
#ifdef _OPENMP
    call MPI_INIT_THREAD(required,provided,ierror)
    if(provided < required) then
        print *,"Not support for MPI_THREAD_MULTIPLE, aborting!"
        call abort()
    end if
    print *,"MPI + OpenMP !!!"
#else
    print *," MPI without threading"
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

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = .true.
    
    ! initial stab at this. Voltron MUST have its own MPI rank. For simplicity, this is always
    !    the last MPI rank for now

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
    if(worldRank == (worldSize-1)) then
        ! voltron rank
        isGamera = .false.
        call MPI_Comm_Split(MPI_COMM_WORLD, MPI_UNDEFINED, worldRank, gamComm, ierror)
        if(ierror /= MPI_Success) then
            call MPI_Error_string( ierror, message, length, ierror)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
        end if
        if(gamComm /= MPI_COMM_NULL) then
            print *,"Voltron rank did not get a null communicator back from mpi_comm_split"
            call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
        endif
    else
        ! gamera rank
        isGamera = .true.
        call MPI_Comm_Split(MPI_COMM_WORLD, 0, worldRank, gamComm, ierror)
        if(ierror /= MPI_Success) then
            call MPI_Error_string( ierror, message, length, ierror)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
        end if
    endif

    if(isGamera) then
        call initGamera_mpi(gApp,userInitFunc,gamComm,doIO=.false.)
        call initGam2Volt(g2vComm,gApp,MPI_COMM_WORLD)

        do while (g2vComm%time < g2vComm%tFin)
            !Start root timer
            call Tic("Omega")
        
            !Advance Gamera MHD
            call stepGamera_mpi(gApp)

            !Do any updates to Voltron
            call performStepVoltron(g2vComm,gApp)
        
            !Coupling    
            deepPerformed = .false.
            call Tic("DeepCoupling")
            if ( (g2vComm%time >= g2vComm%DeepT) .and. g2vComm%doDeep ) then
                call performDeepUpdate(g2vComm, gApp)
                deepPerformed = .true.
            endif
            call Toc("DeepCoupling")

            call Tic("IonCoupling")
            if (g2vComm%time >= g2vComm%ShallowT) then
                !write (*,*) "Gamera performing shallow update at ",g2vComm%time
                call performShallowUpdate(g2vComm, gApp, deepPerformed)
            endif
            call Toc("IonCoupling")
        
            !IO checks
            call Tic("IO")
            !Console output
            if (gApp%Model%IO%doConsole(g2vComm%ts)) then
                !Using console output from Gamera
                if(gApp%Grid%Ri==0 .and. gApp%Grid%Rj==0 .and. gApp%Grid%Rk==0) then
                    call consoleOutput(gApp%Model, gApp%Grid, gApp%State)
                endif
            endif
            !Restart output
            if (gApp%Model%IO%doRestart(gApp%Model%t)) then
                call resOutput(gApp%Model, gApp%Grid, gApp%State)
            endif
            !Data output
            if (gApp%Model%IO%doOutput(gApp%Model%t)) then
                call fOutput(gApp%Model, gApp%Grid, gApp%State)
            endif

            call Toc("IO")

            !Do timing info
            if (gApp%Model%IO%doTimer(g2vComm%ts)) then
                if (gApp%Model%IO%doTimerOut .and. \
                  gApp%Grid%Ri==0 .and. gApp%Grid%Rj==0 .and. gApp%Grid%Rk==0) then
                    call printClocks()
                endif
                call cleanClocks()
            endif

            call Toc("Omega")

        end do

    else
        call initVoltron_mpi(vApp, userInitFunc, MPI_COMM_WORLD)

        ! voltron run loop
        do while (vApp%time < vApp%tFin)
            !Start root timer
            call Tic("Omega")

            !Do any updates to Voltron
            call stepVoltron_mpi(vApp)

            !Coupling
            deepPerformed = .false.
            call Tic("DeepCoupling")
            if ( (vApp%time >= vApp%DeepT) .and. vApp%doDeep ) then
                call DeepUpdate_mpi(vApp, vApp%time)
                deepPerformed = .true.
            endif
            call Toc("DeepCoupling")

            call Tic("IonCoupling")
            if (vApp%time >= vApp%ShallowT) then
                !write (*,*) "Voltron performing shallow update at ",vApp%time
                call ShallowUpdate_mpi(vApp, vApp%time, deepPerformed)
            endif
            call Toc("IonCoupling")

            !IO checks
            call Tic("IO")
            !Console output
            if (vApp%IO%doConsole(vApp%ts)) then
                call consoleOutputVOnly(vApp,vApp%gAppLocal%Model%MJD0)
            endif
            !Restart output
            if (vApp%IO%doRestart(vApp%time)) then
                call resOutputVOnly(vApp)
            endif
            !Data output
            if (vApp%IO%doOutput(vApp%time)) then
                call fOutputVOnly(vApp)
            endif

            call Toc("IO")

            !Do timing info
            if (vApp%IO%doTimer(vApp%ts)) then
                if (vApp%IO%doTimerOut) then
                    call printClocks()
                endif
                call cleanClocks()
            endif

            call Toc("Omega")

        end do

    endif

    call MPI_FINALIZE(ierror)
    print *,"This is the END!"

end program voltron_mpix

