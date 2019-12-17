!Driver for mpi decomposed Gamera coupled with Voltron (remix only)

program voltron_mpix
    use clocks
    use gamapp_mpi
    use voltapp_mpi
    use voltio
    use uservoltic

    implicit none

    type(gamAppMpi_T) :: gApp
    type(voltAppMpi_T) :: vApp
    procedure(StateIC_T), pointer :: userInitFunc => initUser

    integer :: ierror, length
    integer :: required=MPI_THREAD_MULTIPLE
    integer :: provided
    character( len = MPI_MAX_ERROR_STRING) :: message

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
    
    ! initial stab at this. rank 0 will be both a Gamera rank and the Voltron rank

    call initGamera_mpi(gApp,userInitFunc,MPI_COMM_WORLD,doIO=.false.)
    call initVoltron_mpi(vApp, gApp, MPI_COMM_WORLD)

    do while (vApp%time < vApp%tFin)
        !Start root timer
        call Tic("Omega")
        
    !Advance Gamera MHD
        call stepGamera_mpi(gApp)

    !Do any updates to Voltron
        call stepVoltron(vApp,gApp)
        
    !Coupling    
        call Tic("DeepCoupling")
        if ( (vApp%time >= vApp%DeepT) .and. vApp%doDeep ) then
            call DeepUpdate_mpi(vApp, gApp, vApp%time)
        endif
        call Toc("DeepCoupling")

        call Tic("IonCoupling")
        if (vApp%time >= vApp%ShallowT) then
            call ShallowUpdate_mpi(vApp, gApp, vApp%time)
        endif
        call Toc("IonCoupling")
        
    !IO checks
        call Tic("IO")
        !Console output
        if (vApp%IO%doConsole(vApp%ts)) then
            !Using console output from Gamera
            call consoleOutputV(vApp,gApp)
        endif
        !Restart output
        if (vApp%IO%doRestart(vApp%time)) then
            call resOutputV(vApp,gApp)
        endif
        !Data output
        if (vApp%IO%doOutput(vApp%time)) then
            call fOutputV(vApp,gApp)
        endif

        call Toc("IO")

    !Do timing info
        if (vApp%IO%doTimer(vApp%ts)) then
            if (vApp%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif

        call Toc("Omega")

    end do

end program voltron_mpix

