!Driver for Gamera with MPI decomposition

program gamera_mpix
    use clocks
    use gamapp_mpi
    use usergamic
    use mpidefs

    implicit none

    type(gamAppMpi_T) :: gameraAppMpi
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

    !call printConfigStamp()
    call initClocks()

    gameraAppMpi%Model%isLoud = .true.
    
    ! this is a gamera-only application so all MPI ranks are gamera-only processes
    call initGamera_mpi(gameraAppMpi,userInitFunc,MPI_COMM_WORLD)

    do while (gameraAppMpi%Model%t < gameraAppMpi%Model%tFin)
        call Tic("Omega")
        !Start root timer

        call stepGamera_mpi(gameraAppMpi)

        !Output if necessary
        call Tic("IO")
        
        if (gameraAppMpi%Model%IO%doConsole(gameraAppMpi%Model%ts)) then
            call consoleOutput_mpi(gameraAppMpi)
        endif

        if (gameraAppMpi%Model%IO%doOutput(gameraAppMpi%Model%t)) then
            call fOutput(gameraAppMpi%Model,gameraAppMpi%Grid,gameraAppMpi%State)
        endif

        if (gameraAppMpi%Model%IO%doRestart(gameraAppMpi%Model%t)) then
            call resOutput(gameraAppMpi%Model,gameraAppMpi%Grid,gameraAppMpi%oState,gameraAppMpi%State)
        endif

        call Toc("IO")

    !Do timing info
        if (gameraAppMpi%Model%IO%doTimer(gameraAppMpi%Model%ts)) then
            if(gameraAppMpi%Model%IO%doTimerOut .and. &
               gameraAppMpi%Grid%Ri==0 .and. gameraAppMpi%Grid%Rj==0 .and. gameraAppMpi%Grid%Rk==0) then
                call printClocks()
            endif
            call cleanClocks()
        endif

        call Toc("Omega")
    end do

    call MPI_FINALIZE(ierror)
    write(*,*) "Fin"

end program gamera_mpix

