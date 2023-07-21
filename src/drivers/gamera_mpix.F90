!Driver for Gamera with MPI decomposition

program gamera_mpix
    use gamapp_mpi
    use drivertypes
    use usergamic
    use mpidefs

    implicit none

    type(AppDriver_T) :: appDriver
    type(GameraAppMpi_T), allocatable :: gAppMpi

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


    ! create instance of gamera app, and perform any needed configuration
    allocate(gAppMpi)
    gAppMpi%gOptions%userInitFunc => initUser
    gAppMpi%gOptionsMpi%gamComm = MPI_COMM_WORLD

    ! use the single-app helper function to move the gamera app into the driver structure
    ! note that after this point, gamApp will no longer be allocated
    allocate(appDriver%appPointers(1))
    call move_alloc(gAppMpi, appDriver%appPointers(1)%p)

    ! run the gamera app
    call appDriver%RunApps()

    call MPI_FINALIZE(ierror)
    write(*,*) "Fin Mpi"

end program gamera_mpix

