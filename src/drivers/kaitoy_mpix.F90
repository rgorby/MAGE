!Test driver for compilation, OpenMP, and MPI debugging

program kaitoy_mpix
    use kdefs
    use mpidefs
    use clocks
    use mpi_f08

    implicit none

    integer :: ierror, length, provided, worldSize, worldRank
    type(MPI_Comm) :: subComm
    integer :: required=MPI_THREAD_MULTIPLE
    character( len = MPI_MAX_ERROR_STRING) :: message
    integer(KIND=MPI_AN_MYADDR) :: tagMax
    logical :: tagSet

    integer :: i
    real(rp) :: testValue
    real(rp), dimension(:), allocatable :: testArray
    integer, dimension(:), allocatable :: otherRanks
    integer, dimension(:), allocatable :: adjacency
    integer, dimension(:), allocatable :: zeroCounts
    integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: zeroDispls
    type(MPI_Datatype), dimension(:), allocatable :: zeroTypes
    integer, dimension(:), allocatable :: transmitCounts
    integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: transmitDispls
    type(MPI_Datatype), dimension(:), allocatable :: transmitTypes

    ! initialize MPI
    !Set up MPI with or without thread support
#ifdef _OPENMP
    print *,"MPI + OpenMP"
    call MPI_INIT_THREAD(required,provided,ierror)
    if(provided < required) then
        print *,"Not support for MPI_THREAD_MULTIPLE, aborting!"
        call abort()
    end if
#else
    print *,"MPI without threading"
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
    if(worldRank == 0) then
        call printConfigStamp()
#ifdef USEMKL
        write(*,*) 'USEMKL is set ...'
#else
        write(*,*) 'USEMKL is not set ...'
#endif
        print *, 'Tag Upper-Bound = ', tagMax
        print *, 'This was called with ', worldSize, ' MPI ranks'
    endif

    if(worldSize == 1) then
        ! mpi comm only has the 1 rank
        print *, "Can't test MPI comms with only 1 MPI rank. Please launch this test with at least 2 MPI ranks."
        call printClocks()
        call MPI_Finalize(ierror)
        stop
    endif

    ! do a test broadcast
    if(worldRank == 0) then
        print *, 'Testing MPI Bcast'
        testValue = 1
    else
        testValue = 0
    endif
    call mpi_bcast(testValue, 1, MPI_MYFLOAT, 0, MPI_COMM_WORLD, ierror)
    if(testValue /= 1) print *, 'MPI broadcast failed on rank ', worldRank

    ! create a new MPI communicator and do neighborhood comms
    if(worldRank == 0) print *, 'Creating a neighborhood subcommunicator'
    if(worldRank == (worldSize-1)) then
        allocate(testArray(1:10*(worldSize-1)))
        allocate(otherRanks(1:(worldSize-1)))
        allocate(adjacency(1:(worldSize-1)))
        allocate(zeroCounts(1:(worldSize-1)))
        allocate(zeroDispls(1:(worldSize-1)))
        allocate(zeroTypes(1:(worldSize-1)))
        allocate(transmitCounts(1:(worldSize-1)))
        allocate(transmitDispls(1:(worldSize-1)))
        allocate(transmitTypes(1:(worldSize-1)))
        otherRanks = (/ (i,i=0,(worldSize-2)) /)
        adjacency(:) = 1
        transmitCounts(:) = 10
        transmitTypes(:) = MPI_MYFLOAT
        transmitDispls = 80*(/ (i,i=0,(worldSize-2)) /)
    else
        allocate(testArray(1:10))
        allocate(otherRanks(1:1))
        allocate(adjacency(1:1))
        allocate(zeroCounts(1:1))
        allocate(zeroDispls(1:1))
        allocate(zeroTypes(1:1))
        allocate(transmitCounts(1:1))
        allocate(transmitDispls(1:1))
        allocate(transmitTypes(1:1))
        otherRanks(1) = worldSize-1
        adjacency(1) = 1
        transmitCounts(1) = 10
        transmitTypes(1) = MPI_MYFLOAT
        transmitDispls(1) = 0
    endif
    zeroCounts(:) = 0
    zeroDispls(:) = 0
    zeroTypes(:) = MPI_INTEGER

    call mpi_dist_graph_create_adjacent(MPI_COMM_WORLD, &
            size(otherRanks),otherRanks,adjacency, &
            size(otherRanks),otherRanks,adjacency, &
            MPI_INFO_NULL, .false., subComm, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    if(worldRank == 0) print *, 'Gathering data on final rank'
    if(worldRank == (worldSize-1)) then
        testArray = 0
        call mpi_neighbor_alltoallw(0, zeroCounts, &
                            zeroDispls, zeroTypes, &
                            testArray, transmitCounts, &
                            transmitDispls, transmitTypes, &
                            subComm, ierror)
    else
        testArray = (worldRank+1)*(/ (i,i=0,9) /)
        call mpi_neighbor_alltoallw(testArray, transmitCounts, &
                            transmitDispls, transmitTypes, &
                            0, zeroCounts, &
                            zeroDispls, zeroTypes, &
                            subComm, ierror)
    endif
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    if(worldRank == (worldSize-1)) then
        if(sum(testArray) .le. 0) then
            print *, 'Neighborhood collection did not work'
            call abort()
        endif
        testArray = 2*testArray
    endif

    if(worldRank == 0) print *, 'Scattering data back to lower ranks'
    if(worldRank == (worldSize-1)) then
        call mpi_neighbor_alltoallw(testarray, transmitCounts, &
                            transmitDispls, transmitTypes, &
                            0, zeroCounts, &
                            zeroDispls, zeroTypes, &
                            subComm, ierror)
    else
        call mpi_neighbor_alltoallw(0, zeroCounts, &
                            zeroDispls, zeroTypes, &
                            testArray, transmitCounts, &
                            transmitDispls, transmitTypes, &
                            subComm, ierror)
    endif
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    if(worldRank < (worldSize-1)) then
        testArray = testArray - 2*(worldRank+1)*(/ (i,i=0,9) /)
        if(abs(sum(testArray)) > 1e-14) then
            print *, 'Rank ', worldRank, ' did not perform neighborhood comms correctly. Error = ', abs(sum(testArray))
            call abort()
        endif
    endif

    if(worldRank == 0) print *, 'Everything worked!'
    call printClocks()
    call MPI_FINALIZE(ierror)

end program kaitoy_mpix

