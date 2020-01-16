! Collection of data and objects for the voltron middle man
! MPI version

module voltapp_mpi
    use voltapp
    use gamapp_mpi
    use gamapp
    use mpi
    
    implicit none

    type, extends(voltApp_T) :: voltAppMpi_T
        integer :: voltMpiComm = MPI_COMM_NULL
        integer :: myRank
        type(gamApp_T) :: gAppLocal

        ! array of all zeroes to simplify various send/receive calls
        integer, dimension(:), allocatable :: zeroArray

        ! Shallow Gas Data Receive Variables
        integer, dimension(:), allocatable :: recvCountsGasShallow, recvTypesGasShallow
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: recvDisplsGasShallow
        ! Shallow Magnetic Field Data Receive Variables
        integer, dimension(:), allocatable :: recvCountsBxyzShallow, recvTypesBxyzShallow
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: recvDisplsBxyzShallow
        ! Shallow inEijk Data Send Variables
        integer, dimension(:), allocatable :: sendCountsIneijkShallow, sendTypesIneijkShallow
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsIneijkShallow
        ! Shallow inExyz Data Send Variables
        integer, dimension(:), allocatable :: sendCountsInexyzShallow, sendTypesInexyzShallow
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsInexyzShallow

        ! Deep Gas Data Receive Variables
        integer, dimension(:), allocatable :: recvCountsGasDeep, recvTypesGasDeep
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: recvDisplsGasDeep
        ! Deep Magnetic Field Data Receive Variables
        integer, dimension(:), allocatable :: recvCountsBxyzDeep, recvTypesBxyzDeep
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: recvDisplsBxyzDeep
        ! Deep Gas0 Data Send Variables
        integer, dimension(:), allocatable :: sendCountsGas0Deep, sendTypesGas0Deep
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsGas0Deep

    end type voltAppMpi_T

    contains

    !Initialize Voltron (after Gamera has already been initialized)
    subroutine initVoltron_mpi(vApp, userInitFunc, voltComm, optFilename)
        type(voltAppMpi_T), intent(inout) :: vApp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        integer, intent(in) :: voltComm
        character(len=*), optional, intent(in) :: optFilename

        integer :: commSize, ierr, numNeighbors, numCells, length, ic
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder
        integer, allocatable, dimension(:) :: neighborRanks, inData, outData
        integer, allocatable, dimension(:) :: iRanks, jRanks, kRanks

        ! create a new communicator using MPI Topology
        call MPI_Comm_Size(voltComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        allocate(neighborRanks(1:commSize-1))
        allocate(inData(1:commSize-1))
        allocate(outData(1:commSize-1))
        allocate(iRanks(1:commSize-1))
        allocate(jRanks(1:commSize-1))
        allocate(kRanks(1:commSize-1))

        allocate(vApp%zeroArray(1:commSize-1))
        vApp%zeroArray(:) = 0

        ! doing a very very rough approximation of data transferred to help MPI reorder
        ! for deep updates, assume each rank sends data equal to its # physical cells
        ! for shallow updates, i=0 ranks send that much data again

        ! get i/j/k ranks from each Gamera mpi rank
        call mpi_gather(MPI_IN_PLACE, 0, 0, iRanks, commSize-1, MPI_INT, commSize-1, voltComm, ierr)
        call mpi_gather(MPI_IN_PLACE, 0, 0, jRanks, commSize-1, MPI_INT, commSize-1, voltComm, ierr)
        call mpi_gather(MPI_IN_PLACE, 0, 0, kRanks, commSize-1, MPI_INT, commSize-1, voltComm, ierr)

        ! get the number of physical cells from rank 0
        call mpi_recv(numCells, 1, MPI_INT, 0, 97500, voltComm, MPI_STATUS_IGNORE, ierr)

        do ic=1,commSize-1
            neighborRanks(ic) = ic-1
            if(iRanks(ic) == 0) then
                ! i=0 rank, does shallow update
                inData(ic) = numCells
                outData(ic) = numCells
            else
                ! i != 0 rank, does not do shallow update
                inData(ic) = 0
                outData(ic) = 0
            endif
            inData(ic) = inData(ic) + numCells
            outData(ic) = outData(ic) + numCells
        enddo

        reorder = .true. ! allow MPI to reorder the ranks
        call mpi_dist_graph_create_adjacent(voltComm, &
            numNeighbors,neighborRanks,inData, &
            numNeighbors,neighborRanks,outData, &
            MPI_INFO_NULL, reorder, vApp%voltMpiComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Comm_rank(vApp%voltMpiComm, vApp%myRank, ierr)

        ! get i/j/k ranks again in case MPI ranks were reordered in the new communicator
        call mpi_gather(MPI_IN_PLACE, 0, 0, iRanks, commSize-1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_gather(MPI_IN_PLACE, 0, 0, jRanks, commSize-1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_gather(MPI_IN_PLACE, 0, 0, kRanks, commSize-1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)

        ! create the MPI datatypes needed to transfer state data
        call createVoltDataTypes(vApp)

        ! create a local Gamera object which contains the entire domain
        if(present(optFilename)) then
            call initGamera(vApp%gAppLocal,userInitFunc,optFilename,doIO=.false.)
        else
            call initGamera(vApp%gAppLocal,userInitFunc,doIO=.false.)
        endif

        ! use standard voltron with local gamApp object
        if(present(optFilename)) then
            call initVoltron(vApp, vApp%gAppLocal, optFilename)
        else
            call initVoltron(vApp, vApp%gAppLocal)
        endif

        ! send all of the initial voltron parameters to the gamera ranks
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%DeepT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%ShallowT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%MJD, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%ts, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%doDeep, 1, MPI_LOGICAL, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%gAppLocal%Model%MJD0, 1, MPI_MYFLOAT, vAPp%myRank, vApp%voltMpiComm, ierr)

        deallocate(neighborRanks, inData, outData, iRanks, jRanks, kRanks)

    end subroutine initVoltron_mpi

    ! MPI version of updating voltron variables
    subroutine stepVoltron_mpi(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        real(rp) :: g_t
        integer :: g_ts, ierr

        ! get gApp%Model%t,ts from gamera. All ranks have the same, just receive from one of them
        call mpi_recv(g_t, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97600, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)
        call mpi_recv(g_ts, 1, MPI_INT, MPI_ANY_SOURCE, 97700, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)

        vApp%gAppLocal%Model%t = g_t
        vApp%gAppLocal%Model%ts = g_ts

        call stepVoltron(vApp, vApp%gAppLocal)

        ! send vApp%time,MJD,ts to all gamera ranks
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%MJD, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%ts, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)

    end subroutine stepVoltron_mpi

!----------
!Shallow coupling stuff
    subroutine ShallowUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp

        real(rp) :: time
        integer :: ierr

        ! fetch data from Gamera ranks
        ! Receive Shallow Gas Data
        call mpi_neighbor_alltoallw(0, vApp%zeroArray, &
                                    vApp%zeroArray, vApp%zeroArray, &
                                    vApp%gAppLocal%State%Gas, vApp%recvCountsGasShallow, &
                                    vApp%recvDisplsGasShallow, vApp%recvTypesGasShallow, &
                                    vApp%voltMpiComm, ierr)
        ! Receive Shallow Bxyz Data
        call mpi_neighbor_alltoallw(0, vApp%zeroArray, &
                                    vApp%zeroArray, vApp%zeroArray, &
                                    vApp%gAppLocal%State%Bxyz, vApp%recvCountsBxyzShallow, &
                                    vApp%recvDisplsBxyzShallow, vApp%recvTypesBxyzShallow, &
                                    vApp%voltMpiComm, ierr)

        ! call base update function with local data
        call ShallowUpdate(vApp, vApp%gAppLocal, time)

        ! send updated data to Gamera ranks
        ! voltron updates inEijk and inExyz in the IonInnerBC_T
        ! find the remix BC to read data from
        SELECT type(iiBC=>vApp%gAppLocal%Grid%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                ! Send Shallow inEijk Data
                call mpi_neighbor_alltoallw(iiBC%inEijk, vApp%sendCountsIneijkShallow, &
                                            vApp%sendDisplsIneijkShallow, vApp%sendTypesIneijkShallow, &
                                            0, vApp%zeroArray, vApp%zeroArray, vApp%zeroArray, &
                                            vApp%voltMpiComm, ierr)

                ! Send Shallow inExyz Data
                call mpi_neighbor_alltoallw(iiBC%inExyz, vApp%sendCountsInexyzShallow, &
                                            vApp%sendDisplsInexyzShallow, vApp%sendTypesInexyzShallow, &
                                            0, vApp%zeroArray, vApp%zeroArray, vApp%zeroArray, &
                                            vApp%voltMpiComm, ierr)
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in Voltron MPI ShallowUpdate_mpi'
                stop
        END SELECT

        ! send next time for shallow calculation to all gamera ranks
        call mpi_bcast(vApp%ShallowT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)

    end subroutine ShallowUpdate_mpi

!----------
!Deep coupling stuff (time coming from vApp%time, so in seconds)
    subroutine DeepUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv
        integer :: ierr

        if (.not. vApp%doDeep) then
            !Why are you even here?
            return
        endif

        ! fetch data from Gamera ranks
        ! Receive Deep Gas Data
        call mpi_neighbor_alltoallw(0, vApp%zeroArray, vApp%zeroArray, vApp%zeroArray, &
                                    vApp%gAppLocal%State%Gas, vApp%recvCountsGasDeep, &
                                    vApp%recvDisplsGasDeep, vApp%recvTypesGasDeep, &
                                    vApp%voltMpiComm, ierr)
        ! Receive Deep Bxyz Data
        call mpi_neighbor_alltoallw(0, vApp%zeroArray, vApp%zeroArray, vApp%zeroArray, &
                                    vApp%gAppLocal%State%Bxyz, vApp%recvCountsBxyzDeep, &
                                    vApp%recvDisplsBxyzDeep, vApp%recvTypesBxyzDeep, &
                                    vApp%voltMpiComm, ierr)

        ! call base update function with local data
        call DeepUpdate(vApp, vApp%gAppLocal, time)

        ! send updated data to Gamera ranks
        ! Send Deep Gas0 Data
        call mpi_neighbor_alltoallw(vApp%gAppLocal%Grid%Gas0, vApp%sendCountsGas0Deep, &
                                    vApp%sendDisplsGas0Deep, vApp%sendTypesGas0Deep, &
                                    0, vApp%zeroArray, vApp%zeroArray, vApp%zeroArray, &
                                    vApp%voltMpiComm, ierr)

        ! send next time for deep calculation to all gamera ranks
        call mpi_bcast(vApp%DeepT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)

    end subroutine DeepUpdate_mpi

    subroutine createVoltDataTypes(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp
    end subroutine createVoltDataTypes

end module voltapp_mpi

