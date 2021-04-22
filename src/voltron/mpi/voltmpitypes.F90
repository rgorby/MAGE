! define types for mpi specific voltron

module voltmpitypes
    use voltapp
    use mpi
    use mpidefs

    implicit none

    type, extends(voltApp_T) :: voltAppMpi_T
        ! voltron to helpers comms variables
        integer :: vHelpComm = MPI_COMM_NULL
        integer :: vHelpRank
        integer :: vHelpWin = MPI_WIN_NULL
        integer, dimension(:), allocatable :: vHelpIdle
        logical :: amHelper = .false., useHelpers = .false.
        logical :: doSquishHelp = .false., masterSquish = .false.

        ! voltron to gamera comms variables
        integer :: voltMpiComm = MPI_COMM_NULL
        integer :: myRank
        type(gamApp_T) :: gAppLocal
        logical :: doSerialVoltron = .false., doAsyncShallow = .true.
        logical :: firstShallowUpdate = .true., firstDeepUpdate = .true., firstStepUpdate = .true.

        ! array of all zeroes to simplify various send/receive calls
        integer, dimension(:), allocatable :: zeroArrayCounts, zeroArrayTypes
        integer(kind=MPI_MYADDRESS), dimension(:), allocatable ::  zeroArrayDispls

        ! list of gamera ranks to communicate with
        integer, dimension(:), allocatable :: sendRanks, recvRanks

        ! STEP VOLTRON VARIABLES
        integer :: timeReq=MPI_REQUEST_NULL, timeStepReq=MPI_REQUEST_NULL
        real(rp) :: timeBuffer
        integer :: timeStepBuffer

        ! SHALLOW COUPLING VARIABLES
        integer, dimension(:), allocatable :: recvCountsGasShallow, recvTypesGasShallow
        integer(kind=MPI_MYADDRESS), dimension(:), allocatable :: recvDisplsGasShallow
        integer, dimension(:), allocatable :: recvCountsBxyzShallow, recvTypesBxyzShallow
        integer(kind=MPI_MYADDRESS), dimension(:), allocatable :: recvDisplsBxyzShallow
        integer, dimension(:), allocatable :: sendCountsIneijkShallow, sendTypesIneijkShallow
        integer(kind=MPI_MYADDRESS), dimension(:), allocatable :: sendDisplsIneijkShallow
        integer, dimension(:), allocatable :: sendCountsInexyzShallow, sendTypesInexyzShallow
        integer(kind=MPI_MYADDRESS), dimension(:), allocatable :: sendDisplsInexyzShallow
        ! SHALLOW ASYNCHRONOUS VARIABLES
        integer :: shallowIneijkSendReq=MPI_REQUEST_NULL, shallowInexyzSendReq=MPI_REQUEST_NULL
        integer :: asyncShallowBcastReq=MPI_REQUEST_NULL

        ! DEEP COUPLING VARIABLES
        integer, dimension(:), allocatable :: recvCountsGasDeep, recvTypesGasDeep
        integer(kind=MPI_MYADDRESS), dimension(:), allocatable :: recvDisplsGasDeep
        integer, dimension(:), allocatable :: recvCountsBxyzDeep, recvTypesBxyzDeep
        integer(kind=MPI_MYADDRESS), dimension(:), allocatable :: recvDisplsBxyzDeep
        integer, dimension(:), allocatable :: sendCountsGas0Deep, sendTypesGas0Deep
        integer(kind=MPI_MYADDRESS), dimension(:), allocatable :: sendDisplsGas0Deep
        logical :: deepProcessingInProgress = .false.

    end type voltAppMpi_T

end module

