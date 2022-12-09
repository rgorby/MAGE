! define types for mpi specific voltron

module voltmpitypes
    use voltapp
    use mpi_f08
    use mpidefs

    implicit none

    type, extends(voltApp_T) :: voltAppMpi_T
        ! voltron to helpers comms variables
        type(MPI_Comm) :: vHelpComm
        integer :: vHelpRank
        real(rp) :: lastSquishTime
        type(MPI_Win) :: vHelpWin
        integer, dimension(:), allocatable :: vHelpIdle
        logical :: amHelper = .false., useHelpers = .false.
        logical :: doSquishHelp = .false., masterSquish = .false.
        logical :: squishLoadBalance = .true.

        ! voltron to gamera comms variables
        type(MPI_Comm) :: voltMpiComm
        integer :: myRank
        type(gamApp_T) :: gAppLocal
        logical :: doSerialVoltron = .false.
        logical :: firstDeepUpdate = .true., firstStepUpdate = .true.

        ! array of all zeroes to simplify various send/receive calls
        integer, dimension(:), allocatable :: zeroArrayCounts
        type(MPI_Datatype), dimension(:), allocatable :: zeroArrayTypes
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable ::  zeroArrayDispls

        ! list of gamera ranks to communicate with
        integer, dimension(:), allocatable :: sendRanks, recvRanks

        ! STEP VOLTRON VARIABLES
        type(MPI_Request) :: timeReq, timeStepReq
        real(rp) :: timeBuffer
        integer :: timeStepBuffer

        ! DEEP COUPLING VARIABLES
        integer, dimension(:), allocatable :: recvCountsGasDeep
        type(MPI_Datatype), dimension(:), allocatable :: recvTypesGasDeep
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: recvDisplsGasDeep
        integer, dimension(:), allocatable :: recvCountsBxyzDeep
        type(MPI_Datatype), dimension(:), allocatable :: recvTypesBxyzDeep
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: recvDisplsBxyzDeep
        integer, dimension(:), allocatable :: sendCountsGas0Deep
        type(MPI_Datatype), dimension(:), allocatable :: sendTypesGas0Deep
        integer(kind=MPI_AN_MYADDR), dimension(:), allocatable :: sendDisplsGas0Deep
        logical :: deepProcessingInProgress = .false.

    end type voltAppMpi_T

end module

