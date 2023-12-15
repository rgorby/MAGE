module volttypes_mpi
    use volttypes
    use mpi_f08
    use mpidefs

    implicit none

    ! mpi voltron specific options
    type, extends(BaseOptions_T) :: VoltOptionsMpi_T
        type(MPI_Comm) :: allComm
        type(MPI_Comm) :: allVoltComm

        contains
    end type voltOptionsMpi_T

    type, extends(voltApp_T) :: voltAppMpi_T
        ! voltron to helpers comms variables
        type(MPI_Comm) :: vHelpComm
        integer :: vHelpRank
        real(rp) :: lastSquishTime
        type(MPI_Win) :: vHelpWin
        integer, dimension(:), allocatable :: vHelpIdle
        logical :: amHelper = .false., useHelpers = .false.
        logical :: doSquishHelp = .false., masterSquish = .false.
        logical :: squishLoadBalance = .true., deepProcessingInProgress=.false.

        ! mpi voltron specific options
        type(VoltOptionsMpi_T) :: vOptionsMpi

    end type voltAppMpi_T

end module

