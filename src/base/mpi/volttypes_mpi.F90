module volttypes_mpi
    use volttypes
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
        logical :: squishLoadBalance = .true., deepProcessingInProgress=.false.

    end type voltAppMpi_T

end module

