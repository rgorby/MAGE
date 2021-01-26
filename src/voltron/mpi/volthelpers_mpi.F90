! Collection of data and objects for additional voltron helper ranks

module volthelpers_mpi
    use voltapp_mpi
    use mpi
    use ebsquish, only : SquishBlocksRemain, DoSquishBlock
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN, IEEE_QUIET_NAN

    implicit none

    enum, bind(C)
        enumerator :: VHSTEP=1,VHSQUISH
    endenum

    contains

    subroutine vhReqStep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call mpi_bcast(VHSTEP, 1, MPI_INT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
    end subroutine

    subroutine vhHandleStep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
    end subroutine

    subroutine vhReqSquish(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call mpi_bcast(VHSQUISH, 1, MPI_INT, 0, vApp%vHelpComm, ierr)
    end subroutine

    subroutine vhHandleSquish(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp
    end subroutine

end module

