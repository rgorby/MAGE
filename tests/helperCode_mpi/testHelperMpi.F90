module testHelperMpi

    ! limiting mpi namespace pollution from pFUnit
    use testHelper
    use pFUnit, only: MpiTestMethod
    use mpi_f08

    implicit none

    contains

    ! this returns a fortran-2008 syntax version of the pFUnit mpi communicator
    function getMpiF08Communicator(mtm) result(comm)
        class(MpiTestMethod), intent(in) :: mtm
        type(MPI_Comm) :: comm
        comm%MPI_VAL = mtm%getMpiCommunicator()
    end function getMpiF08COmmunicator

end module testHelperMpi

