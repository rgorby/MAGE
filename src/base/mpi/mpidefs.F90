module mpidefs

  use mpi
  use kdefs, ONLY: sp,dp,rp

  implicit none

  integer, public :: MPI_MYFLOAT

contains

  subroutine setMpiReal()

    MPI_MYFLOAT = MPI_DOUBLE_PRECISION
    if(rp ==  sp) MPI_MYFLOAT = MPI_REAL

  end subroutine setMpiReal

end module mpidefs
