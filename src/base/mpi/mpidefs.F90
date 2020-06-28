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

  ! helper function to print info about custom MPI datatypes
  recursive function printDataType(datatype) result(retVal)
      integer, intent(in) :: datatype
      logical :: retVal

      integer :: numInts, numAdds, numDTs, combiner, ierr, i
      integer, dimension(:), allocatable :: arrayInts, arrayDTs
      integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: arrayAdds
      call mpi_type_get_envelope(datatype, numInts, numAdds, numDTs, combiner, ierr)

      SELECT CASE(combiner)
          CASE (MPI_COMBINER_NAMED)
              SELECT CASE (datatype)
                  CASE (MPI_INT)
                      write (*,*) 'Datatype is named: MPI_INT'
                  CASE (MPI_FLOAT)
                      write (*,*) 'Datatype is named: MPI_FLOAT'
                  CASE (MPI_DOUBLE)
                      write (*,*) 'Datatype is named: MPI_DOUBLE'
                  CASE (MPI_DOUBLE_PRECISION)
                      write (*,*) 'Datatype is named: MPI_DOUBLE_PRECISION'
                  CASE DEFAULT
                      write (*,*) 'Unhandled base named datatype in printDataType'
              ENDSELECT
              retVal = .false.
              RETURN
          CASE (MPI_COMBINER_STRUCT)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is struct containing ',arrayInts(1),' datatypes:'
              do i=1,arrayInts(1)
                  write (*,*) 'blocklength ',arrayInts(i+1),', displacement ',arrayAdds(i),' type:'
                  if (printDataType(arrayDTs(i))) then
                      call mpi_type_free(arrayDTs(i),ierr)
                  endif
              enddo
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE (MPI_COMBINER_HINDEXED)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is hindexed containing ',arrayInts(1),' blocks:'
              do i=1,arrayInts(1)
                  write (*,*) 'blocklength ',arrayInts(i+1),' displacement ',arrayAdds(i)
              enddo
              write (*,*) 'hindexed datatype:'
              if (printDataType(arrayDTs(1))) then
                  call mpi_type_free(arrayDTs(1), ierr)
              endif
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE (MPI_COMBINER_HVECTOR)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is hvector containing ',arrayInts(1),' blocks:'
              write (*,*) 'blocklength ',arrayInts(2),' displacement ',arrayAdds(1),' type:'
              if (printDataType(arrayDTs(1))) then
                  call mpi_type_free(arrayDTs(1), ierr)
              endif
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE (MPI_COMBINER_CONTIGUOUS)
              allocate(arrayInts(numInts))
              allocate(arrayAdds(numAdds))
              allocate(arrayDTs(numDTs))
              call mpi_type_get_contents(datatype, numInts, numAdds, numDTs, &
                                         arrayInts, arrayAdds, arrayDTs, ierr)
              write (*,*) 'Datatype is contiguous conaining ',arrayInts(1),' repetitions of datatype:'
              if (printDataType(arrayDTs(1))) then
                  call mpi_type_free(arrayDTs(1), ierr)
              endif
              deallocate(arrayInts)
              deallocate(arrayAdds)
              deallocate(arrayDTs)
          CASE DEFAULT
              write (*,*) 'Unknown combiner type in printDataType'
              retVal = .false.
              RETURN
      ENDSELECT

      retVal = .true.

  end function printDataType

end module mpidefs
