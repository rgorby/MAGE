module couplingHelpers
    use mpi_f08
    use mpidefs

    implicit none

    contains

    ! this helper routine wait its turn for Voltron to couple with this app
    ! participating in splitting until it is its turn
    subroutine appWaitForVoltronSplit(startingComm, appId, key, voltCoupledComm)
        type(MPI_Comm), intent(in) :: startingComm
        integer, intent(in) :: appId, key
        type(MPI_Comm), intent(inout) :: voltCoupledComm

        type(MPI_Comm) :: currentPool
        integer :: ierr, currentApp, voltRank

        currentApp = appId + 1 ! ensure the loop is entered
        currentPool = startingComm
        do while(currentApp .ne. appId)
            call MPI_Allreduce(0, voltRank, 1, MPI_INTEGER, MPI_MAX, currentPool, ierr)
            call MPI_Bcast(currentApp, 1, MPI_INTEGER, voltRank, currentPool, ierr)
            if(currentApp == appId) then
                ! it's my turn, split with volt and then skip making the second comm
                call MPI_comm_split(currentPool, appId, key, voltCoupledComm, ierr)
                call MPI_comm_split(currentPool, MPI_UNDEFINED, key, currentPool, ierr)
            else
                ! it's not my turn, don't split with volt, and then join the second comm
                call MPI_comm_split(currentPool, MPI_UNDEFINED, key, voltCoupledComm, ierr)
                call MPI_comm_split(currentPool, voltId, key, currentPool, ierr)
            endif
        enddo

    end subroutine

    ! this helper routine allows voltron to couple with a specific app
    ! and create a smaller communicator containing the remaining apps
    subroutine voltronSplitWithApp(couplingPool, appId, key, coupledComm)
        type(MPI_Comm), intent(inout) :: couplingPool
        integer, intent(in) :: appId, key
        type(MPI_Comm), intent(inout) :: coupledComm

        integer :: ierr, myRank

        ! tell everyone I'm the broadcasting root
        ! broadcast which app I'm creating a communicator with, split with it, and then
        ! create a smaller pool that excludes that app
        call MPI_comm_rank(couplingPool, myRank, ierr)
        call MPI_Allreduce(myRank, myRank, 1, MPI_INTEGER, MPI_MAX, couplingPool, ierr)
        call MPI_Bcast(appId, 1, MPI_INTEGER, myRank, couplingPool, ierr)
        call MPI_comm_split(couplingPool, appId, key, coupledComm, ierr)
        call MPI_comm_split(couplingPool, voltId, key, couplingPool, ierr)

    end subroutine

end module

