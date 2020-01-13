! Functions and data for MPI Gamera to talk to Voltron running on its own rank

module gam2VoltComm_mpi
    use gamapp_mpi
    use uservoltic

    implicit none

    type :: gam2VoltCommMpi_T
        integer :: voltMpiComm = MPI_COMM_NULL
        integer :: myRank, voltRank

        real(rp) :: time, tFin, DeepT, ShallowT, MJD
        integer :: ts
        logical :: doDeep

        ! array of all zeroes to simplify various send/receive calls
        integer, dimension(1) :: zeroArray = (/ 0 /)

        ! Shallow Gas Data Send Variables
        integer, dimension(1) :: sendCountsGasShallow, sendTypesGasShallow
        integer(MPI_ADDRESS_KIND), dimension(1) :: sendDisplsGasShallow
        ! Shallow Magnetic Field Data Send Variables
        integer, dimension(1) :: sendCountsBxyzShallow, sendTypesBxyzShallow
        integer(MPI_ADDRESS_KIND), dimension(1) :: sendDisplsBxyzShallow
        ! Shallow inEijk Data Receive Variables
        integer, dimension(1) :: recvCountsIneijkShallow, recvTypesIneijkShallow
        integer(MPI_ADDRESS_KIND), dimension(1) :: recvDisplsIneijkShallow
        ! Shallow inExyz Data Receive Variables
        integer, dimension(1) :: recvCountsInexyzShallow, recvTypesInexyzShallow
        integer(MPI_ADDRESS_KIND), dimension(1) :: recvDisplsInexyzShallow
        
        ! Deep Gas Data Send Variables
        integer, dimension(1) :: sendCountsGasDeep, sendTypesGasDeep
        integer(MPI_ADDRESS_KIND), dimension(1) :: sendDisplsGasDeep
        ! Deep Magnetic Field Data Send Variables
        integer, dimension(1) :: sendCountsBxyzDeep, sendTypesBxyzDeep
        integer(MPI_ADDRESS_KIND), dimension(1) :: sendDisplsBxyzDeep
        ! Deep Gas0 Data Receive Variables
        integer, dimension(1) :: recvCountsGas0Deep, recvTypesGas0Deep
        integer(MPI_ADDRESS_KIND), dimension(1) :: recvDisplsGas0Deep

    end type gam2VoltCommMpi_T

    contains

    ! setup the MPI communicator to talk to voltron, and send grid data
    subroutine initGam2Volt(g2vComm, gApp, voltComm)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(in) :: gApp
        integer, intent(in) :: voltComm

        integer :: length, commSize, ierr, numCells, dataCount, numInNeighbors, numOutNeighbors
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder, wasWeighted

        ! create a new communicator using MPI Topology
        call MPI_Comm_Size(voltComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Comm_rank(voltComm, g2vComm%myRank, ierr)

        ! send my i/j/k ranks to the voltron rank
        call mpi_gather(gApp%Grid%Ri, 1, MPI_INT, 0, 0, 0, commSize-1, voltComm, ierr)
        call mpi_gather(gApp%Grid%Rj, 1, MPI_INT, 0, 0, 0, commSize-1, voltComm, ierr)
        call mpi_gather(gApp%Grid%Rk, 1, MPI_INT, 0, 0, 0, commSize-1, voltComm, ierr)

        numCells = gApp%Grid%Nip*gApp%Grid%Njp*gApp%Grid%Nkp
        ! rank 0 send the number of physical cells to voltron rank
        if(g2vComm%myRank == 0) then
            call mpi_send(numCells, 1, MPI_INT, commSize-1, 97500, voltComm, ierr)
        endif

        ! each gamera rank only talks to the voltron rank on this new communicator
        ! appoximate amount of data transfer as num physical cells for both shallow and deep
        dataCount = numCells ! deep data
        if(gApp%Grid%Ri == 0) then
            dataCount = dataCount + numCells ! shallow data
        endif

        reorder = .true. ! allow MPI to reorder the ranks
        call mpi_dist_graph_create_adjacent(voltComm, &
            (/1/),(/commSize-1/),(/dataCount/), &
            (/1/),(/commSize-1/),(/dataCount/), &
            MPI_INFO_NULL, reorder, g2vComm%voltMpiComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Comm_rank(g2vComm%voltMpiComm, g2vComm%myRank, ierr)

        call mpi_dist_graph_neighbors_count(g2vComm%voltMpiComm,numInNeighbors,numOutNeighbors,wasWeighted,ierr)
        if(numInNeighbors /= 1) then
            print *,'Number of in edges was not 1 for rank ', g2vComm%myRank
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif
        if(numOutNeighbors /= 1) then
            print *,'Number of out edges was not 1 for rank ', g2vComm%myRank
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        ! get the rank of voltron in the new communicator
        call mpi_dist_graph_neighbors(g2vComm%voltMpiComm, numInNeighbors, g2vComm%voltRank, dataCount, &
                                      numOutNeighbors, g2vComm%voltRank, dataCount, ierr)

        ! send i/j/k ranks again since my rank may have changed in the new communicator
        call mpi_gather(gApp%Grid%Ri, 1, MPI_INT, 0, 0, 0, g2vComm%voltRank, voltComm, ierr)
        call mpi_gather(gApp%Grid%Rj, 1, MPI_INT, 0, 0, 0, g2vComm%voltRank, voltComm, ierr)
        call mpi_gather(gApp%Grid%Rk, 1, MPI_INT, 0, 0, 0, g2vComm%voltRank, voltComm, ierr)

        ! create the MPI datatypes for communicating state data with voltron
        call createG2VDataTypes(g2vComm)

        ! initialize all of the starting parameters from the voltron rank
        call mpi_bcast(g2vComm%time, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%tFin, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%DeepT, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%ShallowT, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%MJD, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%ts, 1, MPI_INT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%doDeep, 1, MPI_LOGICAL, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

    end subroutine initGam2Volt

    ! transmit data to voltron so that it can keep up with the simulation
    subroutine performStepVoltron(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(in) :: gApp

        integer :: ierr

        ! only the rank with Ri/Rj/Rk==0 should send the time values to voltron
        if(gApp%Grid%Ri==0 .and. gApp%Grid%Rj==0 .and. gApp%Grid%Rk==0) then
            call mpi_send(gApp%Model%t, 1, MPI_MYFLOAT, g2vComm%voltRank, 97600, g2vComm%voltMpiComm, ierr)
            call mpi_send(gApp%Model%ts, 1, MPI_INT, g2vComm%voltRank, 97700, g2vComm%voltMpiComm, ierr)
        endif

        ! all ranks receive the new time data from voltron after it has updated
        call mpi_bcast(g2vComm%time, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%MJD, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%ts, 1, MPI_INT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

    end subroutine performStepVoltron

    ! transmit state data to voltron over MPI and receive new data
    subroutine performShallowUpdate(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp

        integer :: ierr

        ! send state data to voltron
        ! Send Shallow Gas Data
        call mpi_neighbor_alltoallw(gApp%State%Gas, g2vComm%sendCountsGasShallow, &
                                    g2vComm%sendDisplsGasShallow, g2vComm%sendTypesGasShallow, &
                                    0, g2vComm%zeroArray, g2vComm%zeroArray, g2vComm%zeroArray, &
                                    g2vComm%voltMpiComm, ierr)
        ! Send Shallow Bxyz Data
        call mpi_neighbor_alltoallw(gApp%State%Bxyz, g2vComm%sendCountsBxyzShallow, &
                                    g2vComm%sendDisplsBxyzShallow, g2vComm%sendTypesBxyzShallow, &
                                    0, g2vComm%zeroArray, g2vComm%zeroArray, g2vComm%zeroArray, &
                                    g2vComm%voltMpiComm, ierr)

        ! Receive updated data from voltron
        ! The data goes into inEijk and inExyz in the IonInnerBC_T
        ! find the remix BC to write data to
        SELECT type(iiBC=>gApp%Grid%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                ! Recv Shallow inEijk Data
                call mpi_neighbor_alltoallw(0, g2vComm%zeroArray, g2vComm%zeroArray, g2vComm%zeroArray, &
                                            iiBC%inEijk, g2vComm%recvCountsIneijkShallow, &
                                            g2vComm%recvDisplsIneijkShallow, g2vComm%recvTypesIneijkShallow, &
                                            g2vComm%voltMpiComm, ierr)

                ! Recv Shallow inExyz Data
                call mpi_neighbor_alltoallw(0, g2vComm%zeroArray, g2vComm%zeroArray, g2vComm%zeroArray, &
                                            iiBC%inExyz, g2vComm%recvCountsInexyzShallow, &
                                            g2vComm%recvDisplsInexyzShallow, g2vComm%recvTypesInexyzShallow, &
                                            g2vComm%voltMpiComm, ierr)
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in Voltron MPI ShallowUpdate_mpi'
                stop
        END SELECT

        ! receive next time for shallow calculation
        call mpi_bcast(g2vComm%ShallowT, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)


    end subroutine performShallowUpdate

    ! transmit state data to voltron over MPI and receive new data
    subroutine performDeepUpdate(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp

        integer :: ierr

        if (.not. g2vComm%doDeep) then
            !Why are you even here?
            return
        endif

        ! send state data to voltron
        ! Send Deep Gas Data
        call mpi_neighbor_alltoallw(gApp%State%Gas, g2vComm%sendCountsGasDeep, &
                                    g2vComm%sendDisplsGasDeep, g2vComm%sendTypesGasDeep, &
                                    0, g2vComm%zeroArray, g2vComm%zeroArray, g2vComm%zeroArray, &
                                    g2vComm%voltMpiComm, ierr)
        ! Send Deep Bxyz Data
        call mpi_neighbor_alltoallw(gApp%State%Bxyz, g2vComm%sendCountsBxyzDeep, &
                                    g2vComm%sendDisplsBxyzDeep, g2vComm%sendTypesBxyzDeep, &
                                    0, g2vComm%zeroArray, g2vComm%zeroArray, g2vComm%zeroArray, &
                                    g2vComm%voltMpiComm, ierr)

        ! Receive updated data from Voltron
        ! Receive Deep Gas0 Data
        call mpi_neighbor_alltoallw(0, g2vComm%zeroArray, g2vComm%zeroArray, g2vComm%zeroArray, &
                                    gApp%Grid%Gas0, g2vComm%recvCountsGas0Deep, &
                                    g2vComm%recvDisplsGas0Deep, g2vComm%recvTypesGas0Deep, &
                                    g2vComm%voltMpiComm, ierr)

        ! receive next time for deep calculation
        call mpi_bcast(g2vComm%DeepT, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

    end subroutine performDeepUpdate

    subroutine createG2VDataTypes(g2vComm)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
    end subroutine createG2VDataTypes

end module gam2VoltComm_mpi

