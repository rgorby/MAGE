! Functions and data for MPI Gamera to talk to Voltron running on its own rank

module gam2VoltComm_mpi
    use gamapp_mpi
    use uservoltic
    use mpidefs
    use mpi_f08
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN, IEEE_QUIET_NAN

    implicit none

    type :: gam2VoltCommMpi_T
        type(MPI_Comm) :: voltMpiComm
        integer :: myRank, voltRank
        logical :: doSerialVoltron = .false., doAsyncCoupling = .true.
        logical :: firstDeepUpdate = .true.

        real(rp) :: time, tFin, DeepT, MJD
        integer :: ts, JpSt, JpSh, PsiSt, PsiSh
        logical :: doDeep

        ! array of all zeroes to simplify various send/receive calls
        integer, dimension(1) :: zeroArrayCounts = (/ 0 /)
        type(MPI_Datatype), dimension(1) :: zeroArrayTypes
        integer(kind=MPI_AN_MYADDR), dimension(1) :: zeroArrayDispls = (/ 0 /)

        ! SHALLOW COUPLING VARIABLES
        integer, dimension(1) :: recvCountsIneijkShallow
        type(MPI_Datatype), dimension(1) :: recvTypesIneijkShallow
        integer(kind=MPI_AN_MYADDR), dimension(1) :: recvDisplsIneijkShallow
        integer, dimension(1) :: recvCountsInexyzShallow
        type(MPI_Datatype), dimension(1) :: recvTypesInexyzShallow
        integer(kind=MPI_AN_MYADDR), dimension(1) :: recvDisplsInexyzShallow

        ! DEEP COUPLING VARIABLES
        integer, dimension(1) :: sendCountsGasDeep
        type(MPI_Datatype), dimension(1) :: sendTypesGasDeep
        integer(kind=MPI_AN_MYADDR), dimension(1) :: sendDisplsGasDeep
        integer, dimension(1) :: sendCountsBxyzDeep
        type(MPI_Datatype), dimension(1) :: sendTypesBxyzDeep
        integer(kind=MPI_AN_MYADDR), dimension(1) :: sendDisplsBxyzDeep
        integer, dimension(1) :: recvCountsGas0Deep
        type(MPI_Datatype), dimension(1) :: recvTypesGas0Deep
        integer(kind=MPI_AN_MYADDR), dimension(1) :: recvDisplsGas0Deep

    end type gam2VoltCommMpi_T

    contains

    ! setup the MPI communicator to talk to voltron, and send grid data
    subroutine initGam2Volt(g2vComm, gApp, allComm, optFilename, doIO)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp
        type(MPI_Comm), intent(in) :: allComm
        character(len=*), optional, intent(in) :: optFilename
        logical, optional, intent(in) :: doIO

        integer :: length, commSize, ierr, numCells, dataCount, numInNeighbors, numOutNeighbors
        integer :: commId
        type(MPI_Comm) :: voltComm
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder, wasWeighted, doIOX
        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        integer, dimension(1) :: rankArray, weightArray

        commId = gamId + voltId

        ! initialize F08 MPI objects
        g2vComm%voltMpiComm = MPI_COMM_NULL
        g2vComm%zeroArraytypes = (/ MPI_DATATYPE_NULL /)

        ! split voltron helpers off of the communicator
        ! split allComm into a communicator with only the non-helper voltron rank
        call MPI_Comm_rank(allComm, commSize, ierr)
        call MPI_comm_split(allComm, commId, commSize, voltComm, ierr)

        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
        else
            !Find input deck
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")
        !write(*,*) 'Reading input deck from ', trim(inpXML)
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Voltron',.true.)
        if (.not. gApp%Model%isLoud) call xmlInp%BeQuiet()

        call xmlInp%Set_Val(g2vComm%doSerialVoltron,"coupling/doSerial",.false.)
        call xmlInp%Set_Val(g2vComm%doAsyncCoupling,"coupling/doAsyncCoupling",.true.)
        if(g2vComm%doSerialVoltron) then
            ! don't do asynchronous coupling if comms are serial
            g2vComm%doAsyncCoupling = .false.
        endif


        if (present(doIO)) then
            doIOX = doIO
        else
            doIOX = .true.
        endif

        ! create a new communicator using MPI Topology
        call MPI_Comm_Size(voltComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Comm_rank(voltComm, g2vComm%myRank, ierr)

        ! send my i/j/k ranks to the voltron rank
        call mpi_gather(gApp%Grid%Ri, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, commSize-1, voltComm, ierr)
        call mpi_gather(gApp%Grid%Rj, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, commSize-1, voltComm, ierr)
        call mpi_gather(gApp%Grid%Rk, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, commSize-1, voltComm, ierr)

        numCells = gApp%Grid%Nip*gApp%Grid%Njp*gApp%Grid%Nkp
        ! rank 0 send the number of physical cells to voltron rank
        if(g2vComm%myRank == 0) then
            call mpi_send(numCells, 1, MPI_INTEGER, commSize-1, 97500, voltComm, ierr)
        endif

        ! each gamera rank only talks to the voltron rank on this new communicator
        ! appoximate amount of data transfer as num physical cells
        dataCount = numCells ! deep data

        reorder = .true. ! allow MPI to reorder the ranks
        call mpi_dist_graph_create_adjacent(voltComm, &
            1,(/commSize-1/),(/dataCount/), &
            1,(/commSize-1/),(/dataCount/), &
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
        call mpi_dist_graph_neighbors(g2vComm%voltMpiComm, numInNeighbors, rankArray, weightArray, &
                                      numOutNeighbors, rankArray, weightArray, ierr)
        g2vComm%voltRank = rankArray(1)

        ! send i/j/k ranks again since my rank may have changed in the new communicator
        call mpi_gather(gApp%Grid%Ri, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_gather(gApp%Grid%Rj, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_gather(gApp%Grid%Rk, 1, MPI_INTEGER, 0, 0, MPI_DATATYPE_NULL, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

        ! Send restart number so that voltron can ensure they have the same number
        ! only the rank with Ri/Rj/Rk==0 should send the value to voltron
        if(gApp%Grid%Ri==0 .and. gApp%Grid%Rj==0 .and. gApp%Grid%Rk==0) then
            call mpi_send(gApp%Model%IO%nRes, 1, MPI_INTEGER, g2vComm%voltRank, 97520, g2vComm%voltMpiComm, ierr)
        endif

        ! initialize all of the starting parameters from the voltron rank
        call mpi_bcast(g2vComm%time, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%tFin, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%DeepT, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%MJD, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%ts, 1, MPI_INTEGER, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%doDeep, 1, MPI_LOGICAL, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%JpSt, 1, MPI_INTEGER, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%JpSh, 1, MPI_INTEGER, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%PsiSt, 1, MPI_INTEGER, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(g2vComm%PsiSh, 1, MPI_INTEGER, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

        ! get updated Gamera parameters from the voltron rank
        call mpi_bcast(gApp%Model%t, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(gApp%Model%MJD0, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(gApp%Model%tFin, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(gApp%Model%dt, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

        ! send updated gamera parameter to the voltron rank
        ! only the rank with Ri/Rj/Rk==0 should send the values to voltron
        if(gApp%Grid%Ri==0 .and. gApp%Grid%Rj==0 .and. gApp%Grid%Rk==0) then
            call mpi_send(gApp%Model%dt0, 1, MPI_MYFLOAT, g2vComm%voltRank, 97510, g2vComm%voltMpiComm, ierr)
        endif

        ! synchronize IO timing
        call mpi_bcast(gApp%Model%IO%tOut,  1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(gApp%Model%IO%tRes,  1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(gApp%Model%IO%dtOut, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(gApp%Model%IO%dtRes, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)
        call mpi_bcast(gApp%Model%IO%tsOut, 1, MPI_INTEGER,     g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

        ! create the MPI datatypes for communicating state data with voltron
        call createG2VDataTypes(g2vComm, gApp)

        ! perform initial deep update if appropriate
        if (g2vComm%time>=g2vComm%DeepT) then
            ! deep update
            call performDeepUpdate(g2vComm, gApp)
        endif

        if (doIOX) then
            if (.not. gApp%Model%isRestart) then
                ! reset tOut before correcting initial file output
                gApp%Model%IO%tOut = gApp%Model%IO%tOut - gApp%Model%IO%dtOut
                call fOutput(gApp%Model,gApp%Grid,gApp%State)
            endif
            call consoleOutput(gApp%Model,gApp%Grid,gApp%State)
        endif

    end subroutine initGam2Volt

    subroutine endGam2VoltWaits(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp

        logical :: reqStat
        integer :: ierr

        if(g2vComm%doAsyncCoupling) then
            ! voltron sent shallow data which can't be cancelled
            !call recvShallowData(g2vComm, gApp)

            ! voltron sends all data at once now, so it won't send async shallow past the end time
        endif

    end subroutine endGam2VoltWaits

    ! update voltron time units locally while actual voltron is running
    subroutine localStepVoltronTime(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(in) :: gApp

        g2vComm%time = gApp%Model%t*gApp%Model%Units%gT0 !Time in seconds
        g2vComm%MJD = T2MJD(g2vComm%time,gApp%Model%MJD0)
        g2vComm%ts = gApp%Model%ts
    end subroutine localStepVoltronTime

    ! transmit data to voltron so that it can keep up with the simulation
    subroutine performStepVoltron(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(in) :: gApp

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! only the rank with Ri/Rj/Rk==0 should send the time values to voltron
        if(gApp%Grid%Ri==0 .and. gApp%Grid%Rj==0 .and. gApp%Grid%Rk==0) then
            call mpi_send(gApp%Model%t, 1, MPI_MYFLOAT, g2vComm%voltRank, 97600, g2vComm%voltMpiComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            call mpi_send(gApp%Model%ts, 1, MPI_INTEGER, g2vComm%voltRank, 97700, g2vComm%voltMpiComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
        endif

    end subroutine performStepVoltron

    ! transmit state data to voltron over MPI and receive new data
    subroutine performDeepUpdate(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp

        if(g2vComm%doSerialVoltron) then
            call doSerialDeepUpdate(g2vComm, gApp)
        else
            if(g2vComm%firstDeepUpdate) then
                call sendDeepData(g2vComm, gApp)
                g2vComm%firstDeepUpdate = .false.
            endif

            call doConcurrentDeepUpdate(g2vComm, gApp)

        endif

    end subroutine performDeepUpdate

    subroutine doSerialDeepUpdate(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp

        integer :: ierr

        ! send deep data
        call Tic("DeepSend")
        call sendDeepData(g2vComm, gApp)
        call Toc("DeepSend")

        ! receive shallow data
        call Tic("ShallowRecv")
        call recvShallowData(g2vComm, gApp)
        call Toc("ShallowRecv")

        ! receive deep data
        if(g2vComm%doDeep) then
            call Tic("DeepRecv")
            call recvDeepData(g2vComm, gApp)
            call Toc("DeepRecv")
        endif

        ! receive next time for deep calculation
        call mpi_bcast(g2vComm%DeepT, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

    end subroutine

    subroutine doConcurrentDeepUpdate(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp

        integer :: ierr

        ! now reverse process, get (waiting) data
        call Tic("ShallowRecv")
        call recvShallowData(g2vComm, gApp)
        call Toc("ShallowRecv")

        if(g2vComm%doDeep) then
            call Tic("DeepRecv")
            call recvDeepData(g2vComm, gApp)
            call Toc("DeepRecv")
        endif

        ! receive next time for deep calculation
        call mpi_bcast(g2vComm%DeepT, 1, MPI_MYFLOAT, g2vComm%voltRank, g2vComm%voltMpiComm, ierr)

        ! and send data for voltron to work on while gamera is busy
        call Tic("DeepSend")
        call sendDeepData(g2vComm, gApp)
        call Toc("DeepSend")

    end subroutine

    ! send deep state data to voltron over MPI
    subroutine sendDeepData(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(in) :: g2vComm
        type(gamAppMpi_T), intent(in) :: gApp

        integer :: ierr

        ! send state data to voltron
        ! Send Deep Gas Data
        call mpi_neighbor_alltoallw(gApp%State%Gas, g2vComm%sendCountsGasDeep, &
                                    g2vComm%sendDisplsGasDeep, g2vComm%sendTypesGasDeep, &
                                    gApp%State%Gas, g2vComm%zeroArrayCounts, &
                                    g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                    g2vComm%voltMpiComm, ierr)
        ! Send Deep Bxyz Data
        call mpi_neighbor_alltoallw(gApp%State%Bxyz, g2vComm%sendCountsBxyzDeep, &
                                    g2vComm%sendDisplsBxyzDeep, g2vComm%sendTypesBxyzDeep, &
                                    gApp%State%Bxyz, g2vComm%zeroArrayCounts, &
                                    g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                    g2vComm%voltMpiComm, ierr)

    end subroutine sendDeepData

    ! receive deep state data from voltron over MPI
    subroutine recvDeepData(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp

        integer :: ierr

        ! Receive updated data from Voltron
        ! Receive Deep Gas0 Data
        call Tic("VoltSync", .true.)
        call mpi_neighbor_alltoallw(gApp%Grid%Gas0, g2vComm%zeroArrayCounts, &
                                    g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                    gApp%Grid%Gas0, g2vComm%recvCountsGas0Deep, &
                                    g2vComm%recvDisplsGas0Deep, g2vComm%recvTypesGas0Deep, &
                                    g2vComm%voltMpiComm, ierr)
        call Toc("VoltSync", .true.)

    end subroutine recvDeepData

    ! receive shallow state data from voltron over MPI
    subroutine recvShallowData(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(inout) :: gApp

        integer :: ierr
        real(rp) :: nanValue

        ! Receive updated data from voltron
        ! The data goes into inEijk and inExyz in the IonInnerBC_T
        ! find the remix BC to write data to
        if(gApp%Grid%hasLowerBC(IDIR)) then
            SELECT type(iiBC=>gApp%Grid%externalBCs(INI)%p)
                TYPE IS (IonInnerBC_T)

                    iiBC%inEijk(:,:,:,:) = IEEE_VALUE(nanValue, IEEE_SIGNALING_NAN)
                    iiBC%inExyz(:,:,:,:) = IEEE_VALUE(nanValue, IEEE_SIGNALING_NAN)

                    ! Recv Shallow inEijk Data
                    call Tic("VoltSync", .true.)
                    call mpi_neighbor_alltoallw(0, g2vComm%zeroArrayCounts, &
                                                g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                                iiBC%inEijk, g2vComm%recvCountsIneijkShallow, &
                                                g2vComm%recvDisplsIneijkShallow, g2vComm%recvTypesIneijkShallow, &
                                                g2vComm%voltMpiComm, ierr)
                    call Toc("VoltSync", .true.)

                    ! Recv Shallow inExyz Data
                    call mpi_neighbor_alltoallw(0, g2vComm%zeroArrayCounts, &
                                                g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                                iiBC%inExyz, g2vComm%recvCountsInexyzShallow, &
                                                g2vComm%recvDisplsInexyzShallow, g2vComm%recvTypesInexyzShallow, &
                                                g2vComm%voltMpiComm, ierr)
                CLASS DEFAULT
                    write(*,*) 'Could not find Ion Inner BC in Voltron MPI performShallowUpdate'
                    stop
            END SELECT
        else
            ! not a rank with remix BC, but still need to call mpi_neighbor_alltoallw
            ! Recv nothing step 1
            call Tic("VoltSync", .true.)
            call mpi_neighbor_alltoallw(0, g2vComm%zeroArrayCounts, &
                                        g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                        0, g2vComm%zeroArrayCounts, &
                                        g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                        g2vComm%voltMpiComm, ierr)
            call Toc("VoltSync", .true.)

            ! Recv nothing step 2
            call mpi_neighbor_alltoallw(0, g2vComm%zeroArrayCounts, &
                                        g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                        0, g2vComm%zeroArrayCounts, &
                                        g2vComm%zeroArrayDispls, g2vComm%zeroArrayTypes, &
                                        g2vComm%voltMpiComm, ierr)
        endif

    end subroutine recvShallowData

    subroutine createG2VDataTypes(g2vComm, gApp)
        type(gam2VoltCommMpi_T), intent(inout) :: g2vComm
        type(gamAppMpi_T), intent(in) :: gApp

        integer :: ierr, dataSize, sendDataOffset, recvDataOffset
        type(MPI_Datatype) :: iPSI, iPSI1, Eijk2, EIjk3, Eijk4, Exyz2, Exyz3, Exyz4
        type(MPI_Datatype) :: iP,iPjP,iPjPkP,iPjPkP4Gas,iPjPkP4Bxyz,iPjPkP5Gas
        type(MPI_Datatype) :: iPG2,iPG2jPG2,iPG2jPG2kPG2,iPG2jPG2kPG24Gas,iPG2jPG2kPG25Gas

        associate(Grid=>gApp%Grid,Model=>gApp%Model, &
                  JpSt=>g2vComm%JpSt,JpSh=>g2vComm%JpSh, &
                  PsiSt=>g2vComm%PsiSt,PsiSh=>g2vComm%PsiSh)

        ! no need to allocate, all arrays are fixed size (1) since we only talk to voltron

        ! counts always 1
        ! displacements always 0
        ! datatypes to null by default
        g2vComm%recvCountsInexyzShallow = 1
        g2vComm%recvCountsIneijkShallow = 1
        g2vComm%recvDisplsInexyzShallow = 0
        g2vComm%recvDisplsIneijkShallow = 0
        g2vComm%recvTypesInexyzShallow = MPI_DATATYPE_NULL
        g2vComm%recvTypesIneijkShallow = MPI_DATATYPE_NULL

        g2vComm%sendCountsGasDeep = 1
        g2vComm%sendCountsBxyzDeep = 1
        g2vComm%sendDisplsGasDeep = 0
        g2vComm%sendDisplsBxyzDeep = 0
        g2vComm%sendTypesGasDeep = MPI_DATATYPE_NULL
        g2vComm%sendTypesBxyzDeep = MPI_DATATYPE_NULL

        if(g2vComm%doDeep) then
            g2vComm%recvCountsGas0Deep = 1
            g2vComm%recvDisplsGas0Deep = 0
            g2vComm%recvTypesGas0Deep = MPI_DATATYPE_NULL
        endif

        ! assemble datatypes
        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per value

        ! I dimension
        call mpi_type_contiguous(PsiSh, MPI_MYFLOAT, iPSI, ierr) ! PsiSh i
        call mpi_type_contiguous(PsiSh+1, MPI_MYFLOAT, iPSI1, ierr) ! PsiSh+1 i
        call mpi_type_contiguous(Grid%Nip, MPI_MYFLOAT, iP, ierr)
        call mpi_type_contiguous(Grid%Ni, MPI_MYFLOAT, iPG2, ierr)

        ! J dimension
        call mpi_type_hvector(Grid%Nj, 1, PsiSh*dataSize, iPSI, Exyz2, ierr) ! PsiSh i - p+2g j
        call mpi_type_hvector(Grid%Nj+1, 1, (PsiSh+1)*dataSize, iPSI1, Eijk2, ierr) ! PsiSh+1 i - p+2g+1 j
        call mpi_type_hvector(Grid%Njp, 1, Grid%Ni*datasize, iP, iPjP, ierr)
        call mpi_type_hvector(Grid%Nj, 1, Grid%Ni*datasize, iPG2, iPG2jPG2, ierr)

        ! K dimension
        call mpi_type_hvector(Grid%Nk, 1, Grid%Nj*PsiSh*dataSize, &
                              Exyz2, Exyz3, ierr) ! PsiSh i - p+2g j - p+2g k
        call mpi_type_hvector(Grid%Nk+1, 1, (Grid%Nj+1)*(PsiSh+1)*dataSize,&
                              Eijk2, Eijk3, ierr) ! PsiSh+1 i - p+2g+1 j - p+2g+1 k
        call mpi_type_hvector(Grid%Nkp, 1, Grid%Ni*Grid%Nj*datasize, iPjP, iPjPkP, ierr)
        call mpi_type_hvector(Grid%Nk, 1, Grid%Ni*Grid%Nj*datasize, iPG2jPG2, iPG2jPG2kPG2, ierr)

        ! 4th dimension
        call mpi_type_hvector(NDIM, 1, PsiSh*Grid%Nj*Grid%Nk*dataSize, Exyz3, Exyz4, ierr)
        call mpi_type_hvector(NDIM, 1, (PsiSh+1)*(Grid%Nj+1)*(Grid%Nk+1)*dataSize, Eijk3, Eijk4, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*datasize, iPjPkP, iPjPkP4Gas, ierr)
        call mpi_type_hvector(NDIM, 1, Grid%Ni*Grid%Nj*Grid%Nk*datasize, iPjPkP, iPjPkP4Bxyz, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*datasize, iPG2jPG2kPG2, iPG2jPG2kPG24Gas, ierr)

        ! 5th dimension
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,iPjPkP4Gas,iPjPkP5Gas,ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,iPG2jPG2kPG24Gas, &
                              iPG2jPG2kPG25Gas,ierr)

        ! create appropriate MPI Datatypes
        if(.not. Grid%hasLowerBC(IDIR)) then
            ! only gamera ranks with lower I boundary participate in shallow updates
            g2vComm%recvCountsInexyzShallow = 0
            g2vComm%recvCountsIneijkShallow = 0
            ! set these types to non null because MPI complains
            g2vComm%recvTypesInexyzShallow = MPI_INTEGER
            g2vComm%recvTypesIneijkShallow = MPI_INTEGER
        else
            !Inexyz
            recvDataOffset = 0
            call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, Exyz4, &
                                   g2vComm%recvTypesInexyzShallow(1), ierr)

            !Ineijk
            recvDataOffset = 0
            call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, Eijk4, &
                                   g2vComm%recvTypesIneijkShallow(1), ierr)
        endif

        ! Gas
        sendDataOffset = Model%nG*Grid%Nj*Grid%Ni + &
                         Model%nG*Grid%Ni + &
                         Model%nG
        call mpi_type_hindexed(1,(/1/),sendDataOffset*dataSize,iPjPkP5Gas,g2vComm%sendTypesGasDeep(1),ierr)

        ! Bxyz
        call mpi_type_hindexed(1,(/1/),sendDataOffset*dataSize,iPjPkP4Bxyz,g2vComm%sendTypesBxyzDeep(1),ierr)

        if(g2vComm%doDeep) then
            ! Gas0
            recvDataOffset = 0
            call mpi_type_hindexed(1,(/1/),recvDataOffset*dataSize,iPG2jPG2kPG25Gas, &
                                   g2vComm%recvTypesGas0Deep(1),ierr)
        endif

        ! commit new mpi datatypes
        call mpi_type_commit(g2vComm%recvTypesIneijkShallow(1), ierr)
        call mpi_type_commit(g2vComm%recvTypesInexyzShallow(1), ierr)

        call mpi_type_commit(g2vComm%sendTypesGasDeep(1), ierr)
        call mpi_type_commit(g2vComm%sendTypesBxyzDeep(1), ierr)

        if(g2vComm%doDeep) then
            call mpi_type_commit(g2vComm%recvTypesGas0Deep(1), ierr)
        endif

        end associate

    end subroutine createG2VDataTypes

end module gam2VoltComm_mpi

