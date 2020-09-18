! Collection of data and objects for the voltron middle man
! MPI version

module voltapp_mpi
    use voltapp
    use gamapp_mpi
    use gamapp
    use mpi
    use ebsquish, only : SquishBlocksRemain, DoSquishBlock
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN, IEEE_QUIET_NAN
    
    implicit none

    type, extends(voltApp_T) :: voltAppMpi_T
        integer :: voltMpiComm = MPI_COMM_NULL
        integer :: myRank
        type(gamApp_T) :: gAppLocal
        logical :: doSerialVoltron = .false., doAsyncShallow = .true.
        logical :: firstShallowUpdate = .true., firstDeepUpdate = .true.

        ! array of all zeroes to simplify various send/receive calls
        integer, dimension(:), allocatable :: zeroArrayCounts, zeroArrayTypes
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable ::  zeroArrayDispls

        ! list of gamera ranks to communicate with
        integer, dimension(:), allocatable :: sendRanks, recvRanks

        ! STEP VOLTRON VARIABLES
        integer :: timeReq=MPI_REQUEST_NULL, timeStepReq=MPI_REQUEST_NULL
        real(rp) :: timeBuffer
        integer :: timeStepBuffer

        ! SHALLOW COUPLING VARIABLES
        integer, dimension(:), allocatable :: recvCountsGasShallow, recvTypesGasShallow
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: recvDisplsGasShallow
        integer, dimension(:), allocatable :: recvCountsBxyzShallow, recvTypesBxyzShallow
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: recvDisplsBxyzShallow
        integer, dimension(:), allocatable :: sendCountsIneijkShallow, sendTypesIneijkShallow
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsIneijkShallow
        integer, dimension(:), allocatable :: sendCountsInexyzShallow, sendTypesInexyzShallow
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsInexyzShallow
        ! SHALLOW ASYNCHRONOUS VARIABLES
        integer :: shallowIneijkSendReq=MPI_REQUEST_NULL, shallowInexyzSendReq=MPI_REQUEST_NULL

        ! DEEP COUPLING VARIABLES
        integer, dimension(:), allocatable :: recvCountsGasDeep, recvTypesGasDeep
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: recvDisplsGasDeep
        integer, dimension(:), allocatable :: recvCountsBxyzDeep, recvTypesBxyzDeep
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: recvDisplsBxyzDeep
        integer, dimension(:), allocatable :: sendCountsGas0Deep, sendTypesGas0Deep
        integer(MPI_ADDRESS_KIND), dimension(:), allocatable :: sendDisplsGas0Deep
        logical :: deepProcessingInProgress = .false.

    end type voltAppMpi_T

    contains

    !end lingering MPI waits, if there are any
    subroutine endVoltronWaits(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        logical :: reqStat
        integer :: ierr

        if(vApp%timeReq /= MPI_REQUEST_NULL) then
            call MPI_REQUEST_GET_STATUS(vApp%timeReq,reqStat,MPI_STATUS_IGNORE,ierr)
            if(.not. reqStat) then
                call MPI_CANCEL(vApp%timeReq, ierr)
                call MPI_WAIT(vApp%timeReq, MPI_STATUS_IGNORE, ierr)
            endif
        endif

        if(vApp%timeStepReq /= MPI_REQUEST_NULL) then
            call MPI_REQUEST_GET_STATUS(vApp%timeStepReq,reqStat,MPI_STATUS_IGNORE,ierr)
            if(.not. reqStat) then
                call MPI_CANCEL(vApp%timeStepReq, ierr)
                call MPI_WAIT(vApp%timeStepReq, MPI_STATUS_IGNORE, ierr)
            endif
        endif

        call MPI_WAIT(vApp%shallowIneijkSendReq, MPI_STATUS_IGNORE, ierr)

        call MPI_WAIT(vApp%shallowInexyzSendReq, MPI_STATUS_IGNORE, ierr)

    end subroutine endVoltronWaits

    !Initialize Voltron (after Gamera has already been initialized)
    subroutine initVoltron_mpi(vApp, userInitFunc, voltComm, optFilename)
        type(voltAppMpi_T), intent(inout) :: vApp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        integer, intent(in) :: voltComm
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        integer :: commSize, ierr, numCells, length, ic, numInNeighbors, numOutNeighbors
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder, wasWeighted
        integer, allocatable, dimension(:) :: neighborRanks, inData, outData
        integer, allocatable, dimension(:) :: iRanks, jRanks, kRanks

        vApp%isSeparate = .true. ! running on a different process from the actual gamera ranks
        vApp%gAppLocal%Grid%lowMem = .true. ! tell Gamera to limit its memory usage

        ! create a new communicator using MPI Topology
        call MPI_Comm_Size(voltComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        allocate(neighborRanks(1:commSize-1))
        allocate(inData(1:commSize-1))
        allocate(outData(1:commSize-1))
        allocate(iRanks(1:commSize))
        allocate(jRanks(1:commSize))
        allocate(kRanks(1:commSize))

        allocate(vApp%zeroArrayCounts(1:commSize-1))
        allocate(vApp%zeroArrayTypes(1:commSize-1))
        allocate(vAPp%zeroArrayDispls(1:commSize-1))
        vApp%zeroArrayCounts(:) = 0
        vApp%zeroArrayTypes(:) = MPI_INT ! MPI_DATATYPE_NULL
        vApp%zeroArrayDispls(:) = 0

        ! doing a very very rough approximation of data transferred to help MPI reorder
        ! for deep updates, assume each rank sends data equal to its # physical cells
        ! for shallow updates, i=0 ranks send that much data again

        ! get i/j/k ranks from each Gamera mpi rank
        call mpi_gather(-1, 1, MPI_INT, iRanks, 1, MPI_INT, commSize-1, voltComm, ierr)
        call mpi_gather(-1, 1, MPI_INT, jRanks, 1, MPI_INT, commSize-1, voltComm, ierr)
        call mpi_gather(-1, 1, MPI_INT, kRanks, 1, MPI_INT, commSize-1, voltComm, ierr)

        ! get the number of physical cells from rank 0
        call mpi_recv(numCells, 1, MPI_INT, 0, 97500, voltComm, MPI_STATUS_IGNORE, ierr)

        do ic=1,commSize-1
            neighborRanks(ic) = ic-1
            if(iRanks(ic) == 0) then
                ! i=0 rank, does shallow update
                inData(ic) = numCells
                outData(ic) = numCells
            else
                ! i != 0 rank, does not do shallow update
                inData(ic) = 0
                outData(ic) = 0
            endif
            inData(ic) = inData(ic) + numCells
            outData(ic) = outData(ic) + numCells
        enddo

        reorder = .true. ! allow MPI to reorder the ranks
        call mpi_dist_graph_create_adjacent(voltComm, &
            commSize-1,neighborRanks,inData, &
            commSize-1,neighborRanks,outData, &
            MPI_INFO_NULL, reorder, vApp%voltMpiComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call MPI_Comm_rank(vApp%voltMpiComm, vApp%myRank, ierr)

        call mpi_dist_graph_neighbors_count(vApp%voltMpiComm,numInNeighbors,numOutNeighbors,wasWeighted,ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if (numInNeighbors /= numOutNeighbors) then
            print *,'Number of in edges and out edges did not match for voltron mpi'
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        allocate(vApp%sendRanks(numOutNeighbors))
        allocate(vApp%recvRanks(numInNeighbors))

        ! don't care about the weights, dump them into an existing array
        call mpi_dist_graph_neighbors(vApp%voltMpiComm, numInNeighbors, vApp%recvRanks, inData, &
                                      numOutNeighbors, vApp%sendRanks, inData, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! get i/j/k ranks again in case MPI ranks were reordered in the new communicator
        call mpi_gather(-1, 1, MPI_INT, iRanks, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_gather(-1, 1, MPI_INT, jRanks, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_gather(-1, 1, MPI_INT, kRanks, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)

        ! use standard voltron with local gamApp object
        if(present(optFilename)) then
            inpXML = optFilename
        else
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck in initVoltron_mpi, exiting ...")
        xmlInp = New_XML_Input(trim(inpXML),'Gamera',.true.)
        call xmlInp%Set_Val(vApp%doSerialVoltron,"/Voltron/coupling/doSerial",.false.)
        call xmlInp%Set_Val(vApp%doAsyncShallow, "/Voltron/coupling/doAsyncShallow",.true.)
        if(vApp%doSerialVoltron) then
            ! don't do asynchronous shallow if comms are serial
            vApp%doAsyncShallow = .false.
        endif
        vApp%gAppLocal%Grid%ijkShift(1:3) = 0
        call ReadCorners(vApp%gAppLocal%Model,vApp%gAppLocal%Grid,xmlInp,childGameraOpt=.true.)
        call SetRings(vApp%gAppLocal%Model,vApp%gAppLocal%Grid,xmlInp)
        call Corners2Grid(vApp%gAppLocal%Model,vApp%gAppLocal%Grid)
        call DefaultBCs(vApp%gAppLocal%Model,vApp%gAppLocal%Grid)
        call PrepState(vApp%gAppLocal%Model,vApp%gAppLocal%Grid,&
            vApp%gAppLocal%oState,vAPp%gApplocal%State,xmlInp,userInitFunc)

        ! now initialize basic voltron structures from gamera data
        if(present(optFilename)) then
            call initVoltron(vApp, vApp%gAppLocal, optFilename)
        else
            call initVoltron(vApp, vApp%gAppLocal)
        endif

        ! receive current time information in case of a restart
        if(vApp%gAppLocal%Model%isRestart) then
            call mpi_recv(vApp%IO%nOut, 1, MPI_INT, MPI_ANY_SOURCE, 97520, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(vApp%IO%nRes, 1, MPI_INT, MPI_ANY_SOURCE, 97530, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(vApp%gAppLocal%Model%t, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97540, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(vApp%gAppLocal%Model%ts, 1, MPI_INT, MPI_ANY_SOURCE, 97550, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)
        endif

        ! update voltron app time, MJD and ts values
        call stepVoltron(vApp, vApp%gAppLocal)

        ! correct initial shallow time, which always starts immediately
        vApp%ShallowT = vApp%time
        ! Deep start time is user-defined and may not have occurred before restart
        if(vApp%time > vApp%DeepT) vApp%DeepT = vApp%time

        ! send all of the initial voltron parameters to the gamera ranks
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%DeepT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%ShallowT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%MJD, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%ts, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%doDeep, 1, MPI_LOGICAL, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%mhd2mix%JStart, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%mhd2mix%JShells, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%mix2mhd%PsiStart, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%mix2mhd%PsiShells, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)

        ! send updated Gamera parameters to the gamera ranks
        call mpi_bcast(vApp%gAppLocal%Model%t, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%gAppLocal%Model%MJD0, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%gAppLocal%Model%tFin, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%gAppLocal%Model%dt, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)

        ! receive updated gamera parameters from gamera rank
        call mpi_recv(vApp%gAppLocal%Model%dt0, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97510, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)

        ! calculate what the next output and restart timing should be for gamera
        if(vApp%gAppLocal%Model%isRestart) then
            vApp%IO%tOut = floor(vApp%time/vApp%IO%dtOut)*vApp%IO%dtOut
            vApp%IO%tRes = vApp%time + vApp%IO%dtRes
        endif
        ! synchronize IO timing
        call mpi_bcast(vApp%IO%tOut/vApp%gAppLocal%Model%Units%gT0, &
                       1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%IO%tRes/vApp%gAppLocal%Model%Units%gT0, &
                       1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%IO%dtOut/vApp%gAppLocal%Model%Units%gT0, &
                       1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%IO%dtRes/vApp%gAppLocal%Model%Units%gT0, &
                       1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%IO%tsOut, 1, MPI_INT, vApp%myRank, vApp%voltMpiComm, ierr)

        ! create the MPI datatypes needed to transfer state data
        call createVoltDataTypes(vApp, iRanks, jRanks, kRanks)

        deallocate(neighborRanks, inData, outData, iRanks, jRanks, kRanks)

        ! perform initial shallow and deep updates if appropriate
        call ShallowUpdate_mpi(vApp, vApp%time)

        if (vApp%doDeep .and. vApp%time >= vApp%DeepT) then
            call DeepUpdate_mpi(vApp, vApp%time)
        endif

        ! if doing concurrent, start the asynchronous comms
        if(.not. vApp%doSerialVoltron) then
            call mpi_Irecv(vApp%timeBuffer, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97600, vApp%voltMpiComm, vApp%timeReq, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            call mpi_Irecv(vApp%timeStepBuffer, 1, MPI_INT, MPI_ANY_SOURCE, 97700, vApp%voltMpiComm, vApp%timeStepReq, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
        endif

    end subroutine initVoltron_mpi

    function gameraStepReady(vApp)
        type(voltAppMpi_T), intent(in) :: vApp
        logical :: gameraStepReady

        integer :: ierr

        call Tic("GameraSync")
        if(vApp%doSerialVoltron) then
            gameraStepReady = .true.
        else
            call MPI_REQUEST_GET_STATUS(vApp%timeReq,gameraStepReady,MPI_STATUS_IGNORE,ierr)
        endif
        call Toc("GameraSync")

    end function gameraStepReady

    subroutine waitForGameraStep(vApp)
        type(voltAppMpi_T), intent(in) :: vApp
        
        integer :: ierr

        call Tic("GameraSync")
        if(.not. vApp%doSerialVoltron) then
            call MPI_WAIT(vApp%timeReq, MPI_STATUS_IGNORE, ierr)
        endif
        call Toc("GameraSync")

    end subroutine waitForGameraStep

    ! MPI version of updating voltron variables
    subroutine stepVoltron_mpi(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr

        ! get gApp%Model%t,ts from gamera. All ranks have the same, just receive from one of them
        if(vApp%doSerialVoltron) then
            call mpi_recv(vApp%timeBuffer, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97600, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)

            call mpi_recv(vApp%timeStepBuffer, 1, MPI_INT, MPI_ANY_SOURCE, 97700, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)
        else
            call mpi_wait(vApp%timeReq, MPI_STATUS_IGNORE, ierr)
            call mpi_wait(vApp%timeStepReq, MPI_STATUS_IGNORE, ierr)
        endif

        vApp%gAppLocal%Model%t = vApp%timeBuffer
        vApp%gAppLocal%Model%ts = vApp%timeStepBuffer

        call stepVoltron(vApp, vApp%gAppLocal)

        if(.not. vApp%doSerialVoltron) then
            call mpi_Irecv(vApp%timeBuffer, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97600, vApp%voltMpiComm, vApp%timeReq, ierr)

            call mpi_Irecv(vApp%timeStepBuffer, 1, MPI_INT, MPI_ANY_SOURCE, 97700, vApp%voltMpiComm, vApp%timeStepReq, ierr)
        endif

    end subroutine stepVoltron_mpi

    ! special function for when we need to do both updates
    subroutine shallowAndDeepUpdate_Mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        integer :: ierr, asyncShallowBcastReq

        if(vApp%doSerialVoltron) then
            ! deep first
            call deepSerialUpdate_mpi(vApp, time)

            ! then shallow but don't resend data
            call shallowSerialUpdate_mpi(vApp, time, .true.)
        else
            if(vApp%firstDeepUpdate .or. vApp%firstShallowUpdate) then
                call deepSerialUpdate_mpi(vApp, time)
                call shallowSerialUpdate_mpi(vApp, time, .true.)
                vApp%firstDeepUpdate = .false.
                vApp%firstShallowUpdate = .false.
            else
                ! ensure deep update is complete
                do while(deepInProgress(vApp))
                    call doDeepBlock(vApp)
                enddo

                ! send both solutions
                call Tic("ShallowSend")
                call sendShallowData_mpi(vApp)
                call Toc("ShallowSend")
                if(vApp%doAsyncShallow) then
                    ! cannot mix bcast and Ibcast, so this must also by Ibcast even though I want it to be synchronous
                    call mpi_Ibcast(vApp%ShallowT + vApp%ShallowDT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, asyncShallowBcastReq, ierr)
                    call mpi_wait(asyncShallowBcastReq, MPI_STATUS_IGNORE, ierr)
                else
                    ! synchronous
                    call mpi_bcast(vApp%ShallowT + vApp%ShallowDT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
                endif
                call Tic("DeepSend")
                call sendDeepData_mpi(vApp)
                call Toc("DeepSend")
                call mpi_bcast(vApp%DeepT + vApp%DeepDT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)

                ! then receive just deep data
                call Tic("DeepRecv")
                call recvDeepData_mpi(vApp)
                call Toc("DeepRecv")

                ! then calculate shallow
                call Tic("ShallowUpdate")
                call ShallowUpdate(vApp, vApp%gAppLocal, time)
                call Toc("ShallowUpdate")

                ! then start deep
                ! setup squish operation but don't yet perform the computations
                call Tic("DeepUpdate")
                call PreSquishDeep(vApp, vApp%gAppLocal)
                call Toc("DeepUpdate")
                vApp%deepProcessingInProgress = .true.
            endif
        endif
    end subroutine
!----------
!Shallow coupling stuff
    subroutine ShallowUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        if(vApp%doSerialVoltron) then
            call shallowSerialUpdate_mpi(vApp, time)
        else
            ! if this if the first sequence, process initial data
            if(vApp%firstShallowUpdate) then
                call shallowSerialUpdate_mpi(vApp, time)
                vApp%firstShallowUpdate = .false.
            else
                call shallowConcurrentUpdate_mpi(vApp, time)
            endif
        endif

    end subroutine

    subroutine shallowSerialUpdate_mpi(vApp, time, skipUpdateGamera)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time
        logical, optional, intent(in) :: skipUpdateGamera

        logical :: doSkipUpdate
        integer :: ierr

        if(present(skipUpdateGamera)) then
            doSkipUpdate = skipUpdateGamera
        else
            doSkipUpdate = .false.
        endif
        if(doSkipUpdate) then
            ! do nothing here, do not update the incoming gamera data
        else
            ! fetch data from Gamera ranks
            call Tic("ShallowRecv")
            call recvShallowData_mpi(vApp)
            call Toc("ShallowRecv")
        endif

        ! call base update function with local data
        call Tic("ShallowUpdate")
        call ShallowUpdate(vApp, vApp%gAppLocal, time)
        call Toc("ShallowUpdate")

        ! send updated data to Gamera ranks
        call Tic("ShallowSend")
        call sendShallowData_mpi(vApp)
        call Toc("ShallowSend")

        ! send next time for shallow calculation to all gamera ranks
        call mpi_bcast(vApp%ShallowT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)

    end subroutine

    subroutine shallowConcurrentUpdate_mpi(vApp, time, skipUpdateGamera)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time
        logical, optional, intent(in) :: skipUpdateGamera

        integer :: ierr, asyncShallowBcastReq

        ! send updated data to Gamera ranks
        call Tic("ShallowSend")
        call sendShallowData_mpi(vApp)
        call Toc("ShallowSend")

        ! send next time for shallow calculation to all gamera ranks
        if(vApp%doAsyncShallow) then
            ! asynchronous
            ! cannot mix bcast and Ibcast, so this must also by Ibcast even though I want it to be synchronous
            call mpi_Ibcast(vApp%ShallowT + vApp%ShallowDT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, asyncShallowBcastReq, ierr)
            call mpi_wait(asyncShallowBcastReq, MPI_STATUS_IGNORE, ierr)
        else
            ! synchronous
            call mpi_bcast(vApp%ShallowT + vApp%ShallowDT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        endif

        if(present(skipUpdateGamera) .and. skipUpdateGamera) then
            ! do nothing here, do not update the incoming gamera data
        else
            ! fetch data from Gamera ranks
            call Tic("ShallowRecv")
            call recvShallowData_mpi(vApp)
            call Toc("ShallowRecv")
        endif

        ! call base update function with local data
        call Tic("ShallowUpdate")
        call ShallowUpdate(vApp, vApp%gAppLocal, time)
        call Toc("ShallowUpdate")

    end subroutine

    subroutine recvShallowData_mpi(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        real(rp) :: nanValue

        ! Receive Shallow Gas Data
        call mpi_neighbor_alltoallw(0, vApp%zeroArrayCounts, &
                                    vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                    vApp%gAppLocal%State%Gas, vApp%recvCountsGasShallow, &
                                    vApp%recvDisplsGasShallow, vApp%recvTypesGasShallow, &
                                    vApp%voltMpiComm, ierr)

        ! Receive Shallow Bxyz Data
        call mpi_neighbor_alltoallw(0, vApp%zeroArrayCounts, &
                                    vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                    vApp%gAppLocal%State%Bxyz, vApp%recvCountsBxyzShallow, &
                                    vApp%recvDisplsBxyzShallow, vApp%recvTypesBxyzShallow, &
                                    vApp%voltMpiComm, ierr)

    end subroutine recvShallowData_mpi

    subroutine sendShallowData_mpi(vApp)
        type(voltAppMpi_T), intent(in) :: vApp

        integer :: ierr

        ! send updated data to Gamera ranks
        ! voltron updates inEijk and inExyz in the IonInnerBC_T
        ! find the remix BC to read data from
        SELECT type(iiBC=>vApp%gAppLocal%Grid%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                if(vApp%doAsyncShallow) then
                    ! asynchronous
                    call mpi_wait(vApp%shallowIneijkSendReq, MPI_STATUS_IGNORE, ierr)
                    call mpi_Ineighbor_alltoallw(iiBC%inEijk, vApp%sendCountsIneijkShallow, &
                                                 vApp%sendDisplsIneijkShallow, vApp%sendTypesIneijkShallow, &
                                                 0, vApp%zeroArrayCounts, &
                                                 vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                                 vApp%voltMpiComm, vApp%shallowIneijkSendReq, ierr)

                    call mpi_wait(vApp%shallowInexyzSendReq, MPI_STATUS_IGNORE, ierr)
                    call mpi_Ineighbor_alltoallw(iiBC%inExyz, vApp%sendCountsInexyzShallow, &
                                                 vApp%sendDisplsInexyzShallow, vApp%sendTypesInexyzShallow, &
                                                 0, vApp%zeroArrayCounts, &
                                                 vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                                 vApp%voltMpiComm, vApp%shallowInexyzSendReq, ierr)
                else
                    ! synchronous
                    ! Send Shallow inEijk Data
                    call mpi_neighbor_alltoallw(iiBC%inEijk, vApp%sendCountsIneijkShallow, &
                                                vApp%sendDisplsIneijkShallow, vApp%sendTypesIneijkShallow, &
                                                0, vApp%zeroArrayCounts, &
                                                vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                                vApp%voltMpiComm, ierr)

                    ! Send Shallow inExyz Data
                    call mpi_neighbor_alltoallw(iiBC%inExyz, vApp%sendCountsInexyzShallow, &
                                                vApp%sendDisplsInexyzShallow, vApp%sendTypesInexyzShallow, &
                                                0, vApp%zeroArrayCounts, &
                                                vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                                vApp%voltMpiComm, ierr)
                endif
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in Voltron MPI ShallowUpdate_mpi'
                stop
        END SELECT

    end subroutine sendShallowData_mpi

!----------
!Deep coupling stuff (time coming from vApp%time, so in seconds)

    function deepInProgress(vApp)
        type(voltAppMpi_T), intent(in) :: vApp
        logical :: deepInProgress

        deepInProgress = vApp%deepProcessingInProgress

    end function deepInProgress

    subroutine doDeepBlock(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        if(.not. vApp%deepProcessingInProgress) return

        call Tic("Squish")
        call DoSquishBlock(vApp)
        call Toc("Squish")

        if(.not. SquishBlocksRemain()) then
            vApp%deepProcessingInProgress = .false.
        endif

        if(.not. vApp%deepProcessingInProgress) then
            call PostSquishDeep(vApp, vApp%gAppLocal)
        endif

    end subroutine doDeepBlock

    subroutine DeepUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv
        integer :: ierr

        if (.not. vApp%doDeep) then
            !Why are you even here?
            return
        else
            if(vApp%doSerialVoltron) then
                call deepSerialUpdate_mpi(vApp, time)
            else
                if(vApp%firstDeepUpdate) then
                    call deepSerialUpdate_mpi(vApp, time)
                    vApp%firstDeepUpdate = .false.
                else
                    ! ensure deep update is complete
                    do while(deepInProgress(vApp))
                        call doDeepBlock(vApp)
                    enddo

                    call deepConcurrentUpdate_mpi(vApp, time)
                endif
            endif
        endif
    end subroutine

    subroutine deepSerialUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv
        integer :: ierr

        if (.not. vApp%doDeep) then
            !Why are you even here?
            return
        else
            ! fetch data from Gamera ranks
            call Tic("DeepRecv")
            call recvDeepData_mpi(vApp)
            call Toc("DeepRecv")

            ! call base update function with local data
            call Tic("DeepUpdate")
            call DeepUpdate(vApp, vApp%gAppLocal, time)
            call Toc("DeepUpdate")

            ! send updated data to Gamera ranks
            call Tic("DeepSend")
            call sendDeepData_mpi(vApp)
            call Toc("DeepSend")

            ! send next time for deep calculation to all gamera ranks
            call mpi_bcast(vApp%DeepT,1,MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        endif

    end subroutine

    subroutine deepConcurrentUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv
        integer :: ierr

        if (.not. vApp%doDeep) then
            !Why are you even here?
            return
        else
            ! send updated data to Gamera ranks
            call Tic("DeepSend")
            call sendDeepData_mpi(vApp)
            call Toc("DeepSend")

            ! send next time for deep calculation to all gamera ranks
            call mpi_bcast(vApp%DeepT+vApp%DeepDT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)

            ! fetch data from Gamera ranks
            call Tic("DeepRecv")
            call recvDeepData_mpi(vApp)
            call Toc("DeepRecv")

            ! setup squish operation but don't yet perform the computations
            call Tic("DeepUpdate")
            call PreSquishDeep(vApp, vApp%gAppLocal)
            call Toc("DeepUpdate")
            vApp%deepProcessingInProgress = .true.

        endif
    end subroutine

    subroutine recvDeepData_mpi(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr

        ! Receive Deep Gas Data
        call mpi_neighbor_alltoallw(vApp%gAppLocal%State%Gas, vApp%zeroArrayCounts, &
                                    vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                    vApp%gAppLocal%State%Gas, vApp%recvCountsGasDeep, &
                                    vApp%recvDisplsGasDeep, vApp%recvTypesGasDeep, &
                                    vApp%voltMpiComm, ierr)
        ! Receive Deep Bxyz Data
        call mpi_neighbor_alltoallw(vApp%gAppLocal%State%Bxyz, vApp%zeroArrayCounts, &
                                    vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                    vApp%gAppLocal%State%Bxyz, vApp%recvCountsBxyzDeep, &
                                    vApp%recvDisplsBxyzDeep, vApp%recvTypesBxyzDeep, &
                                    vApp%voltMpiComm, ierr)

    end subroutine recvDeepData_mpi

     subroutine sendDeepData_mpi(vApp)
        type(voltAppMpi_T), intent(in) :: vApp

        integer :: ierr

        ! Send Deep Gas0 Data
        call mpi_neighbor_alltoallw(vApp%gAppLocal%Grid%Gas0, vApp%sendCountsGas0Deep, &
                                    vApp%sendDisplsGas0Deep, vApp%sendTypesGas0Deep, &
                                    vApp%gAppLocal%Grid%Gas0, vApp%zeroArrayCounts, &
                                    vApp%zeroArrayDispls, vApp%zeroArrayTypes, &
                                    vApp%voltMpiComm, ierr)

    end subroutine sendDeepData_mpi

    subroutine createVoltDataTypes(vApp, iRanks, jRanks, kRanks)
        type(voltAppMpi_T), intent(inout) :: vApp
        integer, dimension(1:SIZE(vApp%recvRanks)+1), intent(in) :: iRanks, jRanks, kRanks

        integer :: ierr, NiRanks, NjRanks, NkRanks, NipT, NjpT, NkpT, dataSize
        integer :: r, rRank, recvDataOffset, recvDatatype
        integer :: iJP, iJPjP, iJPjPkP, iJPjPkP4Gas, iJPjPkP5Gas, iJPjPkP4Bxyz, iJPjPkP5Bxyz
        integer :: iJP3, iJP3jP, iJP3jPG, iJP3jPG2, iJP3jPkPG2, iJP3jPGkPG2, iJP3jPG2kPG2
        integer :: iJP3jPkPG24Bxyz, iJP3jPGkPG24Bxyz, iJP3jPG2kPG24Bxyz
        integer :: sRank,sendDataOffset,iPSI,iPSI1,Exyz2,Eijk2,Exyz3,Eijk3,Exyz4,Eijk4
        integer :: iP,iPjP,iPjPkP,iPjPkP4Bxyz,iPjPkP4Gas,iPjPkP5Gas
        integer :: iPG2,iPG2jPG2,iPG2jPG2kPG2,iPG2jPG2kPG24Gas,iPG2jPG2kPG25Gas

        associate(Grid=>vApp%gAppLocal%Grid,Model=>vApp%gAppLocal%Model, &
                  JpSt=>vApp%mhd2mix%JStart,JpSh=>vApp%mhd2mix%JShells, &
                  PsiSt=>vApp%mix2mhd%PsiStart,PsiSh=>vApp%mix2mhd%PsiShells)

        allocate(vApp%recvCountsGasShallow(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvDisplsGasShallow(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvTypesGasShallow(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvCountsBxyzShallow(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvDisplsBxyzShallow(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvTypesBxyzShallow(1:SIZE(vApp%recvRanks)))
        allocate(vApp%sendCountsInexyzShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendDisplsInexyzShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendTypesInexyzShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendCountsIneijkShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendDisplsIneijkShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendTypesIneijkShallow(1:SIZE(vApp%sendRanks)))

        if(vApp%doDeep) then
            allocate(vApp%recvCountsGasDeep(1:SIZE(vApp%recvRanks)))
            allocate(vApp%recvDisplsGasDeep(1:SIZE(vApp%recvRanks)))
            allocate(vApp%recvTypesGasDeep(1:SIZE(vApp%recvRanks)))
            allocate(vApp%recvCountsBxyzDeep(1:SIZE(vApp%recvRanks)))
            allocate(vApp%recvDisplsBxyzDeep(1:SIZE(vApp%recvRanks)))
            allocate(vApp%recvTypesBxyzDeep(1:SIZE(vApp%recvRanks)))
            allocate(vApp%sendCountsGas0Deep(1:SIZE(vApp%sendRanks)))
            allocate(vApp%sendDisplsGas0Deep(1:SIZE(vApp%sendRanks)))
            allocate(vApp%sendTypesGas0Deep(1:SIZE(vApp%sendRanks)))
        endif

        ! counts are always 1 because we're sending a single (complicated) mpi datatype
        vApp%recvCountsGasShallow(:) = 1
        vApp%recvCountsBxyzShallow(:) = 1
        vApp%sendCountsInexyzShallow(:) = 1
        vAPp%sendCountsIneijkShallow(:) = 1

        ! displacements are always 0 because the displacements are baked into each mpi datatype
        vApp%recvDisplsGasShallow(:) = 0
        vApp%recvDisplsBxyzShallow(:) = 0
        vApp%sendDisplsInexyzShallow(:) = 0
        vApp%sendDisplsIneijkShallow(:) = 0

        ! set all datatypes to null by default
        vApp%recvTypesGasShallow(:) = MPI_DATATYPE_NULL
        vApp%recvTypesBxyzShallow(:) = MPI_DATATYPE_NULL
        vApp%sendTypesInexyzShallow(:) = MPI_DATATYPE_NULL
        vApp%sendTypesIneijkShallow(:) = MPI_DATATYPE_NULL

        if(vApp%doDeep) then
            vApp%recvCountsGasDeep(:) = 1
            vApp%recvCountsBxyzDeep(:) = 1
            vApp%sendCountsGas0Deep(:) = 1

            vApp%recvDisplsGasDeep(:) = 0
            vApp%recvDisplsBxyzDeep(:) = 0
            vApp%sendDisplsGas0Deep(:) = 0

            vApp%recvTypesGasDeep(:) = MPI_DATATYPE_NULL
            vApp%recvTypesBxyzDeep(:) = MPI_DATATYPE_NULL
            vApp%sendTypesGas0Deep(:) = MPI_DATATYPE_NULL
        endif

        ! assemble the different datatypes
        call mpi_type_extent(MPI_MYFLOAT, dataSize, ierr) ! number of bytes per array entry
        NiRanks = maxval(iRanks)+1 ! gamera i decomposition
        NjRanks = maxval(jRanks)+1 ! gamera j decomposition
        NkRanks = maxval(kRanks)+1 ! gamera k decomposition
        NipT = Grid%Nip / NiRanks ! number of physical cells per gamera rank in the i direction
        NjpT = Grid%Njp / NjRanks ! number of physical cells per gamera rank in the j direction
        NkpT = Grid%Nkp / NkRanks ! number of physical cells per gamera rank in the k direction

        ! I dimension
        call mpi_type_contiguous(JpSh, MPI_MYFLOAT, iJP, ierr) ! JpSh i
        call mpi_type_contiguous(JpSh+3, MPI_MYFLOAT, iJP3, ierr) ! JpSh+3 i
        call mpi_type_contiguous(PsiSh, MPI_MYFLOAT, iPSI, ierr) ! PsiSh i
        call mpi_type_contiguous(PsiSh+1, MPI_MYFLOAT, iPSI1, ierr) ! PsiSh+1 i
        call mpi_type_contiguous(NipT, MPI_MYFLOAT, iP, ierr) ! physical i
        call mpi_type_contiguous(2*Model%nG+NipT, MPI_MYFLOAT, iPG2, ierr) ! physical + 2*ghosts i

        ! J dimension
        call mpi_type_hvector(NjpT, 1, Grid%Ni*dataSize, iJP, iJPjP, ierr) ! JpSh i - physical j
        call mpi_type_hvector(NjpT, 1, Grid%Ni*dataSize, iJP3, iJP3jP, ierr) ! JpSh+3 i - physical j
        call mpi_type_hvector(Model%nG+NjpT, 1, Grid%Ni*dataSize, iJP3, iJP3jPG, ierr) ! JpSh+3 i - p+g j
        call mpi_type_hvector(2*Model%nG+NjpT, 1, Grid%Ni*dataSize, iJP3, iJP3jPG2, ierr) ! JpSh+3 i - p+2g j
        call mpi_type_hvector(2*Model%nG+NjpT, 1, PsiSh*dataSize, iPSI, Exyz2, ierr) ! PsiSh i - p+2g j
        call mpi_type_hvector(2*Model%nG+NjpT+1, 1, (PsiSh+1)*dataSize, iPSI1, Eijk2, ierr) ! PsiSh+1 i - p+2g+1 j
        call mpi_type_hvector(NjpT, 1, Grid%Ni*dataSize, iP, iPjP, ierr) ! physical i - physical j
        call mpi_type_hvector(2*Model%nG+NjpT, 1, Grid%Ni*dataSize, iPG2, iPG2jPG2, ierr) ! p+2g i - p+2g j

        ! K dimension - currently assume NO k decomposition
        call mpi_type_hvector(NkpT, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iJPjP, iJPjPkP, ierr) ! JpSh i - physical j - physical k
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iJP3jP, iJP3jPkPG2, ierr) ! JpSh+3 i - p j - p+2g k
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iJP3jPG, iJP3jPGkPG2, ierr) ! JpSh+3 i - p+g j - p+2g k
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Ni*Grid%Nj*dataSize, &
                              iJP3jPG2, iJP3jPG2kPG2, ierr) ! JpSh+3 i - p+2g j - p+2g k
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Nj*PsiSh*dataSize, &
                              Exyz2, Exyz3, ierr) ! PsiSh i - p+2g j - p+2g k
        call mpi_type_hvector(2*Model%ng+NkpT+1, 1, (Grid%Nj+1)*(PsiSh+1)*dataSize,&
                              Eijk2, Eijk3, ierr) ! PsiSh+1 i - p+2g+1 j - p+2g+1 k
        call mpi_type_hvector(NkpT, 1, Grid%Ni*Grid%Nj*dataSize, iPjP, iPjPkP, ierr)
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Ni*Grid%Nj*dataSize, iPG2jPG2, iPG2jPG2kPG2, ierr)

        ! 4th dimension
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iJPjPkP, iJPjPkP4Gas, ierr)
        call mpi_type_hvector(NDIM, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iJP3jPkPG2,   iJP3jPkPG24Bxyz,   ierr)
        call mpi_type_hvector(NDIM, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iJP3jPGkPG2,  iJP3jPGkPG24Bxyz,  ierr)
        call mpi_type_hvector(NDIM, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iJP3jPG2kPG2, iJP3jPG2kPG24Bxyz, ierr)
        call mpi_type_hvector(NDIM, 1, PsiSh*Grid%Nj*Grid%Nk*dataSize, Exyz3, Exyz4, ierr)
        call mpi_type_hvector(NDIM, 1, (PsiSh+1)*(Grid%Nj+1)*(Grid%Nk+1)*dataSize, Eijk3, Eijk4, ierr)
        call mpi_type_hvector(NDIM, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPjPkP, iPjPkP4Bxyz, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPjPkP, iPjPkP4Gas, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPG2jPG2kPG2, iPG2jPG2kPG24Gas, ierr)

        ! 5th dimension
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,iJPjPkP4Gas,iJPjPkP5Gas,ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,iPjPkP4Gas,iPjPkP5Gas,ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, &
                              iPG2jPG2kPG24Gas, iPG2jPG2kPG25Gas, ierr)

        ! figure out exactly what data needs to be sent to (and received from) each gamera rank
        ! create custom MPI datatypes to perform these transfers
        do r=1,SIZE(vApp%recvRanks)
            rRank = vApp%recvRanks(r)+1

            if(iRanks(rRank) .gt. 0) then
                ! never get shallow data from any rank but minimum i
                vApp%recvCountsGasShallow(r) = 0
                vApp%recvCountsBxyzShallow(r) = 0
                ! set these types to non null because MPI complains
                vApp%recvTypesGasShallow(r) = MPI_INT
                vApp%recvTypesBxyzShallow(r) = MPI_INT
            else
                ! calculate the byte offset to the start of the data

                ! gas
                recvDataOffset = (Model%nG + kRanks(rRank)*NkpT)*Grid%Nj*Grid%Ni + &
                                 (Model%nG + jRanks(rRank)*NjpT)*Grid%Ni + &
                                 (Model%nG + JpSt - 1) ! turning JpSt into an offset

                call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, iJPjPkP5Gas, &
                                       vApp%recvTypesGasShallow(r), ierr)

                !Bxyz
                recvDatatype = MPI_DATATYPE_NULL
                recvDataOffset  = Model%nG + JpSt - 2 ! JpSt-1
                if(NjRanks == 1) then
                    ! only one J rank, ghosts on both sides
                    ! no change to recvDataOffset
                    recvDatatype = iJP3jPG2kPG24Bxyz
                elseif(jRanks(rRank) == 0) then
                    ! min J rank, ghosts on low side
                    ! no change to recvDataOffset
                    recvDatatype = iJP3jPGkPG24Bxyz
                elseif(jRanks(rRank) == (NjRanks-1)) then
                    ! max J rank, ghosts on high side
                    recvDataOffset = recvDataOffset + (Model%nG + jRanks(rRank)*NjpT)*Grid%Ni
                    recvDatatype = iJP3jPGkPG24Bxyz
                else
                    ! no ghosts
                    recvDataOffset = recvDataOffset + (Model%nG + jRanks(rRank)*NjpT)*Grid%Ni
                    recvDatatype = iJP3jPkPG24Bxyz
                endif

                call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, recvDatatype, &
                                       vApp%recvTypesBxyzShallow(r), ierr)
            endif

            if(vApp%doDeep) then
                ! gas
                recvDataOffset = (Model%nG + kRanks(rRank)*NkpT)*Grid%Nj*Grid%Ni + &
                                 (Model%nG + jRanks(rRank)*NjpT)*Grid%Ni + &
                                 (Model%nG + iRanks(rRank)*NipT)
                call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, iPjPkP5Gas, &
                                       vApp%recvTypesGasDeep(r), ierr)
                ! Bxyz
                call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, iPjPkP4Bxyz, &
                                       vApp%recvTypesBxyzDeep(r), ierr)
            endif
        enddo

        do r=1,SIZE(vApp%sendRanks)
            sRank = vApp%sendRanks(r)+1

            if(iRanks(sRank) .gt. 0) then
                ! never send shallow data to any rank but minimum i
                vApp%sendCountsInexyzShallow(r) = 0
                vApp%sendCountsIneijkShallow(r) = 0
                ! set these types to non null because MPI complains
                vApp%sendTypesInexyzShallow(r) = MPI_INT
                vApp%sendTypesIneijkShallow(r) = MPI_INT
            else
                ! calculate the byte offset to the start of the data

                ! Inexyz
                sendDataOffset = kRanks(sRank)*NkpT*Grid%Nj*PsiSh + &
                                 jRanks(sRank)*NjpT*PsiSh

                call mpi_type_hindexed(1, (/1/), sendDataOffset*dataSize, Exyz4, &
                                       vApp%sendTypesInexyzShallow(r), ierr)

                !Ineijk
                sendDataOffset = kRanks(sRank)*NkpT*(Grid%Nj+1)*(PsiSh+1) + &
                                 jRanks(sRank)*NjpT*(PsiSh+1)

                call mpi_type_hindexed(1, (/1/), sendDataOffset*dataSize, Eijk4, &
                                       vApp%sendTypesIneijkShallow(r), ierr)
            endif

            if(vApp%doDeep) then
                ! gas0
                sendDataOffset = kRanks(sRank)*NkpT*Grid%Nj*Grid%Ni + &
                                 jRanks(sRank)*NjpT*Grid%Ni + &
                                 iRanks(sRank)*NipT
                call mpi_type_hindexed(1, (/1/), sendDataOffset*dataSize, iPG2jPG2kPG25Gas, &
                                       vApp%sendTypesGas0Deep(r), ierr)
            endif
        enddo

        do r=1,size(vApp%recvTypesGasShallow)
            call mpi_type_commit(vApp%recvTypesGasShallow(r), ierr)
            call mpi_type_commit(vApp%recvTypesBxyzShallow(r), ierr)
            call mpi_type_commit(vApp%sendTypesInexyzShallow(r), ierr)
            call mpi_type_commit(vApp%sendTypesIneijkShallow(r), ierr)

            if(vApp%doDeep) then
                call mpi_type_commit(vApp%recvTypesGasDeep(r), ierr)
                call mpi_type_commit(vApp%recvTypesBxyzDeep(r), ierr)
                call mpi_type_commit(vApp%sendTypesGas0Deep(r), ierr)
            endif
        enddo

        end associate

    end subroutine createVoltDataTypes

end module voltapp_mpi

