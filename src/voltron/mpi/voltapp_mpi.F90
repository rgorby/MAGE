! Collection of data and objects for the voltron middle man
! MPI version

module voltapp_mpi
    use voltmpitypes
    use gamapp_mpi
    use gamapp
    use mpi_f08
    use ebsquish, only : SquishBlocksRemain, DoSquishBlock
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN, IEEE_QUIET_NAN
    use volthelpers_mpi
    
    implicit none

    contains

    !end lingering MPI waits, if there are any
    subroutine endVoltronWaits(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        logical :: reqStat
        integer :: ierr

        call MPI_TEST(vApp%timeReq,reqStat,MPI_STATUS_IGNORE,ierr)
        if(.not. reqStat) then
            call MPI_CANCEL(vApp%timeReq, ierr)
            call MPI_WAIT(vApp%timeReq, MPI_STATUS_IGNORE, ierr)
        endif

        call MPI_TEST(vApp%timeStepReq,reqStat,MPI_STATUS_IGNORE,ierr)
        if(.not. reqStat) then
            call MPI_CANCEL(vApp%timeStepReq, ierr)
            call MPI_WAIT(vApp%timeStepReq, MPI_STATUS_IGNORE, ierr)
        endif

        call MPI_TEST(vApp%shallowIneijkSendReq,reqStat,MPI_STATUS_IGNORE,ierr)
        if(.not. reqStat) then
            ! async neighborhood ops don't support cancel
            call MPI_WAIT(vApp%shallowIneijkSendReq, MPI_STATUS_IGNORE, ierr)
        endif

        call MPI_TEST(vApp%shallowInexyzSendReq,reqStat,MPI_STATUS_IGNORE,ierr)
        if(.not. reqStat) then
            ! async neighborhood ops don't support cancel
            call MPI_WAIT(vApp%shallowInexyzSendReq, MPI_STATUS_IGNORE, ierr)
        endif

        if(vApp%vHelpWin /= MPI_WIN_NULL) then
            call MPI_WIN_FREE(vApp%vHelpWin, ierr)
        endif

    end subroutine endVoltronWaits

    !Initialize Voltron (after Gamera has already been initialized)
    subroutine initVoltron_mpi(vApp, userInitFunc, helperComm, allComm, optFilename)
        type(voltAppMpi_T), intent(inout) :: vApp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        type(MPI_Comm), intent(in) :: helperComm
        type(MPI_Comm), intent(in) :: allComm
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        integer :: commSize, ierr, numCells, length, ic, numInNeighbors, numOutNeighbors
        type(MPI_Comm) :: voltComm
        integer :: nHelpers, gamNRES, commId
        character( len = MPI_MAX_ERROR_STRING) :: message
        logical :: reorder, wasWeighted
        integer, allocatable, dimension(:) :: neighborRanks, inData, outData
        integer, allocatable, dimension(:) :: iRanks, jRanks, kRanks
        integer(KIND=MPI_BASE_MYADDR) :: winsize

        ! initialize F08 MPI objects
        vApp%vHelpComm = MPI_COMM_NULL
        vApp%vHelpWin = MPI_WIN_NULL
        vApp%voltMpiComm = MPI_COMM_NULL
        vApp%timeReq = MPI_REQUEST_NULL
        vApp%timeStepReq = MPI_REQUEST_NULL
        vApp%shallowIneijkSendReq = MPI_REQUEST_NULL
        vApp%shallowInexyzSendReq = MPI_REQUEST_NULL


        vApp%isSeparate = .true. ! running on a different process from the actual gamera ranks
        vApp%gAppLocal%Grid%lowMem = .true. ! tell Gamera to limit its memory usage

        ! get info about voltron-only mpi communicator
        vApp%vHelpComm = helperComm
        call MPI_Comm_rank(vApp%vHelpComm, vApp%vHelpRank, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(vApp%vHelpRank > 0) vApp%amHelper = .true.

        ! split allComm into a communicator with only the non-helper voltron rank
        call MPI_Comm_rank(allComm, commSize, ierr)
        commId = gamId+voltId
        if(vApp%amHelper) then
            call MPI_comm_split(allComm, MPI_UNDEFINED, commSize, voltComm, ierr)
        else
            call MPI_comm_split(allComm, commId, commSize, voltComm, ierr)
        endif

        ! helpers don't do full voltron initialization
        if(vApp%amHelper) then
            vApp%isLoud = .false.
            vApp%writeFiles = .false.

            ! read helper XML options
            if(present(optFilename)) then
                inpXML = optFilename
            else
                call getIDeckStr(inpXML)
            endif
            call CheckFileOrDie(inpXML,"Error opening input deck in initVoltron_mpi, exiting ...")
            xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Gamera',.false.)

            call xmlInp%Set_Val(vApp%useHelpers,"/Kaiju/Voltron/Helpers/useHelpers",.false.)
            call xmlInp%Set_Val(vApp%doSquishHelp,"/Kaiju/Voltron/Helpers/doSquishHelp",.true.)
            call xmlInp%Set_Val(vApp%masterSquish,"/Kaiju/Voltron/Helpers/masterSquish",.false.)
            call xmlInp%Set_Val(vApp%squishLoadBalance,"/Kaiju/Voltron/Helpers/squishLoadBalance",.true.)
            call xmlInp%Set_Val(nHelpers,"/Kaiju/Voltron/Helpers/numHelpers",0)
            call MPI_Comm_Size(vApp%vHelpComm, commSize, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            if(.not. vApp%useHelpers) then
                print *,"Voltron helpers were created, but the helping option is disabled."
                print *,"Please either turn the /Voltron/Helpers/useHelpers option on, or "
                print *,"  remove the unnecessary Voltron helper ranks."
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif
            if(nHelpers .ne. commSize-1) then
                print *,"The number of voltron helpers is not correct."
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            endif

	    ! add topology to the helper communicator to permit neighborhood operations
	    reorder = .false. ! don't allow MPI to reorder the ranks, master must remain master
	    call mpi_dist_graph_create_adjacent(helperComm, &
                1,(/0/),(/1/), &
                1,(/0/),(/1/), &
                MPI_INFO_NULL, reorder, vApp%vHelpComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            ! initialize a full baby voltron
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

            ! get starting time values from master voltron rank
            call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
            call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)

            ! accept the mpi window on the voltron master rank, expose none of my own memory
            winsize = 0 ! have to use this because it's an odd size integer
            call mpi_win_create(MPI_BOTTOM, winsize, 1, MPI_INFO_NULL, vApp%vHelpComm, vApp%vHelpWin, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            return
        endif

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
        vApp%zeroArrayTypes(:) = MPI_INTEGER ! MPI_DATATYPE_NULL
        vApp%zeroArrayDispls(:) = 0

        ! doing a very very rough approximation of data transferred to help MPI reorder
        ! for deep updates, assume each rank sends data equal to its # physical cells

        ! get i/j/k ranks from each Gamera mpi rank
        call mpi_gather(-1, 1, MPI_INTEGER, iRanks, 1, MPI_INTEGER, commSize-1, voltComm, ierr)
        call mpi_gather(-1, 1, MPI_INTEGER, jRanks, 1, MPI_INTEGER, commSize-1, voltComm, ierr)
        call mpi_gather(-1, 1, MPI_INTEGER, kRanks, 1, MPI_INTEGER, commSize-1, voltComm, ierr)

        ! get the number of physical cells from rank 0
        call mpi_recv(numCells, 1, MPI_INTEGER, 0, 97500, voltComm, MPI_STATUS_IGNORE, ierr)

        do ic=1,commSize-1
            neighborRanks(ic) = ic-1
            inData(ic) = numCells
            outData(ic) = numCells
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
        call mpi_gather(-1, 1, MPI_INTEGER, iRanks, 1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_gather(-1, 1, MPI_INTEGER, jRanks, 1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_gather(-1, 1, MPI_INTEGER, kRanks, 1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)

        ! use standard voltron with local gamApp object
        if(present(optFilename)) then
            inpXML = optFilename
        else
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck in initVoltron_mpi, exiting ...")
        if (vApp%amHelper) then
            xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Gamera',.false.)
        else
            xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Gamera',.true.)
        endif
        call xmlInp%Set_Val(vApp%doSerialVoltron,"/Kaiju/Voltron/coupling/doSerial",.false.)
        call xmlInp%Set_Val(vApp%doAsyncCoupling, "/Kaiju/Voltron/coupling/doAsyncCoupling",.true.)
        call xmlInp%Set_Val(vApp%useHelpers,"/Kaiju/Voltron/Helpers/useHelpers",.false.)
        call xmlInp%Set_Val(vApp%doSquishHelp,"/Kaiju/Voltron/Helpers/doSquishHelp",.true.)
        call xmlInp%Set_Val(vApp%masterSquish,"/Kaiju/Voltron/Helpers/masterSquish",.false.)
        call xmlInp%Set_Val(vApp%squishLoadBalance,"/Kaiju/Voltron/Helpers/squishLoadBalance",.true.)
        call xmlInp%Set_Val(nHelpers,"/Kaiju/Voltron/Helpers/numHelpers",0)

        if(vApp%doSerialVoltron) then
            ! don't do asynchronous coupling if comms are serial
            vApp%doAsyncCoupling = .false.
        endif


        if(vApp%masterSquish .and. vApp%squishLoadBalance) then
            print *,"Dynamic load balancing of squish helpers is not supported if the"
            print *,"     voltron master rank is also calculating squish."
            print *,"     Please set either:"
            print *,"        /Kaiju/Voltron/Helpers/masterSquish       to False"
            print *,"     or"
            print *,"        /Kaiju/Voltron/Helpers/squishLoadBalance  to False"
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        call MPI_Comm_Size(vApp%vHelpComm, commSize, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if(.not. vApp%useHelpers) nHelpers = 0
        if(nHelpers .eq. 0) vApp%useHelpers = .false.
        if(nHelpers .ne. commSize-1) then
            print *,"The number of voltron helpers is not correct."
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        endif

        vApp%gAppLocal%Grid%ijkShift(1:3) = 0
        call ReadCorners(vApp%gAppLocal%Model,vApp%gAppLocal%Grid,xmlInp,childGameraOpt=.true.)
        call SetRings(vApp%gAppLocal%Model,vApp%gAppLocal%Grid,xmlInp)
        call Corners2Grid(vApp%gAppLocal%Model,vApp%gAppLocal%Grid)
        call DefaultBCs(vApp%gAppLocal%Model,vApp%gAppLocal%Grid)
        call PrepState(vApp%gAppLocal%Model,vApp%gAppLocal%Grid,&
            vApp%gAppLocal%oState,vAPp%gApplocal%State,xmlInp,userInitFunc)

        ! now initialize basic voltron structures from gamera data
        if(vApp%gcmCplComm /= MPI_COMM_NULL) then
          call init_gcm_mix_mpi(vApp%gcm,vApp%gcmCplComm,vApp%gcmCplRank)
        endif
        if(present(optFilename)) then
            call initVoltron(vApp, vApp%gAppLocal, optFilename)
        else
            call initVoltron(vApp, vApp%gAppLocal)
        endif

        ! Receive Gamera's restart number and ensure Voltron has the same restart number
        call mpi_recv(gamNRES, 1, MPI_INTEGER, MPI_ANY_SOURCE, 97520, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)
        if (vApp%gAppLocal%Model%isRestart .and. vApp%IO%nRes /= gamNRES) then
            write(*,*) "Gamera and Voltron disagree on restart number, you should sort that out."
            write(*,*) "Error code: A house divided cannot stand"
            write(*,*) "   Voltron nRes = ", vApp%IO%nRes
            write(*,*) "   Gamera  nRes = ", gamNRES
            stop
        endif

        ! send all of the initial voltron parameters to the gamera ranks
        call mpi_bcast(vApp%time,     1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%tFin,     1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%DeepT,    1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%MJD,      1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%ts,       1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%doDeep,   1, MPI_LOGICAL, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%mhd2mix%JStart,    1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%mhd2mix%JShells,   1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%mix2mhd%PsiStart,  1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%mix2mhd%PsiShells, 1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)

        ! send updated Gamera parameters to the gamera ranks
        call mpi_bcast(vApp%gAppLocal%Model%t,    1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%gAppLocal%Model%MJD0, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%gAppLocal%Model%tFin, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%gAppLocal%Model%dt,   1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)

        ! receive updated gamera parameters from gamera rank
        call mpi_recv(vApp%gAppLocal%Model%dt0, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97510, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)

        ! synchronize IO timing
        call mpi_bcast(vApp%IO%tOut/vApp%gAppLocal%Model%Units%gT0, &
                       1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%IO%tRes/vApp%gAppLocal%Model%Units%gT0, &
                       1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%IO%dtOut/vApp%gAppLocal%Model%Units%gT0, &
                       1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%IO%dtRes/vApp%gAppLocal%Model%Units%gT0, &
                       1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
        call mpi_bcast(vApp%IO%tsOut, 1, MPI_INTEGER, vApp%myRank, vApp%voltMpiComm, ierr)

	! create the MPI datatypes needed to transfer state data
        call createVoltDataTypes(vApp, iRanks, jRanks, kRanks)

        deallocate(neighborRanks, inData, outData, iRanks, jRanks, kRanks)

        if(vApp%useHelpers) then
            if(vApp%doSquishHelp) then
                ! over-ride the number of squish blocks
                if(vApp%masterSquish) then
                    vApp%ebTrcApp%ebSquish%numSquishBlocks = nHelpers + 1
                else
                    vApp%ebTrcApp%ebSquish%numSquishBlocks = nHelpers
                endif

                if(vApp%ebTrcApp%ebSquish%numSquishBlocks /= size(vApp%ebTrcApp%ebSquish%blockStartIndices)) then
                    ! number of squish blocks changed
                    deallocate(vApp%ebTrcApp%ebSquish%blockStartIndices)
                    allocate(vApp%ebTrcApp%ebSquish%blockStartIndices(vApp%ebTrcApp%ebSquish%numSquishBlocks))
                    call LoadBalanceBlocks(vApp) ! redo the default load balancing for this number of blocks
                endif

            endif

            call MPI_Comm_Size(vApp%vHelpComm, commSize, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            ! add topology to the helper communicator to permit neighborhood operations
            reorder = .false. ! don't allow MPI to reorder the ranks, master must remain master
            allocate(neighborRanks(1:commSize-1))
            allocate(inData(1:commSize-1))
            allocate(outData(1:commSize-1))
            do ic=1,commSize-1
                neighborRanks(ic) = ic
            enddo
            ! currently all helpers only send and receive data from rank 0
            inData(:) = 1
            outData(:) = 1
            call mpi_dist_graph_create_adjacent(helperComm, &
                commSize-1,neighborRanks,inData, &
                commSize-1,neighborRanks,outData, &
                MPI_INFO_NULL, reorder, vApp%vHelpComm, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            deallocate(neighborRanks, inData, outData)

            ! send initial timing data to helper voltron ranks
            call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
            call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)

            ! allocate data for helpers to set status, and share with them
            allocate(vApp%vHelpIdle(nHelpers))
            vApp%vHelpIdle = 0
            call mpi_type_extent(MPI_INTEGER, length, ierr) ! size of integer
            winsize = nHelpers*length ! have to use this because it's an odd size integer
            call mpi_win_create(vApp%vHelpIdle, winsize, length, MPI_INFO_NULL, vApp%vHelpComm, vApp%vHelpWin, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
        endif

        ! perform initial deep update if appropriate
        call Tic("Coupling")
        if (vApp%time >= vApp%DeepT) then
            ! do deep
            call DeepUpdate_mpi(vApp, vApp%time)
        endif
        call Toc("Coupling")

    end subroutine initVoltron_mpi

    function gameraStepReady(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp
        logical :: gameraStepReady

        integer :: ierr

        call Tic("GameraSync")
        if(vApp%doSerialVoltron) then
            gameraStepReady = .true.
        else
            call MPI_TEST(vApp%timeReq,gameraStepReady,MPI_STATUS_IGNORE,ierr)
        endif
        call Toc("GameraSync")

    end function gameraStepReady

    subroutine waitForGameraStep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp
        
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
        if(vApp%doSerialVoltron .or. vApp%firstStepUpdate) then
            call mpi_recv(vApp%timeBuffer, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97600, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)

            call mpi_recv(vApp%timeStepBuffer, 1, MPI_INTEGER, MPI_ANY_SOURCE, 97700, vApp%voltMpiComm, MPI_STATUS_IGNORE, ierr)
            vApp%firstStepUpdate = .false.
        else
            call mpi_wait(vApp%timeReq, MPI_STATUS_IGNORE, ierr)
            call mpi_wait(vApp%timeStepReq, MPI_STATUS_IGNORE, ierr)
        endif

        vApp%gAppLocal%Model%t = vApp%timeBuffer
        vApp%gAppLocal%Model%ts = vApp%timeStepBuffer

        call stepVoltron(vApp, vApp%gAppLocal)

        if(.not. vApp%doSerialVoltron) then
            call mpi_Irecv(vApp%timeBuffer, 1, MPI_MYFLOAT, MPI_ANY_SOURCE, 97600, vApp%voltMpiComm, vApp%timeReq, ierr)

            call mpi_Irecv(vApp%timeStepBuffer, 1, MPI_INTEGER, MPI_ANY_SOURCE, 97700, vApp%voltMpiComm, vApp%timeStepReq, ierr)
        endif

    end subroutine stepVoltron_mpi

!----------
!Deep coupling stuff (time coming from vApp%time, so in seconds)

    subroutine startDeep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        ! convert gamera data to mixInput
        call Tic("G2R")
        call convertGameraToRemix(vApp%mhd2mix, vApp%gAppLocal, vApp%remixApp)
        call Toc("G2R")

        if (vApp%doGCM .and. vApp%time >=0) then
            call Tic("GCM2MIX")
            call coupleGCM2MIX(vApp%gcm,vApp%remixApp%ion,vApp%MJD,vApp%time,vApp%gcmCplComm,vApp%gcmCplRank)
            call Toc("GCM2MIX")
        end if

        ! run remix
        call Tic("ReMIX")
        call runRemix(vApp)
        call Toc("ReMIX")

        ! convert mixOutput to gamera data
        call Tic("R2G")
        call convertRemixToGamera(vApp%mix2mhd, vApp%remixApp, vApp%gAppLocal)
        call Toc("R2G")

        !Update coupling time now so that voltron knows what to expect
        vApp%DeepT = vApp%DeepT + vApp%DeepDT

        ! only do imag after spinup with deep enabled
        if(vApp%doDeep .and. vApp%time > 0) then
            call Tic("DeepUpdate")
            call PreDeep(vApp, vApp%gAppLocal)
            call SquishStart(vApp)

            if(vApp%useHelpers .and. vApp%doSquishHelp) then
                call Tic("VoltHelpers")
                call vhReqStep(vApp)
                call vhReqSquishStart(vApp)
                call Toc("VoltHelpers")
            endif

            ! moving this to after voltron helpers are started
            call DoImag(vApp)

            vApp%deepProcessingInProgress = .true.
            call Toc("DeepUpdate")
        else
            vApp%gAppLocal%Grid%Gas0 = 0
        endif

    end subroutine startDeep

    subroutine endDeep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        vApp%deepProcessingInProgress = .false.

        ! only do imag after spinup with deep enabled
        if(vApp%doDeep .and. vApp%time >= 0) then
            call Tic("DeepUpdate")

            if(vApp%useHelpers .and. vApp%doSquishHelp) then
                call Tic("VoltHelpers")
                call vhReqSquishEnd(vApp)
                call Toc("VoltHelpers")
            endif

            call SquishEnd(vApp)
            call PostDeep(vApp, vApp%gAppLocal)
            call Toc("DeepUpdate")
        endif

    end subroutine endDeep

    function deepInProgress(vApp)
        type(voltAppMpi_T), intent(in) :: vApp
        logical :: deepInProgress

        deepInProgress = vApp%deepProcessingInProgress

    end function deepInProgress

    subroutine doDeepBlock(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        if(.not. vApp%deepProcessingInProgress) return

        if(SquishBlocksRemain(vApp)) then
            call Tic("DeepUpdate")
            call Tic("Squish")
            call DoSquishBlock(vApp)
            call Toc("Squish")
            call Toc("DeepUpdate")
        endif

        if(.not. SquishBlocksRemain(vApp)) then
            if((.not. vApp%useHelpers) .or. &
               (vApp%useHelpers .and. .not. vApp%doSquishHelp) .or. &
               (vApp%useHelpers .and. vApp%doSquishHelp .and. allHelpersIdle(vApp)) ) then
                call endDeep(vApp)
            endif
        endif

    end subroutine doDeepBlock

    subroutine firstDeep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        real(rp) :: saveDeepT

        vApp%firstDeepUpdate = .false.

        if(vApp%doSerialVoltron) then
            ! serial needs no prep
            return
        endif

        call Tic("DeepRecv")
        call recvDeepData_mpi(vApp)
        call Toc("DeepRecv")
        saveDeepT = vApp%DeepT ! save the DeepT so it won't change during initial update
        call startDeep(vApp)
        vApp%firstDeepUpdate = .false.
        do while(deepInProgress(vApp))
            call doDeepBlock(vApp)
        enddo
        vApp%DeepT = saveDeepT ! restore DeepT

    end subroutine

    subroutine DeepUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv
        integer :: ierr

        if(vApp%firstDeepUpdate) call firstDeep(vApp)

        if(vApp%doSerialVoltron) then
            call deepSerialUpdate_mpi(vApp, time)
        else
            ! ensure deep update is complete
            do while(deepInProgress(vApp))
                call doDeepBlock(vApp)
            enddo
            call deepConcurrentUpdate_mpi(vApp, time)
        endif
    end subroutine

    subroutine deepSerialUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv
        integer :: ierr

            ! Update coupling DT
            vApp%DeepDT = vApp%TargetDeepDT

            ! fetch data from Gamera ranks
            call Tic("DeepRecv")
            call recvDeepData_mpi(vApp)
            call Toc("DeepRecv")

            ! call base update function with local data
            call Tic("DeepUpdate")
            call DeepUpdate(vApp, vApp%gAppLocal)
            call Toc("DeepUpdate")

            ! send updated data to Gamera ranks
            call Tic("ShallowSend")
            call sendShallowData_mpi(vApp)
            call Toc("ShallowSend")

            ! send updated data to Gamera ranks
            call Tic("DeepSend")
            if(vApp%doDeep) then
                call sendDeepData_mpi(vApp)
            endif

            ! send next time for deep calculation to all gamera ranks
            call mpi_bcast(vApp%DeepT,1,MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
            call Toc("DeepSend")

    end subroutine

    subroutine deepConcurrentUpdate_mpi(vApp, time)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv
        integer :: ierr

            ! send updated data to Gamera ranks
            call Tic("ShallowSend")
            call sendShallowData_mpi(vApp)
            call Toc("ShallowSend")

            ! send updated data to Gamera ranks
            call Tic("DeepSend")
            if(vApp%doDeep) then
                call sendDeepData_mpi(vApp)
            endif

            ! Update coupling DT
            vApp%DeepDT = vApp%TargetDeepDT

            ! send next time for deep calculation to all gamera ranks
            call mpi_bcast(vApp%DeepT+vApp%DeepDT, 1, MPI_MYFLOAT, vApp%myRank, vApp%voltMpiComm, ierr)
            call Toc("DeepSend")

            ! fetch data from Gamera ranks
            call Tic("DeepRecv")
            call recvDeepData_mpi(vApp)
            call Toc("DeepRecv")

            ! setup squish operation but don't yet perform the computations
            call startDeep(vApp)

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

    subroutine sendShallowData_mpi(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr

        ! send updated data to Gamera ranks
        ! voltron updates inEijk and inExyz in the IonInnerBC_T
        ! find the remix BC to read data from
        SELECT type(iiBC=>vApp%gAppLocal%Grid%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                if(vApp%doAsyncCoupling) then
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


    subroutine createVoltDataTypes(vApp, iRanks, jRanks, kRanks)
        type(voltAppMpi_T), intent(inout) :: vApp
        integer, dimension(1:SIZE(vApp%recvRanks)+1), intent(in) :: iRanks, jRanks, kRanks

        integer :: ierr, NiRanks, NjRanks, NkRanks, NipT, NjpT, NkpT, dataSize
        integer :: r, rRank, recvDataOffset, sRank, sendDataOffset
        type(MPI_Datatype) :: recvDatatype
        type(MPI_Datatype) :: iPSI,iPSI1,Exyz2,Eijk2,Exyz3,Eijk3,Exyz4,Eijk4
        type(MPI_Datatype) :: iP,iPjP,iPjPkP,iPjPkP4Bxyz,iPjPkP4Gas,iPjPkP5Gas
        type(MPI_Datatype) :: iPG2,iPG2jPG2,iPG2jPG2kPG2,iPG2jPG2kPG24Gas,iPG2jPG2kPG25Gas

        associate(Grid=>vApp%gAppLocal%Grid,Model=>vApp%gAppLocal%Model, &
                  JpSt=>vApp%mhd2mix%JStart,JpSh=>vApp%mhd2mix%JShells, &
                  PsiSt=>vApp%mix2mhd%PsiStart,PsiSh=>vApp%mix2mhd%PsiShells)

        allocate(vApp%sendCountsInexyzShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendDisplsInexyzShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendTypesInexyzShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendCountsIneijkShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendDisplsIneijkShallow(1:SIZE(vApp%sendRanks)))
        allocate(vApp%sendTypesIneijkShallow(1:SIZE(vApp%sendRanks)))

        allocate(vApp%recvCountsGasDeep(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvDisplsGasDeep(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvTypesGasDeep(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvCountsBxyzDeep(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvDisplsBxyzDeep(1:SIZE(vApp%recvRanks)))
        allocate(vApp%recvTypesBxyzDeep(1:SIZE(vApp%recvRanks)))

        if(vApp%doDeep) then
            allocate(vApp%sendCountsGas0Deep(1:SIZE(vApp%sendRanks)))
            allocate(vApp%sendDisplsGas0Deep(1:SIZE(vApp%sendRanks)))
            allocate(vApp%sendTypesGas0Deep(1:SIZE(vApp%sendRanks)))
        endif

        ! counts are always 1 because we're sending a single (complicated) mpi datatype
        ! displacements are always 0 because the displacements are baked into each mpi datatype
        ! set all datatypes to null by default
        vApp%sendCountsInexyzShallow(:) = 1
        vApp%sendCountsIneijkShallow(:) = 1
        vApp%sendDisplsInexyzShallow(:) = 0
        vApp%sendDisplsIneijkShallow(:) = 0
        vApp%sendTypesInexyzShallow(:) = MPI_DATATYPE_NULL
        vApp%sendTypesIneijkShallow(:) = MPI_DATATYPE_NULL

        vApp%recvCountsGasDeep(:) = 1
        vApp%recvCountsBxyzDeep(:) = 1
        vApp%recvDisplsGasDeep(:) = 0
        vApp%recvDisplsBxyzDeep(:) = 0
        vApp%recvTypesGasDeep(:) = MPI_DATATYPE_NULL
        vApp%recvTypesBxyzDeep(:) = MPI_DATATYPE_NULL

        if(vApp%doDeep) then
            vApp%sendCountsGas0Deep(:) = 1
            vApp%sendDisplsGas0Deep(:) = 0
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
        call mpi_type_contiguous(PsiSh, MPI_MYFLOAT, iPSI, ierr) ! PsiSh i
        call mpi_type_contiguous(PsiSh+1, MPI_MYFLOAT, iPSI1, ierr) ! PsiSh+1 i
        call mpi_type_contiguous(NipT, MPI_MYFLOAT, iP, ierr) ! physical i
        call mpi_type_contiguous(2*Model%nG+NipT, MPI_MYFLOAT, iPG2, ierr) ! physical + 2*ghosts i

        ! J dimension
        call mpi_type_hvector(2*Model%nG+NjpT, 1, PsiSh*dataSize, iPSI, Exyz2, ierr) ! PsiSh i - p+2g j
        call mpi_type_hvector(2*Model%nG+NjpT+1, 1, (PsiSh+1)*dataSize, iPSI1, Eijk2, ierr) ! PsiSh+1 i - p+2g+1 j
        call mpi_type_hvector(NjpT, 1, Grid%Ni*dataSize, iP, iPjP, ierr) ! physical i - physical j
        call mpi_type_hvector(2*Model%nG+NjpT, 1, Grid%Ni*dataSize, iPG2, iPG2jPG2, ierr) ! p+2g i - p+2g j

        ! K dimension - currently assume NO k decomposition
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Nj*PsiSh*dataSize, &
                              Exyz2, Exyz3, ierr) ! PsiSh i - p+2g j - p+2g k
        call mpi_type_hvector(2*Model%ng+NkpT+1, 1, (Grid%Nj+1)*(PsiSh+1)*dataSize,&
                              Eijk2, Eijk3, ierr) ! PsiSh+1 i - p+2g+1 j - p+2g+1 k
        call mpi_type_hvector(NkpT, 1, Grid%Ni*Grid%Nj*dataSize, iPjP, iPjPkP, ierr)
        call mpi_type_hvector(2*Model%nG+NkpT, 1, Grid%Ni*Grid%Nj*dataSize, iPG2jPG2, iPG2jPG2kPG2, ierr)

        ! 4th dimension
        call mpi_type_hvector(NDIM, 1, PsiSh*Grid%Nj*Grid%Nk*dataSize, Exyz3, Exyz4, ierr)
        call mpi_type_hvector(NDIM, 1, (PsiSh+1)*(Grid%Nj+1)*(Grid%Nk+1)*dataSize, Eijk3, Eijk4, ierr)
        call mpi_type_hvector(NDIM, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPjPkP, iPjPkP4Bxyz, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPjPkP, iPjPkP4Gas, ierr)
        call mpi_type_hvector(NVAR, 1, Grid%Ni*Grid%Nj*Grid%Nk*dataSize, iPG2jPG2kPG2, iPG2jPG2kPG24Gas, ierr)

        ! 5th dimension
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize,iPjPkP4Gas,iPjPkP5Gas,ierr)
        call mpi_type_hvector(Model%nSpc+1,1,NVAR*Grid%Ni*Grid%Nj*Grid%Nk*dataSize, &
                              iPG2jPG2kPG24Gas, iPG2jPG2kPG25Gas, ierr)

        ! figure out exactly what data needs to be sent to (and received from) each gamera rank
        ! create custom MPI datatypes to perform these transfers
        do r=1,SIZE(vApp%recvRanks)
            rRank = vApp%recvRanks(r)+1

            ! gas
            recvDataOffset = (Model%nG + kRanks(rRank)*NkpT)*Grid%Nj*Grid%Ni + &
                             (Model%nG + jRanks(rRank)*NjpT)*Grid%Ni + &
                             (Model%nG + iRanks(rRank)*NipT)
            call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, iPjPkP5Gas, &
                             vApp%recvTypesGasDeep(r), ierr)
            ! Bxyz
            call mpi_type_hindexed(1, (/1/), recvDataOffset*dataSize, iPjPkP4Bxyz, &
                             vApp%recvTypesBxyzDeep(r), ierr)
        enddo

        do r=1,SIZE(vApp%sendRanks)
            sRank = vApp%sendRanks(r)+1

            if(iRanks(sRank) .gt. 0) then
                ! never send shallow data to any rank but minimum i
                vApp%sendCountsInexyzShallow(r) = 0
                vApp%sendCountsIneijkShallow(r) = 0
                ! set these types to non null because MPI complains
                vApp%sendTypesInexyzShallow(r) = MPI_INTEGER
                vApp%sendTypesIneijkShallow(r) = MPI_INTEGER
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

        do r=1,size(vApp%sendTypesInexyzShallow)
            call mpi_type_commit(vApp%sendTypesInexyzShallow(r), ierr)
            call mpi_type_commit(vApp%sendTypesIneijkShallow(r), ierr)

            call mpi_type_commit(vApp%recvTypesGasDeep(r), ierr)
            call mpi_type_commit(vApp%recvTypesBxyzDeep(r), ierr)

            if(vApp%doDeep) then
                call mpi_type_commit(vApp%sendTypesGas0Deep(r), ierr)
            endif
        enddo

        end associate

    end subroutine createVoltDataTypes

    subroutine helpVoltron(vApp, helperQuit)
        type(voltAppMpi_T), intent(inout) :: vApp
        logical, intent(out) :: helperQuit ! should the helper quit

        integer :: ierr, helpType
        type(MPI_Request) :: helpReq

        helperQuit = .false. ! don't quit normally

        ! assumed to only be in this function if helpers are enabled

        ! async wait for command type
        call mpi_Ibcast(helpType, 1, MPI_INTEGER, 0, vApp%vHelpComm, helpReq, ierr)
        call mpi_wait(helpReq, MPI_STATUS_IGNORE, ierr)

        !write (*,*) 'Helper got request type: ', helpType

        ! not idle now
        call setIdleStatus(vApp, helpType)

        ! then call appropriate function to deal with command
        select case(helpType)
            CASE (VHSTEP)
                call vhHandleStep(vApp)
            CASE (VHSQUISHSTART)
                call vhHandleSquishStart(vApp)
            CASE (VHSQUISHEND)
                call vhHandleSquishEnd(vApp)
            case (VHQUIT)
                helperQuit = .true.
            CASE DEFAULT
                print *,"Unknown voltron helper request type received."
                stop
        end select

        ! idle now
        call setIdleStatus(vApp, 0)

    end subroutine

end module voltapp_mpi

