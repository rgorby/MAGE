! Collection of data and objects for additional voltron helper ranks

module volthelpers_mpi
    use voltmpitypes
    use mpi_f08
    use ebsquish, only : SquishBlocksRemain, DoSquishBlock
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN, IEEE_QUIET_NAN

    implicit none

    enum, bind(C)
        enumerator :: VHSTEP=1,VHSQUISHSTART,VHSQUISHEND,VHQUIT
    endenum

    contains

    ! helpers idle status functions

    function allHelpersIdle(vApp)
        type(voltAppMpi_T), intent(in) :: vApp
        logical :: allHelpersIdle

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        if (.not. vApp%useHelpers) then
            ! shouldn't be here
            allHelpersIdle = .true.
            return
        endif

        call Tic("VoltHelpers")
        call Tic("VHWait")

        ! lock the data window to read status
        call mpi_win_lock(MPI_LOCK_SHARED, 0, 0, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        allHelpersIdle = all(vApp%vHelpIdle .le. 0)

        ! unlock the data window
        call mpi_win_unlock(0, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call Toc("VHWait")
        call Toc("VoltHelpers")

    end function

    subroutine setIdleStatus(vApp, newStatus)
        type(voltAppMpi_T), intent(inout) :: vApp
        integer, intent(in) :: newStatus

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message
        integer(KIND=MPI_BASE_MYADDR) :: disp

        ! lock the data window to set my status
        call mpi_win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        disp = vApp%vHelpRank-1 ! my slot in the target's array
        call mpi_put(newStatus, 1, MPI_INTEGER, 0, disp, 1, MPI_INTEGER, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! if desired, could try to use mpi_win_flush_local here to allow this rank to move on

        ! unlock the data window
        call mpi_win_unlock(0, vApp%vHelpWin, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    ! chimp data update functions

    subroutine sendChimpStateData(ebState, vHelpComm)
        type(ebState_T), intent(in) :: ebState
        type(MPI_Comm), intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! Fields
        !  eb1
        call mpi_bcast(ebState%eb1%dB, size(ebState%eb1%dB), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb1%E, size(ebState%eb1%E), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb1%W, size(ebState%eb1%W), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb1%time, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        !  eb2
        call mpi_bcast(ebState%eb2%dB, size(ebState%eb2%dB), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb2%E, size(ebState%eb2%E), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb2%W, size(ebState%eb2%W), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb2%time, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine recvChimpStateData(ebState, vHelpComm)
        type(ebState_T), intent(inout) :: ebState
        type(MPI_Comm), intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! Fields
        !  eb1
        call mpi_bcast(ebState%eb1%dB, size(ebState%eb1%dB), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb1%E, size(ebState%eb1%E), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb1%W, size(ebState%eb1%W), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb1%time, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        !  eb2
        call mpi_bcast(ebState%eb2%dB, size(ebState%eb2%dB), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb2%E, size(ebState%eb2%E), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb2%W, size(ebState%eb2%W), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%eb2%time, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine sendChimpUpdate(vApp)
        type(voltAppMpi_T), intent(in) :: vApp

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call sendChimpStateData(vApp%ebTrcApp%ebState, vApp%vHelpComm)

        call mpi_bcast(vApp%iDeep, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%qkSquishStride, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%doQkSquish, 1, MPI_LOGICAL, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%rTrc, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%nTrc, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call mpi_bcast(vApp%ebTrcApp%ebSquish%Rinner, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine recvChimpUpdate(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call recvChimpStateData(vApp%ebTrcApp%ebState, vApp%vHelpComm)

        call mpi_bcast(vApp%iDeep, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%qkSquishStride, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%doQkSquish, 1, MPI_LOGICAL, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%rTrc, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%nTrc, 1, MPI_INTEGER, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if        
        call mpi_bcast(vApp%ebTrcApp%ebSquish%Rinner, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

! ----------------------------------------------------------

! Functions below for sending and receiving help requests

    subroutine vhRequestType(vApp, rType)
        type(voltAppMpi_T), intent(in) :: vApp
        integer, intent(in) :: rType

        integer :: ierr
        type(MPI_Request) :: helpReq

        ! async to match waiting helper nodes
        call mpi_Ibcast(rType, 1, MPI_INT, 0, vApp%vHelpComm, helpReq, ierr)
        call mpi_wait(helpReq, MPI_STATUS_IGNORE, ierr)

    end subroutine

    subroutine vhReqStep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call Tic("VHReqStep")
        call vhRequestType(vApp, VHSTEP)
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call Toc("VHReqStep")
    end subroutine

    subroutine vhHandleStep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
    end subroutine

    subroutine vhReqSquishStart(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call Tic("VHReqSquishS")
        call vhRequestType(vApp, VHSQUISHSTART)

        ! send eb data
        call sendChimpUpdate(vApp)

        ! send which sections of squish data to solve
        if(vApp%masterSquish) then
             ! 1 block per volton rank, including master and helpers
            vApp%ebTrcApp%ebSquish%myNumBlocks = 1
        else
            ! 1 block per helper, none for master
            vApp%ebTrcApp%ebSquish%myNumBlocks = 0
        endif
        call mpi_bcast(vApp%ebTrcApp%ebSquish%numSquishBlocks, 1, MPI_INT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%ebTrcApp%ebSquish%blockStartIndices, vApp%ebTrcApp%ebSquish%numSquishBlocks, MPI_INT, 0, vApp%vHelpComm, ierr)
        vApp%ebTrcApp%ebSquish%myFirstBlock = 1 ! master is always the first block, if it has one
        vApp%ebTrcApp%ebSquish%curSquishBlock = vApp%ebTrcApp%ebSquish%myFirstBlock
        call Toc("VHReqSquishS")

    end subroutine

    subroutine vhHandleSquishStart(vApp, squishTime)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), optional, intent(out) :: squishTime

        integer :: ierr
        integer :: clockStart, clockEnd, clockRate, clockMax

        call recvChimpUpdate(vApp)

        call mpi_bcast(vApp%ebTrcApp%ebSquish%numSquishBlocks, 1, MPI_INT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%ebTrcApp%ebSquish%blockStartIndices, vApp%ebTrcApp%ebSquish%numSquishBlocks, MPI_INT, 0, vApp%vHelpComm, ierr)
        vApp%ebTrcApp%ebSquish%myNumBlocks = 1 ! helpers always have 1 block
        if(vApp%masterSquish) then
            ! master voltron rank is helping, rank 1 works on block 2, etc...
            vApp%ebTrcApp%ebSquish%myFirstBlock = vApp%vHelpRank + 1
        else
            ! master voltron rank is not helping, rank 1 works on block 1, etc...
            vApp%ebTrcApp%ebSquish%myFirstBlock = vApp%vHelpRank
        endif
        
        call SquishStart(vApp)

        call system_clock(clockStart,clockRate,clockMax)
        do while(SquishBlocksRemain(vApp))
            call DoSquishBlock(vApp)
        end do
        call system_clock(count=clockEnd)
        if(present(squishTime)) then
            if(clockEnd < clockStart) then ! wrap
                squishTime = real(clockMax-clockStart+clockEnd, RP)/real(clockRate, RP)
            else
                squishTime = real(clockEnd-clockStart, RP)/real(clockRate, RP)
            endif
        endif

    end subroutine

    subroutine vhReqSquishEnd(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr,length,i,firstBlock,ks,ke, oldSizes(4), newSizes(4), offsets(4)
        type(MPI_Datatype) :: newtype
        character( len = MPI_MAX_ERROR_STRING) :: message
        real(rp), dimension(:), allocatable, save :: helperTimes

        call Tic("VHReqSquishE")
        call vhRequestType(vApp, VHSQUISHEND)

        oldSizes = shape(vApp%chmp2mhd%xyzSquish)
        newSizes = (/vApp%iDeep+2-vApp%ebTrcApp%ebState%ebGr%is, &
                     vApp%ebTrcApp%ebState%ebGr%je+2-vApp%ebTrcApp%ebState%ebGr%js, &
                     0, 2/)
        offsets = (/0, 0, 0, 0/)

        ! collect solved squish data
        if(vApp%masterSquish) then
            firstBlock = 1 ! master rank handles block 0
        else
            firstBlock = 0 ! collect all data
        endif
        do i = firstBlock,vApp%ebTrcApp%ebSquish%numSquishBlocks-1
            call GetSquishBds(vApp, ks, ke, i)
            offsets(3) = ks-1
            newSizes(3) = ke+1-ks

            call mpi_type_create_subarray(4, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, MPI_MYFLOAT, newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_type_commit(newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_recv(vApp%chmp2mhd%xyzSquish, 1, newtype, i-firstBlock+1, 78111, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_type_free(newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if

            ! isGood is the same size as xyzsquish, but only 3 dimensions
            call mpi_type_create_subarray(3, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, MPI_LOGICAL, newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_type_commit(newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_recv(vApp%chmp2mhd%isGood, 1, newtype, i-firstBlock+1, 78112, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if
            call mpi_type_free(newtype, ierr)
            if(ierr /= MPI_Success) then
                call MPI_Error_string( ierr, message, length, ierr)
                print *,message(1:length)
                call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
            end if 

        enddo

        if(vApp%squishLoadBalance) then
            if(.not. allocated(helperTimes)) then
                ! empty space in front for master
                allocate(helperTimes(1+vApp%ebTrcApp%ebSquish%numSquishBlocks))
            endif
            call mpi_gather(MPI_IN_PLACE, 0, MPI_MYFLOAT, helperTimes, 1, MPI_MYFLOAT, 0, vApp%vhelpComm, ierr)

            ! perform load balancing of squish block sizes
            call LoadBalanceBlocks(vApp, helperTimes(2:))
        endif
        call Toc("VHReqSquishE")

    end subroutine

    subroutine vhHandleSquishEnd(vApp, squishTime)
        type(voltAppMpi_T), intent(inout) :: vApp
        real(rp), optional, intent(in) :: squishTime

        integer :: ierr,length,ks,ke, oldSizes(4), newSizes(4), offsets(4)
        type(MPI_Datatype) :: newtype
        character( len = MPI_MAX_ERROR_STRING) :: message
        real(rp) :: squishTimeSendVal

        oldSizes = shape(vApp%chmp2mhd%xyzSquish)
        newSizes = (/vApp%iDeep+2-vApp%ebTrcApp%ebState%ebGr%is, &
                     vApp%ebTrcApp%ebState%ebGr%je+2-vApp%ebTrcApp%ebState%ebGr%js, &
                     0, 2/)
        offsets = (/0, 0, 0, 0/)

        ! send solved squish data
        call GetSquishBds(vApp, ks, ke, vApp%ebTrcApp%ebSquish%myFirstBlock)
        offsets(3) = ks-1
        newSizes(3) = ke+1-ks

        call mpi_type_create_subarray(4, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, MPI_MYFLOAT, newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_type_commit(newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_send(vApp%chmp2mhd%xyzSquish, 1, newtype, 0, 78111, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_type_free(newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! isGood is the same size as xyzsquish, but only 3 dimensions
        call mpi_type_create_subarray(3, oldSizes, newSizes, offsets, MPI_ORDER_FORTRAN, MPI_LOGICAL, newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_type_commit(newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_send(vApp%chmp2mhd%isGood, 1, newtype, 0, 78112, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_type_free(newtype, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        if(vApp%squishLoadBalance) then
            if(present(squishTime)) then
                squishTimeSendVal = squishTime
            else
                squishTimeSendVal = -1.0_rp ! default
            endif
            call mpi_gather([squishTimeSendval], 1, MPI_MYFLOAT, [squishTimeSendVal], 1, MPI_MYFLOAT, 0, vApp%vhelpComm, ierr)
        endif

    end subroutine

    subroutine vhReqHelperQuit(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call vhRequestType(vApp, VHQUIT)

    end subroutine

end module

