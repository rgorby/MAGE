! Collection of data and objects for additional voltron helper ranks

module volthelpers_mpi
    use voltmpitypes
    use mpi
    use ebsquish, only : SquishBlocksRemain, DoSquishBlock
    use, intrinsic :: ieee_arithmetic, only: IEEE_VALUE, IEEE_SIGNALING_NAN, IEEE_QUIET_NAN

    implicit none

    enum, bind(C)
        enumerator :: VHSTEP=1,VHSQUISHSTART,VHSQUISHEND
    endenum

    contains

! Chimp specific transfer functions

    subroutine sendChimpModelData(ebModel, vHelpComm)
        type(chmpModel_T), intent(in) :: ebModel
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call mpi_bcast(ebModel, sizeof(ebModel), MPI_BYTE, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
    end subroutine

    subroutine recvChimpModelData(ebModel, vHelpComm)
        type(chmpModel_T), intent(inout) :: ebModel
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call mpi_bcast(ebModel, sizeof(ebModel), MPI_BYTE, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
    end subroutine

    subroutine sendChimpStateSizeAndData(ebState, vHelpComm)
        type(ebState_T), intent(in) :: ebState
        integer, intent(in) :: vHelpComm

        call sendChimpStateSize(ebState, vHelpComm)
        call sendChimpStateData(ebState, vHelpComm)
    end subroutine

    subroutine sendChimpStateSize(ebState, vHelpComm)
        type(ebState_T), intent(in) :: ebState
        integer, intent(in) :: vHelpComm

        integer :: ierr, length, sizes(10)
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! Grid
        sizes(1:4) = lbound(ebState%ebGr%xyz)
        sizes(5:8) = ubound(ebState%ebGr%xyz)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:4) = lbound(ebState%ebGr%xyzcc)
        sizes(5:8) = ubound(ebState%ebGr%xyzcc)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:4) = lbound(ebState%ebGr%B0cc)
        sizes(5:8) = ubound(ebState%ebGr%B0cc)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:5) = lbound(ebState%ebGr%Txi)
        sizes(6:10) = ubound(ebState%ebGr%Txi)
        call mpi_bcast(sizes, 10, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:5) = lbound(ebState%ebGr%Tix)
        sizes(6:10) = ubound(ebState%ebGr%Tix)
        call mpi_bcast(sizes, 10, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Fields
        !  eb1
        sizes(1:4) = lbound(ebState%eb1%dB)
        sizes(5:8) = ubound(ebState%eb1%dB)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:4) = lbound(ebState%eb1%E)
        sizes(5:8) = ubound(ebState%eb1%E)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:4) = lbound(ebState%eb1%W)
        sizes(5:8) = ubound(ebState%eb1%w)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        !  eb2
        sizes(1:4) = lbound(ebState%eb2%dB)
        sizes(5:8) = ubound(ebState%eb2%dB)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:4) = lbound(ebState%eb2%E)
        sizes(5:8) = ubound(ebState%eb2%E)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:4) = lbound(ebState%eb2%W)
        sizes(5:8) = ubound(ebState%eb2%W)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! Tab
        call mpi_bcast(size(ebState%ebTab%gStrs,1), 1, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:1) = lbound(ebState%ebTab%times)
        sizes(2:2) = ubound(ebState%ebTab%times)
        call mpi_bcast(sizes, 2, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:1) = lbound(ebState%ebTab%MJDs)
        sizes(2:2) = ubound(ebState%ebTab%MJDs)
        call mpi_bcast(sizes, 2, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine sendChimpStateData(ebState, vHelpComm)
        type(ebState_T), intent(in) :: ebState
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! Grid
        call mpi_bcast(ebState%ebGr, 19, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%xyz, size(ebState%ebGr%xyz), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%xyzcc, size(ebState%ebGr%xyzcc), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%B0cc, size(ebState%ebGr%B0cc), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%Txi, size(ebState%ebGr%Txi), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%Tix, size(ebState%ebGr%Tix), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

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

        ! Tab
        call mpi_bcast(ebState%ebTab%N, 1, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%bStr, strLen, MPI_CHAR, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%gStrs, size(ebState%ebTab%gStrs), MPI_CHAR, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%times, size(ebState%ebTab%times), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%MJDs, size(ebState%ebTab%MJDs), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%isMPI, 1, MPI_LOGICAL, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%hasMJD, 1, MPI_LOGICAL, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%Ri, 6, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! one leftover logical
        call mpi_bcast(ebState%doStatic, 1, MPI_LOGICAL, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine recvChimpStateSizeAndData(ebState, vHelpComm)
        type(ebState_T), intent(inout) :: ebState
        integer, intent(in) :: vHelpComm

        call recvChimpStateSize(ebState, vHelpComm)
        call recvChimpStateData(ebState, vHelpComm)
    end subroutine

    subroutine recvChimpStateSize(ebState, vHelpComm)
        type(ebState_T), intent(inout) :: ebState
        integer, intent(in) :: vHelpComm

        integer :: ierr, length, sizes(10)
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! Grid
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%ebGr%xyz(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%ebGr%xyzcc(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%ebGr%B0cc(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        call mpi_bcast(sizes, 10, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%ebGr%Txi(sizes(1):sizes(6),sizes(2):sizes(7),sizes(3):sizes(8),sizes(4):sizes(9),sizes(5):sizes(10)))
        call mpi_bcast(sizes, 10, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%ebGr%Tix(sizes(1):sizes(6),sizes(2):sizes(7),sizes(3):sizes(8),sizes(4):sizes(9),sizes(5):sizes(10)))

        ! Fields
        !  eb1
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%eb1%dB(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%eb1%E(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%eb1%W(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        !  eb2
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%eb2%dB(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%eb2%E(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%eb2%W(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))

        ! Tab
        call mpi_bcast(sizes, 1, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%ebTab%gStrs(sizes(1)))
        call mpi_bcast(sizes, 2, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%ebTab%times(sizes(1):sizes(2)))
        call mpi_bcast(sizes, 2, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(ebState%ebTab%MJDs(sizes(1):sizes(2)))

    end subroutine

    subroutine recvChimpStateData(ebState, vHelpComm)
        type(ebState_T), intent(inout) :: ebState
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! Grid
        call mpi_bcast(ebState%ebGr, 19, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%xyz, size(ebState%ebGr%xyz), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%xyzcc, size(ebState%ebGr%xyzcc), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%B0cc, size(ebState%ebGr%B0cc), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%Txi, size(ebState%ebGr%Txi), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebGr%Tix, size(ebState%ebGr%Tix), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

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

        ! Tab
        call mpi_bcast(ebState%ebTab%N, 1, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%bStr, strLen, MPI_CHAR, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%gStrs, size(ebState%ebTab%gStrs), MPI_CHAR, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%times, size(ebState%ebTab%times), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%MJDs, size(ebState%ebTab%MJDs), MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%isMPI, 1, MPI_LOGICAL, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%hasMJD, 1, MPI_LOGICAL, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(ebState%ebTab%Ri, 6, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        ! one leftover logical
        call mpi_bcast(ebState%doStatic, 1, MPI_LOGICAL, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine sendChimpSquishData(ebSquish, vHelpComm)
        type(ebSquish_T), intent(in) :: ebSquish
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call mpi_bcast(ebSquish, sizeof(ebSquish), MPI_BYTE, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
    end subroutine

    subroutine recvChimpSquishData(ebSquish, vHelpComm)
        type(ebSquish_T), intent(inout) :: ebSquish
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call mpi_bcast(ebSquish, sizeof(ebSquish), MPI_BYTE, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
    end subroutine

! Mhd Chimp Interface Transfer Functions

    subroutine sendMhd2ChmpData(mhd2Chmp, vHelpComm)
        type(mhd2Chmp_T), intent(in) :: mhd2Chmp
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call mpi_bcast(mhd2Chmp, sizeof(mhd2Chmp), MPI_BYTE, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine recvMhd2ChmpData(mhd2Chmp, vHelpComm)
        type(mhd2Chmp_T), intent(inout) :: mhd2Chmp
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call mpi_bcast(mhd2Chmp, sizeof(mhd2Chmp), MPI_BYTE, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine sendChmp2MhdSizeAndData(chmp2Mhd, vHelpComm)
        type(chmp2Mhd_T), intent(in) :: chmp2Mhd
        integer, intent(in) :: vHelpComm

        call sendChmp2MhdSize(chmp2Mhd, vHelpComm)
        call sendChmp2MhdData(chmp2Mhd, vHelpComm)
    end subroutine

    subroutine sendChmp2MhdSize(chmp2Mhd, vHelpComm)
        type(chmp2Mhd_T), intent(in) :: chmp2Mhd
        integer, intent(in) :: vHelpComm

        integer :: ierr, length, sizes(8)
        character( len = MPI_MAX_ERROR_STRING) :: message

        sizes(1:4) = lbound(chmp2Mhd%xyzSquish)
        sizes(5:8) = ubound(chmp2Mhd%xyzSquish)
        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:3) = lbound(chmp2Mhd%isGood)
        sizes(4:6) = ubound(chmp2Mhd%isGood)
        call mpi_bcast(sizes, 6, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        sizes(1:3) = lbound(chmp2Mhd%isEdible)
        sizes(4:6) = ubound(chmp2Mhd%isEdible)
        call mpi_bcast(sizes, 6, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine sendChmp2MhdData(chmp2Mhd, vHelpComm)
        type(chmp2Mhd_T), intent(in) :: chmp2Mhd
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! don't send fields, recv side will fill with default values

        call mpi_bcast(chmp2Mhd%iMax, 1, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call mpi_bcast(chmp2Mhd%epsSquish, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call mpi_bcast(chmp2Mhd%epsds0, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine recvChmp2MhdSizeAndData(chmp2Mhd, vHelpComm)
        type(chmp2Mhd_T), intent(inout) :: chmp2Mhd
        integer, intent(in) :: vHelpComm

        call recvChmp2MhdSize(chmp2Mhd, vHelpComm)
        call recvChmp2MhdData(chmp2Mhd, vHelpComm)
    end subroutine

    subroutine recvChmp2MhdSize(chmp2Mhd, vHelpComm)
        type(chmp2Mhd_T), intent(inout) :: chmp2Mhd
        integer, intent(in) :: vHelpComm

        integer :: ierr, length, sizes(8)
        character( len = MPI_MAX_ERROR_STRING) :: message

        call mpi_bcast(sizes, 8, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(chmp2Mhd%xyzSquish(sizes(1):sizes(5),sizes(2):sizes(6),sizes(3):sizes(7),sizes(4):sizes(8)))
        call mpi_bcast(sizes, 6, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(chmp2Mhd%isGood(sizes(1):sizes(4),sizes(2):sizes(5),sizes(3):sizes(6)))
        call mpi_bcast(sizes, 6, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        allocate(chmp2Mhd%isEdible(sizes(1):sizes(4),sizes(2):sizes(5),sizes(3):sizes(6)))

    end subroutine

    subroutine recvChmp2MhdData(chmp2Mhd, vHelpComm)
        type(chmp2Mhd_T), intent(inout) :: chmp2Mhd
        integer, intent(in) :: vHelpComm

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! fill fields with default values
        chmp2Mhd%xyzSquish = 0.0_rp
        chmp2Mhd%isGood = .false.
        chmp2Mhd%isEdible = .false.

        call mpi_bcast(chmp2Mhd%iMax, 1, MPI_INTEGER, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call mpi_bcast(chmp2Mhd%epsSquish, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

        call mpi_bcast(chmp2Mhd%epsds0, 1, MPI_MYFLOAT, 0, vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

! Data initialization and data update functions

    subroutine sendChimpInitialization(vApp)
        type(voltAppMpi_T), intent(in) :: vApp

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        ! model is static size, only send data
        call sendChimpModelData(vApp%ebTrcApp%ebModel, vApp%vHelpComm)
        ! state is dynamic size
        call sendChimpStateSizeAndData(vApp%ebTrcApp%ebState, vApp%vHelpComm)
        ! squish is static size, only send data
        call sendChimpSquishData(vApp%ebTrcApp%ebSquish, vApp%vHelpComm)

        ! mhd2Chmp is static size, only send data
        call sendMhd2ChmpData(vApp%mhd2Chmp, vApp%vHelpComm)
        ! chmp2Mhd is dynamic size
        call sendChmp2MhdSizeAndData(vApp%chmp2Mhd, vApp%vHelpComm)

        call mpi_bcast(vApp%imType, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%prType, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if

    end subroutine

    subroutine initializeAndReceiveChimp(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr, length
        character( len = MPI_MAX_ERROR_STRING) :: message

        call recvChimpModelData(vApp%ebTrcApp%ebModel, vApp%vHelpComm)
        call recvChimpStateSizeAndData(vApp%ebTrcApp%ebState, vApp%vHelpComm)
        call recvChimpSquishData(vApp%ebTrcApp%ebSquish, vApp%vHelpComm)
        call recvMhd2ChmpData(vApp%mhd2Chmp, vApp%vHelpComm)
        call recvChmp2MhdSizeAndData(vApp%chmp2Mhd, vApp%vHelpComm)

        call mpi_bcast(vApp%imType, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        if(ierr /= MPI_Success) then
            call MPI_Error_string( ierr, message, length, ierr)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        call mpi_bcast(vApp%prType, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
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

    end subroutine

! ----------------------------------------------------------

! Functions below for sending and receiving help requests

    subroutine vhRequestType(vApp, rType)
        type(voltAppMpi_T), intent(in) :: vApp
        integer, intent(in) :: rType

        integer :: ierr, helpReq = MPI_REQUEST_NULL

        ! async to match waiting helper nodes
        call mpi_Ibcast(rType, 1, MPI_INT, 0, vApp%vHelpComm, helpReq, ierr)
        call mpi_wait(helpReq, MPI_STATUS_IGNORE, ierr)

    end subroutine

    subroutine vhReqStep(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr
        call vhRequestType(vApp, VHSTEP)
        call mpi_bcast(vApp%tFin, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%time, 1, MPI_MYFLOAT, 0, vApp%vHelpComm, ierr)
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
        call vhRequestType(vApp, VHSQUISHSTART)

        ! send eb data
        call sendChimpUpdate(vApp)

        ! send which sections of squish data to solve
        call MPI_Comm_Size(vApp%vHelpComm, vApp%ebTrcApp%ebSquish%numSquishBlocks, ierr)
        vApp%ebTrcApp%ebSquish%myNumBlocks = 1 ! 1 block per volton rank, including master and helpers
        call mpi_bcast(vApp%ebTrcApp%ebSquish%numSquishBlocks, 1, MPI_INT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%ebTrcApp%ebSquish%myNumBlocks, 1, MPI_INT, 0, vApp%vHelpComm, ierr)
        vApp%ebTrcApp%ebSquish%myFirstBlock = vApp%vHelpRank ! work on my own rank block
        vApp%ebTrcApp%ebSquish%curSquishBlock = vApp%ebTrcApp%ebSquish%myFirstBlock

    end subroutine

    subroutine vhHandleSquishStart(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr

        call recvChimpUpdate(vApp)

        call mpi_bcast(vApp%ebTrcApp%ebSquish%numSquishBlocks, 1, MPI_INT, 0, vApp%vHelpComm, ierr)
        call mpi_bcast(vApp%ebTrcApp%ebSquish%myNumBlocks, 1, MPI_INT, 0, vApp%vHelpComm, ierr)
        vApp%ebTrcApp%ebSquish%myFirstBlock = vApp%vHelpRank ! work on my own rank block
        vApp%ebTrcApp%ebSquish%curSquishBlock = vApp%ebTrcApp%ebSquish%myFirstBlock

        do while(SquishBlocksRemain(vApp))
            call DoSquishBlock(vApp)
        end do

    end subroutine

    subroutine vhReqSquishEnd(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr,length,i,ks,ke, oldSizes(4), newSizes(4), offsets(4), newtype
        character( len = MPI_MAX_ERROR_STRING) :: message
        call vhRequestType(vApp, VHSQUISHEND)

        oldSizes = shape(vApp%chmp2mhd%xyzSquish)
        newSizes = (/vApp%iDeep+2-vApp%ebTrcApp%ebState%ebGr%is, &
                     vApp%ebTrcApp%ebState%ebGr%je+2-vApp%ebTrcApp%ebState%ebGr%js, &
                     0, 2/)
        offsets = (/0, 0, 0, 0/)

        ! collect solved squish data
        do i = 1,vApp%ebTrcApp%ebSquish%numSquishBlocks-1
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
            call mpi_recv(vApp%chmp2mhd%xyzSquish, 1, newtype, i, 78111, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
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
            call mpi_recv(vApp%chmp2mhd%isGood, 1, newtype, i, 78112, vApp%vHelpComm, MPI_STATUS_IGNORE, ierr)
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

    end subroutine

    subroutine vhHandleSquishEnd(vApp)
        type(voltAppMpi_T), intent(inout) :: vApp

        integer :: ierr,length,ks,ke, oldSizes(4), newSizes(4), offsets(4), newtype
        character( len = MPI_MAX_ERROR_STRING) :: message

        oldSizes = shape(vApp%chmp2mhd%xyzSquish)
        newSizes = (/vApp%iDeep+2-vApp%ebTrcApp%ebState%ebGr%is, &
                     vApp%ebTrcApp%ebState%ebGr%je+2-vApp%ebTrcApp%ebState%ebGr%js, &
                     0, 2/)
        offsets = (/0, 0, 0, 0/)

        ! send solved squish data
        call GetSquishBds(vApp, ks, ke, vApp%vHelpRank)
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

    end subroutine

end module

