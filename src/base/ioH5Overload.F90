!Only used by ioH5, holds various repeated routines for overloading
!b/c Fortran doesn't do templating

module ioH5Overload
    use kdefs
    use ioH5types
    use strings
    
    implicit none

    !Overloader for function to copy N-rank data into 1D buffer
    interface LoadIO
        module procedure LoadIO_5D,LoadIO_4D,LoadIO_3D,LoadIO_2D,LoadIO_1D,LoadIO_0D
    end interface

	contains

!-------------------------------------------
!Simple helper functions for navigating IOVars_T
    !Finds element of IO chain matching input id-string
    function FindIO(IOVars,inStr,doFailO) result(nOut)
        type(IOVAR_T), dimension(:), intent(in) :: IOVars
        character(len=*), intent(in) :: inStr
        logical, intent(in), optional :: doFailO

        integer :: nOut, Nv,i
        logical :: doFail

        if (present(doFailO)) then
            doFail = doFailO
        else
            doFail = .false.
        endif
        nOut = -1
        Nv = size(IOVars)
        do i=1,Nv
            if (trim(toUpper(IOVars(i)%idStr)) == trim(toUpper(inStr))) then
                nOut = i
            endif
        enddo

        if ( (nOut == -1) .and. doFail ) then
            !Failed to find variable, blow this thing up
            Nv = NextIO(IOVars)
            write(*,*) 'Failed to find variable: ', trim(inStr)
            write(*,*) 'Variables present: '
            do i=1,Nv-1
                write(*,*) '-- ', trim(toUpper(IOVars(i)%idStr))
            enddo
            write(*,*) 'Exiting ...'
            stop
        endif
    end function FindIO

    !Finds first unused element of IO chain
    !FIXME: Probably a cleaner way with minloc (but maybe not standard?)
    function NextIO(IOVars) result(nOut)
        type(IOVAR_T), dimension(:), intent(in) :: IOVars

        integer :: nOut, Nv,i
        nOut = -1
        Nv = size(IOVars)
        do i=1,Nv
            if ( (.not. IOVars(i)%toRead) .and. (.not. IOVars(i)%toWrite) ) then
                nOut = i
                exit !Break out of loop
            endif
        enddo

        if (nOut == Nv) then
            !No more room left
            write(*,*) 'ERROR: Overflow on IO Chain, exiting ...'
            stop
        endif
    end function NextIO

!-------------------------------------------
!These routines reshape data of various ranks/types into array
    !Fill 1D array
    subroutine IOArray1DFill(IOVars,vID,Q)
        type(IOVAR_T), dimension(:), intent(in) :: IOVars
        character(len=*), intent(in) :: vID
        real(rp), dimension(:), intent(inout) :: Q

        integer :: nvar
        nvar = FindIO(IOVars,vID,.true.)
        if (.not. IOVars(nvar)%isDone) call FailArrayFill(vID)

        Q = IOVars(nvar)%data
    end subroutine IOArray1DFill

    !Fill 2D array
    subroutine IOArray2DFill(IOVars,vID,Q)
        type(IOVAR_T), dimension(:), intent(in) :: IOVars
        character(len=*), intent(in) :: vID
        real(rp), dimension(:,:), intent(inout) :: Q

        integer :: nvar
        nvar = FindIO(IOVars,vID,.true.)
        if (.not. IOVars(nvar)%isDone) call FailArrayFill(vID)

        Q = reshape(IOVars(nvar)%data,[IOVars(nvar)%dims(1),IOVars(nvar)%dims(2)])
    end subroutine IOArray2DFill

    !Fill 3D array
    subroutine IOArray3DFill(IOVars,vID,Q)
        type(IOVAR_T), dimension(:), intent(in) :: IOVars
        character(len=*), intent(in) :: vID
        real(rp), dimension(:,:,:), intent(inout) :: Q

        integer :: nvar
        nvar = FindIO(IOVars,vID,.true.)
        if (.not. IOVars(nvar)%isDone) call FailArrayFill(vID)

        Q = reshape(IOVars(nvar)%data,[IOVars(nvar)%dims(1),IOVars(nvar)%dims(2),IOVars(nvar)%dims(3)])
    end subroutine IOArray3DFill

    subroutine FailArrayFill(vID)
        character(len=*), intent(in) :: vID

        write(*,*) 'Failed to fill array w/ variable: ', trim(vID)
        stop
        
    end subroutine FailArrayFill
!-------------------------------------------
!These routines add data for output to IO chain
!All routines find first unused link and add there
    subroutine AddOut_5D(IOVars,idStr,Q,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:,:,:,:,:) :: Q
        character(len=*), intent(in), optional :: uStr

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)

        call LoadIO(IOVars(n),Q)

    end subroutine AddOut_5D

    subroutine AddOut_4D(IOVars,idStr,Q,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:,:,:,:) :: Q
        character(len=*), intent(in), optional :: uStr

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)

        call LoadIO(IOVars(n),Q)

    end subroutine AddOut_4D

    subroutine AddOut_3D(IOVars,idStr,Q,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:,:,:) :: Q
        character(len=*), intent(in), optional :: uStr

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)

        call LoadIO(IOVars(n),Q)

        !Catch for fake 3D (ie nk=1)
        if (IOVars(n)%dims(3) == 1) then
            IOVars(n)%Nr = 2
            IOVars(n)%dims(3) = 0
        endif

    end subroutine AddOut_3D

    subroutine AddOut_2D(IOVars,idStr,Q,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:,:) :: Q
        character(len=*), intent(in), optional :: uStr

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)

        call LoadIO(IOVars(n),Q)
    end subroutine AddOut_2D

    subroutine AddOut_1D(IOVars,idStr,Q,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:) :: Q
        character(len=*), intent(in), optional :: uStr

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)

        call LoadIO(IOVars(n),Q)
    end subroutine AddOut_1D

    subroutine AddOut_DP(IOVars,idStr,X,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(dp), intent(in) :: X
        character(len=*), intent(in), optional :: uStr

        integer :: n
        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        IOVars(n)%vType = IOREAL
        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)

        call LoadIO(IOVars(n),real(X,rp))
    end subroutine AddOut_DP

    subroutine AddOut_SP(IOVars,idStr,X,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(sp), intent(in) :: X
        character(len=*), intent(in), optional :: uStr
        
        integer :: n
        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        IOVars(n)%vType = IOREAL
        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)

        call LoadIO(IOVars(n),real(X,rp))
    end subroutine AddOut_SP

    subroutine AddOut_Int(IOVars,idStr,X,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        integer, intent(in) :: X
        character(len=*), intent(in), optional :: uStr

        integer :: n
        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        IOVars(n)%vType = IOINT
        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)

        call LoadIO(IOVars(n),real(X,rp))
    end subroutine AddOut_Int

    subroutine AddOut_Str(IOVars,idStr,X,uStr)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr, X
        character(len=*), intent(in), optional :: uStr

        integer :: n
        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        IOVars(n)%vType = IOSTR
        IOVars(n)%dStr = trim(X)
        if (present(uStr)) IOVars(n)%unitStr = trim(uStr)
        
    end subroutine AddOut_Str
!---------------------------
!Loads structured data into buffers

    subroutine LoadIO_5D(IOVar,Q)
        type(IOVAR_T), intent(inout) :: IOVar
        real(rp), intent(in), dimension(:,:,:,:,:) :: Q
        
        IOVar%Nr = 5
        IOVar%dims = 0
        
        !Find size/shape of Q
        IOVar%dims(1:IOVar%Nr) = shape(Q)
        IOVar%N = product(IOVar%dims(1:IOVar%Nr))

        !Store data into 1D array in IOVar type
        call CopyIO(IOVar%data,Q,IOVar%N)

    end subroutine LoadIO_5D

    subroutine LoadIO_4D(IOVar,Q)
        type(IOVAR_T), intent(inout) :: IOVar
        real(rp), intent(in), dimension(:,:,:,:) :: Q
        
        IOVar%Nr = 4
        IOVar%dims = 0
        
        !Find size/shape of Q
        IOVar%dims(1:IOVar%Nr) = shape(Q)
        IOVar%N = product(IOVar%dims(1:IOVar%Nr))

        !Store data into 1D array in IOVar type
        call CopyIO(IOVar%data,Q,IOVar%N)

    end subroutine LoadIO_4D

    subroutine LoadIO_3D(IOVar,Q)
        type(IOVAR_T), intent(inout) :: IOVar
        real(rp), intent(in), dimension(:,:,:) :: Q
        
        IOVar%Nr = 3
        IOVar%dims = 0
        
        !Find size/shape of Q
        IOVar%dims(1:IOVar%Nr) = shape(Q)
        IOVar%N = product(IOVar%dims(1:IOVar%Nr))

        !Store data into 1D array in IOVar type
        call CopyIO(IOVar%data,Q,IOVar%N)

    end subroutine LoadIO_3D

    subroutine LoadIO_2D(IOVar,Q)
        type(IOVAR_T), intent(inout) :: IOVar
        real(rp), intent(in), dimension(:,:) :: Q
        
        IOVar%Nr = 2
        IOVar%dims = 0
        
        !Find size/shape of Q
        IOVar%dims(1:IOVar%Nr) = shape(Q)
        IOVar%N = product(IOVar%dims(1:IOVar%Nr))

        !Store data into 1D array in IOVar type
        call CopyIO(IOVar%data,Q,IOVar%N)

    end subroutine LoadIO_2D

    subroutine LoadIO_1D(IOVar,Q)
        type(IOVAR_T), intent(inout) :: IOVar
        real(rp), intent(in), dimension(:) :: Q
        
        IOVar%Nr = 1
        IOVar%dims = 0
        
        !Find size/shape of Q
        IOVar%dims(1:IOVar%Nr) = shape(Q)
        IOVar%N = product(IOVar%dims(1:IOVar%Nr))

        !Store data into 1D array in IOVar type
        call CopyIO(IOVar%data,Q,IOVar%N)

    end subroutine LoadIO_1D

    !Handle scalar data, always use rp precision real
    subroutine LoadIO_0D(IOVar,X)
        type(IOVAR_T), intent(inout) :: IOVar
        real(rp), intent(in) :: X

        IOVar%Nr = 0
        IOVar%dims = 0
        IOVar%N = 1

        call CopyIO(IOVar%data,[X],IOVar%N)

    end subroutine LoadIO_0D

    !Workhorse routine to do copies from LoadIO
    !N is size (1D) of data
    subroutine CopyIO(IOBuffer,IOSource,N)
        integer, intent(in) :: N
        real(rp), dimension(:), allocatable, intent(inout) :: IOBuffer
        real(rp), intent(in) :: IOSource(N)

        !For now, always realloc
        !FIXME: Replace with check of current data and realloc only if necessary
        if (allocated(IOBuffer)) then
            if ( size(IOBuffer) /= N ) then
                !Allocated but wrong size
                deallocate(IOBuffer)
                allocate(IOBuffer(N))
            endif
        else !Not allocated
            allocate(IOBuffer(N))
        endif

        !Copy data
        IOBuffer = IOSource
    end subroutine CopyIO


end module ioH5Overload