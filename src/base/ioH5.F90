!Various routines to read/write HDF5 files

module ioH5
    use kdefs
    use ISO_C_BINDING
    use strings
    use hdf5
    use h5lt
    

    !Useful user routines
    !AddOutVar, AddInVar, WriteVars,ReadVars
    !FindIO: Example, call FindIO(IOVars,"Toy",n), IOVars(n) = Toy data/metadata
    implicit none

    integer, parameter :: MAXIODIM = 5 !Max rank for an IO array
    !Simple IO types
    enum, bind(C)
        enumerator :: IONULL=0,IOREAL,IOINT,IOLOGICAL,IOSTR
    endenum
    !Silly hard-coded value of how far to look for first Step#XXX before giving up 
    integer, parameter :: MAXSTEP0 = 10000
    logical :: doSkipG = .false. !Whether to skip (vs overwrite) existing groups

    !--------------
    !IO array type
    !Holds information for heavy input/output
    !idStr/unitStr: String description of variable and its units
    !Nr is number of ranks, and dims is size in each dimension (unused are 0)
    !data holds the information as a 1D array
    !scale is multiplier prior to output
    !vType: holds info about variable type, only partially implemented

    type IOVAR_T
        character(len=strLen) :: idStr="NONE",unitStr="CODE"
        integer :: dims(MAXIODIM) = 0 !Dimension information
        integer :: Nr = 0 !Number of ranks
        integer :: N = 0 !Total number of elements
        real(rp), dimension(:), allocatable :: data !1D holder for data
        real(rp) :: scale=1.0, renorm=0.0
        logical :: toWrite=.false.,toRead=.false. !Read or write this variable
        logical :: isDone=.false. !Whether or not variable has been successfully read/written
        integer :: vType
        character(len=strLen) :: dStr !Optional string data
    end type IOVAR_T

    !Overloader to add data (array or scalar/string) to output chain
    interface AddOutVar
        module procedure AddOut_5D,AddOut_4D,AddOut_3D,AddOut_2D,AddOut_1D,AddOut_Int,AddOut_DP,AddOut_SP,AddOut_Str
    end interface
    !Overloader for function to copy N-rank data into 1D buffer
    interface LoadIO
        module procedure LoadIO_5D,LoadIO_4D,LoadIO_3D,LoadIO_2D,LoadIO_1D,LoadIO_0D
    end interface

  !Necessary for HDF5 routines
  integer :: herror
  
contains

    !Get number of groups of form "Step#XXX" and start/end
    subroutine StepInfo(fStr,s0,sE,Nstp)
        character(len=strLen), intent(in) :: fStr
        integer, intent(out) :: s0,sE,Nstp
        integer :: herr,i
        logical :: gExist,fExist,isEnd,idFirst
        integer(HID_T) :: h5fId
        character(len=strLen) :: gStr

        call h5open_f(herr) !Setup Fortran interface
        !Open file
        call h5fopen_f(trim(fStr), H5F_ACC_RDONLY_F, h5fId, herr)

        isEnd = .false.
        idFirst = .false. !Found first step
        i = 0
        do while (.not. isEnd)
            write(gStr,'(A,I0)') "Step#", i
            !write(*,*) 'Checking ', trim(gStr)
            call h5lexists_f(h5fId,gStr,gExist,herr)
            if (.not. gExist) then
                !Step doesn't exist, either haven't found first or passed last
                if (idFirst) then
                    !Already found first, so we passed the last
                    isEnd = .true.
                else
                    !Haven't yet found first, keep going
                    i = i+1
                endif
            else
                !Group exists, check if first
                if (.not. idFirst) then
                    !write(*,*) 'Found first step at ', i
                    s0 = i
                    idFirst = .true.
                endif
                !Increment no matter what
                i = i+1
            endif

            !Give up after a while
            if (i>=MAXSTEP0) then
                write(*,*) "Error, can't read file ",fStr
                write(*,*) "Bailing out"
                stop
            endif
        enddo
        call h5fclose_f(h5fId,herr) !Close file
        call h5close_f(herr) !Close intereface
        sE = i-1
        Nstp = sE-s0+1
        !write(*,*) 'Start/Stop/Num = ', s0,sE,NStp
        
    end subroutine StepInfo

    !Find step times between s0,sE and store in pre-allocated array
    subroutine StepTimes(fStr,s0,sE,Ts)
        character(len=strLen), intent(in) :: fStr
        integer, intent(in) :: s0,sE
        real(rp), intent(inout) :: Ts(1:sE-s0+1)

        integer :: n, Nstp,herr
        character(len=strLen) :: aStr
        integer(HID_T) :: h5fId
        real(sp) :: t(1)
        Ts = 0.0
        Nstp = sE-s0+1

        !Do HDF5 prep
        call h5open_f(herr) !Setup Fortran interface
        !Open file
        call h5fopen_f(trim(fStr), H5F_ACC_RDONLY_F, h5fId, herr)

        do n=1,Nstp
            write(aStr,('(A,I0)')) "/Step#",s0+n-1
            call h5ltget_attribute_float_f(h5fId,trim(aStr),"time",t,herr)
            Ts(n) = t(1)
        enddo
        call h5fclose_f(h5fId,herr) !Close file
        call h5close_f(herr) !Close intereface

    end subroutine StepTimes

    !Count number of groups of form "Step#XXX"
    function NumSteps(fStr) result(Nstp)
        character(len=strLen), intent(in) :: fStr
        integer :: Nstp
        integer :: s0,sE

        call StepInfo(fStr,s0,sE,Nstp)
    end function NumSteps
    
    !Checks if object exists in file fStr
    !fStr:/vStr (if no gStrO), otherwise
    !fStr:/gStrO/vStr
    function ioExist(fStr,vStr,gStrO)
        character(len=*), intent(in) :: fStr,vStr
        character(len=*), intent(in), optional :: gStrO
        logical :: ioExist

        logical :: fExist, gExist
        integer(HID_T) :: h5fId, gId
        integer :: herr

        !Start w/ file
        inquire(file=fStr,exist=fExist)
        if (.not. fExist) then
            ioExist = .false.
            return
        endif

        call h5open_f(herr) !Setup Fortran interface
        !Open file
        call h5fopen_f(trim(fStr), H5F_ACC_RDONLY_F, h5fId, herr)

        !Set outId to be either gId/fId depending on group/root reading
        if (present(gStrO)) then
            !Check if group exists
            call h5lexists_f(h5fId,trim(gStrO),gExist,herr)
            if (.not. gExist) then
                ioExist = .false.
                call h5fclose_f(h5fId,herr)
                return
            endif
            !Open group
            call h5gopen_f(h5fId,trim(gStrO),gId,herr)
            call h5lexists_f(gId,trim(vStr),ioExist,herr)
            !Close up group
            call h5gclose_f(gId,herr)
        else
            !Read from root
            call h5lexists_f(h5fId,trim(vStr),ioExist,herr)
        endif !gStrO
        !Close rest
        call h5fclose_f(h5fId,herr)

        call h5close_f(herr) !Close intereface
    end function ioExist

    !Delete file if it already exists
    subroutine CheckAndKill(fStr)
        character(len=strLen), intent(in) :: fStr

        logical :: fExist

        inquire(file=trim(fStr),exist=fExist)
        if (fExist) then
            write(*,'(3a)') '<',trim(fStr),' already exists, deleting ...>'
            call EXECUTE_COMMAND_LINE('rm ' // trim(fStr) , wait=.true.)
            !write(*,*) ''
        endif
    end subroutine CheckAndKill
       
    !Read into array of IOVar from file fOut/optional group ID gStr
    subroutine ReadVars(IOVars,doIOp,baseStr,gStrO)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        logical, intent(in) :: doIOp
        character(len=strLen), intent(in) :: baseStr
        character(len=strLen), intent(in), optional :: gStrO
        character(len=strLen) :: outStr

        
        integer :: n, Nv
        ! integer(HID_T) :: Nr
        logical :: fExist,gExist,dsExist,aExist
        integer :: herr, dsTest
        integer(HID_T) :: h5fId, gId, inId
        character(len=strLen) :: h5File

        !Set filename to baseStr
        !FIXME: Correct to add .h5 to baseStr
        h5File = baseStr
        inquire(file=h5File,exist=fExist)
        if (.not. fExist) then
            write(*,*) 'Unable to open file: ', h5File
            stop
        endif

        call h5open_f(herr) !Setup Fortran interface
        !Open file
        call h5fopen_f(trim(h5File), H5F_ACC_RDONLY_F, h5fId, herr)

        !Set outId to be either gId/fId depending on group/root reading
        if (present(gStrO)) then
            !Check if group exists
            call h5lexists_f(h5fId,trim(gStrO),gExist,herr)
            if (.not. gExist) then
                write(*,*) 'Unable to find file/group = ', trim(h5File),trim(gStrO)
                stop
            endif
            !Open group
            call h5gopen_f(h5fId,trim(gStrO),gId,herr)
            inId = gId
            
        else
            !Read from root
            inId = h5fId
        endif !gStrO

        Nv = size(IOVars)
        !-----------------------
        !Do input loop
        do n=1,Nv
            if (IOVars(n)%toRead) then !Read this variable
                
                !Check if this variable is a dataset or attribute
                call h5aexists_f(inId,trim(IOVars(n)%idStr),aExist,herr)

                dsTest = h5ltfind_dataset_f(inId,trim(IOVars(n)%idStr))
                if (dsTest == 1) then
                    dsExist = .true.
                else
                    dsExist = .false.
                endif

                if (aExist) then
                    !Attribute
                    call ReadHDFAtt(IOVars(n),inId)
                endif

                if (dsExist) then
                    !Dataset
                    call ReadHDFVar(IOVars(n),inId)
                endif
  
            endif
        enddo        

        !-----------------------
        !Now done, close up shop
        if (present(gStrO)) then
            call h5gclose_f(gId,herr) !Close group
        endif
        call h5fclose_f(h5fId,herr) !Close file
        call h5close_f(herr) !Close intereface

    end subroutine ReadVars

    !Write array of IOVar to file fOut (from baseStr), under (optional) group gStrO
    !If gStrO unspecified, written to root of HDF5
    !doIOp is whether to use IOP (ie output slice) or rp (ie restart)
    subroutine WriteVars(IOVars,doIOp,baseStr,gStrO,gStrOO)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        logical, intent(in) :: doIOp
        character(len=*), intent(in) :: baseStr
        character(len=*), intent(in), optional :: gStrO,gStrOO

        logical :: fExist,gExist
        integer :: herr
        integer(HID_T) :: h5fId, gId, ggId,outId
        character(len=strLen) :: h5File

        !Set filename to baseStr
        !FIXME: Correct to add .h5 to baseStr
        h5File = baseStr

        call h5open_f(herr) !Setup Fortran interface

        !Start by opening file, create if necessary
        inquire(file=h5File,exist=fExist)
        if (fExist) then
            !Open file
            call h5fopen_f(trim(h5File), H5F_ACC_RDWR_F, h5fId, herr)
        else
            !Create file
            call h5fcreate_f(trim(h5File),H5F_ACC_TRUNC_F, h5fId, herr)            
        endif

        !Figure out output location (outId) and create groups as necessary
        if (present(gStrO)) then
            !Write to group
            !Check if group already exists
            call h5lexists_f(h5fId,trim(gStrO),gExist,herr)
            if (gExist .and. .not. present(gStrOO)) then
                !Group exists (and not writing to subgroup)
                !Can either skip it, or kill and replace
                if (doSkipG) then
                    !Group exists, close up and skip it
                    write(*,*) 'Skipping due to pre-existing group ', trim(gStrO)
                    
                    !Close up
                    call h5fclose_f(h5fId, herr)
                    call h5close_f(herr)
                    return
                else
                    !Kill group
                    write(*,*) 'Overwriting group ', trim(gStrO)
                    call h5ldelete_f(h5fId,trim(gStrO),herr)
                    !Reset gExist and let next block of code recreate it
                    gExist = .false.
                endif
            endif
            if (.not. gExist) then
                !Create group
                call h5gcreate_f(h5fId,trim(gStrO),gId,herr)
            else
                !Open group
                call h5gopen_f(h5fId,trim(gStrO),gId,herr)
            endif

            if (present(gStrOO)) then
                !Create subgroup
                call h5gcreate_f(gId,trim(gStrOO),ggId,herr)
                outId = ggId
            else
                outId = gId 
            endif
        else
            !Write to root
            outId = h5fId
            
        endif !gStrO

        !Do writing
        call WriteVars2ID(IOVars,outId,doIOp)

        !Now done, close up shop
        if (present(gStrOO)) call h5gclose_f(ggId,herr)
        if (present(gStrO )) call h5gclose_f( gId,herr)

        call h5fclose_f(h5fId,herr) !Close file
        call h5close_f(herr) !Close intereface

    end subroutine WriteVars

    subroutine WriteVars2ID(IOVars,outId,doIOp)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        logical, intent(in) :: doIOp
        integer(HID_T), intent(in) :: outId

        integer :: n, Nv
        integer(HID_T) :: Nr

        Nv = size(IOVars)

        do n=1,Nv
            if (IOVars(n)%toWrite) then
                !Variable is ready, write it

                !Treat as scalar attribute if Nr = 0
                Nr = IOVars(n)%Nr

                if (Nr == 0) then
                !Scalar attribute
                    call WriteHDFAtt(IOVars(n),outId)
                else
                !N-rank array
                    !Create data space, use rank/dim info from IOVar
                    call WriteHDFVar(IOVars(n),outId,doIOP)

                endif !Nr=0
            endif !isSet
        enddo

    end subroutine WriteVars2ID
!-------------------------------------------
!Routines to read/write individual variables or attributes
    !FIXME: Assuming double precision for input binary data
    subroutine ReadHDFVar(IOVar,gId)
        type(IOVAR_T), intent(inout) :: IOVar
        integer(HID_T), intent(in) :: gId

        integer(HSIZE_T), allocatable, dimension(:) :: dims
        integer :: Nr,herr, N
        integer :: typeClass
        integer(SIZE_T) :: typeSize

        !Start by getting rank, dimensions and total size
        call h5ltget_dataset_ndims_f(gId,trim(IOVar%idStr),Nr,herr)
        allocate(dims(Nr))
        call h5ltget_dataset_info_f(gId,trim(IOVar%idStr), dims, typeClass, typeSize, herr)
        N = product(dims)

        !Set relevant values in input chain and allocate memory
        IOVar%Nr = Nr
        IOVar%N  = N
        IOVar%dims(1:Nr) = dims

        if (allocated(IOVar%data)) then
            deallocate(IOVar%data)
        endif
        allocate(IOVar%data(N))

        !Read based on data type
        select case(IOVar%vType)
        case(IONULL,IOREAL)
            call h5ltread_dataset_f(gId,trim(IOVar%idStr),H5T_NATIVE_DOUBLE,IOVar%data,dims,herr)
        case default
            write(*,*) 'Unknown HDF data type, bailing ...'
            stop
        end select
        IOVar%isDone = .true.
    end subroutine ReadHDFVar

    !FIXME: Add scaling to attributes
    subroutine ReadHDFAtt(IOVar,gId)
        type(IOVAR_T), intent(inout) :: IOVar
        integer(HID_T), intent(in) :: gId

        integer :: herr
        real(rp) :: X

        !Make space
        if (allocated(IOVar%data)) then
            deallocate(IOVar%data)
        endif
        allocate(IOVar%data(1))

        !Write based on data type
        select case(IOVar%vType)
        case(IONULL,IOREAL)
            IOVar%data(1) = readRealHDF(gId,trim(IOVar%idStr))
        case(IOINT)
            IOVar%data(1) = 1.0*readIntHDF(gId,trim(IOVar%idStr))
        case(IOSTR)
            write(*,*) 'Read HDF string not implemented'
            stop
        case default
            write(*,*) 'Unknown HDF data type, bailing ...'
            stop
        end select
        !write(*,*) 'Read attribute ', IOVar%data(1)
        IOVar%isDone = .true.
    end subroutine ReadHDFAtt

    !FIXME: Add scaling to attributes
    subroutine WriteHDFAtt(IOVar,gId)
        type(IOVAR_T), intent(inout) :: IOVar
        integer(HID_T), intent(in) :: gId

        integer :: herr
        real(rp) :: X

        !Write based on data type
        select case(IOVar%vType)
        case(IONULL,IOREAL)
            X = IOVar%data(1)
            call writeReal2HDF(gId,trim(IOVar%idStr),X)
        case(IOINT)
            X = IOVar%data(1)
            call writeInt2HDF(gId,trim(IOVar%idStr),int(X))
        case(IOSTR)
            call writeString2HDF(gId,trim(IOVar%idStr),trim(IOVar%dStr))
            !call h5ltmake_dataset_string_f(gID,trim(IOVar%idStr),trim(IOVar%dStr),herr)
        case default
            write(*,*) 'Unknown HDF data type, bailing ...'
            stop
        end select
        IOVar%isDone = .true.
    end subroutine WriteHDFAtt

    !FIXME: Assuming IOP is single and double otherwise
    !FIXME: Add variable attributes (units, scaling, etc)
    subroutine WriteHDFVar(IOVar,gId,doIOPO)
        type(IOVAR_T), intent(inout) :: IOVar
        integer(HID_T), intent(in) :: gId
        logical, intent(in), optional :: doIOPO

        logical :: doIOP !Do IO precision for reals
        integer :: Nr
        integer :: herr
        integer(HSIZE_T) :: h5dims(MAXIODIM)
        real(rp) :: vScl

        Nr = IOVar%Nr
        h5dims(1:Nr) = IOVar%dims(1:Nr)
        vScl = IOVar%scale

        if (present(doIOPO)) then
            doIOP = doIOPO
        else
            doIOP = .true.
        endif

        !Write based on data type
        select case(IOVar%vType)
        case(IONULL,IOREAL)
            !Assume real by default
            if (doIOP) then
                call h5ltmake_dataset_float_f(gId,trim(IOVar%idStr),Nr,h5dims(1:Nr),real(vScl*IOVar%data,sp),herr)
            else
                call h5ltmake_dataset_double_f(gId,trim(IOVar%idStr),Nr,h5dims(1:Nr),real(vScl*IOVar%data,dp),herr)
            endif
        case(IOINT)
            call h5ltmake_dataset_int_f(gId,trim(IOVar%idStr),Nr,h5dims(1:Nr),int(vScl*IOVar%data),herr)
        case default
            write(*,*) 'Unknown HDF data type, bailing ...'
            stop
        end select

        !Set attached values
        call h5ltset_attribute_string_f(gId,trim(IOVar%idStr),"Units",trim(IOVar%unitStr),herr)
        IOVar%isDone = .true.
    end subroutine WriteHDFVar

!-------------------------------------------
!These routines add data for input to IO chain
    subroutine AddInVar(IOVars,idStr,vTypeO,vSclO)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        integer, intent(in), optional :: vTypeO
        real(rp), intent(in), optional :: vSclO
        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toRead = .true.
        IOVars(n)%idStr = trim(idStr)

        if (present(vTypeO)) then
            IOVars(n)%vType = vTypeO
        else
            IOVars(n)%vType = IONULL
        endif
        if (present(vSclO)) then
            IOVars(n)%scale = vSclO
        else
            IOVars(n)%scale = 1.0
        endif

    end subroutine AddInVar
!-------------------------------------------
!These routines add data for output to IO chain
!All routines find first unused link and add there

    subroutine AddOut_5D(IOVars,idStr,Q)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:,:,:,:,:) :: Q

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        call LoadIO(IOVars(n),Q)

    end subroutine AddOut_5D

    subroutine AddOut_4D(IOVars,idStr,Q)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:,:,:,:) :: Q

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        call LoadIO(IOVars(n),Q)

    end subroutine AddOut_4D

    subroutine AddOut_3D(IOVars,idStr,Q)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:,:,:) :: Q

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        call LoadIO(IOVars(n),Q)

        !Catch for fake 3D (ie nk=1)
        if (IOVars(n)%dims(3) == 1) then
            IOVars(n)%Nr = 2
            IOVars(n)%dims(3) = 0
        endif

    end subroutine AddOut_3D

    subroutine AddOut_2D(IOVars,idStr,Q)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:,:) :: Q

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        call LoadIO(IOVars(n),Q)
    end subroutine AddOut_2D

    subroutine AddOut_1D(IOVars,idStr,Q)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(rp), intent(in), dimension(:) :: Q

        integer :: n

        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)

        call LoadIO(IOVars(n),Q)
    end subroutine AddOut_1D

    subroutine AddOut_DP(IOVars,idStr,X)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(dp), intent(in) :: X
        
        integer :: n
        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        IOVars(n)%vType = IOREAL
        call LoadIO(IOVars(n),real(X,rp))
    end subroutine AddOut_DP

    subroutine AddOut_SP(IOVars,idStr,X)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        real(sp), intent(in) :: X
        
        integer :: n
        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        IOVars(n)%vType = IOREAL
        call LoadIO(IOVars(n),real(X,rp))
    end subroutine AddOut_SP

    subroutine AddOut_Int(IOVars,idStr,X)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr
        integer, intent(in) :: X
        
        integer :: n
        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        IOVars(n)%vType = IOINT
        call LoadIO(IOVars(n),real(X,rp))
    end subroutine AddOut_Int

    subroutine AddOut_Str(IOVars,idStr,X)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        character(len=*), intent(in) :: idStr, X
        
        integer :: n
        !Find first unused
        n = NextIO(IOVars)

        IOVars(n)%toWrite = .true.
        IOVars(n)%idStr = trim(idStr)
        IOVars(n)%vType = IOSTR
        IOVars(n)%dStr = trim(X)

    end subroutine AddOut_Str

!Clears info/memory from an IO chain
    subroutine ClearIO(IOVars)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars

        integer :: i,Nv

        Nv = size(IOVars)
        IOVars(:)%idStr   = "NONE"
        IOVars(:)%unitStr = "CODE"
        IOVars(:)%dStr    = "UNSET"
        IOVars(:)%Nr      = 0
        IOVars(:)%N       = 0
        IOVars(:)%vType   = IONULL
        IOVars(:)%scale   = 1.0
        IOVars(:)%renorm  = 0.0
        IOVars(:)%toWrite = .false.
        IOVars(:)%toRead  = .false.
        IOVars(:)%isDone  = .false.
        
        

        do i=1,Nv
            if (allocated(IOVars(i)%data)) then
                deallocate(IOVars(i)%data)
            endif
        enddo
    end subroutine ClearIO

!Finds element of IO chain matching input id-string
    function FindIO(IOVars,inStr) result(nOut)
        type(IOVAR_T), dimension(:), intent(in) :: IOVars
        character(len=*), intent(in) :: inStr

        integer :: nOut, Nv,i
        nOut = -1
        Nv = size(IOVars)
        do i=1,Nv
            if (trim(toUpper(IOVars(i)%idStr)) == trim(toUpper(inStr))) then
                nOut = i
            endif
        enddo
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

    end function NextIO
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
        if (allocated(IOBuffer)) deallocate(IOBuffer)

        allocate(IOBuffer(N))

        !Copy data
        IOBuffer = IOSource
        !write(*,*) 'Copying N elements, N = ', N
    end subroutine CopyIO

!-----------------------------
!HDF 5 helper routines
    !Converts Fortran real kind to HDF5 precision
    function H5Precision(h5p) result(h5gReal)
        integer, intent(in) :: h5p
        integer(HID_T):: h5gReal
        
        !Set precision
        if (h5p == dp) then
            h5gReal = H5T_NATIVE_DOUBLE
        else
            h5gReal = H5T_NATIVE_REAL
        endif

    end function H5Precision

!-----------------------------
!Write attributes
 
    !Write integer scalar to attribute
    subroutine writeInt2HDF(gId,vId,datIn)
        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        integer, intent(in) :: datIn
        
        integer(HID_T) :: dspcId, atId
        integer(HSIZE_T), dimension(1) :: dims

        dims(1) = 0
        !Create data space
        call h5screate_f(H5S_SCALAR_F, dspcId, herror)
    
        !Create attribute and write data
        !Weird optional argument issue here
        call h5acreate_f(loc_id=gId, name=vId, type_id=H5T_NATIVE_INTEGER, space_id=dspcId, attr_id=atId,hdferr=herror )
        call h5awrite_f(atId, H5T_NATIVE_INTEGER, datIn, dims, herror)
    
        !Close up shop
        call h5aclose_f(atId, herror)
        call h5sclose_f(dspcId,herror)

    end subroutine writeInt2HDF

    !Writes a single scalar as attribute to a group/root
    !Always uses RP precision
    subroutine writeReal2HDF(gId,vId,datIn)

        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        real(rp), intent(in) :: datIn
    
        integer(HID_T) :: dspcId, atId
        integer(HID_T):: h5gReal
        integer(HSIZE_T), dimension(1) :: dims
    

        h5gReal = H5Precision(rp)
        dims(1) = 0
    
        !Create data space
        call h5screate_f(H5S_SCALAR_F, dspcId, herror)
    
        !Create attribute and write data
        !Weird optional argument issue here
        call h5acreate_f(loc_id=gId, name=vId, type_id=h5gReal, space_id=dspcId, attr_id=atId,hdferr=herror )
        call h5awrite_f(atId, h5gReal, datIn, dims, herror)
    
        !Close up shop
        call h5aclose_f(atId, herror)
        call h5sclose_f(dspcId,herror)

    end subroutine writeReal2HDF

    subroutine writeString2HDF(gId,vId,data)
        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        character(len=*), intent(in) :: data        
        integer(HID_T) :: dspcId, atId, strId
        integer(HSIZE_T), dimension(1) :: dims

        dims(1) = len(trim(data))

        !Create data space/string type
        call h5screate_f(H5S_SCALAR_F, dspcId, herror)
        call h5tcopy_f(H5T_FORTRAN_S1,strId,herror)
        call h5tset_size_f(strId,dims(1),herror)
    
        !Create attribute and write data
        call h5acreate_f(loc_id=gId, name=vId, type_id=strId, space_id=dspcId, attr_id=atId,hdferr=herror )
    
        call h5awrite_f(atId, strId, trim(data), dims, herror)
    
        !Close up shop
        call h5aclose_f(atId, herror)
        call h5sclose_f(dspcId,herror)

    end subroutine writeString2HDF

!-----------------------------
!Read attributes

    !Read integer from HDF5 attribute
    function readIntHDF(gId,vId,vDefOpt) result(vOut)
        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        integer, intent(in), optional :: vDefOpt

        integer :: vOut
        integer :: vDef
        integer(HSIZE_T), dimension(1) :: dims
        integer(HID_T) :: attrId

        if (present(vDefOpt)) then
            vDef = vDefOpt
        else
            vDef = 0
        endif

        dims(1) = 1
        call h5aopen_f(gId,vId,attrId,herror)
        if (herror /= 0) then
            vOut = vDef
        else
            call h5aread_f(attrId,H5T_NATIVE_INTEGER,vOut,dims,herror)
        endif
        call h5aclose_f(attrId,herror)
    end function readIntHDF
    
    !Read real (rp) from HDF5 attribute
    function readRealHDF(gId,vId,vDefOpt) result(vOut)
        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        real(rp), intent(in), optional :: vDefOpt

        real(rp) :: vOut,vDef
        integer(HSIZE_T), dimension(1) :: dims
        integer(HID_T) :: attrId, h5gReal

        if (present(vDefOpt)) then
            vDef = vDefOpt
        else
            vDef = 0
        endif

        dims(1) = 1
        h5gReal = H5Precision(rp)

        call h5aopen_f(gId,vId,attrId,herror)
        if (herror /= 0) then
            vOut = vDef
        else
            call h5aread_f(attrId,h5gReal,vOut,dims,herror)
        endif
        call h5aclose_f(attrId,herror)
    end function readRealHDF

end module ioH5
