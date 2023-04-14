!> Various routines to read/write HDF5 files
!> @notes     
!>  Useful user routines
!>  AddOutVar, AddInVar, WriteVars,ReadVars
!>  FindIO: Example, call FindIO(IOVars,"Toy",n), IOVars(n) = Toy data/metadata
module ioH5
    use kdefs
    use ISO_C_BINDING
    use strings
    use hdf5
    use h5lt
    use ioH5Types
    use ioH5Overload
    use files
    use dates
#ifdef __ENABLE_ZFP
    use h5zzfp_props_f
#endif
    implicit none

#ifdef __ENABLE_COMPRESS
    !> Compression Mode variables
    enum, bind(C)
        enumerator :: ZFP,ZSTD,ZLIB,SZIP,BLOSC
    end enum

    integer, private :: cacheSize = 0
    integer(HSIZE_T), parameter, private ::  CACHE_CHUNK_SIZE = 256
    logical, parameter, private :: ENABLE_COMPRESS = .TRUE.
    integer, parameter, private :: COMPRESS_ZSTD = 32015
    integer, parameter, private :: COMPRESS_ZFP = 32013
    integer, parameter, private :: COMPRESS_BLOSC = 32026
    integer, parameter, private :: H5Z_FLAG_MANDATORY = INT(Z'00000000')
    ! Compression Algorithm Selection
#ifdef __ENABLE_SZIP
    integer, parameter, private :: Z_ALG = SZIP
#endif
#ifdef __ENABLE_ZLIB
    integer, parameter, private :: Z_ALG = ZLIB
#endif
#ifdef __ENABLE_ZFP
    integer, parameter, private :: Z_ALG = ZFP
    ! ZFP compression parameters (defaults taken from ZFP header)
    integer(C_INT) :: zfpmode = 5 !1=rate, 2=prec, 3=acc, 4=expert, 5=lossless
    real(dp) :: rate = 4.0_dp
    real(dp) :: acc = 0.00001_dp
    integer(C_INT) :: prec = 11
    integer(C_INT) :: dim = 0
    ! Used for Expert Mode
    integer(C_INT) :: minbits = 1
    integer(C_INT) :: maxbits = 32768
    integer(C_INT) :: maxprec = 64
    integer(C_INT) :: minexp = -1074
#endif
#ifdef __ENABLE_ZSTD
    integer, parameter, private :: Z_ALG = ZSTD
#endif
#ifdef __ENABLE_BLOCSC
    integer, parameter, private :: Z_ALG = BLOCSC
#endif

#else
    logical, parameter, private :: ENABLE_COMPRESS = .FALSE.
#endif

    !Overloader to add data (array or scalar/string) to output chain
    interface AddOutVar
        module procedure AddOut_5D,AddOut_4D,AddOut_3D,AddOut_2D,AddOut_1D,AddOut_Int,AddOut_DP,AddOut_SP,AddOut_Str
    end interface

    !Overloader to fill array from already read IO chain
    !TODO: Add more rank/datatype options
    !TODO: Add bounds checking
    interface IOArrayFill
        module procedure IOArray1DFill,IOArray2DFill,IOArray3DFill,IOArray4DFill
    end interface IOArrayFill

contains

    !> Routine to stamp output file with various information
    subroutine StampIO(fIn)
        character(len=*), intent(in) :: fIn
        character(len=strLen) :: gStr,dtStr,bStr
        type(IOVAR_T), dimension(10) :: IOVars

        !> Check if this file has already been stamped (has githash)
        !> @NOTE: This creates ambiguity when doing restarts using binaries w/ different hashes
        !> But what're you gonna do?
        if ( ioExist(fIn,"GITHASH") ) then
            !Already stamped, let's get outta here
            return
        endif

        call GitHash(gStr)
        call GitBranch(bStr)
        call DateTimeStr(dtStr)

        call ClearIO(IOVars)
        call AddOutVar(IOVars,"GITHASH",gStr)
        call AddOutVar(IOVars,"GITBRANCH",bStr)
#ifdef __INTEL_COMPILER_OLD
        call AddOutVar(IOVars,"COMPILER",gStr)
        call AddOutVar(IOVars,"COMPILEROPTS",gStr)
#else
        call AddOutVar(IOVars,"COMPILER",compiler_version())
        call AddOutVar(IOVars,"COMPILEROPTS",compiler_options())
#endif   
        call AddOutVar(IOVars,"DATETIME",dtStr)
        call WriteVars(IOVars,.true.,fIn,doStampCheckO=.false.)
    end subroutine StampIO

    !> Calculate cdims for a given dataset size
    subroutine chunk_size(elsize, dims, cdims)
        !> Size of elements in bytes
        integer, intent(in)                          :: elsize
        !> Dimensions of dataset
        integer(HSIZE_T), dimension(*), intent(in)   :: dims
        !>  Maximum allowed chunk size/dims
        integer(HSIZE_T), dimension(:), intent(out)  :: cdims
           
        real, parameter :: max_bytes = 1048576
        !1048576 bytes {64,64,32} double
        !2097152 bytes {64,64,64} double 
        !1048576 bytes {64,64,64} float
        !2097152 bytes {64,64,128} float
    
        integer :: d
        real(dp) :: nbytes

        d = size(cdims)
        cdims(1:d) = dims(1:d)
        nbytes = elsize*product(1.d0*cdims)
        do while ((nbytes.gt.max_bytes).and.(d.gt.0))
            cdims(d) = max(floor(cdims(d)*max_bytes/nbytes,Int32),1)
            nbytes = elsize*product(1.d0*cdims)
            d = d - 1
        enddo
    
    end subroutine chunk_size

    ! -------------------------------------------   
    ! Various routines to quickly pull scalars from IOVar_T
    !
    !> Helper funciton to pull INT from an IOVar data
    function GetIOInt(IOVars,vID) result(vOut)
        type(IOVAR_T), dimension(:), intent(in) :: IOVars
        character(len=*), intent(in) :: vID
        integer :: nvar
        integer :: vOut
        nvar = FindIO(IOVars,vID,.true.)
        vOut = IOVars(nvar)%data(1)
    end function GetIOInt

    !> Helper funciton to pull REAL from an IOVar data
    function GetIOReal(IOVars,vID) result(vOut)
        type(IOVAR_T), dimension(:), intent(in) :: IOVars
        character(len=*), intent(in) :: vID
        integer :: nvar

        real(rp) :: vOut
        nvar = FindIO(IOVars,vID,.true.)
        vOut = IOVars(nvar)%data(1)
    end function GetIOReal

    !-------------------------------------------   
    !Various routines to get information about step structure of H5 file

    !> Determine grid size from HDF5 file
    function GridSizeH5(fIn) result(Nijk)
        !> File Name
        character(len=*), intent(in) :: fIn
        integer, dimension(NDIM) :: Nijk

        integer :: herr, Nr
        integer(HID_T) :: h5fId
        integer(HSIZE_T), allocatable, dimension(:) :: dims
        integer :: typeClass
        integer(SIZE_T) :: typeSize

        Nijk = 0

        call CheckFileOrDie(fIn,"Grid file not found ...")

        if ( ioExist(fIn,"X") ) then
            !Found grid variable, get size from it

            call h5open_f(herr) !Setup Fortran interface
            !Open file
            call h5fopen_f(trim(fIn), H5F_ACC_RDONLY_F, h5fId, herr)

            !Start by getting rank, dimensions and total size
            call h5ltget_dataset_ndims_f(h5fId,"X",Nr,herr)
            if (Nr > NDIM) then
                write(*,*) 'Grid X too big, rank > 3! What kinda nonsense are you trying to pull?'
                stop
            endif
            allocate(dims(Nr))
            call h5ltget_dataset_info_f(h5fId,"X", dims, typeClass, typeSize, herr)
            Nijk(1:Nr) = dims(1:Nr)
            
            !Close up shop
            call h5fclose_f(h5fId,herr) !Close file
            call h5close_f(herr) !Close intereface
        else
            write(*,*) 'Unable to find X in grid file, ', trim(fIn)
            stop
        endif
    end function GridSizeH5

    !> Get number of groups of form "Step#XXX" and start/end
    subroutine StepInfo(fStr,s0,sE,Nstp)
        !> File Name
        character(len=*), intent(in) :: fStr
        !> Step Number
        integer, intent(out) :: s0
        !> Step End
        integer, intent(out) :: sE
        !> Number of Steps
        integer, intent(out) :: Nstp
        integer :: herr,i
        logical :: gExist,fExist,isEnd,idFirst
        integer(HID_T) :: h5fId
        character(len=strLen) :: gStr

        call CheckFileOrDie(fStr,"Unable to open file")
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

    !> 
    subroutine StepAttributes_Int(fStr,attr,s0,sE,attrs)
        character(len=*), intent(in) :: fStr
        character(len=*), intent(in) :: attr
        integer, intent(in) :: s0,sE
        real(rp), intent(inout) :: attrs(1:sE-s0+1)

        integer :: n, Nstp,herr
        character(len=strLen) :: aStr
        integer(HID_T) :: h5fId
        integer:: att(1)
        attrs = 0
        Nstp = sE-s0+1

        call CheckFileOrDie(fStr,"Unable to open file")
        !Do HDF5 prep
        call h5open_f(herr) !Setup Fortran interface
        !Open file
        call h5fopen_f(trim(fStr), H5F_ACC_RDONLY_F, h5fId, herr)

        do n=1,Nstp
            write(aStr,('(A,I0)')) "/Step#",s0+n-1
            call h5ltget_attribute_int_f(h5fId,trim(aStr),trim(attr),att,herr)
            attrs(n) = att(1)
        enddo
        call h5fclose_f(h5fId,herr) !Close file
        call h5close_f(herr) !Close intereface

    end subroutine StepAttributes_Int

    !> 
    subroutine StepAttributes_Real(fStr,attr,s0,sE,attrs)
        character(len=*), intent(in) :: fStr
        character(len=*), intent(in) :: attr
        integer, intent(in) :: s0,sE
        real(rp), intent(inout) :: attrs(1:sE-s0+1)

        integer :: n, Nstp,herr
        character(len=strLen) :: aStr
        integer(HID_T) :: h5fId
        real(sp):: att(1)
        attrs = 0
        Nstp = sE-s0+1

        call CheckFileOrDie(fStr,"Unable to open file")
        !Do HDF5 prep
        call h5open_f(herr) !Setup Fortran interface
        !Open file
        call h5fopen_f(trim(fStr), H5F_ACC_RDONLY_F, h5fId, herr)

        do n=1,Nstp
            write(aStr,('(A,I0)')) "/Step#",s0+n-1
            call h5ltget_attribute_float_f(h5fId,trim(aStr),trim(attr),att,herr)
            attrs(n) = att(1)
        enddo
        call h5fclose_f(h5fId,herr) !Close file
        call h5close_f(herr) !Close intereface

    end subroutine StepAttributes_Real

    !> Find step times between s0,sE and store in pre-allocated array
    subroutine StepTimes(fStr,s0,sE,Ts)
        character(len=*), intent(in) :: fStr
        integer, intent(in) :: s0,sE
        real(rp), intent(inout) :: Ts(1:sE-s0+1)

        integer :: n, Nstp,herr
        character(len=strLen) :: aStr
        integer(HID_T) :: h5fId
        real(sp) :: t(1)
        Ts = 0.0
        Nstp = sE-s0+1

        call CheckFileOrDie(fStr,"Unable to open file")
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

    !> Same as StepTimes but for MJDs
    subroutine StepMJDs(fStr,s0,sE,MJDs)
        character(len=*), intent(in) :: fStr
        integer, intent(in) :: s0,sE
        real(rp), intent(inout) :: MJDs(1:sE-s0+1)
        integer :: n, Nstp,herr
        character(len=strLen) :: gStr
        integer(HID_T) :: h5fId,gId
        logical :: aX
        real(rp) :: M(1)

        MJDs = 0.0
        
        Nstp = sE-s0+1

        call CheckFileOrDie(fStr,"Unable to open file")
        !Do HDF5 prep
        call h5open_f(herr) !Setup Fortran interface
        !Open file
        call h5fopen_f(trim(fStr), H5F_ACC_RDONLY_F, h5fId, herr)

        do n=1,Nstp
            write(gStr,('(A,I0)')) "/Step#",s0+n-1
            call h5gopen_f(h5fId,trim(gStr),gId,herr)
            call h5aexists_f(gId,"MJD",aX,herr)
            if (aX) then
                call h5ltget_attribute_double_f(h5fId,trim(gStr),"MJD",M,herr)
            else
                M = -TINY
            endif
            
            MJDs(n) = M(1)
        enddo

        call h5gclose_f(gId,herr)
        call h5fclose_f(h5fId,herr) !Close file
        call h5close_f(herr) !Close intereface

    end subroutine StepMJDs

    !> Count number of groups of form "Step#XXX"
    function NumSteps(fStr) result(Nstp)
        character(len=*), intent(in) :: fStr
        integer :: Nstp
        integer :: s0,sE

        call CheckFileOrDie(fStr,"Unable to open file")
        call StepInfo(fStr,s0,sE,Nstp)
    end function NumSteps
    
    !> Checks if object exists in file fStr
    !> fStr:/vStr (if no gStrO), otherwise
    !> fStr:/gStrO/vStr
    function ioExist(fStr,vStr,gStrO)
        character(len=*), intent(in) :: fStr,vStr
        character(len=*), intent(in), optional :: gStrO
        logical :: ioExist

        logical :: fExist, gExist, dsExist,atExist
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
                call h5close_f(herr) !Close intereface
                return
            endif
            !Open group
            call h5gopen_f(h5fId,trim(gStrO),gId,herr)
            !Check both dataset/attribute
            call h5lexists_f(gId,trim(vStr),dsExist,herr)
            call h5aexists_f(gId,trim(vStr),atExist,herr)
            !Close up group
            call h5gclose_f(gId,herr)
        else
            !Read from root
            call h5lexists_f(h5fId,trim(vStr),dsExist,herr)
            call h5aexists_f(h5fId,trim(vStr),atExist,herr)
        endif !gStrO
        
        !Close rest
        call h5fclose_f(h5fId,herr)
        call h5close_f(herr) !Close intereface
        ioExist = dsExist .or. atExist
    end function ioExist
       
    !> Read into array of IOVar from file fOut/optional group ID gStr
    subroutine ReadVars(IOVars,doIOp,baseStr,gStrO)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        logical, intent(in) :: doIOp
        character(len=*), intent(in) :: baseStr
        character(len=*), intent(in), optional :: gStrO
        character(len=strLen) :: outStr
        
        integer :: n, Nv
        ! integer(HID_T) :: Nr
        logical :: fExist, gExist, dsExist = .false., aExist = .false.
        integer :: herr, dsTest
        integer(HID_T) :: h5fId, gId, inId
        character(len=strLen) :: h5File

        !Set filename to baseStr
        !FIXME: Correct to add .h5 to baseStr
        h5File = baseStr
        inquire(file=h5File,exist=fExist)
        if (.not. fExist) then
            write(*,*) 'Unable to open file: ', trim(h5File)
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
                write(*,*) 'Unable to find file/group = ', trim(h5File),'/',trim(gStrO)
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
                    if(IOVars(n)%useHyperslab) then
                        write(*,*) 'Unable to read attribute "',trim(IOVars(n)%idStr),'" as a hyperslab'
                        stop
                    else
                        call ReadHDFAtt(IOVars(n),inId)
                    endif
                endif

                if (dsExist) then
                    !Dataset
                    if(IOVars(n)%useHyperslab) then
                        call ReadHDFVarHyper(IOVars(n),inId)
                    else
                        call ReadHDFVar(IOVars(n),inId)
                    endif
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
    !-------------------------------------------
    !Routines to read/write individual variables or attributes

    !FIXME: Assuming double precision for input binary data
    !> Read HDF dataset from group
    subroutine ReadHDFVar(IOVar,gId)
        type(IOVAR_T), intent(inout) :: IOVar
        integer(HID_T), intent(in) :: gId

        integer(HSIZE_T), allocatable, dimension(:) :: dims
        integer :: Nr,herr, N
        integer :: typeClass
        integer(SIZE_T) :: typeSize
        character(len=strLen) :: inStr
        logical :: aExists

        !Start by getting rank, dimensions and total size
        call h5ltget_dataset_ndims_f(gId,trim(IOVar%idStr),Nr,herr)
        allocate(dims(Nr))
        call h5ltget_dataset_info_f(gId,trim(IOVar%idStr), dims, typeClass, typeSize, herr)
        N = product(dims)

        ! TODO: Check if chunked/compressed, adjust chunk cach based on ranks (size of dataset for each rank)
        !h5pset_chunk_cache_f()

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

        !Check for attribute strings
        !Unit
        call h5aexists_by_name_f(gID,trim(IOVar%idStr),"Units",aExists,herr)
        if (aExists) then
            call h5ltget_attribute_string_f(gID,trim(IOVar%idStr),"Units",inStr,herr)
        else
            inStr = "NULL"
        endif
        IOVar%unitStr = inStr

        !Description
        call h5aexists_by_name_f(gID,trim(IOVar%idStr),"Description",aExists,herr)
        if (aExists) then
            call h5ltget_attribute_string_f(gID,trim(IOVar%idStr),"Description",inStr,herr)
        else
            inStr = "NULL"
        endif
        IOVar%descStr = inStr
        !write(*,"(A,A,A,I)") 'Read dataset ', trim(IOVar%idStr), ' from ', gId

        IOVar%isDone = .true.
    end subroutine ReadHDFVar

    !> Read a dataset specified as a hyperslab. This assumes a stride of 1 in all dimensions
    subroutine ReadHDFVarHyper(IOVar,gId)
        type(IOVAR_T), intent(inout) :: IOVar
        integer(HID_T), intent(in) :: gId

        integer(HSIZE_T), allocatable, dimension(:) :: dims
        integer(HSIZE_T) :: N
        integer :: Nr,herr
        integer :: typeClass
        integer(SIZE_T) :: typeSize
        integer(HID_T) :: dsetId,dataspace,memspace

        !Start by getting rank and dimensions
        call h5ltget_dataset_ndims_f(gId,trim(IOVar%idStr),Nr,herr)
        allocate(dims(Nr))
        call h5ltget_dataset_info_f(gId,trim(IOVar%idStr), dims, typeClass, typeSize, herr)

        !Ensure the requested hyperslab is valid within the variable's dimensions

        if (allocated(IOVar%data)) then
            deallocate(IOVar%data)
        endif
        allocate(IOVar%data(IOVar%N))
        N = IOVar%N !Convert to HSIZE_T

        !Read based on data type
        select case(IOVar%vType)
        case(IONULL,IOREAL)
            call h5dopen_f(gId,trim(IOVar%idStr),dsetId,herr)
            call h5dget_space_f(dsetId,dataspace,herr)
            call h5sselect_hyperslab_f(dataspace,H5S_SELECT_SET_F,IOVar%offsets,IOVar%dims,herr)
            call h5screate_simple_f(1,(/N/),memspace,herr)
            call h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,(/integer(HSIZE_T)::0/),(/N/),herr)
            call h5dread_f(dsetId,H5T_NATIVE_DOUBLE,IOVar%data,(/N/),herr,memspace,dataspace)
            call h5sclose_f(dataspace,herr)
            call h5sclose_f(memspace,herr)
            call h5dclose_f(dsetId,herr)
        case default
            write(*,*) 'Unknown HDF data type, bailing ...'
            stop
        end select
        IOVar%isDone = .true.
    end subroutine ReadHDFVarHyper

    !FIXME: Add scaling to attributes
    !> Read an HDF attribute from a group
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
        write(*,"(A,A,F,A,I)") 'Read attribute ', trim(IOVar%idStr), IOVar%data(1), " from ", gId
        IOVar%isDone = .true.
    end subroutine ReadHDFAtt

    !> Helper function to get integer step from
    !> string name of Step#XXX format
    function GetStepInt(stepStr) result(nStep)
        character(len=*), intent(in) :: stepStr
        integer :: nStep, status
        read(stepStr(6:),*,iostat=status) nStep
        if(status /= 0) then
            write(*,"(A,I,A,A)"), "Conversion of step ", stepStr, " failed.", & 
            " Update to timeAttributeCache dataset size failed."
            stop
        endif
    end function 

    !> Helper function to check for resizing 
    !> of the timeAttributeCache datasets' size
    !> 
    subroutine CheckAttCacheSize(stepStr, cacheExist, cacheCreated)
        character(len=*), intent(in) :: stepStr
        logical, intent(in) :: cacheExist
        logical, intent(in) :: cacheCreated
        integer :: nStep

        nStep = GetStepInt(stepStr)

        if(cacheExist .and. .not. cacheCreated) then
            !write(*,"(A,1x,I)") "timeAttributeCache exists, and not just created, set cacheSize to ", nstep + 1
            cacheSize = nStep + 1
        else 
            !write(*,"(A,1x,I)") "timeAttributeCache exists, and just created, set cachesize to ", cacheSize + 1
            cacheSize = cacheSize + 1
        end if
        
    end subroutine

    !> Writes an Attribute that is a float or integer to
    !> timeAttributeCache group for each attribute written
    !> to a Step# group
    subroutine WriteCacheAtt(IOVar,gId)
        !> IOVar to write
        type(IOVAR_T), intent(inout) :: IOVar
        !> Group ID of attribute cache
        integer(HID_T), intent(in) :: gId

        integer(HID_T) :: sId, dId, pId, memId
        integer :: herr
        integer :: Nr = 1, memRank = 1
        logical :: dSetExists = .False.
        real(rp) :: X
        integer(HSIZE_T) :: dSize = 1, rank = 1
        integer(HSIZE_T), dimension(1) :: cdims(1), maxdims(1), dataDim(1), memDim(1)
        integer(HSIZE_T), dimension(1,1) :: coord
        dataDim = (/1/)
        memDim = (/1/)
        maxdims(1) = H5S_UNLIMITED_F

        call h5ltpath_valid_f(gId, trim(IOVar%idStr), .True., dSetExists, herr)
        ! Create dataspace and dataset initially in group
        if (.not. dSetExists) then
            !write(*,"(A,1x,A)") "Create var ", trim(IOVar%idStr)
            cdims(1) = 1
            call h5screate_simple_f(Nr, cdims, sId, herr, maxdims=maxdims)
            call h5pcreate_f(H5P_DATASET_CREATE_F, pId, herr)
            cdims(1) = CACHE_CHUNK_SIZE
            call h5pset_chunk_f(pId, Nr, cdims, herr)
            select case(IOVar%vType)
                case(IONULL,IOREAL)
                    call h5dcreate_f(gId, trim(IOVar%idStr), H5T_NATIVE_DOUBLE, sId, &
                                    dId, herr, dcpl_id=pId)
                case(IOINT)
                    call h5dcreate_f(gId, trim(IOVar%idStr), H5T_NATIVE_INTEGER, sId, &
                                    dId, herr, dcpl_id=pId)
            end select

            coord(1,1) = 1
            call h5sselect_elements_f(sId, H5S_SELECT_SET_F, Nr, 1, coord, herr)
            call h5pclose_f(pId,herr)
        else
            !write(*,"(A,1x,A)") "Found var ", trim(IOVar%idStr)
            ! Open dataset on subsequent time steps
            call h5dopen_f(gId, trim(IOVar%idStr), dId, herr)
            ! Get the proper dataspae to select the element position 
            ! of the next attribute element to add and write to the dataset
            call h5dget_space_f(dId, sId, herr)
            ! Resize dataset
            dSize = cacheSize
            call h5dset_extent_f(dId, (/dSize/), herr)
            ! Create memory space for (extended) dataset addition
            cdims(1) = 1
            call h5screate_simple_f(Nr, cdims, memId, herr)
            ! Select the single element at coord = offset in the file space
            coord(1,1) = cacheSize
            call h5sselect_elements_f(sId, H5S_SELECT_SET_F, Nr, 1, coord, herr)
        end if

        select case(IOVar%vType)
            case(IONULL,IOREAL)
                X = IOVar%data(1)
                if (dSetExists) then
                    call h5dwrite_f(dId, H5T_NATIVE_DOUBLE, X, dataDim, herr, &
                                    mem_space_id=memId, file_space_id=sId)
                    call h5sclose_f(memId,herr)
                else
                    call h5dwrite_f(dId, H5T_NATIVE_DOUBLE, X, dataDim, herr)
                end if
            case(IOINT)
                X = IOVar%data(1)
                if (dSetExists) then
                    call h5dwrite_f(dId, H5T_NATIVE_INTEGER, int(X), dataDim, herr, &
                                    mem_space_id=memId, file_space_id=sId)
                    call h5sclose_f(memId,herr)
                else
                    call h5dwrite_f(dId, H5T_NATIVE_INTEGER, int(X), dataDim, herr)
                end if
        end select
        
        ! Cleanup 
        call h5dclose_f(dId,herr)
        call h5sclose_f(sId,herr)
    end subroutine WriteCacheAtt

    !FIXME: Add scaling to attributes
    !> Write a variable as an attribue for the
    !> specified HDF group
    subroutine WriteHDFAtt(IOVar,gId)
        !> Variable to write to group
        type(IOVAR_T), intent(inout) :: IOVar
        !> HDF ID for Group
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
    !> Write a variable to an HDF dataset and add
    !> attributes to the dataset
    subroutine WriteHDFVar(IOVar,gId,doIOPO,doCompress)
        !> IO Var to write
        type(IOVAR_T), intent(inout)    :: IOVar
        !> Group ID
        integer(HID_T), intent(in)      :: gId
        !> Flag to do IO Precision
        logical, intent(in), optional   :: doIOPO
        !> Flag to compress
        logical, intent(in), optional   :: doCompress

        logical :: doIOP !Do IO precision for reals
        integer :: Nr
        integer :: herr
        integer(HSIZE_T) :: h5dims(MAXIODIM)
        integer(HID_T) :: dId, sId, pId
        real(rp) :: vScl
        integer(HSIZE_T), dimension(:), allocatable :: cdims
        !real(dp), dimension(:), allocatable, target:: data
#ifdef __ENABLE_COMPRESS
        logical :: avail
        !type(C_PTR) :: f_ptr
        integer :: status
        ! Used for ZFP HDF5 plugin
        integer(C_INT), dimension(:), allocatable :: cd_values
        integer(C_SIZE_T) :: cd_nelmts = 1
        integer :: szip_options_mask
        integer :: szip_pixels_per_block
#endif
        Nr = IOVar%Nr
        allocate(cdims(Nr))
        !allocate(data(size(IOVar%data)))
        !data = 0.0_dp
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
            if(.not. doCompress) then  
                !Assume real by default
                if (doIOP) then
                    call h5ltmake_dataset_float_f(gId,trim(IOVar%idStr), Nr, h5dims(1:Nr), & 
                            real(vScl*IOVar%data,sp),herr)
                else
                    call h5ltmake_dataset_double_f(gId,trim(IOVar%idStr), Nr, h5dims(1:Nr), &
                            real(vScl*IOVar%data,dp),herr)
                endif
#ifdef __ENABLE_COMPRESS
            else
                call h5screate_simple_f(Nr, h5dims(1:Nr), sid, herr)
                call h5pcreate_f(H5P_DATASET_CREATE_F, pId, herr)
                if (doIOP) then
                    call chunk_size(sp, h5dims(1:Nr), cdims)
                else
                    call chunk_size(dp, h5dims(1:Nr), cdims)
                endif
                call h5pset_chunk_f(pId, Nr, cdims, herr)

                if(Z_ALG == ZLIB) then
                    call h5pset_shuffle(pId, herr)
                    call h5pset_deflate_f(pId, 6, herr)
                elseif(Z_ALG == SZIP) then
                    szip_options_mask = H5_SZIP_NN_OM_F
                    if (doIOP) then
                        szip_pixels_per_block = 16
                    else
                        szip_pixels_per_block = 8
                    endif
                    call H5pset_szip_f(pId, szip_options_mask, szip_pixels_per_block, herr)
                elseif(Z_ALG == ZSTD) then
                    call H5zfilter_avail_f(COMPRESS_ZSTD, avail, status)
                    if(avail) then
                        allocate(cd_values(1))
                        cd_nelmts = 1
                        cd_values(1) = 20
                        call h5pset_filter_f(pId, COMPRESS_ZSTD, H5Z_FLAG_MANDATORY, &
                        cd_nelmts, cd_values, herr)
                    else
                        write(*,*) 'ZSTD filter not initailized, please ensure the ZSTD HDF5 plugin is loaded. \n'
                        write(*,*) 'You may also use the default compression SZIP without needing plugins.'
                        stop
                    endif
#ifdef __ENABLE_ZFP
                elseif(Z_ALG == ZFP) then
                    ! Only necessary for using H5Z_zfp properties calls
                    ! status = H5Z_zfp_initialize()
                    
                    call H5zfilter_avail_f(COMPRESS_ZFP, avail, status)

                    if (avail) then
                        ! Setup ZFP
                        !------------
                        ! Use Plug-in
                        !------------
                        allocate(cd_values(1:H5Z_ZFP_CD_NELMTS_MEM))
                        cd_values = 0
                        cd_nelmts = H5Z_ZFP_CD_NELMTS_MEM
                        
                        if (zfpmode .EQ. H5Z_ZFP_MODE_RATE)  then
                            call H5pset_zfp_rate_cdata(rate, cd_nelmts, cd_values)
                            if(cd_values(1).NE.1 .OR. cd_nelmts.NE.4) then
                                print*,'H5Pset_zfp_rate_cdata failed'
                                stop 1
                            endif
                        else if (zfpmode .EQ. H5Z_ZFP_MODE_PRECISION)  then
                            call H5pset_zfp_precision_cdata(prec, cd_nelmts, cd_values)
                            if(cd_values(1).NE.2 .OR. cd_nelmts.NE.3) then
                                print*,'H5Pset_zfp_precision_cdata failed'
                                stop 1
                            endif
                        else if (zfpmode .EQ. H5Z_ZFP_MODE_ACCURACY) then
                            call H5pset_zfp_accuracy_cdata(0._dp, cd_nelmts, cd_values)
                            if(cd_values(1).NE.3 .OR. cd_nelmts.NE.4) then
                                print*,'H5Pset_zfp_accuracy_cdata failed'
                                stop 1
                            endif
                        else if (zfpmode .EQ. H5Z_ZFP_MODE_EXPERT)  then
                            call H5pset_zfp_expert_cdata(minbits, maxbits, maxprec, minexp, cd_nelmts, cd_values)
                            if(cd_values(1).NE.4 .OR. cd_nelmts.NE.6) then
                                print*,'H5Pset_zfp_expert_cdata failed'
                                stop 1
                            endif
                        else if (zfpmode .EQ. H5Z_ZFP_MODE_REVERSIBLE)  then
                            call H5pset_zfp_reversible_cdata(cd_nelmts, cd_values)
                            if(cd_values(1).NE.5 .OR. cd_nelmts.NE.1) then
                                print*,'H5Pset_zfp_reversible_cdata failed'
                                stop 1
                            endif
                        endif
                
                        call h5pset_filter_f(pId, COMPRESS_ZFP, H5Z_FLAG_MANDATORY, &
                                cd_nelmts, cd_values, herr)

                        !---------------
                        ! Use Properties (this sets the filter)
                        !---------------
                        ! if (zfpmode == H5Z_ZFP_MODE_RATE)  then
                        !     status = H5Pset_zfp_rate(pId, rate)
                        !     !call check("H5Pset_zfp_rate", status, nerr)
                        ! else if (zfpmode == H5Z_ZFP_MODE_PRECISION)  then
                        !     status = H5Pset_zfp_precision(pId, prec)
                        !     !call check("H5Pset_zfp_precision", status, nerr)
                        ! else if (zfpmode == H5Z_ZFP_MODE_ACCURACY) then
                        !     status = H5Pset_zfp_accuracy(pId, acc)
                        !     !call check("H5Pset_zfp_accuracy", status, nerr)
                        ! else if (zfpmode == H5Z_ZFP_MODE_EXPERT)  then
                        !     status = H5Pset_zfp_expert(pId, minbits, maxbits, maxprec, minexp)
                        !     !call check("H5Pset_zfp_expert", status, nerr)
                        ! else if (zfpmode == H5Z_ZFP_MODE_REVERSIBLE)  then
                        !     status = H5Pset_zfp_reversible(pId)
                        !     !call check("H5Pset_zfp_reversible", status, nerr)
                        ! endif
                        ! status = H5Z_zfp_finalize()
                    else
                        write(*,*) 'ZFP filter not initailized, please ensure the ZFP HDF5 plugin is loaded.'
                        write(*,*) 'You may also use the default compression SZIP without needing plugins.'
                        stop
                    endif
#endif              
                endif

                if (doIOP) then
                    call h5dcreate_f(gId, trim(IOVar%idStr), H5T_NATIVE_REAL, sId, &
                    dId, herr, dcpl_id=pId)
                    !data = vScl*IOVar%data
                    !f_ptr = C_LOC(data(1))
                    call h5dwrite_f(dId, H5T_NATIVE_REAL, real(vScl*IOVar%data,sp), h5dims(1:Nr), herr)
                else
                    call h5dcreate_f(gId, trim(IOVar%idStr), H5T_NATIVE_DOUBLE, sId, &
                    dId, herr, dcpl_id=pId)
                    !data = vScl*IOVar%data
                    !f_ptr = C_LOC(data(1))
                    call h5dwrite_f(dId, H5T_NATIVE_DOUBLE, real(vScl*IOVar%data,dp), h5dims(1:Nr), herr)
                endif

                call h5pclose_f(pId, herr)
                call h5dclose_f(dId, herr)
                call h5sclose_f(sId, herr)
#endif           
            endif
        case(IOINT)
            call h5ltmake_dataset_int_f(gId,trim(IOVar%idStr), Nr, h5dims(1:Nr), &
                    int(vScl*IOVar%data),herr)
        case default
            write(*,*) 'Unknown HDF data type, bailing ...'
            stop
        end select

        !Set attached values
        call h5ltset_attribute_string_f(gId,trim(IOVar%idStr),"Units"      ,trim(IOVar%unitStr),herr)
        call h5ltset_attribute_string_f(gId,trim(IOVar%idStr),"Description",trim(IOVar%descStr),herr)

        IOVar%isDone = .true.
    end subroutine WriteHDFVar

    !> Write array of IOVar to file fOut (from baseStr), under (optional) group gStrO
    !> @note If gStrO unspecified, written to root of HDF5
    !> doIOp is whether to use IOP (ie output slice) or rp (ie restart)
    subroutine WriteVars(IOVars,doIOp,baseStr,gStrO,gStrOO,doStampCheckO)
        !> Array of IOVars
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        !> Do IO Precision writes
        logical, intent(in) :: doIOp
        !> Base String
        character(len=*), intent(in) :: baseStr
        !> Group Name
        character(len=*), intent(in), optional :: gStrO
        !> Subgroup Name
        character(len=*), intent(in), optional :: gStrOO
        !> Check if output file has been stamped
        logical         , intent(in), optional :: doStampCheckO

        logical :: fExist, gExist, doStampCheck
        logical :: writeCache, cacheExist, cacheCreate
        logical :: doCompress
        integer :: herr
        integer(HID_T) :: h5fId, gId, ggId, outId, cacheId
        character(len=strLen) :: h5File
        character(len=strLen), parameter :: attrGrpName = "timeAttributeCache"
        type(IOVAR_T) :: stepVar
        !Set filename to baseStr
        !FIXME: Correct to add .h5 to baseStr
        h5File = baseStr
        writeCache = .false.
        cacheCreate = .false.
        doCompress = .false.

        if (present(doStampCheckO)) then
            doStampCheck = doStampCheckO
        else
            doStampCheck = .true.
        endif

        !If we're writing to root of this file, then stamp (will ignore if already stamped)
        if (.not. present(gStrO) .and. doStampCheck) then
            call StampIO(h5File)
        endif      

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
                    write(*,*) 'Overwriting group ', trim(h5File), '/', trim(gStrO)
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

            if(trim(toUpper(gStrO(1:5))) == "STEP#") then
                writeCache = .true.
                !Check if cache group exists
                call h5lexists_f(h5fId,trim(attrGrpName),cacheExist,herr)
                if (.not. cacheExist) then
                    !write(*,*) "Cache group does not exist, creating it."
                    !Create cache group
                    call h5gcreate_f(h5fId,trim(attrGrpName),cacheId,herr)     
                    cacheCreate = .true.    
                endif 
                ! Open attribute cache group
                call h5gopen_f(h5fId,trim(attrGrpName),cacheId,herr)  

                ! Check attribute cache size and resize
                call CheckAttCacheSize(trim(gStrO), cacheExist, cacheCreate)
                ! Write Step# to cache
                stepVar%Nr = 0
                stepVar%idStr = "step"
                stepVar%vType = IOINT
                stepVar%data = [GetStepInt(trim(gStrO))]
                call WriteCacheAtt(stepVar,cacheId)
            endif 

            if (present(gStrOO)) then
                !Create subgroup
                call h5gcreate_f(gId,trim(gStrOO),ggId,herr)
                outId = ggId
            else
                outId = gId 
            endif
            doCompress = ENABLE_COMPRESS
        else
            !Write to root
            outId = h5fId
            doCompress = .False.
        endif !gStrO

        !Do writing
        if(present(doStampCheckO)) then
            call WriteVars2ID(IOVars,outId,doIOp, &
                                doCompress=doCompress,isRoot=.true.)
        else
            ! Write step vars and attributes to group outId
            if(writeCache) then  
                call WriteVars2ID(IOVars,outId,doIOp, &
                    cacheId=cacheId,doCompress=doCompress)
            else
                call WriteVars2ID(IOVars,outId,doIOp,doCompress=doCompress)
            end if
        end if

        !Now done, close up shop
        if (present(gStrOO)) call h5gclose_f(ggId,herr)
        if (present(gStrO )) call h5gclose_f( gId,herr)
        if (writeCache) call h5gclose_f(cacheId,herr) !Close cache group

        call h5fclose_f(h5fId,herr) !Close file
        call h5close_f(herr) !Close intereface

    end subroutine WriteVars

    !> Write out all IOVars to their respective IDs
    subroutine WriteVars2ID(IOVars,outId,doIOp,cacheId,doCompress,isRoot)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        logical, intent(in) :: doIOp
        integer(HID_T), intent(in) :: outId
        integer(HID_T), intent(in), optional :: cacheId
        logical, intent(in), optional :: doCompress
        logical, intent(in), optional :: isRoot
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
                    if(IOVars(n)%useHyperslab) then
                        write(*,*) 'Unable to write attribute "',trim(IOVars(n)%idStr),'" as a hyperslab'
                        stop
                    else
                        call WriteHDFAtt(IOVars(n),outId)
                        if(present(cacheId)) then
                            call WriteCacheAtt(IOVars(n),cacheId)
                        endif
                    endif
                else
                !N-rank array
                    !Create data space, use rank/dim info from IOVar
                    if(IOVars(n)%useHyperslab) then
                        write(*,*) 'Writing dataset "',trim(IOVars(n)%idStr),'" as a hyperslab not yet supported'
                        stop
                    else
                        call WriteHDFVar(IOVars(n),outId,doIOP,doCompress)
                    endif
                endif !Nr=0
            endif !isSet
        enddo

    end subroutine WriteVars2ID

    !-------------------------------------------
    !These routines add data for input to IO chain

    !> Add Input Variable to IOVars pool to read
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

    !> Alternate subroutine to read in a hyperslab from a dataset
    subroutine AddInVarHyper(IOVars,idStr,offsets,counts,vTypeO,vSclO)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars
        integer, dimension(:), intent(in) :: offsets
        integer, dimension(:), intent(in) :: counts
        character(len=*), intent(in) :: idStr
        integer, intent(in), optional :: vTypeO
        real(rp), intent(in), optional :: vSclO
        integer :: n,nr

        nr = SIZE(offsets)

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

        IOVars(n)%Nr = nr
        IOVars(n)%offsets(1:nr) = offsets
        IOVars(n)%dims(1:nr) = counts
        IOVars%N = product(counts)
        IOVars%useHyperslab = .true.

    end subroutine AddInVarHyper

    !> Clears info/memory from an IO chain
    subroutine ClearIO(IOVars)
        type(IOVAR_T), dimension(:), intent(inout) :: IOVars

        integer :: i,Nv

        Nv = size(IOVars)
        IOVars(:)%idStr   = "NONE"
        IOVars(:)%unitStr = "CODE"
        IOVars(:)%descStr = "DESCRIPTION"
        IOVars(:)%dStr    = "UNSET"
        IOVars(:)%Nr      = 0
        IOVars(:)%N       = 0
        IOVars(:)%vType   = IONULL
        IOVars(:)%scale   = 1.0
        IOVars(:)%renorm  = 0.0
        IOVars(:)%toWrite = .false.
        IOVars(:)%toRead  = .false.
        IOVars(:)%isDone  = .false.
        IOVars(:)%useHyperslab = .false.
        
        do i=1,Nv
            IOVars(i)%dims(:)    = 0
            IOVars(i)%offsets(:) = 0

            if (allocated(IOVars(i)%data)) then
                deallocate(IOVars(i)%data)
            endif
        enddo
    end subroutine ClearIO

    !-----------------------------
    !HDF 5 helper routines
    !> Converts Fortran real kind to HDF5 precision
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
    
    !> Write integer scalar to attribute
    subroutine writeInt2HDF(gId,vId,datIn)
        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        integer, intent(in) :: datIn
        
        integer(HID_T) :: dspcId, atId
        integer(HSIZE_T), dimension(1) :: dims
        integer :: herror

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

    !> Writes a single scalar as attribute to a group/root
    !> Always uses RP precision
    subroutine writeReal2HDF(gId,vId,datIn)

        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        real(rp), intent(in) :: datIn
    
        integer(HID_T) :: dspcId, atId
        integer(HID_T):: h5gReal
        integer(HSIZE_T), dimension(1) :: dims
        integer :: herror

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

    !> Write a string to HDF group 
    !> as an attribute
    subroutine writeString2HDF(gId,vId,data)
        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        character(len=*), intent(in) :: data        
        integer(HID_T) :: dspcId, atId, strId
        integer(HSIZE_T), dimension(1) :: dims
        integer :: herror

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

    !> Read integer from HDF5 attribute
    function readIntHDF(gId,vId,vDefOpt) result(vOut)
        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        integer, intent(in), optional :: vDefOpt

        integer :: vOut
        integer :: vDef
        integer(HSIZE_T), dimension(1) :: dims
        integer(HID_T) :: attrId
        integer :: herror

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
    
    !> Read real (rp) from HDF5 attribute
    function readRealHDF(gId,vId,vDefOpt) result(vOut)
        integer(HID_T), intent(in) :: gId
        character(len=*),intent(in) :: vId
        real(rp), intent(in), optional :: vDefOpt

        real(rp) :: vOut,vDef
        integer(HSIZE_T), dimension(1) :: dims
        integer(HID_T) :: attrId, h5gReal
        integer :: herror

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
