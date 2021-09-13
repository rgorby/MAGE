!Base types and definitions for HDF5 IO
module ioH5types
    use kdefs
    use hdf5, ONLY : HSIZE_T
    
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
    !idStr/unitStr: String name of variable, its units, and its short description
    !Nr is number of ranks, and dims is size in each dimension (unused are 0)
    !data holds the information as a 1D array
    !scale is multiplier prior to output
    !vType: holds info about variable type, only partially implemented

    type IOVAR_T
        character(len=strLen) :: idStr="NONE",unitStr="CODE",descStr="DESCRIPTION"
        integer(HSIZE_T) :: dims(MAXIODIM) = 0 !Dimension information
        integer(HSIZE_T) :: offsets(MAXIODIM) = 0 !Offset for reading/writing a hyperslab, optional
        logical :: useHyperslab=.false.
        integer :: Nr = 0 !Number of ranks
        integer :: N = 0 !Total number of elements
        real(rp), dimension(:), allocatable :: data !1D holder for data
        real(rp) :: scale=1.0, renorm=0.0
        logical :: toWrite=.false.,toRead=.false. !Read or write this variable
        logical :: isDone=.false. !Whether or not variable has been successfully read/written
        integer :: vType
        character(len=strLen) :: dStr !Optional string data
    end type IOVAR_T


end module ioH5types
