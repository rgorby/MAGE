module strings
    use kdefs
    use, intrinsic :: iso_fortran_env
    use dates, ONLY : DateTimeStr
    implicit none

    contains

    !Get filename for XML input deck
    !Assuming first command line argument
    subroutine getIDeckStr(IDeckStr)
        character(len=strLen), intent(out) :: IDeckStr

        integer :: Narg
        logical :: fExist
        
        !Find input deck
        Narg = command_argument_count()
        if (Narg .eq. 0) then
            write(*,*) 'No input deck specified, defaulting to Input.xml'
            IDeckStr = "Input.xml"
        else
            call get_command_argument(1,IDeckStr)
        endif

        write(*,*) 'Reading input deck from ', trim(IDeckStr)
        inquire(file=IDeckStr,exist=fExist)
        if (.not. fExist) then
            write(*,*) 'Error opening input deck, exiting ...'
            write(*,*) ''
            stop
        endif

    end subroutine getIDeckStr

    ! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html 
    ! Original author: Clive Page
    function toUpper(strIn) result(strOut)
        character(len=*), intent(in) :: strIn
        character(len=len(strIn)) :: strOut

        integer :: i,j

        do i=1,len(strIn)
            j = iachar(strIn(i:i))
            if ( j>=ichar("a") .and. j<=iachar("z") ) then
                strOut(i:i) = achar(iachar(strIn(i:i))-32)
            else
                strOut(i:i) = strIn(i:i)
            endif
        enddo
    end function toUpper
    
    !Print out basic configuration info
    subroutine printConfigStamp()
        character(len=strLen) :: dtStr,gStr,bStr

        call DateTimeStr(dtStr)
        call GitHash(gStr)
        call GitBranch(bStr)
        write(*,*) ANSIGREEN
        write(*,*) '---------------'
        write(*,*) 'Kaiju configuration'
        write(*,'(2a)') 'Git branch = ', trim(bStr)
        write(*,'(2a)') 'Git hash   = ', trim(gStr)
#ifdef __INTEL_COMPILER_OLD
        write(*,'(2a)') 'Compiler   = ', trim(gStr)
#else
        write(*,'(2a)') 'Compiler   = ', compiler_version()
#endif
        write(*,'(2a)') 'Run starting on: ', trim(dtStr)
        
        !write(*,'(2a)') 'Compiler flags = ', compiler_options()
        write(*,*) '---------------'
        write(*,'(a)',advance="no") ANSIRESET!, ''

    end subroutine printConfigStamp
    
    !Create string with git hash if possible
    subroutine GitHash(gStr)

        character(len=*), intent(inout) :: gStr
        character(len=strLen) :: cOpts
        integer :: n,nOff,nH

        nOff = 16
        nH = 7
#ifdef __INTEL_COMPILER_OLD
        gStr = "XXXXXXX" !Avoid unavailable compiler_options
#else
        cOpts = compiler_options()
        n = index(cOpts,"-DGITCOMMITHASH=")
        gStr = cOpts(n+nOff:n+nOff+nH)        
#endif
    end subroutine GitHash

    !Create string with git hash if possible
    subroutine GitBranch(gStr)

        character(len=*), intent(inout) :: gStr
        character(len=strLen) :: cOpts
        integer :: n,nOff

        nOff = 12
#ifdef __INTEL_COMPILER_OLD
        
        gStr = "XXXXXXX" !Avoid unavailable compiler_options
#else
        cOpts = compiler_options()
        n = index(cOpts,"-DGITBRANCH=")
        !Don't know length of branch name, so need to find next space
        gStr = cOpts(n+nOff:)
        n = index(gStr," ")
        gStr = gStr(1:n-1)
#endif
    end subroutine GitBranch


end module strings
