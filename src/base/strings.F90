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
        character(len=strLen) :: dtStr,gStr

        call DateTimeStr(dtStr)
        call GitHash(gStr)
        write(*,*) ANSIGREEN
        write(*,*) '---------------'
        write(*,*) 'Kaiju configuration'
        write(*,'(2a)') 'Git hash = ', trim(gStr)
        write(*,'(2a)') 'Run starting on: ', trim(dtStr)
        write(*,'(2a)') 'Compiler = ', compiler_version()
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
        cOpts = compiler_options()
        n = index(cOpts,"-DGITCOMMITHASH=")
        gStr = cOpts(n+nOff:n+nOff+nH)        

    end subroutine GitHash
end module strings
