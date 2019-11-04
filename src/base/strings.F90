module strings
    use kdefs, ONLY: strLen

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
    

end module strings
