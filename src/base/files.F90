module files
    use kdefs

    implicit none

    contains

    subroutine CheckFileOrDie(fIn,errStr)
        character(len=*), intent(in) :: fIn,errStr

        logical :: fExist
        inquire(file=fIn,exist=fExist)
        if (.not. fExist) then
            write(*,*) trim(errStr)
            write(*,*) ''
            stop
        endif
        
    end subroutine CheckFileOrDie
    
    !Delete file if it already exists
    subroutine CheckAndKill(fStr)
        character(len=*), intent(in) :: fStr

        logical :: fExist

        inquire(file=trim(fStr),exist=fExist)
        if (fExist) then
            write(*,'(3a)') '<',trim(fStr),' already exists, deleting ...>'
            call EXECUTE_COMMAND_LINE('rm ' // trim(fStr) , wait=.true.)
            !write(*,*) ''
        endif
    end subroutine CheckAndKill

end module files
