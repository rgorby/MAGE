module files
    use kdefs

    implicit none

    contains

    function CheckFile(fIn)
        character(len=*), intent(in) :: fIn
        logical :: CheckFile

        inquire(file=fIn,exist=CheckFile)
    end function CheckFile

    function CheckDir(fIn)
        character(len=*), intent(in) :: fIn
        logical :: CheckDir
#IFDEF __INTEL_COMPILER
        inquire(directory=fIn,exist=CheckDir) !Intel only
#ELSE
        !Might work for gfortran
        inquire(file=trim(fIn)//'/.',exist=CheckDir)
#ENDIF
    end function CheckDir

    subroutine CheckDirOrMake(fIn)
        character(len=*), intent(in) :: fIn
        logical :: isDir

        isDir = CheckDir(fIn)
        if (.not. isDir) then
            write(*,'(5a)') ANSIRED,'<',trim(fIn),' does not exist, creating ...>',ANSIRESET
            call EXECUTE_COMMAND_LINE('mkdir ' // trim(fIn) , wait=.true.)
        endif
        
    end subroutine CheckDirOrMake

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
            write(*,'(5a)') ANSIRED,'<',trim(fStr),' already exists, deleting ...>',ANSIRESET
            call EXECUTE_COMMAND_LINE('rm ' // trim(fStr) , wait=.true.)
        endif
    end subroutine CheckAndKill

end module files
