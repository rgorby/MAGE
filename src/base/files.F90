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

#ifdef __INTEL_COMPILER
        inquire(directory=fIn,exist=CheckDir) !Intel only
#else
        !Might work for gfortran
        inquire(file=trim(fIn)//'/.',exist=CheckDir)
#endif
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
        inquire(file=trim(fIn),exist=fExist)
        if (.not. fExist) then
            write(*,*) trim(errStr)
            write(*,*) "File: ", trim(fIn)
            write(*,*) ''
            stop
        endif
        
    end subroutine CheckFileOrDie
    
    !Delete file if it already exists
    subroutine CheckAndKill(fStr,doVerbO)
        character(len=*), intent(in) :: fStr
        logical, intent(in), optional :: doVerbO

        logical :: fExist,doVerb
        if (present(doVerbO)) then
            doVerb = doVerbO
        else
            doVerb = .true.
        endif

        inquire(file=trim(fStr),exist=fExist)
        if (fExist) then
            if (doVerb) then
                write(*,'(5a)') ANSIRED,'<',trim(fStr),' already exists, deleting ...>',ANSIRESET
            endif
            call EXECUTE_COMMAND_LINE('rm ' // trim(fStr) , wait=.true.)
        endif
    end subroutine CheckAndKill

    !Create symlink connecting fReal and fSym
    subroutine MapSymLink(fReal,fSym)
        character(len=*), intent(in) :: fReal,fSym

        !Kill sym link if it exists
        call CheckAndKill(fSym,.false.)

        !Create sym link
        call EXECUTE_COMMAND_LINE('ln -sf '//trim(fReal)//' '//trim(fSym), wait=.true.)

    end subroutine MapSymLink
    
end module files
