module files
    use kdefs

    implicit none

    !Whether to do sym links, can be disabled due to issue w/ virtual mem blowing up
    logical, private :: doSymLinks = .true. 

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
        integer :: eStat,cStat
        character(len=strLen) :: cmd,cMsg

        isDir = CheckDir(fIn)
        if (.not. isDir) then
            write(*,'(5a)') ANSIRED,'<',trim(fIn),' does not exist, creating ...>',ANSIRESET
            write (cmd,'(A,A)') 'mkdir ',trim(fIn)
            call EXECUTE_COMMAND_LINE(trim(cmd), wait=.true., exitstat=eStat, cmdstat=cStat, cmdmsg=cMsg)
            if(cStat > 0) then
                write (*,*) 'Command line "',trim(cmd),'" failed with error: ',trim(cMsg)
                call printProcessInfo()
                stop
            elseif(cStat < 0) then
                write (*,*) 'EXECUTE_COMMAND_LINE not supported'
                stop
            endif
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
        integer :: eStat,cStat
        character(len=strLen) :: cmd,cMsg

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
            write (cmd,'(A,A)') 'rm ',trim(fStr)
            call EXECUTE_COMMAND_LINE(trim(cmd), wait=.true., exitstat=eStat, cmdstat=cStat, cmdmsg=cMsg)
            if(cStat > 0) then
                write (*,*) 'Command line "',trim(cmd),'" failed with error: ',trim(cMsg)
                call printProcessInfo()
                stop
            elseif(cStat < 0) then
                write (*,*) 'EXECUTE_COMMAND_LINE not supported'
                stop
            endif
        endif
    end subroutine CheckAndKill

    subroutine DisableSymLinks()
        doSymLinks = .false.
    end subroutine
    
    !Create symlink connecting fReal and fSym
    subroutine MapSymLink(fReal,fSym)
        character(len=*), intent(in) :: fReal,fSym

        integer :: eStat,cStat
        character(len=strLen) :: cmd,cMsg

        if (.not. doSymLinks) return

        !Kill sym link if it exists
        call CheckAndKill(fSym,.false.)

        !Create sym link
        write (cmd,'(A,A,A,A)') 'ln -sf ',trim(fReal),' ',trim(fSym)
        call EXECUTE_COMMAND_LINE(trim(cmd), wait=.true., exitstat=eStat, cmdstat=cStat, cmdmsg=cMsg)
        if(cStat > 0) then
            write (*,*) 'Command line "',trim(cmd),'" failed with error: ',trim(cMsg)
            call printProcessInfo()
            stop
        elseif(cStat < 0) then
            write (*,*) 'EXECUTE_COMMAND_LINE not supported'
            stop
        endif

    end subroutine MapSymLink
    
end module files
