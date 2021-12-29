module files
    use kdefs
#ifdef __INTEL_COMPILER
    use ifport
#endif
    use, intrinsic :: iso_c_binding

    implicit none

    !Whether to do sym links, can be disabled due to issue w/ virtual mem blowing up
    logical, private :: doSymLinks = .true. 

    interface
        integer(kind=c_int) function symlink_c(filename, linkname) bind(C, name='symlink')
            import c_int, c_char
            character(kind=c_char), intent(in) :: filename(*)
            character(kind=c_char), intent(in) :: linkname(*)
        end function symlink_c
    end interface

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

        logical :: isDir,mdirStat
        integer :: eStat,cStat
        character(len=strLen) :: cmd,cMsg

        isDir = CheckDir(fIn)
        if (.not. isDir) then
            write(*,'(5a)') ANSIRED,'<',trim(fIn),' does not exist, creating ...>',ANSIRESET
#ifdef __INTEL_COMPILER
            mdirStat = MAKEDIRQQ(trim(fIn))
            if(.not. mdirStat) then
                write (*,*) 'Failed to create directory "',trim(fIn),'"'
                stop
            endif
#else
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
#endif
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
#ifdef __INTEL_COMPILER
            eStat = unlink(fStr)
            if (eStat) then
                write(*,*) 'Unlink error: ', eStat
                stop
            endif
#else
            !System call method
            call EXECUTE_COMMAND_LINE(trim(cmd), wait=.true., exitstat=eStat, cmdstat=cStat, cmdmsg=cMsg)
            if(cStat > 0) then
                write (*,*) 'Command line "',trim(cmd),'" failed with error: ',trim(cMsg)
                call printProcessInfo()
                stop
            elseif(cStat < 0) then
                write (*,*) 'EXECUTE_COMMAND_LINE not supported'
                stop
            endif
#endif

        endif
    end subroutine CheckAndKill

    subroutine DisableSymLinks()
        doSymLinks = .false.
    end subroutine
    
    !Create symlink connecting fReal and fSym
    subroutine MapSymLink(fReal,fSym)
        character(len=*), intent(in) :: fReal,fSym

        integer :: cStat

        if (.not. doSymLinks) return

        !Kill sym link if it exists
        call CheckAndKill(fSym,.false.)

        !Create sym link
        cStat = symlink_c(trim(fReal)//c_null_char, trim(fSym)//c_null_char)
        if(cStat < 0) then
            write (*,*) 'LIBC symlink() of "',trim(fSym),'" to "',trim(fReal),'" failed with error: ',cStat
            stop
        endif

    end subroutine MapSymLink
    
end module files
