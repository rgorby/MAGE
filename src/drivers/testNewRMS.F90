! Simple driver to test remix reader into ShellGrid
program testNewRMS

    use kdefs
    use XML_Input

    use remixReader
    use shellInterp

    use ioh5

    implicit none

    character(len=strLen) :: xmlStr = 'testRM.xml'
    character(len=strLen) :: ftag = 'msphere'
    character(len=100) :: fOutname = "rms.h5"
    type(XML_Input_T) :: inpXML
    type(rmState_T) :: rmState

    integer :: i
    real(rp) :: dt, t

    inpXML = New_XML_Input(trim(xmlStr),"Kaiju/REMIX",.true.)

    call initRM(ftag, inpXML, rmState)
    
    call dump(fOutname, rmState%shGr)

    dt = 0.5 * (rmState%rmTab%times(rmState%rmTab%N) - rmState%rmTab%times(rmState%rmTab%N-1))
    write(*,*)"dt=",dt
    t = rmState%rmTab%times(2)
    write(*,*)t
    !do while (t < dt*20)
    do while (t < rmState%rmTab%times(rmState%rmTab%N) + 4*dt)
        call updateRM(rmState, t)
        call testShellInterp(rmState)
        t = t + dt
    enddo


    contains

    subroutine testShellInterp(rmS)
        type(rmState_T), intent(in) :: rmS

        type(ShellGridVar_T) :: tmpVarSame, tmpVarCC
        integer :: i,j

        integer, dimension(:), allocatable :: i0_arr, j0_arr

        call initShellVar(rmS%shGr, rmS%nsFac(NORTH)%loc, tmpVarSame)
        tmpVarSame%mask = rmS%nsFac(NORTH)%mask
        call initShellVar(rmS%shGr, SHCC, tmpVarCC)
        tmpVarCC%mask = .true.

        allocate(i0_arr(tmpVarSame%isv:tmpVarSame%iev))
        allocate(j0_arr(tmpVarSame%jsv:tmpVarSame%jev))

        call InterpShellVar_TSC_SG(rmS%shGr, rmS%nsFac(NORTH), rmS%shGr, tmpVarSame)
        call InterpShellVar_TSC_SG(rmS%shGr, rmS%nsFac(NORTH), rmS%shGr, tmpVarCC)

        !do j=tmpVarSame%jsv,tmpVarSame%jev
        !    do i=tmpVarSame%isv,tmpVarSame%iev
        !        write(*,*) i,j
        !        write(*,*) "  ",tmpVarSame%data(i,j), tmpVarSame%mask(i,j)
        !        if (j == tmpVarSame%jev+10) then
        !            write(*,*) "  ",rmS%nsFAC(NORTH)%data(i,tmpvarSame%isv+1),rmS%nsFAC(NORTH)%mask(i,tmpvarSame%isv+1)
        !            write(*,*) "  ",tmpVarSame%data(i,j)/rmS%nsFAC(NORTH)%data(i,tmpvarSame%isv+1)
        !        else
        !            write(*,*) "  ",rmS%nsFAC(NORTH)%data(i,j),rmS%nsFAC(NORTH)%mask(i,j)
        !            write(*,*) "  ",tmpVarSame%data(i,j)/rmS%nsFAC(NORTH)%data(i,j)
        !        endif
        !    enddo
        !enddo

        ! Do our own calculation of index mappings
        !do j=varOut%jsv,varOut%jev
        !    do i=varOut%isv,varOut%iev
        !        !if (.not. varOut%mask(i,j)) cycle
        !        varOut%data(i,j) = InterpShellVar_TSC_pnt( \
        !                            sgVar, sgSource,\
        !                            dTheta, dPhi,\
        !                            sgDest%th(i), sgDest%ph(j))
        !    enddo
        !enddo
        do j=tmpVarSame%jsv, tmpVarSame%jev
            call getShellJLoc(rmS%shGr, rmS%nsFac(NORTH)%loc, rmS%shGr%ph(j), j0_arr(j))
        enddo
        do i=tmpVarSame%isv, tmpVarSame%iev
            call getShellILoc(rmS%shGr, rmS%nsFac(NORTH)%loc, rmS%shGr%th(i), i0_arr(i))
        enddo

        call dump(fOutname, rmS%shGr, rmS%nsFac(NORTH), "nsFac")
        call dump(fOutname, rmS%shGr, tmpVarSame, "tmpVarSame")
        call dump(fOutname, rmS%shGr, tmpVarCC, "tmpVarCC")
        call dump1Dint(fOutname, i0_arr, "i0")
        call dump1Dint(fOutname, j0_arr, "j0")


        ! Now interp back to corners
        call InterpShellVar_TSC_SG(rmS%shGr, tmpVarCC, rmS%shGr, tmpVarSame)
        call dump(fOutname, rmS%shGr, tmpVarSame, "tmpVarC2CC2C")

        stop

    end subroutine testShellInterp


    subroutine dump(fname, shGr, varO, vNameO)
        character(len=100), intent(in) :: fname
        type(ShellGrid_T), intent(in) :: shGr
        type(ShellGridVar_T), optional, intent(in) :: varO
        character(len=*), optional, intent(in) :: vNameO

        character(len=100) :: tmp
        type(IOVAR_T), dimension(5) :: IOVars        
        ! If varO not present, we assume we are starting fresh
        ! Wipe anything there and write shellGrid info
        ! If varO present, assume we are writing it out

        if (.not. present(varO)) then
            call CheckAndKill(fname, .true.)

            call ClearIO(IOVars)
            call AddOutVar(IOVars,"sh_th",shGr%th)
            call AddOutVar(IOVars,"sh_ph",shGr%ph)
            call WriteVars(IOVars,.true.,fname)
            return
        endif
        if (.not. present(vNameO)) then
            write(*,*)"vnameO"
            stop
        endif

        ! If still here, varO and vNameO present
        call ClearIO(IOVars)
        write(tmp,'(A,A)')trim(vNameO), "_data"
        call AddOutVar(IOVars,tmp,varO%data)
        write(tmp,'(A,A)')trim(vNameO), "_mask"
        call AddOutVar(IOVars,tmp,varO%mask*1.0_rp)
        call WriteVars(IOVars,.true.,fname)

    end subroutine dump


    subroutine dump1Dint(fname, var, vName)
        character(len=100), intent(in) :: fname
        integer, dimension(:), intent(in) :: var
        character(len=*), intent(in) :: vName

        character(len=100) :: tmp
        type(IOVAR_T), dimension(5) :: IOVars        
        ! If varO not present, we assume we are starting fresh
        ! Wipe anything there and write shellGrid info
        ! If varO present, assume we are writing it out

        
        ! If still here, varO and vNameO present
        call ClearIO(IOVars)
        call AddOutVar(IOVars,vName, var*1.0_rp)
        call WriteVars(IOVars,.true.,fname)


    end subroutine dump1Dint

end program testNewRMS
