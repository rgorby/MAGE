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

        type(ShellGridVar_T) :: tmpVar
        integer :: i,j

        call initShellVar(rmS%shGr, rmS%nsFac(NORTH)%loc, tmpVar)
        tmpVar%mask = rmS%nsFac(NORTH)%mask

        call InterpShellVar_TSC_SG(rmS%nsFac(NORTH), rmS%shGr, rmS%shGr, tmpVar)

        !do j=tmpVar%jsv,tmpVar%jev
        !    do i=tmpVar%isv,tmpVar%iev
        !        write(*,*) i,j
        !        write(*,*) "  ",tmpVar%data(i,j), tmpVar%mask(i,j)
        !        write(*,*) "  ",rmS%nsFAC(NORTH)%data(i,j),rmS%nsFAC(NORTH)%mask(i,j)
        !        write(*,*) "  ",tmpVar%data(i,j)/rmS%nsFAC(NORTH)%data(i,j)
        !    enddo
        !enddo

        call dump(fOutname, rmS%shGr, rmS%nsFac(NORTH), "nsFac")
        call dump(fOutname, rmS%shGr, tmpVar, "tmpVar")
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

end program testNewRMS
