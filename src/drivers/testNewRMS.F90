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
    type(rmReader_T) :: rmReader

    integer :: i
    real(rp) :: dt, t

    inpXML = New_XML_Input(trim(xmlStr),"Kaiju/REMIX",.true.)

    call initRM(ftag, inpXML, rmReader)
    
    call dump(fOutname, rmReader%shGr)

    dt = 0.5 * (rmReader%rmTab%times(rmReader%rmTab%N) - rmReader%rmTab%times(rmReader%rmTab%N-1))
    write(*,*)"dt=",dt
    t = rmReader%rmTab%times(2)
    write(*,*)t
    !do while (t < dt*20)
    do while (t < rmReader%rmTab%times(rmReader%rmTab%N) + 4*dt)
        call updateRM(rmReader, t)
        call testShellInterp(rmReader)
        t = t + dt
    enddo


    contains

    subroutine testShellInterp(rmS)
        type(rmReader_T), intent(in) :: rmS

        type(ShellGridVar_T) :: tmpVarSame, tmpVarCC
        integer :: i,j

        integer, dimension(:), allocatable :: i0_arr, j0_arr

        call initShellVar(rmS%shGr, rmS%nsFac(NORTH)%loc, tmpVarSame)
        tmpVarSame%mask = rmS%nsFac(NORTH)%mask
        call initShellVar(rmS%shGr, SHCC, tmpVarCC)
        tmpVarCC%mask = .true.

        call InterpShellVar_TSC_SG(rmS%shGr, rmS%nsFac(NORTH), rmS%shGr, tmpVarSame) ! RM fac to same location
        call InterpShellVar_TSC_SG(rmS%shGr, rmS%nsFac(NORTH), rmS%shGr, tmpVarCC)   ! RM fac to cell centers
        call dump(fOutname, rmS%shGr, rmS%nsFac(NORTH), "nsFac")
        call dump(fOutname, rmS%shGr, tmpVarSame, "tmpVarSame")
        call dump(fOutname, rmS%shGr, tmpVarCC, "tmpVarCC")

        call InterpShellVar_TSC_SG(rmS%shGr, tmpVarCC, rmS%shGr, tmpVarSame) ! CCs back to corners
        call dump(fOutname, rmS%shGr, tmpVarSame, "tmpVarC2CC2C")

        ! Try out letting InterpShellVar_TSC_pnt calculate its own dx
        tmpVarSame%data = 0.0
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j)
        do j=rmS%shGr%jsg,rmS%shGr%jeg+1
            do i=rmS%shGr%isg,rmS%shGr%ieg+1
                if (.not. tmpVarSame%mask(i,j)) cycle
                call InterpShellVar_TSC_pnt( \
                        rmS%shGr, rmS%nsFac(NORTH),\
                        rmS%shGr%th(i), rmS%shGr%ph(j),\
                        tmpVarSame%data(i,j) )
            enddo
        enddo
        call wrapJ_SGV(rmS%shGr, tmpVarSame)
        call dump(fOutname, rmS%shGr, tmpVarSame, "tmpVar_noDCell") ! Should be equal to "tmpVarSame"

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
