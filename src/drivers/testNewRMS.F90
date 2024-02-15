program testNewRMS

    use kdefs
    use XML_Input

    use remixReader
    use shellInterp

    implicit none

    character(len=strLen) :: xmlStr = 'testRM.xml'
    character(len=strLen) :: ftag = 'msphere'
    type(XML_Input_T) :: inpXML
    type(rmState_T) :: rmState

    integer :: i
    real(rp) :: dt, t

    inpXML = New_XML_Input(trim(xmlStr),"Kaiju/REMIX",.true.)

    call initRM(ftag, inpXML, rmState)

    dt = 0.5 * (rmState%rmTab%times(rmState%rmTab%N) - rmState%rmTab%times(rmState%rmTab%N-1))
    write(*,*)"dt=",dt
    t = 0
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

        call initShellVar(rmS%shGr, rmS%nsFac(NORTH)%loc, tmpVar)

        call InterpShellVar_TSC_SG(rmS%nsFac(NORTH), rmS%shGr, rmS%shGr, tmpVar)

    end subroutine testShellInterp

end program testNewRMS