program testNewRMS

    use kdefs
    use XML_Input

    use remixReader

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
        t = t + dt
    enddo

end program testNewRMS