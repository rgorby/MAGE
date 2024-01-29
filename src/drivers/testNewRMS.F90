program testNewRMS

    use kdefs
    use XML_Input

    use remixReader

    implicit none

    character(len=strLen) :: xmlStr = 'testRM.xml'
    character(len=strLen) :: ftag = 'msphere'
    type(XML_Input_T) :: inpXML
    type(rmState_T) :: rmState


    inpXML = New_XML_Input(trim(xmlStr),"Kaiju/REMIX",.true.)

    call initRM(ftag, inpXML, rmState)

    write(*,*) "th1D:"
    write(*,*) rmState%shGrid%th
    write(*,*) "ph1D:"
    write(*,*) rmState%shGrid%ph


end program testNewRMS