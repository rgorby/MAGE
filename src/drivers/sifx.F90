
program sifx
    use kdefs
    use planethelper
    use xml_input

    use sifdefs
    use sifstarter

    implicit none

    !type(sifModel_T) :: Model
    !type(sifGrid_T ) :: Grid
    !type(sifState_T) :: State

    type(sifApp_T) :: sifApp

    character(len=strLen) :: XMLStr
    type(XML_Input_T) :: inpXML

    
    ! Init xml
    call getIDeckStr(XMLStr)
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/SIF",.true.)

    call sifInitCore(sifApp, inpXML)

    ! Save first state to file
    call sifOutput(sifApp%Model,sifApp%Grid,sifApp%State)

    do while ( sifApp%State%t < (sifApp%Model%tFin + 0.5) )
        if (sifApp%State%IO%doOutput(sifApp%State%t)) then
            call sifOutput(sifApp%Model,sifApp%Grid,sifApp%State)
        endif

        ! No advance yet, just increase time by 1
        sifApp%State%t = sifApp%State%t + 1
        sifApp%State%ts = sifApp%State%ts + 1

    enddo


end program sifx