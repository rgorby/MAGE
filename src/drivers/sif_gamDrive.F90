program sif_gamdrive

    use kdefs
    use planethelper
    use xml_input

    ! Sif stuff
    use sifdefs
    use sifstarter

    ! Chimp stuff
    use ebtypes
    use starter
    use chmpfields

    implicit none

    type(sifApp_T   ) :: sifApp
    type(chmpModel_T) :: chmpModel
    type(ebState_T  ) :: ebState

    character(len=strLen) :: XMLStr
    type(XML_Input_T) :: inpXML    
    
    
    integer :: a = 0


    ! Init xml
    call getIDeckStr(XMLStr)

    ! Init SIF
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/SIF",.true.)
    call sifInit(sifApp, inpXML)

    ! Init CHIMP
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/CHIMP",.true.)


end program sif_gamdrive