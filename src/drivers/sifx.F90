
program sifx
    use kdefs
    use siftypes
    use sifstarter
    use planethelper
    use xml_input

    implicit none


    character(len=strLen) :: XMLStr
    type(XML_Input_T) :: inpXML

    ! Temporary, will be mored to app later
    type(sifModel_T) :: Model
    type(sifGrid_T)  :: Grid
    type(planet_T)   :: planet

    call getIDeckStr(XMLStr)
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/SIF",.true.)

    call sifInitModel(Model, Grid, planet, inpXML)


end program sifx