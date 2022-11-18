
program sifx
    use kdefs
    use siftypes
    use sifstarter
    use planethelper
    use xml_input

    use shellgrid

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

    associate(sh=>Grid%shGrid)

    write(*,*)"Grid:"
    write(*,*)" ShGrid:"
    write(*,*)"  Nt,Np=",sh%Nt,sh%Np
    write(*,*)"  min/max theta=",sh%minTheta,sh%maxTheta
    write(*,*)"  min/max phi  =",sh%minPhi  ,sh%maxPhi
    write(*,*)"  theta  =",sh%th(sh%is:sh%ie)
    write(*,*)"  phi    =",sh%ph(sh%js:sh%je)


    end associate

end program sifx