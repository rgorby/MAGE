
program sifx
    use kdefs
    use planethelper
    use xml_input

    use sifdefs
    use sifstarter

    use shellgrid

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



    call sifInit(sifApp, inpXML)

    !associate(sh=>Grid%shGrid)
    !write(*,*)"Grid:"
    !write(*,*)" ShGrid:"
    !write(*,*)"  Nt,Np=",sh%Nt,sh%Np
    !write(*,*)"  min/max theta=",sh%minTheta,sh%maxTheta
    !write(*,*)"  min/max phi  =",sh%minPhi  ,sh%maxPhi
    !write(*,*)"  theta  =",sh%th(sh%is:sh%ie)
    !write(*,*)"  phi    =",sh%ph(sh%js:sh%je)
    !end associate


end program sifx