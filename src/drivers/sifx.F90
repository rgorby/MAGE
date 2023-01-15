
program sifx
    use kdefs
    use planethelper
    use xml_input

    use sifdefs
    use sifstarter
    use sifio

    use shellgrid

    implicit none

    ! Could use sifApp, but that depends on voltlib which is annoying to build
    type(sifModel_T) :: Model
    type(sifGrid_T) :: Grid
    type(sifState_T) :: State

    type(planet_T) :: planet

    character(len=strLen) :: XMLStr
    type(XML_Input_T) :: inpXML

    
    ! Init xml
    call getIDeckStr(XMLStr)
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/SIF",.true.)

    ! Init planet specs, default Earth
    call getPlanetParams(planet, inpXML)

    ! Init model, grid, state
    call sifInitModel(Model, planet, inpXML)
    call sifInitGrid(Model, Grid, inpXML)

    ! Init output file
    call sifInitIO(Model, Grid)

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