! Main data objects and functions to perform a Gibson-Low CME simulation

module giblowapp
    use gltypes
    use init
    use clocks
    use glsolution

    implicit none

    contains

    !------------------------------------------------------------------------------------
    ! Outline: coordchange -> fieldcalc -> getoutside -> getinside -> add background pressure
    ! Accept array (x,y,z) position from Gamera i.e. Grid%xyzcc(i,j,k,:) -> xyz(3)
    !  or stand-alone xyz(Nx,Ny,Nz) vector
    ! run through giblow model function for position at time t
    !   construct r,theta,phi location for xyz(3) position x -> xyz(1),  y -> xyz(2), z -> xyz(3)
    !   TODO:  solve eqn. 6 in gibson low 1998 > alpha/eta = 0.49 case for accelrating
    !           Need to implement accelerating, decelerating case -> not in IDL
    !           for now just do the alpha = 0 case, or alpha = alpha = 6.762d-8*velmult*velmult?
    !   calculate fields for r,theta,phi 
    !   calculate outside solution
    !   calculate inside solution
    !   return values of solution at xyz position for time
    !       Dens(3,nt)
    ! Init 
    !     
    !        Vals: Pressure, Dens, Temp, BR, BTH, BPH, VR, VTH, VPH, JR, JTH, JPH
    !        From IDL/OLDF - dens, vr, 0, 0, pres, brphys, bthphys, bphphys, inside_mask
    !    Set Defaults? Expect XML to define all parameters?
    !       Giblow precheck - which components necessary? Just defining params or actually
    !           additional information, i.e. volume, velocity? 
    ! Define coords
    ! ! Calc model: accepts shell tile for Gamera-Helio, or full 3D grid for standalone?
    !   Calc fields
    !   Calc outside Bc
    !   Calc inside Bc
    ! Convert (r,theta,phi) -> Gamera/Helio (x,y,z) h5? 
    ! Utilities:
    !   GL accept tile/export tile
    !   Rotate tile/grid from Gamera to r,theta,phi
    !------------------------------------------------------------------------------------


    subroutine initGL(glApp, optFilename, doIO)
        type(glApp_T), intent(inout) :: glApp
        !procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        character(len=*), optional, intent(in) :: optFilename
        logical, optional, intent(in) :: doIO

        character(len=strLen) :: inpXML,kaijuRoot
        type(XML_Input_T) :: xmlInp
        logical :: doIOX

        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
        else
            !Find input deck
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

        if (present(doIO)) then
            doIOX = doIO
        else
            doIOX = .true.
        endif

        !Create XML reader
        write(*,*) 'Reading input deck from ', trim(inpXML)
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju/GL',.true.)

           ! try to verify that the XML file has "Kaiju" as a root element
        kaijuRoot = ""
        call xmlInp%Get_Key_Val("/GL/sim/H5Grid",kaijuRoot,.false.)
        if(len(trim(kaijuRoot)) /= 0) then
            write(*,*) "The input XML appears to be of an old style."
            write(*,*) "As of June 12th, 2021 it needs a root element of <Kaiju>."
            write(*,*) "Please modify your XML config file by adding this line at the top:"
            write(*,*) "<Kaiju>"
            write(*,*) "and this line at the bottom:"
            write(*,*) "</Kaiju>"
            write(*,*) "OR (preferred) convert your configuration to an INI file and use"
            write(*,*) " the XMLGenerator.py script to create conforming XML files."
            write(*,*) "Please refer to the python script or"
            write(*,*) " the [Generating XML Files] wiki page for additional info."
            stop
        endif

        ! read debug flags
        !call xmlInp%Set_Val(writeGhosts,"debug/writeGhosts",.false.)
        !call xmlInp%Set_Val(writeMagFlux,"debug/writeMagFlux",.false.)

        !Initialize Grid/State/Model (Hatch Gamera)
        !Will enforce 1st BCs, caculate 1st timestep, set oldState
        call Tic("Init")
        call initGLStandalone(glApp%Model, glApp%State, glApp%Solution, xmlInp)
        call Toc("Init")
        call cleanClocks()

        ! TODO: Output
        ! if (doIOX) then
        !     call output(glApp%Model,glApp%State,glApp%Solution)
        !     call consoleOutput(glApp%Model,glApp%Grid,glApp%State)
        ! endif
    end subroutine initGL

    subroutine step(glApp)
        type(glApp_T), intent(inout) :: glApp
        !Calculate new timestep
        call Tic("Solution")
        call generateGLSolution(glApp%Model,glApp%State,glApp%Solution)
        call Toc("Solution")
    end subroutine step
end module