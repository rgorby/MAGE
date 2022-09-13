! Main data objects and functions to perform a Gibson-Low CME simulation

module giblowapp
    use clocks
    use gltypes
    use glinit
    use glsolution

    implicit none

    type :: glApp_T
        type(glSolution_T) :: Solution
        type(glModel_T) :: Model
        type(glState_T) :: State
    end type glApp_T

    contains

    subroutine initGL(glApp, optFilename, doIO)
        type(glApp_T), intent(inout) :: glApp
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

        !Initialize Model, State, Solution
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
        call generateCMESolution(glApp%Solution,glApp%Model,glApp%State)
        call Toc("Solution")
    end subroutine step
end module