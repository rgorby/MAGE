! Main data objects and functions to perform a gamera simulation

module gamapp_mpi
    use gamtypes
    use step
    use init
    use mhdgroup
    use gamapp
    use mpi

    implicit none

    type, extends(GamApp_T) :: gamAppMpi_T
        integer :: NumRi=1,NumRj=1,NumRk=1
        integer :: Ri=1,Rj=1,Rk=1
        type(MPI_Comm) :: gamMpiComm
    end type gamAppMpi_T

    contains

    subroutine initGamera_mpi(gamAppMpi, userInitFunc, optFilename)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        
        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
            call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")
        else
            !Find input deck
            call getIDeckStr(inpXML)

        endif

        write(*,*) 'Reading input deck from ', trim(inpXML)
        inquire(file=inpXML,exist=fExist)
        if (.not. fExist) then
            write(*,*) 'Error opening input deck, exiting ...'
            write(*,*) ''
            stop
        endif
        
        !Create XML reader
        xmlInp = New_XML_Input(trim(inpXML),'Gamera',.true.)

        !Initialize Grid/State/Model (Hatch Gamera)
        !Will enforce 1st BCs, caculate 1st timestep, set oldState
        call Hatch_mpi(gamAppMpi,xmlInp,userInitFunc)
        call cleanClocks()

        if (.not. gamAppMpi%Model%isRestart) call fOutput(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)
        call consoleOutput(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State)

    end subroutine initGamera_mpi

    subroutine Hatch_mpi(gamAppMpi,xmlInp,userInitFunc,endTime)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
        type(XML_Input_T), intent(inout) :: xmlInp
        procedure(StateIC_T), pointer, intent(in) :: userInitFunc
        real(rp), optional, intent(in) :: endTime

        ! call appropriate subroutines to read corner info and mesh size data
        call ReadCorners(gamAppMpi%Model,gamAppMpi%Grid,xmlInp,endTime)

        ! adjust corner information to reflect this individual node's grid data
        gamAppMpi%Grid%ijkShift(1:3) = 0

        ! call appropriate subroutines to calculate all appropriate grid data from the corner data
        call CalcGridInfo(gamAppMpi%Model,gamAppMpi%Grid,gamAppMpi%State,gamAppMpi%oState,gamAppMpi%Solver,xmlInp,userInitFunc)

        ! check MPI values to assign responsibility for ring averaging
        if (gamAppMpi%Rj ==               1) gamAppMpi%Model%Ring%doS = .true.
        if (gamAppMpi%Rj == gamAppMpi%NumRj) gamAppMpi%Model%Ring%doE = .true.

    end subroutine Hatch_mpi

    subroutine stepGamera_mpi(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi

        !update the state variables to the next timestep
        call UpdateStateData(gamAppMpi)

        !Update ghost cells
        call Tic("Halos")
        call HaloUpdate(gamAppMpi)
        call Toc("Halos")

        ! calculate new DT, update BCs, write output data, etc...
        call FinishStep(gamAppMpi)

    end subroutine stepGamera_mpi

    subroutine haloUpdate(gamAppMpi)
        type(gamAppMpi_T), intent(inout) :: gamAppMpi
    end subroutine haloUpdate

end module gamapp_mpi

