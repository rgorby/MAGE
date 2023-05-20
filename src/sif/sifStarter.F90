! Initialize kaimag
module sifstarter
    use shellgrid
    use xml_input
    use planethelper

    use sifdefs
    use siftypes
    use sifgrids
    use sifetautils
    use sifout
    use sific
    use sifuseric

    implicit none

    contains


!------
! Main Initialization Routines
!------

    subroutine sifInit(app, iXML)
        type(sifApp_T), intent(inout) :: app
        type(XML_Input_T), intent(in) :: iXML

        ! Init model, grid, state
        call sifInitModel(app%Model, iXML)
        call sifInitGrid(app%Model, app%Grid, iXML)

        ! TODO: Handle restart here. For now, assuming no restart

        ! Init output file
        call sifInitIO(app%Model, app%Grid)

        call sifInitState(app%Model,app%Grid,app%State,iXML)

        ! Initialize IOCLOCK
        call app%State%IO%init(iXML,app%State%t,app%State%ts)

    end subroutine sifInit

    ! Sets up Model, but Grid and State must be set up separately
    subroutine sifInitModel(Model, iXML)
        type(sifModel_T) , intent(inout) :: Model
        type(XML_Input_T), intent(in)    :: iXML
         
        character(len=strLen) :: tmpStr

        write(*,*) "sifInitModel is starting"

        ! Assuming that if being controlled by e.g. Voltron, someone else will set RunID accordingly
        ! If getting here without RunID being set, assume we're running in stand-alone
        ! (Currently no decisions are made based on being SA, just means we set the RunID ourselves)
        if (trim(Model%RunID) .eq. "") then
            write(*,*) "Setting default RunID to be sifSA"
            call iXML%Set_Val(Model%RunID, "prob/RunID","sifSA")  ! sif stand-alone
        endif

        ! Timing info, if provided
        call iXML%Set_Val(Model%t0  ,'time/T0',0.0)
        call iXML%Set_Val(Model%tFin,'time/tFin',60.0)
        call iXML%Set_Val(Model%dt  ,'time/dt',1.0)

        ! Config file
        call iXML%Set_Val(Model%configFName, "config/fname","sifconfig.h5")
        call CheckFileOrDie(Model%configFName,"Unable to open file")

        call iXML%Set_Val(Model%isMPI, "mpi/isMPI",.false.)
        if (Model%isMPI) then
            write(*,*) "MPI not implemented for SIF yet, dying."
            stop
        endif

        call iXML%Set_Val(Model%isRestart, "Kaiju/Gamera/restart/doRes",.false.)
            ! TODO: Is this the right place to look for restart? Could be given as an extra arg
        call iXML%Set_Val(Model%isLoud, "debug/isLoud",.false.)
        call iXML%Set_Val(Model%writeGhosts, "debug/writeGhosts",.false.)

        ! Plasmasphere settings
        call iXML%Set_Val(Model%doPlasmasphere, "plasmasphere/doPsphere",.false.)
        ! Determine number of species. First set default, then read from xml to overwrite if present
        if (Model%doPlasmasphere) then
            Model%nSpc = 3
        else
            Model%nSpc = 2
        endif
        call iXML%Set_Val(Model%nSpc, "prob/nSpc",Model%nSpc)


        ! Lambda channel settings
        call iXML%Set_Val(Model%doDynamicLambdaRanges, "lambdas/dynamicRanges",.false.)

        ! Determine which DP2eta mapping should be used (e.g. Maxwellian, Kappa, user-defined)
        call iXML%Set_Val(tmpStr, "moments/DP2EtaMap","MAXWELL")
        select case(tmpStr)
            case("MAXWELL")
                Model%dp2etaMap => Maxwell2Eta
            case("KAPPA")
                Model%dp2etaMap => Kappa2Eta
                call iXML%Set_Val(Model%kappa, "moments/kappa",6.0)
            ! TODO: Add user option which will point to whatever is in the chosen IC file
            case DEFAULT
                write(*,*) "SIF received unavailable DP2Eta mapping: ",tmpStr
                write(*,*) " Dying."
                stop
        end select

        ! Set planet params
        !! This should only be kept for as long as planet_T doesn't contain pointers
        !! In this current case, there should be a full copy to our own planet params
        call getPlanetParams(Model%planet, iXML)
        
        ! Set up timing
        !!TODO
        ! If we are running stand-alone, look for timing info inside SIF XML block
        ! If voltron is present, look for timing information there

        !TODO: Add flags to output certain data, like coupling information

    end subroutine sifInitModel


    subroutine sifInitGrid(Model, Grid, iXML, shGridO)
        type(sifModel_T) , intent(inout) :: Model
        type(sifGrid_T)  , intent(inout) :: Grid
        type(XML_Input_T), intent(in)   :: iXML
        type(ShellGrid_T), optional, intent(in)    :: shGridO

        character(len=strLen) :: tmpStr
        type(ShellGrid_T) :: shGrid

        ! Set grid params
        ! Grid settings
        call iXML%Set_Val(tmpStr, "grid/gType","UNISPH")

        ! Fill out Grid object depending on chosen method
        select case(tmpStr)
            case("UNISPH")
                Grid%gType = G_UNISPH
                ! Generate our own grid from scratch
                call sifGenUniSphGrid(Grid, iXML)
            case("SHGRID")
                Grid%gType = G_SHGRID
                ! Then we should be receiving a predefined ShellGrid that Voltron has set up
                if(present(shGridO)) then
                    shGrid = shGridO
                    call sifGenGridFromShGrid(Grid, shGrid)
                else
                    write(*,*) "SIF expecting a ShellGrid_T but didn't receive one. Dying."
                endif
            case DEFAULT
                write(*,*) "SIF Received invalid grid definition: ",Grid%gType
                write(*,*) " Dying."
                stop
        end select

        ! Finalize the spatial part of the grid
        call finalizeLLGrid(Grid)
        

        ! Now handle lambda grid
        Grid%nSpc = Model%nSpc  ! Make a copy of nSpc
        call iXML%Set_Val(Grid%ignoreConfigMismatch, "config/ignoreMismatch",.false.)
        call initLambdaGrid(Model, Grid, Model%configFName)

    end subroutine sifInitGrid


    subroutine sifInitState(Model, Grid, State, iXML)
        type(sifModel_T), intent(inout) :: Model
        type(sifGrid_T) , intent(in)    :: Grid
        type(sifState_T), intent(inout) :: State
        type(XML_Input_T), intent(in)   :: iXML

        character(len=strLen) :: icStr

    ! Allocate arrays
        ! Where we keep all our stuff
        allocate( State%eta (Grid%shGrid%Nt, Grid%shGrid%Np, Grid%Nk) )
        ! Interface and cell velocities
        allocate( State%iVel(Grid%shGrid%Nt+1, Grid%shGrid%Np+1, Grid%Nk, 2) )
        allocate( State%cVel(Grid%shGrid%Nt  , Grid%shGrid%Np  , Grid%Nk, 2) )
        ! Coupling input moments
        allocate( State%Pavg(Grid%shGrid%Nt, Grid%shGrid%Np, Grid%nSpc) )
        allocate( State%Davg(Grid%shGrid%Nt, Grid%shGrid%Np, Grid%nSpc) )
        ! Bmin surface
        allocate( State%Bmin  (Grid%shGrid%Nt  , Grid%shGrid%Np  , 3 ) )
        allocate( State%xyzMin(Grid%shGrid%Nt+1, Grid%shGrid%Np+1, 3 ) )
        ! 2D corner quantities
        allocate( State%topo  (Grid%shGrid%Nt+1, Grid%shGrid%Np+1) )
        ! 2D cell-centered quantities
        allocate( State%active(Grid%shGrid%Nt, Grid%shGrid%Np) )
        allocate( State%thc   (Grid%shGrid%Nt, Grid%shGrid%Np) )
        allocate( State%phc   (Grid%shGrid%Nt, Grid%shGrid%Np) )
        allocate( State%espot (Grid%shGrid%Nt, Grid%shGrid%Np) )
        allocate( State%bvol  (Grid%shGrid%Nt, Grid%shGrid%Np) )
        ! Coupling output data
        allocate( State%precipFlux(Grid%shGrid%Nt, Grid%shGrid%Np, Grid%Nk) )
        allocate( State%precipAvgE(Grid%shGrid%Nt, Grid%shGrid%Np, Grid%Nk) )
        allocate( State%Den  (Grid%shGrid%Nt, Grid%shGrid%Np, Grid%nSpc+1) )
        allocate( State%Press(Grid%shGrid%Nt, Grid%shGrid%Np, Grid%nSpc+1) )
        allocate( State%vAvg (Grid%shGrid%Nt, Grid%shGrid%Np, Grid%nSpc+1) )

        !! TODO: Handle restart somewhere around here
        ! For now, just set t to tStart and ts to 0
        State%t = Model%t0
        State%ts = 0

    ! Init state
        call iXML%Set_Val(icStr, "prob/IC","DIP")
        select case(trim(icStr))
            case("DIP")
                !! Simple Maxwellian in a dipole field
                Model%initState => initSifIC_DIP
            case("USER")
                ! Call the IC in the module sifuseric
                ! This module is set in cmake via the SIFIC variable
                call SIFinitStateUserIC(Model, Grid, State, iXML)
            case DEFAULT
                write(*,*)"Invalid IC name to SIF, see sifStarter.F90:sifInitState. Bye."
                stop
        end select

        call Model%initState(Model, Grid, State, iXML)

    end subroutine sifInitState


!------
! Defaults
!------




end module sifstarter