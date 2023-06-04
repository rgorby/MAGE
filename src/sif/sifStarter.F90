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
    use sifICHelpers

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

        ! Certain physical constants that shouldn't be constants
        call iXML%Set_Val(Model%tiote, "prob/tiote",4.0)
        call iXML%Set_Val(Model%worthyFrac, "prob/worthyFrac",fracWorthyDef)

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

        call iXML%Set_Val(Model%doFatOutput, "output/doFat",.false.)

        ! Set planet params
        !! This should only be kept for as long as planet_T doesn't contain pointers
        !! In this current case, there should be a full copy to our own planet params
        call getPlanetParams(Model%planet, iXML)


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
        call iXML%Set_Val(Grid%nB, "grid/Nbnd", 4  )  ! Number of cells between open boundary and active domain
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
        call finalizeLLGrid(Grid, Model%planet)
        

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

        ! Allocate arrays

        associate(sh=>Grid%shGrid)

            ! dt for every lambda channel
            allocate( State%dtk (Grid%Nk) )
            ! Where we keep all our stuff
            allocate( State%eta (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            ! Effective potential
            allocate( State%pEff(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            ! I shells shat should be evolved for each k
            allocate( State%activeShells (sh%isg:sh%ieg, Grid%Nk) )
            ! Interface and cell velocities
            allocate( State%iVel(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
            allocate( State%cVel(sh%isg:sh%ieg  , sh%jsg:sh%jeg  , Grid%Nk, 2) )
            ! Coupling input moments
            allocate( State%Pavg(sh%isg:sh%ieg  , sh%jsg:sh%jeg, Grid%nSpc) )
            allocate( State%Davg(sh%isg:sh%ieg  , sh%jsg:sh%jeg, Grid%nSpc) )
            ! Bmin surface
            allocate( State%Bmin    (sh%isg:sh%ieg  , sh%jsg:sh%jeg  , 3 ) )
            allocate( State%xyzMin  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 3 ) )
            allocate( State%xyzMincc(sh%isg:sh%ieg  , sh%jsg:sh%jeg  , 3 ) )
            allocate( State%thcon   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1    ) )
            allocate( State%phcon   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1    ) )
            ! 2D corner quantities
            allocate( State%topo  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            ! 2D cell-centered quantities
            allocate( State%active (sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            allocate( State%OCBDist(sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            allocate( State%espot  (sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            allocate( State%bvol   (sh%isg:sh%ieg, sh%jsg:sh%jeg) )

            ! Coupling output data
            allocate( State%precipFlux(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%precipAvgE(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%Den  (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%nSpc+1) )
            allocate( State%Press(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%nSpc+1) )
            allocate( State%vAvg (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%nSpc+1) )
        
        end associate
        !! TODO: Handle restart somewhere around here
        ! For now, just set t to tStart and ts to 0
        State%t = Model%t0
        State%ts = 0

    ! Init state
        call iXML%Set_Val(Model%icStr, "prob/IC","DIP")
        select case(trim(Model%icStr))
            case("DIP")
                !! Simple Maxwellian in a dipole field
                Model%initState => initSifIC_DIP
            case("USER")
                ! Call the IC in the module sifuseric
                ! This module is set in cmake via the SIFIC variable
                !Model%initState => userInitStateFunc
                write(*,*)"User initState not yet implemented"
                stop
            case DEFAULT
                write(*,*)"Invalid IC name to SIF, see sifStarter.F90:sifInitState. Bye."
                stop
        end select

        call Model%initState(Model, Grid, State, iXML)

    end subroutine sifInitState


end module sifstarter