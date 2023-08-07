! Initialize kaimag
module raijustarter
    use shellgrid
    use xml_input
    use planethelper

    use raijudefs
    use raijutypes
    use raijugrids
    use raijuetautils
    use raijuout
    use raijuICHelpers

    implicit none

    contains


!------
! Main Initialization Routines
!------

    subroutine raijuInit(app, iXML)
        type(raijuApp_T), intent(inout) :: app
        type(XML_Input_T), intent(in) :: iXML
        ! Init model, grid, state
        call raijuInitModel(app%Model, iXML)
        call raijuInitGrid(app%Model, app%Grid, iXML)

        ! TODO: Handle restart here. For now, assuming no restart

        ! Init output file
        call raijuInitIO(app%Model, app%Grid)

        call raijuInitState(app%Model,app%Grid,app%State,iXML)

        ! Initialize IOCLOCK
        call app%State%IO%init(iXML,app%State%t,app%State%ts)

    end subroutine raijuInit

    ! Sets up Model, but Grid and State must be set up separately
    subroutine raijuInitModel(Model, iXML)
        type(raijuModel_T) , intent(inout) :: Model
        type(XML_Input_T), intent(in)    :: iXML
         
        character(len=strLen) :: tmpStr

        write(*,*) "raijuInitModel is starting"

        ! Assuming that if being controlled by e.g. Voltron, someone else will set RunID accordingly
        ! If getting here without RunID being set, assume we're running in stand-alone
        ! (Currently no decisions are made based on being SA, just means we set the RunID ourselves)
        if (trim(Model%RunID) .eq. "") then
            call iXML%Set_Val(Model%RunID, "prob/RunID","raijuSA")  ! raiju stand-alone
        endif

        ! Timing info, if provided
        call iXML%Set_Val(Model%t0  ,'time/T0',0.0)
        call iXML%Set_Val(Model%tFin,'time/tFin',60.0)
        call iXML%Set_Val(Model%dt  ,'time/dt',1.0)

        ! Config file
        call iXML%Set_Val(Model%configFName, "config/fname","raijuconfig.h5")
        call CheckFileOrDie(Model%configFName,"Unable to open file")

        call iXML%Set_Val(Model%isMPI, "mpi/isMPI",.false.)
        if (Model%isMPI) then
            write(*,*) "MPI not implemented for RAIJU yet, dying."
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
        call iXML%Set_Val(Model%doLosses, "prob/doLosses",.true.)
        call iXML%Set_Val(Model%doGeoCorot, "prob/doGeoCorot",.false.)

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
                write(*,*) "RAIJU received unavailable DP2Eta mapping: ",tmpStr
                write(*,*) " Dying."
                stop
        end select


        ! Losses
        call iXML%Set_Val(Model%doSS , "losses/doSS" ,.true. )
        call iXML%Set_Val(Model%doCC , "losses/doCC" ,.true.)
        call iXML%Set_Val(Model%doCX , "losses/doCX" ,.false.)
        call iXML%Set_Val(Model%doFLC, "losses/doFLC",.false.)


        call iXML%Set_Val(Model%doFatOutput, "output/doFat",.false.)

        ! Set planet params
        !! This should only be kept for as long as planet_T doesn't contain pointers
        !! In this current case, there should be a full copy to our own planet params
        call getPlanetParams(Model%planet, iXML)


        !TODO: Add flags to output certain data, like coupling information

    end subroutine raijuInitModel


    subroutine raijuInitGrid(Model, Grid, iXML, shGridO)
        type(raijuModel_T) , intent(inout) :: Model
        type(raijuGrid_T)  , intent(inout) :: Grid
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
                call raijuGenUniSphGrid(Grid, iXML)
            case("SHGRID")
                Grid%gType = G_SHGRID
                ! Then we should be receiving a predefined ShellGrid that Voltron has set up
                if(present(shGridO)) then
                    shGrid = shGridO
                    call raijuGenGridFromShGrid(Grid, shGrid)
                else
                    write(*,*) "RAIJU expecting a ShellGrid_T but didn't receive one. Dying."
                endif
            case DEFAULT
                write(*,*) "RAIJU Received invalid grid definition: ",Grid%gType
                write(*,*) " Dying."
                stop
        end select

        ! Finalize the spatial part of the grid
        call finalizeLLGrid(Grid, Model%planet)
        

        ! Now handle lambda grid
        Grid%nSpc = Model%nSpc  ! Make a copy of nSpc
        call iXML%Set_Val(Grid%ignoreConfigMismatch, "config/ignoreMismatch",.false.)
        call initLambdaGrid(Model, Grid, Model%configFName)

    end subroutine raijuInitGrid


    subroutine raijuInitState(Model, Grid, State, iXML)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T) , intent(in)    :: Grid
        type(raijuState_T), intent(inout) :: State
        type(XML_Input_T), intent(in)   :: iXML

        ! Allocate arrays

        associate(sh=>Grid%shGrid)

            ! dt for every lambda channel
            allocate( State%dtk (Grid%Nk) )
            ! Where we keep all our stuff
            allocate( State%eta (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            ! I shells shat should be evolved for each k
            allocate( State%activeShells (sh%isg:sh%ieg, Grid%Nk) )
            ! Effective potential (used for output only)
            allocate( State%pEff(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            ! Gradient of ionspheric potential
            allocate( State%gradPotE     (sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
            ! Gradient of corotation potential
            allocate( State%gradPotCorot (sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
            ! Gradient of (flux tube volume ^ -2/3)
            allocate( State%gradVM      (sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
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
            allocate( State%precipNFlux(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%precipEFlux(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
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
                Model%initState => initRaijuIC_DIP
            case("USER")
                ! Call the IC in the module raijuuseric
                ! This module is set in cmake via the RAIJUIC variable
                !Model%initState => userInitStateFunc
                write(*,*)"User initState not yet implemented"
                stop
            case DEFAULT
                write(*,*)"Invalid IC name to RAIJU, see raijuStarter.F90:raijuInitState. Bye."
                stop
        end select

        call Model%initState(Model, Grid, State, iXML)

    end subroutine raijuInitState


end module raijustarter