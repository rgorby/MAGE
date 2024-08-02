! Initialize kaimag
module raijustarter
    
    ! Base
    use shellgrid
    use xml_input
    use planethelper

    ! Raiju
    use raijudefs
    use raijutypes
    use raijugrids
    use raijuetautils
    use raijuRecon, only : maxOrderSupported
    use raijuout
    use raijuICHelpers
    use raijuELossWM

    ! Cmake points to this
    use raijuuseric

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
        call raijuInitIO(app%Model, app%Grid, app%Model%writeGhosts)

        call raijuInitState(app%Model,app%Grid,app%State,iXML)

        ! Initialize IOCLOCK
        call app%State%IO%init(iXML,app%State%t,app%State%ts)

        ! Some sub-models need RAIJU to make up its mind before they can do their full init
        if (app%Model%doLosses .and. app%Model%eLossModel .eq. RaiELOSS_WM) then
            call initEWM(app%Model%eLossWM, app%Model%configFName, iXML, app%Grid%shGrid)
        endif


        ! TODO: Add some final checks before returning
        ! e.g., if we are mapping fluids, are we mapping something to HOTP?
        !  (totally just thought of it randomly and not because I spent 4 hours realizing this was my problem)

    end subroutine raijuInit


    ! Sets up Model; Grid and State must be set up separately
    subroutine raijuInitModel(Model, iXML)
        type(raijuModel_T), intent(inout) :: Model
        type(XML_Input_T) , intent(in)    :: iXML
         
        character(len=strLen) :: tmpStr
        real(rp) :: cfl0

        write(*,*) "raijuInitModel is starting"

        ! Assume we are running with others, but allow for solo run
        call iXML%Set_Val(Model%isSA, "driver/isSA", .false.)

        ! Assuming that if being controlled by e.g. Voltron, someone else will set RunID accordingly
        ! If we are in SA mode, need to set it ourselves
        if (Model%isSA) then
            call iXML%Set_Val(Model%RunID, "prob/RunID","raijuSA")  ! raiju stand-alone
        endif

        ! Timing info, if provided
        call iXML%Set_Val(Model%t0  ,'time/T0',0.0)
        call iXML%Set_Val(Model%tFin,'time/tFin',60.0)
        call iXML%Set_Val(Model%dt  ,'time/dt',1.0)

        ! Config file
        call iXML%Set_Val(Model%configFName, "config/fname","raijuconfig.h5")
        call CheckFileOrDie(Model%configFName,"RAIJU unable to open config file")

        call iXML%Set_Val(Model%isMPI, "mpi/isMPI",.false.)
        if (Model%isMPI) then
            write(*,*) "(RAIJU) MPI not implemented yet, dying."
            stop
        endif

        ! Restart time
        if (Model%isSA) then
            tmpStr = "restart/doRes"
        else
            tmpStr = "Kaiju/Gamera/restart/doRes"
        endif
        call iXML%Set_Val(Model%isRestart, trim(tmpStr),.false.)
        if (Model%isRestart) then
            if (Model%isSA) then
                tmpStr = "restart/nRes"
            else
                tmpStr = "Kaiju/Gamera/restart/nRes"
            endif
            call iXML%Set_Val(Model%nResIn, trim(tmpStr), Model%nResIn)
            call genResInFname(Model, Model%ResF)  ! Determine filename to read from
        endif

        call iXML%Set_Val(Model%isLoud, "debug/isLoud",.false.)
        call iXML%Set_Val(Model%writeGhosts, "debug/writeGhosts",.false.)
        call iXML%Set_Val(Model%doDebugOutput, "debug/debugOutput",.false.)
        
        ! Plasmasphere settings
        call iXML%Set_Val(Model%doPlasmasphere, "plasmasphere/doPsphere",.false.)
        ! Determine number of species. First set default, then read from xml to overwrite if present
        if (Model%doPlasmasphere) then
            Model%nSpc = 3
        else
            Model%nSpc = 2
        endif
        call iXML%Set_Val(Model%doExcessToPsph, "prob/doExcessMap",.true.)
            !! Allow mapping of excess H+ to plasmasphere channel
        if (Model%doExcesstoPsph .and. .not. Model%doPlasmasphere) then
            write(*,*)"Error: Can't map excess to psph because psph doesn't exist!"
        endif

        call iXML%Set_Val(Model%nSpc, "prob/nSpc",Model%nSpc)

        ! Domain constraints
        call iXML%Set_Val(Model%maxTail_buffer, "domain/tail_buffer", def_maxTail_buffer)
        call iXML%Set_Val(Model%maxSun_buffer , "domain/sun_buffer" , def_maxSun_buffer)
        call iXML%Set_Val(Model%maxTail_active, "domain/tail_active", def_maxTail_active)
        call iXML%Set_Val(Model%maxSun_active , "domain/sun_active" , def_maxSun_active)
        ! Store all distances as positive values, we'll add signs as needed later
        Model%maxTail_buffer = abs(Model%maxTail_buffer)
        Model%maxSun_buffer  = abs(Model%maxSun_buffer)
        Model%maxTail_active = abs(Model%maxTail_active)
        Model%maxSun_active  = abs(Model%maxSun_active)

        ! Solver params
        call iXML%Set_Val(Model%doUseVelLRs,'sim/useVelLRs',def_doUseVelLRs)
        call iXML%Set_Val(Model%maxItersPerSec,'sim/maxIter',def_maxItersPerSec)
        call iXML%Set_Val(Model%maxOrder,'sim/maxOrder',7)
        if (Model%maxOrder > maxOrderSupported) then
            write(*,*) "(RAIJU) Too much order, allowable orders <= ",maxOrderSupported
            stop
        endif
        call iXML%Set_Val(Model%PDMB,'sim/pdmb',def_pdmb)
        cfl0 = min(0.5/(Model%PDMB+0.5), cflMax) !Set CFL based on PDM

        !Set CFL from XML
        call iXML%Set_Val(Model%CFL ,'sim/cfl'  ,cfl0)
        if (Model%CFL > cfl0) then
            write(*,*) '-------------------------------------'
            write(*,*) '(RAIJU) WARNING, CFL is above critical value!'
            write(*,*) 'CFL/Critical/PDMB = ', Model%CFL,cfl0,Model%PDMB
            write(*,*) '-------------------------------------'
        else
            write(*,*) 'CFL # = ', Model%CFL
        endif

        ! Certain physical constants that shouldn't be constants
        call iXML%Set_Val(Model%tiote, "prob/tiote",4.0)

        ! Active shell settings
        call iXML%Set_Val(Model%doActiveShell, "activeShell/doAS",.true.)
        call iXML%Set_Val(Model%worthyFrac, "prob/worthyFrac",fracWorthyDef)

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
                write(*,*) "(RAIJU) Received unavailable DP2Eta mapping: ",tmpStr
                write(*,*) " Dying."
                stop
        end select


        ! Losses
        call iXML%Set_Val(Model%doLosses, "losses/doLosses",.true.)
        call iXML%Set_Val(Model%doSS    , "losses/doSS" ,.true. )
        call iXML%Set_Val(Model%doCC    , "losses/doCC" ,.true. )
        call iXML%Set_Val(Model%doCX    , "losses/doCX" ,.false.)
        call iXML%Set_Val(Model%doFLC   , "losses/doFLC",.false.)

        ! Electron loss model
        call iXML%Set_Val(tmpStr, "losses/eLossModel","WM")
        select case (tmpStr)
            case ("WM")
                write(*,*) "(RAIJU) Using Wang-Bao electron wave model"
                Model%eLossModel = RaiELOSS_WM
                Model%eLossRateFn => calcELossRate_WM
                ! We init later now, up in raijuInit, since initEWM needs Grid info for diagnostic output
            case default
                write(*,*) "(RAIJU) Did not get a valid electron loss model, goodbye"
                stop
        end select

        call iXML%Set_Val(Model%doFatOutput, "output/doFat",.false.)
        ! TODO: Add flags to output certain data, like coupling information

        ! Model Hax?
        call iXML%Set_Val(Model%doHack_rcmEtaPush, "hax/rcmEtaPush",.false.)
        

        ! Determine coupling information
        call setMHD2RaiInfo(Model, iXML)

        ! Set planet params
        call getPlanetParams(Model%planet, iXML, doLoudO=.true.)

    end subroutine raijuInitModel


    subroutine setMHD2RaiInfo(Model, iXML)
        type(raijuModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: iXML

        integer :: nFluids
        integer :: f
        character(len=strLen) :: xmlFluidTag


        ! Domain determination
        call iXML%Set_Val(Model%vaFracThresh, "cpl/vaFracThresh", def_vaFracThresh)
        call iXML%Set_Val(Model%bminThresh, "cpl/bminThresh", def_bminThresh)

        ! Fluid mapping
        call iXML%Set_Val(nFluids, "cpl/nFluidsIn",0)
        Model%nFluidIn = nFluids
        
        if (nFluids > 0) then
            allocate(Model%fluidInMaps(nFluids))
            do f=1,nFluids
                write(xmlFluidTag,'(A,I0,A)')"fluidIn",f,"/imhd"
                call iXML%Set_Val(Model%fluidInMaps(f)%idx_mhd, xmlFluidTag, 0)
                write(xmlFluidTag,'(A,I0,A)')"fluidIn",f,"/flav"
                call iXML%Set_Val(Model%fluidInMaps(f)%flav, xmlFluidTag, -1)
                write(xmlFluidTag,'(A,I0,A)')"fluidIn",f,"/excessToPsph"
                call iXML%Set_Val(Model%fluidInMaps(f)%doExcessToPsph, xmlFluidTag, .false.)
            enddo
        endif

    end subroutine setMHD2RaiInfo


    subroutine raijuInitGrid(Model, Grid, iXML, shGridO)
        type(raijuModel_T) , intent(inout) :: Model
        type(raijuGrid_T)  , intent(inout) :: Grid
        type(XML_Input_T), intent(in)   :: iXML
        type(ShellGrid_T), optional, intent(in)    :: shGridO

        character(len=strLen) :: tmpStr
        type(ShellGrid_T) :: shGrid
        logical :: doHack_constB

        ! Set grid params
        call iXML%Set_Val(Grid%nB, "grid/Nbnd", 4      )  ! Number of cells between open boundary and active domain
        call iXML%Set_Val(tmpStr, "grid/gType","UNISPH")

        if (.not. Model%isRestart) then
            ! Fill out Grid object depending on chosen method
            select case(tmpStr)
                case("UNISPH")
                    Grid%gType = RAI_G_UNISPH
                    ! Generate our own grid from scratch
                    call raijuGenUniSphGrid(Model, Grid, iXML)
                case("WARPSPH")
                    Grid%gType = RAI_G_WARPSPH
                    ! Generate our own grid from scratch
                    !call raijuGenWarpSphGrid_Fok2021(Model, Grid, iXML)
                    call raijuGenWarpSphGrid_Shafee2008(Model, Grid, iXML)
                case("SHGRID")
                    Grid%gType = RAI_G_SHGRID
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
        else
            call GenShellGridFromFile(Grid%shGrid, RAI_SG_NAME, Model%ResF)
        endif

        ! Finalize the spatial part of the grid
        call finalizeLLGrid(Grid, Model%planet)

        ! Hax?
        call iXML%Set_Val(doHack_constB, "hax/doConstB",.false.)
        if (doHack_constB) then
            Grid%Bmag = Model%planet%magMoment*G2nT
            Grid%Brcc = Model%planet%magMoment*G2nT
            Grid%BrFace = Model%planet%magMoment*G2nT
        endif

        ! Now handle lambda grid
        Grid%nSpc = Model%nSpc  ! Make a copy of nSpc
        Grid%nFluidIn = Model%nFluidIn  ! Make a copy of nFluidIn
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
            ! nSteps for each channel
            allocate( State%nStepk(Grid%Nk) )
            ! Where we keep all our stuff
            allocate( State%eta      (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            ! Where we keep all our stuff but a half-step ahead of now
            allocate( State%eta_half (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            ! Where we kept all our stuff one step ago
            allocate( State%eta_last (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            ! I shells shat should be evolved for each k
            allocate( State%activeShells (sh%isg:sh%ieg, Grid%Nk) )
            ! Effective potential (used for output only)
            allocate( State%pEff(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk) )
            ! Gradient of ionspheric potential
            allocate( State%gradPotE     (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 2) )
            ! Gradient of corotation potential
            allocate( State%gradPotCorot (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 2) )
            ! Gradient of (flux tube volume ^ -2/3)
            allocate( State%gradVM      (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 2) )
            ! Interface and cell velocities
            allocate( State%gradPotE_cc    (sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
            allocate( State%gradPotCorot_cc(sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
            allocate( State%gradVM_cc      (sh%isg:sh%ieg, sh%jsg:sh%jeg, 2) )
            allocate( State%iVel (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
            allocate( State%iVelL(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
            allocate( State%iVelR(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
            allocate( State%cVel (sh%isg:sh%ieg  , sh%jsg:sh%jeg  , Grid%Nk, 2) )
            
            ! Coupling input moments
            allocate( State%Pavg(sh%isg:sh%ieg  , sh%jsg:sh%jeg, 0:Grid%nFluidIn) )
            allocate( State%Davg(sh%isg:sh%ieg  , sh%jsg:sh%jeg, 0:Grid%nFluidIn) )
            ! Bmin surface
            allocate( State%Bmin    (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 3 ) )
            allocate( State%xyzMin  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, 3 ) )
            allocate( State%xyzMincc(sh%isg:sh%ieg  , sh%jsg:sh%jeg  , 3 ) )
            allocate( State%thcon   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1    ) )
            allocate( State%phcon   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1    ) )
            ! 2D corner quantities
            allocate( State%topo   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            allocate( State%espot  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            allocate( State%bvol   (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            allocate( State%bvol_cc(sh%isg:sh%ieg  , sh%jsg:sh%jeg  ) )
            allocate( State%vaFrac(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1) )
            ! 2D cell-centered quantities
            allocate( State%active      (sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            allocate( State%active_last (sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            allocate( State%OCBDist(sh%isg:sh%ieg, sh%jsg:sh%jeg) )

            ! Coupling output data
            allocate( State%lossRates      (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%precipType_ele (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%lossRatesPrecip(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%precipNFlux    (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%precipEFlux    (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( State%Den  (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%nSpc+1) )
            allocate( State%Press(sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%nSpc+1) )
            allocate( State%vAvg (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%nSpc+1) )

            ! Only bother allocating persistent versions of debug stuff if we ned them
            if (Model%doDebugOutput) then
                allocate( State%etaFaceReconL(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                allocate( State%etaFaceReconR(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                allocate( State%etaFacePDML  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                allocate( State%etaFacePDMR  (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
                allocate( State%etaFlux      (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1, Grid%Nk, 2) )
            endif
        
        end associate
        
        if (Model%isRestart) then
            call raijuResInput(Model, Grid, State)
            return
        endif
        
        ! For now, just set t to tStart and ts to 0
        State%t = Model%t0
        State%ts = 0
        State%nStepk = 0

        ! Set all activeShells to true. There should be doActiveShell checks anywhere it is used, but just in case, make sure default has all shells on
        State%activeShells = .true.
        ! Similarly, set vaFrac to safe value in case stand-alone never writes to it
        State%vaFrac = 1.0

    ! Init state
        call iXML%Set_Val(Model%icStr, "prob/IC","DIP")
        select case(trim(Model%icStr))
            case("DIP")
                !! Simple Maxwellian in a dipole field
                Model%initState => initRaijuIC_DIP
            case("TTA")
                !! Test Theta Advection
                Model%initState => initRaijuIC_testThetaAdvection
            case("USER")
                ! Call the IC in the module raijuuseric
                ! This module is set in cmake via the RAIJUIC variable
                !Model%initState => userInitStateFunc
                call raijuInitState_useric(Model, Grid, State, iXML)
            case DEFAULT
                write(*,*)"Invalid IC name to RAIJU, see raijuStarter.F90:raijuInitState. Bye."
                stop
        end select

        call Model%initState(Model, Grid, State, iXML)

    end subroutine raijuInitState


end module raijustarter
