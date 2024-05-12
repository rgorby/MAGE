module raijuCpl
    !! Contains functions used to pass data between other models and RAIJU

    ! Base
    use math
    use planethelper
    use shellGrid
    use voltCplTypes

    ! Raiju
    use raijudefs
    use raijuTypes
    use raijuCplTypes
    use raijuCplHelpers
    
    ! Cmake points to this
    use raijuuseric

    implicit none

    contains

    subroutine raijuCpl_init(vApp, raiApp, cplBase)
        type(voltApp_T), intent(in) :: vApp
        type(raijuApp_T), intent(in) :: raiApp
        type(raiju_cplBase_T), intent(inout) :: cplBase

        integer, dimension(4) :: shGhosts

        associate(fromV=>cplBase%fromV, toV=>cplBase%toV, &
            !sh=>raiApp%Grid%shGrid, nFluidIn=>vApp%ebTrcApp%ebModel%nSpc)
            sh=>raiApp%Grid%shGrid, nFluidIn=>raiApp%Model%nFluidIn)

        ! Init fromV first
            ! Allocations
            allocate(fromV%magLines (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))
            allocate(fromV%ijTubes(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))

            ! Shell Grid inits
            !shGhosts(NORTH) = sh%Ngn
            !shGhosts(SOUTH) = sh%Ngs
            !shGhosts(EAST) = sh%Nge
            !shGhosts(WEST) = sh%Ngw
            shGhosts(1) = sh%Ngn
            shGhosts(2) = sh%Ngs
            shGhosts(3) = sh%Nge
            shGhosts(4) = sh%Ngw
            call GenChildShellGrid(sh, fromV%shGr, "raijuCpl", nGhosts=shGhosts)
            call initShellVar(fromV%shGr, SHGR_CORNER, fromV%pot)

            ! Initial values
            fromV%tLastUpdate = -1.0*HUGE
            fromV%pot%data = 0.0
            fromV%pot%mask = .true.

        ! Init toV next
            ! Allocations
            !Initial values
            toV%tLastUpdate = -1.0*HUGE

        end associate

        ! If using user IC, let user determine coupling
        !  (this assumes icStr was already set by raijuInitState)
        if (trim(raiApp%Model%icStr) .eq. "USER") then
            ! Set defaults, let user override if they want to
            cplBase%convertToRAIJU => raijuCpl_Volt2RAIJU
            cplBase%convertToVoltron => raijuCpl_RAIJU2Volt
            cplBase%fromV%mhd2spcMap => defaultMHD2SpcMap
            call raijuCpl_init_useric(vApp, raiApp, cplBase)
        else
            ! Point to default coupling functions
            cplBase%convertToRAIJU => raijuCpl_Volt2RAIJU
            cplBase%convertToVoltron => raijuCpl_RAIJU2Volt
            cplBase%fromV%mhd2spcMap => defaultMHD2SpcMap
        endif

    end subroutine raijuCpl_init

    
    subroutine raijuCpl_Volt2RAIJU(cplBase, vApp, raiApp)
        !! Take info from cplBase%fromV and put it into RAIJU state
        class(raiju_cplBase_T), intent(in) :: cplBase
        type(voltApp_T)  , intent(in   ) :: vApp
        type(raijuApp_T ), intent(inout) :: raiApp

        ! Populate raiju state with coupling info
        ! Tubes
        call imagTubes2RAIJU(raiApp%Model, raiApp%Grid, raiApp%State, &
                cplBase%fromV%ijTubes, &
                cplBase%fromV%mhd2spcMap)
        ! Potential
        raiApp%State%espot(:,:) = cplBase%fromV%pot%data(:,:)

    end subroutine raijuCpl_Volt2RAIJU


    subroutine raijuCpl_RAIJU2Volt(cplBase, vApp, raiApp)
        class(raiju_cplBase_T), intent(inout) :: cplBase
        type(voltApp_T), intent(inout) :: vApp
        type(raijuApp_T) , intent(in   ) :: raiApp
    end subroutine raijuCpl_RAIJU2Volt



end module raijuCpl