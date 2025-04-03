submodule (volttypes) raijuCplTypesSub

    use raijutypes
    use raijustarter
    use raijuCplHelper
    use raijuColdStartHelper
    use raijuDomain

    use shellInterp
    use imaghelper
    use dstutils
    use math
    use ieee_arithmetic

    implicit none

    contains

    module subroutine raiCplInitModel(App, xml)
        class(raijuCoupler_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        ! Allocate our contained raiju app
        allocate(raijuApp_T :: App%raiApp)
        ! Init raiju app itself
        !call App%raiApp%InitModel(xml)
        ! Note: we are bypassing raiApp%InitModel because we can't pass the voltron grid as an argument
        call raijuInit(App%raiApp, xml, App%opt%voltGrid)
        ! Update MJD with whatever voltron handed us
        ! If we are restarting, this will get replaced with whatever's in file later
        App%raiApp%State%mjd = App%opt%mjd0
        write(*,*)"MJD0=",App%opt%mjd0
        if (App%opt%doColdStart) then
            ! We are gonna cold start, so ignore plasma ingestion rules for first coupling
            App%raiApp%State%isFirstCpl = .false.
        endif
        ! Then allocate and initialize coupling variables based on raiju app
        call raijuCpl_init(App, xml)

    end subroutine raiCplInitModel


    module subroutine volt2RAIJU(App, vApp)
        class(raijuCoupler_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: vApp

        logical :: doFirstColdStart
        logical :: doUpdateColdStart
        real(rp) :: BSDst

        doFirstColdStart = .false.
        doUpdateColdStart = .false.

        associate(raiApp=>App%raiApp)

            ! If we are running realtime, its our job to get everything we need from vApp into raiCpl
            if (.not. App%raiApp%Model%isSA) then
                ! First, determine if we should cold start, i.e. Completely reset raiju's eta's to match some target conditions
                ! Determine if we should cold start before packing coupler because it will set tLastUpdate to vApp%time and then we can't do the checks we want
                ! But actually do cold start after coupler packing completes so we can use real field line info

                ! Do we do our very first coldstart ever
                if (App%opt%doColdStart .and. App%tLastUpdate < 0.0 .and. vApp%time >= 0.0) then
                    doFirstColdStart = .true.
                endif
                ! Do we do "updates" to our coldstart during pre-conditioning period
                if(App%opt%doColdStart .and. App%tLastUpdate > 0.0 .and. vApp%time < App%startup_blendTscl) then
                    doUpdateColdStart = .true.
                endif

                call packRaijuCoupler_RT(App, vApp)
            endif

            ! Someone updated raiCpl's coupling variables by now, stuff it into RAIJU proper
            call raiCpl2RAIJU(App)

            if (doFirstColdStart) then
                ! Its happening, everybody stay calm
                write(*,*) "RAIJU Doing first cold start..."
                ! NOTE: By this point we have put coupling info into raiju (e.g. bVol, xyzmin, MHD moments)
                ! But haven't calculated active domain yet because that happens in preadvancer
                ! So we jump in and do it here so we have it for cold starting
                call setActiveDomain(raiApp%Model, raiApp%Grid, raiApp%State)
                ! Calc voltron dst ourselves since vApp%BSDst is only set on console output
                call EstDST(vApp%gApp%Model,vApp%gApp%Grid,vApp%gApp%State,BSDst0=BSDst)
                call raijuGeoColdStart(raiApp%Model, raiApp%Grid, raiApp%State, vApp%time, BSDst, doCXO=App%doColdstartCX,doPsphO=.true.)
            endif
            if (doUpdateColdStart) then
                write(*,*)"RAIJU doing update cold start at t=",vApp%time
                write(*,*)" (calculating model BSDst,)",vApp%time
                call setActiveDomain(raiApp%Model, raiApp%Grid, raiApp%State)
                call EstDST(vApp%gApp%Model,vApp%gApp%Grid,vApp%gApp%State,BSDst0=BSDst)
                call raijuGeoColdStart(raiApp%Model, raiApp%Grid, raiApp%State, vApp%time, BSDst, doCXO=App%doColdstartCX,doPsphO=.false.)
            endif
        end associate
    end subroutine volt2RAIJU


    module subroutine getMomentsRAIJU(App,th,ph,t,imW,isEdible)
        !! Get fluid moments from RAIJU, formatted by enum in volttypes.F90 for a given theta,phi location within the RAIJU domain
        class(raijuCoupler_T), intent(inout) :: App
        real(rp), intent(in) :: th
            !! Theta [rad]
        real(rp), intent(in) :: ph
            !! Phi [rad]
        real(rp), intent(in) :: t
            !! Time since run start [s]
        real(rp), intent(out) :: imW(IM_D_RING:IM_TSCL)
        logical , intent(out) :: isEdible

        integer :: s  ! Iterators
        integer :: i0, j0  ! i,j cell that provided th,ph are in
        real(rp) :: active_interp
        real(rp) :: d_cold, t_cold, d_hot, p_hot
        real(rp) :: tScl, rampC

        associate(Model=>App%raiApp%Model, State=>App%raiApp%State, sh=>App%raiApp%Grid%shGrid, spcList=>App%raiApp%Grid%spc)

        ! Default
        imW = 0.0
        isEdible = .false.

        d_cold = 0
        t_cold = TINY
        d_hot = 0
        p_hot = 0

        ! Is this a good point?
        if (th < sh%minTheta .or. th > sh%maxTheta) then
            return ! Off grid, return default
        endif

        ! Active check
        call getSGCellILoc(sh, th, i0)
        call getSGCellJLoc(sh, ph, j0)
        if (State%active(i0,j0) .ne. RAIJUACTIVE) then
            return
        endif

        ! Otherwise we are good, gonna return stuff
        isEdible = .true.

        do s=1, Model%nSpc
            if (spcList(s)%flav == F_PSPH) then
                call InterpShellVar_TSC_pnt(sh, State%Den(s), th, ph, d_cold)
                imW(IM_D_COLD) = d_cold  ! [#/cc]
                t_cold = PsphTemp_Genestreti(d_cold)  ! [keV]
                imW(IM_P_COLD) = DkT2P(d_cold, t_cold)  ! [nPa]
            elseif (spcList(s)%spcType == RAIJUHPLUS) then
                call InterpShellVar_TSC_pnt(sh, State%Den(s)  , th, ph, d_hot)
                call InterpShellVar_TSC_pnt(sh, State%Press(s), th, ph, p_hot)
                imW(IM_D_RING) = imW(IM_D_RING) + d_hot  ! [nPa]
                imW(IM_P_RING) = imW(IM_P_RING) + p_hot  ! [#/cc]
            elseif (spcList(s)%spcType == RAIJUELE) then
                ! Don't add number density
                call InterpShellVar_TSC_pnt(sh, State%Press(s), th, ph, p_hot)
                imW(IM_P_RING) = imW(IM_P_RING) + p_hot  ! [#/cc]
            endif
        enddo

        
        !call InterpShellVar_TSC_pnt(sh, State%Tb, th, ph, imW(IM_TSCL))
        !imW(IM_TSCL) = Model%nBounce*imW(IM_TSCL)  ! [s]
        tScl = 10.0_rp  ! [s]
        !tScl = 10.0_rp/App%vaFrac(i0,j0)  ! [s]

        ! Adjust IM_TSCL if we wanna ramp up over time
        if (t < App%startup_blendTscl) then
            rampC = RampDown(t, 0.0_rp, App%startup_blendTscl)
            !tScl = sqrt(tScl*App%startup_blendTscl)*rampC + (1-rampC)*tScl  ! idk
            tScl = rampC*50.0_rp*tScl + (1-rampC)*tScl  ! No good reason for 50 except for wanting starting tScl to be ~8-10 minutes
        endif
        
        imW(IM_TSCL) = tScl

        end associate
    end subroutine getMomentsRAIJU


    module subroutine getMomentsPrecipRAIJU(App,rai_fluxes)
        class(raijuCoupler_T), intent(inout) :: App
        real(rp), dimension(:,:,:), allocatable, intent(out) :: rai_fluxes
        integer :: i,j,k,s
        integer :: is, ie, js, je, ks, ke
        real(rp) :: valPrecip  ! threshold value of grid status for accumulating precip.

        associate(Model=>App%raiApp%Model, State=>App%raiApp%State, sh=>App%raiApp%Grid%shGrid, spcList=>App%raiApp%Grid%spc)
        is = sh%is
        ie = sh%ie
        js = sh%js
        je = sh%je

        ! Calculate accumulated precipitation for each species
!        allocate(gtype        (is:ie,js:je))

        allocate(rai_fluxes (is:ie,js:je,nVars_imag2mix))
        rai_fluxes = 0.0_rp

        rai_fluxes(:,:,RAI_GTYPE) = (State%active(is:ie,js:je)*1.0_rp+1.0)/2.0  ! "-1=Inactive, 0=Buffer, 1=Active" -> 0, 0.5, and 1.0
        valPrecip = (RAIJUBUFFER*1.0_rp + 1.0)/2.0
        rai_fluxes(:,:,RAI_THCON) = State%thcon(is:ie,js:je) ! conjugate co-lat in radians, 0-pi
        rai_fluxes(:,:,RAI_PHCON) = State%phcon(is:ie,js:je) ! conjugate long in radians, 0-2pi
        do s=1, Model%nSpc
            ks = spcList(s)%kStart
            ke = spcList(s)%kEnd
            if (spcList(s)%flav == F_PSPH) then
                rai_fluxes(:,:,RAI_NPSP ) = rai_fluxes(:,:,RAI_NPSP ) + State%Den(s)%data*1.0e6 ! uStr="#/cc" -> /m^3 , dStr="Density from RAIJU flavors"
            elseif (spcList(s)%spcType == RAIJUHPLUS) then
                ! add proton precipitation later.
            elseif (spcList(s)%spcType == RAIJUELE) then
                rai_fluxes(:,:,RAI_ENFLX) = rai_fluxes(:,:,RAI_ENFLX) + sum(State%precipNFlux(is:ie,js:je,ks:ke), dim=3) ! uStr="#/cm^2/s"
                rai_fluxes(:,:,RAI_EFLUX) = rai_fluxes(:,:,RAI_EFLUX) + sum(State%precipEFlux(is:ie,js:je,ks:ke), dim=3) ! uStr="erg/cm^2/s"
                rai_fluxes(:,:,RAI_CCHF ) = rai_fluxes(:,:,RAI_CCHF ) + sum(State%CCHeatFlux (is:ie,js:je,ks:ke), dim=3) ! uStr="erg/cm^2/s"
                rai_fluxes(:,:,RAI_EPRE ) = rai_fluxes(:,:,RAI_EPRE ) + State%Press(s)%data*1.0e-9 ! uStr="nPa" -> Pa , dStr="Pressure from RAIJU flavors"
                rai_fluxes(:,:,RAI_EDEN ) = rai_fluxes(:,:,RAI_EDEN ) + State%Den  (s)%data*1.0e6  ! uStr="#/cc" -> /m^3 , dStr="Density from RAIJU flavors"
            endif
        enddo

        where (rai_fluxes(:,:,RAI_GTYPE) < valPrecip)
            ! zero out fluxes on Inactive grids.
            rai_fluxes(:,:,RAI_NPSP ) = 0.0_rp
            rai_fluxes(:,:,RAI_ENFLX) = 0.0_rp
            rai_fluxes(:,:,RAI_EFLUX) = 0.0_rp
            rai_fluxes(:,:,RAI_CCHF ) = 0.0_rp
            rai_fluxes(:,:,RAI_EPRE ) = 0.0_rp
            rai_fluxes(:,:,RAI_EDEN ) = 0.0_rp
            rai_fluxes(:,:,RAI_EAVG ) = 0.0_rp
        end where

        where (rai_fluxes(:,:,RAI_ENFLX) > TINY)
            ! derive mean energy where nflux is non-trivial.
            rai_fluxes(:,:,RAI_EAVG) = rai_fluxes(:,:,RAI_EFLUX)/rai_fluxes(:,:,RAI_ENFLX) * erg2kev  ! Avg E [keV]
        end where
        end associate
    end subroutine getMomentsPrecipRAIJU


!------
! raijuCoupler_T prodecure => raijuApp_T procedure
!------
    module subroutine raiCplInitIO(App, xml)
        class(raijuCoupler_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        call App%raiApp%InitIO(xml)

    end subroutine raiCplInitIO

    module subroutine raiCplWriteRestart(App, nRes)
        class(raijuCoupler_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! Write raiApp info into runid.raiju.Res.h5
        call App%raiApp%WriteRestart(nRes)
        ! And now the coupler's info
        call writeRaiCplRes(App, nRes)

    end subroutine

    module subroutine raiCplReadRestart(App, resId, nRes)
        class(raijuCoupler_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        ! Read runid.raiju.Res.h5 into raiApp
        call App%raiApp%ReadRestart(resId, nRes)
        ! Now coupler
        call readRaiCplRes(App,resId,nRes)

    end subroutine

    module subroutine raiCplWriteConsoleOutput(App)
        class(raijuCoupler_T), intent(inout) :: App

        if (App%tLastUpdate < 0) then
            ! We haven't even started yet, nothing to report
            return
        endif
        call App%raiApp%WriteConsoleOutput()

    end subroutine

    module subroutine raiCplWriteFileOutput(App, nStep)
        class(raijuCoupler_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call App%raiApp%WriteFileOutput(nStep)

    end subroutine

    module subroutine raiCplWriteSlimFileOutput(App, nStep)
        class(raijuCoupler_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call App%raiApp%WriteSlimFileOutput(nStep)

    end subroutine

    module subroutine raiCplAdvanceModel(App, dt)
        class(raijuCoupler_T), intent(inout) :: App
        real(rp), intent(in) :: dt

        call App%raiApp%AdvanceModel(dt)
        App%raiApp%State%t = App%raiApp%State%t + dt
        App%raiApp%State%ts = App%raiApp%State%ts + 1
        App%raiApp%State%mjd = T2MJD(dt,App%raiApp%State%mjd)

    end subroutine

    module subroutine raiCplCleanup(App)
        class(raijuCoupler_T), intent(inout) :: App

        call App%raiApp%Cleanup()

    end subroutine
    
end submodule raijuCplTypesSub