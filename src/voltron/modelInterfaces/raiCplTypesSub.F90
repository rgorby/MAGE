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

    use planethelper

    implicit none

    contains

    module subroutine raiCplInitModel(App, xml)
        class(raijuCoupler_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        real(rp), dimension(NDIM) :: rxyz
        real(rp) :: lowLatBnd, thetaU

        ! Allocate our contained raiju app
        allocate(raijuApp_T :: App%raiApp)


        ! Set option overrides
        !  ThetaU - calculate where the upper theta boundary should be such that the active domain covers all of the gamera ghost cells
        rxyz(XDIR) = max(App%opt%mhdRin - 0.25, 1.0_rp)
        rxyz(YDIR) = 0.0_rp
        rxyz(ZDIR) = 0.0_rp
        lowLatBnd = InvLatitude(rxyz)*rad2deg
        thetaU = (90 - lowLatBnd)
        call App%raiApp%opt%thetaU%set(thetaU)
        ! Init raiju app itself
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

            if (.not. raiApp%State%coldStarter%doneFirstCS .or. vApp%time < raiApp%State%coldStarter%tEnd) then
                !! Make sure we run at least once
                call setActiveDomain(raiApp%Model, raiApp%Grid, raiApp%State)
                ! Calc voltron dst ourselves since vApp%BSDst is only set on console output
                call EstDST(vApp%gApp%Model,vApp%gApp%Grid,vApp%gApp%State,BSDst0=BSDst)
                call raijuGeoColdStart(raiApp%Model, raiApp%Grid, raiApp%State, vApp%time, BSDst)
            endif
            !if (doFirstColdStart) then
            !    ! Its happening, everybody stay calm
            !    write(*,*) "RAIJU Doing first cold start..."
            !    ! NOTE: By this point we have put coupling info into raiju (e.g. bVol, xyzmin, MHD moments)
            !    ! But haven't calculated active domain yet because that happens in preadvancer
            !    ! So we jump in and do it here so we have it for cold starting
            !    call setActiveDomain(raiApp%Model, raiApp%Grid, raiApp%State)
            !    ! Calc voltron dst ourselves since vApp%BSDst is only set on console output
            !    call EstDST(vApp%gApp%Model,vApp%gApp%Grid,vApp%gApp%State,BSDst0=BSDst)
            !    call raijuGeoColdStart(raiApp%Model, raiApp%Grid, raiApp%State, vApp%time, BSDst, doCXO=App%doColdstartCX,doPsphO=.true.)
            !endif
            !if (doUpdateColdStart) then
            !    write(*,*)"RAIJU doing update cold start at t=",vApp%time
            !    write(*,*)" (calculating model BSDst,)",vApp%time
            !    call setActiveDomain(raiApp%Model, raiApp%Grid, raiApp%State)
            !    call EstDST(vApp%gApp%Model,vApp%gApp%Grid,vApp%gApp%State,BSDst0=BSDst)
            !    call raijuGeoColdStart(raiApp%Model, raiApp%Grid, raiApp%State, vApp%time, BSDst, doCXO=App%doColdstartCX,doPsphO=.false.)
            !endif
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
        imW = 0.0_rp
        isEdible = .false.

        d_cold = 0.0_rp
        t_cold = TINY
        d_hot = 0.0_rp
        p_hot = 0.0_rp

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
                call InterpShellVar_TSC_pnt(sh, State%Den_avg(s), th, ph, d_cold)
                imW(IM_D_COLD) = d_cold  ! [#/cc]
                t_cold = PsphTemp_Genestreti(d_cold)  ! [keV]
                imW(IM_P_COLD) = DkT2P(d_cold, t_cold)  ! [nPa]
            elseif (spcList(s)%spcType == RAIJUHPLUS) then
                call InterpShellVar_TSC_pnt(sh, State%Den_avg(s)  , th, ph, d_hot)
                call InterpShellVar_TSC_pnt(sh, State%Press_avg(s), th, ph, p_hot)
                imW(IM_D_RING) = imW(IM_D_RING) + d_hot  ! [nPa]
                imW(IM_P_RING) = imW(IM_P_RING) + p_hot  ! [#/cc]
            elseif (spcList(s)%spcType == RAIJUELE) then
                ! Don't add number density
                call InterpShellVar_TSC_pnt(sh, State%Press_avg(s), th, ph, p_hot)
                imW(IM_P_RING) = imW(IM_P_RING) + p_hot  ! [#/cc]
            endif
        enddo

        
        call InterpShellVar_TSC_pnt(sh, State%Tb, th, ph, tScl)
        tScl = Model%nBounce*tScl  ! [s]
        !tScl = 10.0_rp  ! [s]

        ! 1/(x)^4 for x from 1 to 0.5 goes from 1 to 16. Higher exponents means stronger ramp-up
        !tScl = 15.0_rp/(App%vaFrac%data(i0,j0))**4  ! [s]

        ! Adjust IM_TSCL if we wanna ramp up over time
        if (t < App%startup_blendTscl) then
            rampC = RampDown(t, 0.0_rp, App%startup_blendTscl)
            !tScl = sqrt(tScl*App%startup_blendTscl)*rampC + (1-rampC)*tScl  ! idk
            tScl = rampC*30.0_rp*tScl + (1-rampC)*tScl  ! No good reason for 30 except for wanting starting tScl to be ~8-10 minutes
            !if (th > 50*deg2rad .and. th < 55*deg2rad .and. ph > 130*deg2rad .and. ph < 150*deg2rad) then
            !    write(*,*)"--",t,App%startup_blendTscl,rampC,tScl
            !endif
        endif
        
        imW(IM_TSCL) = tScl

        end associate
    end subroutine getMomentsRAIJU


    module subroutine getMomentsPrecipRAIJU(App,th,ph,imP,isEdible)
        !! Get precipitation quantities from RAIJU, formatted by enum in volttypes.F90 for a given theta,phi location within the RAIJU domain
        class(raijuCoupler_T), intent(inout) :: App
        real(rp), intent(in) :: th
            !! Theta [rad]
        real(rp), intent(in) :: ph
            !! Phi [rad]
        real(rp), intent(out) :: imP(nVars_imag2mix)
        logical , intent(out) :: isEdible

        integer :: s  ! Iterators
        integer :: i0, j0  ! i,j cell that provided th,ph are in
        real(rp) :: tScl, rampC

        integer :: i,j,k
        integer :: is, ie, js, je, ks, ke
        integer :: idx_proton
        real(rp) :: d_cold, d_hot, p_hot, dn_flux, de_flux
        real(rp) :: d_ion, p_ion, pie_frac, pe_kT, pe_eflux, pe_nflux

        associate(Model=>App%raiApp%Model, State=>App%raiApp%State, sh=>App%raiApp%Grid%shGrid, spcList=>App%raiApp%Grid%spc)
            ! Default
            imP = 0.0_rp
            isEdible = .false.

            idx_proton = spcIdx(App%raiApp%Grid, F_HOTP)

            d_cold = 0.0_rp
            d_hot = 0.0_rp
            p_hot = 0.0_rp
            dn_flux = 0.0_rp
            de_flux = 0.0_rp
    
            ! Is this a good point?
            if (th < sh%minTheta .or. th > sh%maxTheta) then
                return ! Off grid, return default
            endif
    
            ! Active check
            call getSGCellILoc(sh, th, i0)
            call getSGCellJLoc(sh, ph, j0)

            ! record the grid type info
            ! gtype is a proxy of raiju active, buffer, inactive
            imP(RAI_GTYPE) = (State%active(i0,j0)*1.0_rp+1.0)/2.0  ! "-1=Inactive, 0=Buffer, 1=Active" -> 0, 0.5, and 1.0        
            imP(RAI_THCON) = State%thcon(i0,j0) ! conjugate co-lat in radians, 0-pi
            imP(RAI_PHCON) = State%phcon(i0,j0) ! conjugate long in radians, 0-2pi
            if (State%active(i0,j0) .eq. RAIJUINACTIVE) then
                return
            endif
    
            ! Otherwise we are good, gonna return stuff
            isEdible = .true.

            do s=1, Model%nSpc
                ks = spcList(s)%kStart
                ke = spcList(s)%kEnd
                if (spcList(s)%flav == F_PSPH) then
                    !call InterpShellVar_TSC_pnt(sh, State%Den(s), th, ph, d_cold)
                    !imW(IM_D_COLD) = d_cold  ! [#/cc]
                    call InterpShellVar_TSC_pnt(sh, State%Den_avg(s), th, ph, d_cold)
                    imP(RAI_NPSP) = imP(RAI_NPSP) + d_cold*1.0e6 ! uStr="#/cc" -> /m^3 , dStr="Density from RAIJU flavors"
                elseif (spcList(s)%spcType == RAIJUHPLUS) then
                    ! add proton precipitation later.
                elseif (spcList(s)%spcType == RAIJUELE) then

                    call InterpShellVar_TSC_pnt(sh, State%Den_avg(s)  , th, ph, d_hot)  ! electron density (will ignore)
                    call InterpShellVar_TSC_pnt(sh, State%Press_avg(s), th, ph, p_hot)  ! Actual electron pressure

                    
                    imP(RAI_EDEN ) = imP(RAI_EDEN ) + d_hot*1.0e6 ! uStr="#/cc" -> /m^3 , dStr="Density from RAIJU flavors"
                    imP(RAI_EPRE ) = imP(RAI_EPRE ) + p_hot*1.0e-9 ! uStr="nPa" -> Pa , dStr="Pressure from RAIJU flavors"
                    do k=ks,ke
                        !if(.not. isnan(State%precipNFlux(i0,j0,k))) then
                        if(.not. isnan(State%precipNFlux(k)%data(i0,j0))) then ! use mask?
                            call InterpShellVar_TSC_pnt(sh, State%precipNFlux(k), th, ph, dn_flux)
                            imP(RAI_ENFLX) = imP(RAI_ENFLX) + dn_flux ! uStr="#/cm^2/s"
                            call InterpShellVar_TSC_pnt(sh, State%precipEFlux(k), th, ph, de_flux)
                            imP(RAI_EFLUX) = imP(RAI_EFLUX) + de_flux ! uStr="erg/cm^2/s"
                        endif
                    enddo


                    ! Mock up cold electrons balancing hot protons and see if it produces meaningful flux
                    call InterpShellVar_TSC_pnt(sh, State%Den_avg(idx_proton)    , th, ph, d_ion)
                    call InterpShellVar_TSC_pnt(sh, State%Press_avg(idx_proton)  , th, ph, p_ion)
                    pie_frac = 0.05  ! Fraction of ion pressure contained by these neutralizing electrons
                    pe_kT = DP2kT(d_ion, p_ion*pie_frac)  ! [keV]
                    pe_nflux = imP(RAI_ENFLX)*d_ion/d_hot  ! Scale number flux to same loss processes except there were d_ion density instead of d_electron density
                    pe_eflux = (pe_kT*kev2erg)*pe_nflux  ! [erg/cm^2/s]
                    if (pe_eflux > imP(RAI_EFLUX)) then  ! Use in place of normal flux only if energy flux for these are greater than hot electron channel fluxes
                        imP(RAI_EFLUX) = pe_eflux
                        imP(RAI_ENFLX) = pe_nflux
                        imP(RAI_EDEN ) = d_ion*1.0e6 ! [#/m^3]
                        imP(RAI_EPRE ) = p_ion*pie_frac*1.0e-9  ! [Pa]
                    endif
                endif
            enddo
            ! derive mean energy where nflux is non-trivial.
            if (imP(RAI_ENFLX) > TINY) imP(RAI_EAVG) = imP(RAI_EFLUX)/imP(RAI_ENFLX) * erg2kev  ! Avg E [keV]
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