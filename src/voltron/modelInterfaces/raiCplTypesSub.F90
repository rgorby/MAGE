submodule (volttypes) raijuCplTypesSub

    use raijutypes
    use raijustarter
    use raijuCplHelper
    use raijuColdStartHelper

    use shellInterp
    use imaghelper

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
        call raijuCpl_init(App)

    end subroutine raiCplInitModel


    module subroutine volt2RAIJU(App, vApp)
        class(raijuCoupler_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: vApp

        logical :: doColdStart
        doColdStart = .false.

        associate(raiApp=>App%raiApp)

            ! If we are running realtime, its our job to do tracing and get all other stuff from vApp into raiCpl
            if (.not. App%raiApp%Model%isSA) then
                ! Determine if we should cold start
                ! e.g. Completely reset raiju's eta's to match some target conditions
                ! Determine if we should cold start before packing coupler because it will set tLastUpdate to vApp%time and then we can't do the check we want
                ! But actually do cold start after normal coupling completes so we can use real field line info
                if (App%opt%doColdStart .and. App%tLastUpdate < 0.0 .and. vApp%time >= 0.0) then
                    doColdStart = .true.
                endif

                call packRaijuCoupler_RT(App, vApp)
            endif

            ! Someone updated raiCpl's coupling variables by now, stuff it into RAIJU proper
            !call imagTubes2RAIJU(raiApp%Model, raiApp%Grid, raiApp%State, App%ijTubes)
            !raiApp%State%espot(:,:) = App%pot%data(:,:) ! They live on the same grid so this is okay
            call raiCpl2RAIJU(App)

            if (doColdStart) then
                ! Its happening, everybody stay calm
                write(*,*) "RAIJU Cold starting..."
                call raijuGeoColdStart(raiApp%Model, raiApp%Grid, raiApp%State, vApp%time, vApp%BSDst)
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

        ! IM_D_RING=1,IM_P_RING,IM_D_COLD, IM_P_COLD, IM_TSCL
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
        call getSGCellILoc(sh, ph, j0)
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

        call InterpShellVar_TSC_pnt(sh, State%Tb, th, ph, imW(IM_TSCL))
        imW(IM_TSCL) = Model%nBounce*imW(IM_TSCL)  ! [s]

        end associate
    end subroutine getMomentsRAIJU


    module subroutine getMomentsPrecipRAIJU(App,th,ph,t,imW,isEdible)
        class(raijuCoupler_T), intent(inout) :: App
        real(rp), intent(in) :: th,ph,t
        real(rp), intent(out) :: imW(NVARIMAG0)
        logical, intent(out) :: isEdible

        imW = 0.0
        isEdible = .false.
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

        call App%raiApp%WriteRestart(nRes)

    end subroutine

    module subroutine raiCplReadRestart(App, resId, nRes)
        class(raijuCoupler_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        call App%raiApp%ReadRestart(resId, nRes)

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