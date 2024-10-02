submodule (volttypes) raijuCplTypesSub

    use raijutypes
    use raijuCplHelper
    use raijuColdStartHelper

    implicit none

    contains

    module subroutine raiCplInitModel(App, xml)
        class(raijuCoupler_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        ! Allocate our contained raiju app
        allocate(raijuApp_T :: App%raiApp)
        ! Init raiju app itself
        call App%raiApp%InitModel(xml)
        ! Update MJD with whatever voltron handed us
        ! If we are restarting, this will get replaced with whatever's in file later
        App%raiApp%State%mjd = App%opt%mjd0
        write(*,*)"MJD0=",App%opt%mjd0
        if (App%opt%doColdStart) then
            ! We are gonna cold start, so ignore plasma ingestion rules for first coupling
            Ap%raiApp%State%isFirstCpl = .false.
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
            call imagTubes2RAIJU(raiApp%Model, raiApp%Grid, raiApp%State, App%ijTubes)
            raiApp%State%espot(:,:) = App%pot%data(:,:) ! They live on the same grid so this is okay

            if (doColdStart) then
                ! Its happening, everybody stay calm
                write(*,*) "RAIJU Cold starting..."
                call raijuGeoColdStart(raiApp%Model, raiApp%Grid, raiApp%State, vApp%time, vApp%BSDst)
            endif
        end associate
    end subroutine volt2RAIJU


    module subroutine evalRAIJU(App,x1,x2,t,imW,isEdible)
        class(raijuCoupler_T), intent(inout) :: App
        real(rp), intent(in) :: x1,x2,t
        real(rp), intent(out) :: imW(NVARIMAG)
        logical, intent(out) :: isEdible

        imW = 0.0
        isEdible = .false.
    end subroutine evalRAIJU


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