submodule (volttypes) volttypessub
    use volttypes
    use gamCouple

    implicit none

    contains

    ! implementations for functions from volt types, which would otherwise cause cyclic inheritance issues

    ! procedures for gamCoupler_T
    module subroutine gamCplInitIO(App, Xml)
        class(gamCoupler_T), intent(inout) :: App
        type(XML_Input_T), intent(inout) :: Xml

        ! initialize parent's IO
        call gamInitIO(App, Xml)

        ! then my own
        App%vh5File = trim(App%Model%RunID) // ".gamCpl.h5"

    end subroutine

    module subroutine gamCplWriteRestart(App)
        class(gamCoupler_T), intent(inout) :: App

        ! write parent's restart data
        call gamWriteRestart(App)

        ! then my own
        call writeGamCouplerRestart(App)

    end subroutine

    module subroutine gamCplReadRestart(App, resId, nRes)
        class(gamCoupler_T), intent(inout) :: App
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        ! read parent's restart data
        call gamReadRestart(App, resId, nRes)

        ! then my own
        call readGamCouplerRestart(App, resId, nRes)

    end subroutine

    module subroutine gamCplWriteConsoleOutput(App)
        class(gamCoupler_T), intent(inout) :: App

        real(rp) :: cpcp(2) = 0.0

        ! write parent's console info
        call gamWriteConsoleOutput(App)

        ! then my own
        if(App%Model%isLoud) then
            call getCPCP(App%mixOutput,cpcp)
             write (*, '(a,2f8.3,a)')             '      CPCP = ' , cpcp(NORTH), cpcp(SOUTH), ' [kV, N/S]'
        endif

    end subroutine

    module subroutine gamCplWriteFileOutput(App)
        class(gamCoupler_T), intent(inout) :: App

        ! write parent's file
        call gamWriteFileOutput(App)

        ! then my own
        call writeCouplerFileOutput(App)

    end subroutine

    module subroutine gamCplWriteSlimFileOutput(App)
        class(gamCoupler_T), intent(inout) :: App

        call gamCplWriteFileOutput(App)

    end subroutine

    module subroutine gamInitMhdCoupler(App, voltApp)
        class(gamCoupler_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: voltApp

        type(TimeSeries_T) :: tsMJD

        allocate(App%SrcNC(App%Grid%is:voltApp%chmp2mhd%iMax+1,App%Grid%js:App%Grid%je+1,App%Grid%ks:App%Grid%ke+1,1:NVARIMAG))

        ! over-ride some Gamera parameters with Voltron values
        App%Model%t = voltApp%time / App%Model%Units%gT0
        App%Model%tFin = voltApp%tFin / App%Model%Units%gT0

        tsMJD%wID = voltApp%tilt%wID
        call tsMJD%initTS("MJD",doLoudO=.false.)
        App%Model%MJD0 = tsMJD%evalAt(0.0_rp) !Evaluate at T=0

    end subroutine

    module subroutine gamUpdateMhdData(App, voltApp)
        class(gamCoupler_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: voltApp

        real(rp) :: stepDT

        ! update to DeepT time
        stepDT = voltApp%DeepT - voltApp%time

        call App%AdvanceModel(stepDT)

    end subroutine


end submodule
