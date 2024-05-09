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

    module subroutine gamCplWriteRestart(App, nRes)
        class(gamCoupler_T), intent(inout) :: App
        integer, intent(in) :: nRes

        ! write parent's restart data
        call gamWriteRestart(App, nRes)

        ! then my own
        call writeGamCouplerRestart(App, nRes)

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

    module subroutine gamCplWriteFileOutput(App, nStep)
        class(gamCoupler_T), intent(inout) :: App
        integer, intent(in) :: nStep

        ! write parent's file
        call gamWriteFileOutput(App, nStep)

        ! then my own
        call writeCouplerFileOutput(App, nStep)

    end subroutine

    module subroutine gamCplWriteSlimFileOutput(App, nStep)
        class(gamCoupler_T), intent(inout) :: App
        integer, intent(in) :: nStep

        call gamCplWriteFileOutput(App, nStep)

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

        call ioSync(voltApp%IO, App%Model%IO, 1.0_rp/App%Model%Units%gT0)
        App%Model%IO%nRes = voltApp%IO%nRes
        App%Model%IO%nOut = voltApp%IO%nOut

    end subroutine

    module subroutine gamStartUpdateMhdData(App, voltApp)
        class(gamCoupler_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: voltApp

        ! for local coupling this function doesn't do anything
        ! all of the work is in gamFinishUpdateMhdData

    end subroutine

    module subroutine gamFinishUpdateMhdData(App, voltApp)
        class(gamCoupler_T), intent(inout) :: App
        class(voltApp_T), intent(inout) :: voltApp

        real(rp) :: stepDT

        ! update to DeepT time
        stepDT = voltApp%DeepT - voltApp%time

        call App%AdvanceModel(stepDT)

    end subroutine

end submodule
