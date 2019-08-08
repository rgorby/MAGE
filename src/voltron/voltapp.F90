! Collection of data and objects for the voltron middle man

module voltapp
    use mixapp
    use mhd2mix_interface
    use mix2mhd_interface

    implicit none

    type voltApp_T
        type(mixApp_T) :: remixApp

        real(rp) :: fastShallowTime
        real(rp) :: fastShallowDt

    end type voltApp_T

    contains

    subroutine initVoltron(voltronApp, gameraApp,optFilename)
        type(gamApp_T), intent(inout) :: gameraApp
        type(voltApp_T), intent(inout) :: voltronApp
        character(len=*), optional, intent(in) :: optFilename

        voltronApp%fastShallowTime = 0.0_rp
        voltronApp%fastShallowDt = 0.1_rp

        if(present(optFilename)) then
            ! read from the prescribed file
            call initializeRemixFromGamera(voltronApp%remixApp, gameraApp, optFilename)
        else
            call initializeRemixFromGamera(voltronApp%remixApp, gameraApp)
        endif

    end subroutine initVoltron

    subroutine fastShallowUpdate(voltronApp, gameraApp, time)
        type(gamApp_T), intent(inout) :: gameraApp
        type(voltApp_T), intent(inout) :: voltronApp
        real(rp) :: time

        ! convert gamera data to mixInput
        call convertGameraToRemix(gameraApp, voltronApp%remixApp)

        ! run remix
        call runRemix(voltronApp%remixApp, time)

        ! convert mixOutput to gamera data
        call convertRemixToGamera(gameraApp, voltronApp%remixApp)

        voltronApp%fastShallowTime = time + voltronApp%fastShallowDt

    end subroutine fastShallowUpdate

    subroutine initializeRemixFromGamera(remixApp, gameraApp, optFilename)
        type(mixApp_T), intent(inout) :: remixApp
        type(gamApp_T), intent(inout) :: gameraApp
        character(len=*), optional, intent(in) :: optFilename

        real(rp), allocatable, dimension(:,:,:,:,:) :: mhdJGrid,mhdPsiGrid

        ! initialize REMIX
        call init_remix_grids(remixApp, gameraApp, mhdJGrid, mhdPsiGrid)

        if(present(optFilename)) then
            ! read from the prescribed file
            call init_mix_mhd_interface(remixApp,mhdJGrid,mhdPsiGrid,optFilename)
        else
            call init_mix_mhd_interface(remixApp,mhdJGrid,mhdPsiGrid)
        endif

        deallocate(mhdJGrid,mhdPsiGrid)
    end subroutine initializeRemixFromGamera

    subroutine runRemix(remixApp, time)
        type(mixApp_T), intent(inout) :: remixApp
        real(rp), intent(in) :: time

        ! convert gamera inputs to remix
        call mhd2mix(remixApp)

        ! solve for remix output
        call run_mix(remixApp%ion,remixApp%tilt,remixApp%conductance)

        ! get stuff from mix to gamera
        call mix2mhd(remixApp)

        ! output remix info
        call mix_mhd_output(remixApp%ion,remixApp%mixOutput,hmsphrs,time)

    end subroutine runRemix

end module voltapp

