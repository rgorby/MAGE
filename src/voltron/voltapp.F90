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

    subroutine initVoltron(voltronApp, gameraApp)
        type(gamApp_T), intent(inout) :: gameraApp
        type(voltApp_T), intent(inout) :: voltronApp

        voltronApp%fastShallowTime = 0.0_rp
        voltronApp%fastShallowDt = 0.1_rp

        call initializeRemixFromGamera(voltronApp%remixApp, gameraApp)

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

    subroutine initializeRemixFromGamera(remixApp, gameraApp)
        type(mixApp_T), intent(inout) :: remixApp
        type(gamApp_T), intent(inout) :: gameraApp

        real(rp), allocatable, dimension(:,:,:,:,:) :: mhdJGrid,mhdPsiGrid

        ! initialize REMIX
        call init_remix_grids(remixApp, gameraApp, mhdJGrid, mhdPsiGrid)

        call init_mix_mhd_interface(remixApp,mhdJGrid,mhdPsiGrid)

        deallocate(mhdJGrid,mhdPsiGrid)
    end subroutine initializeRemixFromGamera

    subroutine runRemix(remixApp, time)
        type(mixApp_T), intent(inout) :: remixApp
        real(rp), intent(in) :: time

        ! calls run_mix inside
        call mhd2mix(remixApp)

        ! get stuff from mix to gamera
        call mix2mhd(remixApp)

        ! output remix info
        call mix_mhd_output(remixApp%ion,remixApp%mixOutput,hmsphrs,time)

    end subroutine runRemix

end module voltapp

