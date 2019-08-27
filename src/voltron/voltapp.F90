! Collection of data and objects for the voltron middle man

module voltapp
    use mixtypes
    use mhd2mix_interface
    use mix2mhd_interface

    implicit none

    type voltApp_T
        type(mixApp_T) :: remixApp
        type(mhd2Mix_T) :: mhd2mix
        type(mix2Mhd_T) :: mix2mhd

        real(rp) :: fastShallowTime
        real(rp) :: fastShallowDt

        real(rp) :: tilt

    end type voltApp_T

    contains

    subroutine initVoltron(voltronApp, gameraApp,optFilename)
        type(gamApp_T), intent(inout) :: gameraApp
        type(voltApp_T), intent(inout) :: voltronApp
        character(len=*), optional, intent(in) :: optFilename

        voltronApp%fastShallowTime = 0.0_rp
        voltronApp%fastShallowDt = 0.1_rp

        voltronApp%tilt = 0.0_rp

        if(present(optFilename)) then
            ! read from the prescribed file
            call initializeRemixFromGamera(voltronApp, gameraApp, optFilename)
        else
            call initializeRemixFromGamera(voltronApp, gameraApp)
        endif

    end subroutine initVoltron

    subroutine fastShallowUpdate(voltronApp, gameraApp, time)
        type(gamApp_T), intent(inout) :: gameraApp
        type(voltApp_T), intent(inout) :: voltronApp
        real(rp) :: time

        ! convert gamera data to mixInput
        call Tic("G2R")
        call convertGameraToRemix(voltronApp%mhd2mix, gameraApp, voltronApp%remixApp)
        call Toc("G2R")

        ! run remix
        call Tic("ReMIX")
        call runRemix(voltronApp, time)
        call Toc("ReMIX")

        ! convert mixOutput to gamera data
        call Tic("R2G")
        call convertRemixToGamera(voltronApp%mix2mhd, voltronApp%remixApp, gameraApp)
        call Toc("R2G")

        voltronApp%fastShallowTime = time + voltronApp%fastShallowDt

    end subroutine fastShallowUpdate

    subroutine initializeRemixFromGamera(voltronApp, gameraApp, optFilename)
        type(voltApp_T), intent(inout) :: voltronApp
        type(gamApp_T), intent(inout) :: gameraApp
        character(len=*), optional, intent(in) :: optFilename

        if(present(optFilename)) then
            ! read from the prescribed file
            call init_mix(voltronApp%remixApp%ion,[NORTH, SOUTH],optFilename)
        else
            call init_mix(voltronApp%remixApp%ion,[NORTH, SOUTH])
        endif

        call init_mhd2Mix(voltronApp%mhd2mix, gameraApp, voltronApp%remixApp)
        call init_mix2Mhd(voltronApp%mix2mhd, voltronApp%remixApp, gameraApp)

    end subroutine initializeRemixFromGamera

    subroutine runRemix(voltronApp, time)
        type(voltApp_T), intent(inout) :: voltronApp
        real(rp), intent(in) :: time

        ! convert gamera inputs to remix
        call mapGameraToRemix(voltronApp%mhd2mix, voltronApp%remixApp)

        ! solve for remix output
        call run_mix(voltronApp%remixApp%ion,voltronApp%tilt)

        ! get stuff from mix to gamera
        call mapRemixToGamera(voltronApp%mix2mhd, voltronApp%remixApp)

        ! output remix info
        call mix_mhd_output(voltronApp%remixApp%ion,voltronApp%mix2mhd%mixOutput,time)

    end subroutine runRemix

end module voltapp

