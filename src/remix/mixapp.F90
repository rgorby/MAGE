! Main data objects and functions to perform a remix simulation

module remixapp
    use mixtypes
    use gamapp

    implicit none

    integer,parameter :: hmsphrs(2) = [NORTH, SOUTH]

    type remixApp_T
        type(mixIon_T),dimension(2) :: ion  ! two ionosphere objects: north(1) and south(2)
        real(rp), dimension(:,:,:,:,:) :: mixInput,mixOutput
        real(rp) :: tilt
    end type remixApp_T

    contains

    subroutine initRemix(remixApp, gameraApp)
        type(remixApp_T), intent(inout) :: remixApp
        type(gameraApp_T), intent(in) :: gameraApp
    end subroutine initRemix

    subroutine runRemix(remixApp)
        type(remixApp_T), intent(inout) :: remixApp

        ! calls run_mix inside
        call mhd2mix(remixApp%ion,remixApp%remixInput,remixApp%tilt,remixApp%conductance)

        ! get stuff from mix to gamera
        call mix2mhd(remixApp%ion,remixApp%remixOutput)

    end subroutine stepRemix

end module remixapp

