! Main data objects and functions to perform a remix simulation

module remixapp
    use mixtypes
    use gamapp
    use mixconductance

    implicit none

    integer,parameter :: hmsphrs(2) = [NORTH, SOUTH]

    type remixApp_T
        type(mixIon_T),dimension(2) :: ion  ! two ionosphere objects: north(1) and south(2)
        real(rp), dimension(:,:,:,:,:), allocatable :: mixInput,mixOutput
        real(rp) :: tilt
        type(mixConductance_T) :: conductance
    end type remixApp_T

    contains

    subroutine initRemix(remixApp, mhdJGrid, mhdPsiGrid)
        type(remixApp_T), intent(inout) :: remixApp
        real(rp), dimension(:,:,:,:,:), intent(in) :: mhdJGrid,mhdPsiGrid

        call init_mix_mhd_interface(remixApp%ion,hmsphrs,mhdJGrid,mhdPsiGrid)

    end subroutine

    subroutine runRemix(remixApp, time)
        type(remixApp_T), intent(inout) :: remixApp
        real(rp), intent(in) :: time

        ! calls run_mix inside
        call mhd2mix(remixApp%ion,remixApp%mixInput,remixApp%tilt,remixApp%conductance)

        ! get stuff from mix to gamera
        call mix2mhd(remixApp%ion,remixApp%mixOutput)

        ! output remix info
        call mix_mhd_output(remixApp%ion,remixApp%mixOutput,hmsphrs,time)

    end subroutine runRemix

end module remixapp

