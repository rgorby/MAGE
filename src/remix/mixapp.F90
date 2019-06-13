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

        real(rp), allocatable, dimension(:,:,:,:,:) :: mhdJGrid,mhdPsiGrid

        ! initialize REMIX
        call init_remix_grids(mhdJGrid, mhdPsiGrid, gamGridXyz, Rion)

        call init_mix_mhd_interface(ion,hmsphrs,mhdJGrid,mhdPsiGrid)

        deallocate(mhdJGrid,mhdPsiGrid)

    end subroutine initRemix

    subroutine runRemix(remixApp)
        type(remixApp_T), intent(inout) :: remixApp

        ! calls run_mix inside
        call mhd2mix(remixApp%ion,remixApp%remixInput,remixApp%tilt,remixApp%conductance)

        ! get stuff from mix to gamera
        call mix2mhd(remixApp%ion,remixApp%remixOutput)

        ! output remix info
        call mix_mhd_output(ion,remixOutputs,hmsphrs,time)

    end subroutine stepRemix

end module remixapp

