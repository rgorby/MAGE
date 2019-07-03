! Main data objects and functions to perform a remix simulation

module mixapp
    use mixtypes
    use mixconductance
    use msphutils, ONLY: Rion ! Rion should be moved into gamapp

    implicit none

    integer,parameter :: hmsphrs(2) = [NORTH, SOUTH]

    integer, parameter :: mhd2mix_varn = 3  ! how many variables are we sending (should be consistent with the enumerator in mixdefs.F90)
    integer, parameter :: mix2mhd_varn = 1  ! for now just the potential is sent back

    integer, parameter :: MAXMIXIOVAR = 10

    type mixApp_T
        type(mixIon_T),dimension(2) :: ion  ! two ionosphere objects: north(1) and south(2)
        real(rp), dimension(:,:,:,:,:), allocatable :: mixInput,mixOutput
        real(rp) :: tilt
        type(mixConductance_T) :: conductance

        real(rp), dimension(:,:,:,:), allocatable :: gJ ! temporary storage

        real(rp), dimension(:,:,:), allocatable :: gPsi ! output from REmix

        ! Gamera normalization
        ! Scaling factor for remix potential [kV]
        real(rp) :: rm2g

        type(Map_T), allocatable, dimension(:) :: PsiMaps,Jmaps
        type(mixGrid_T) :: mixGfpd

        integer :: PsiShells = 5
        integer :: JShells = 1

        integer :: PsiStart = -3
        integer :: JStart = 2

    end type mixApp_T

    contains

end module mixapp

