!Voltron types

module volttypes
    use kdefs
    use kronos
    use ioclock
    use mixtypes
    use ebtypes

    implicit none

    enum, bind(C)
        enumerator :: IMAGRCM=1,IMAGSST
    endenum

    !Projection types
    enum, bind(C) 
        !L-phi (equatorial), Lat-Lon (northern ionospheric)
        enumerator :: LPPROJ=1,LLPROJ
    endenum

    !Data for inner mag => gamera variables
    enum, bind(C)
        enumerator :: IMDEN=1,IMLSCL,IMTSCL,IMVAR1,IMPR
    endenum
    integer, parameter :: NVARIMAG = 5

    ! data for remix -> gamera conversion
    type mix2Mhd_T
        real(rp), dimension(:,:,:,:,:), allocatable :: mixOutput
        real(rp), dimension(:,:,:), allocatable :: gPsi
        type(Map_T), allocatable, dimension(:) :: PsiMaps
        integer :: PsiStart = -3, PsiShells = 5
        real(rp) :: rm2g
    end type mix2Mhd_T

    ! data for gamera -> remix conversion
    type mhd2Mix_T
        real(rp), dimension(:,:,:,:,:), allocatable :: mixInput
        real(rp), dimension(:,:,:,:), allocatable :: gJ
        type(Map_T), allocatable, dimension(:) :: Jmaps
        integer :: JStart = 2, JShells = 1
        type(mixGrid_T) :: mixGfpd
    end type mhd2mix_T

    ! data for chimp -> gamera conversion
    type chmp2Mhd_T
        !Ni,Nj,Nk,2 array
        !Holds mapping from cell-centered xyz => x1,x1 (projection coordinates)
        !Projection coordinates can be R,phi (cylindrical) or lat,lon (ionospheric)
        real(rp), dimension(:,:,:,:), allocatable :: xyzSquish
        integer :: iMax !Possibly changing i-boundary of squish mapping
    end type chmp2Mhd_T

    ! data for gamera -> chimp conversion
    type mhd2Chmp_T
        logical :: isLonely = .true.
    end type mhd2Chmp_T

    type voltApp_T

        !Voltron state information
        type(TimeSeries_T) :: tilt
        real(rp) :: time, MJD,tFin
        integer :: ts

        !Voltron output/restart info
        type (IOClock_T) :: IO

        !Apps
        type(mixApp_T) :: remixApp
        type(mhd2Mix_T) :: mhd2mix
        type(mix2Mhd_T) :: mix2mhd

        type(ebTrcApp_T)  :: ebTrcApp
        type(mhd2Chmp_T)  :: mhd2chmp
        type(chmp2Mhd_T)  :: chmp2mhd

        !Shallow coupling information
        real(rp) :: ShallowT
        real(rp) :: ShallowDT

        !Deep coupling information
        real(rp) :: DeepT
        real(rp) :: DeepDT
        logical  :: doDeep = .false. !Whether to do deep coupling
        real(rp) :: rDeep !Radius (in code units) to do deep coupling
        integer  :: iDeep  = 0 !Index of max i shell containing deep coupling radius
        integer  :: imType = 0 !Type of inner magnetosphere model (0 = None)
        integer  :: prType = 0 !Type of projection for coupling   (0 = None)
    end type voltApp_T

end module volttypes