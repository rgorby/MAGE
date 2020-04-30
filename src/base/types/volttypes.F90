!Voltron types

module volttypes
    use kdefs
    use cmidefs
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
        enumerator :: IMDEN=1,IMX1,IMX2,IMTSCL,IMPR
    endenum
    integer, parameter :: NVARIMAG = 5

    ! data for remix -> gamera conversion
    type mix2Mhd_T
        real(rp), dimension(:,:,:,:,:), allocatable :: mixOutput
        real(rp), dimension(:,:,:), allocatable :: gPsi
        type(Map_T), allocatable, dimension(:) :: PsiMaps
        integer :: PsiStart = PsiSt, PsiShells = PsiSh !Coming from cmidefs
        real(rp) :: rm2g
    end type mix2Mhd_T

    ! data for imag => remix for conductance
    type imag2Mix_T
        !Assuming IMag data is coming on northern hemisphere
        logical :: isFresh = .false.
        logical :: isInit = .false.
        !Working on logically Cartesian grid w/ lat, lon 
        real(rp), dimension(:), allocatable :: gcolat,glong
        real(rp), dimension(:,:), allocatable :: eflux,eavg
        real(rp), dimension(:,:), allocatable :: iflux,iavg
        !latc/lonc are the mappings to the southern hemisphere
        real(rp), dimension(:,:), allocatable :: latc,lonc
        real(rp), dimension(:,:), allocatable :: fac
        
        logical, dimension(:,:), allocatable :: isClosed
    end type imag2Mix_T

    ! data for gamera -> remix conversion
    type mhd2Mix_T
        real(rp), dimension(:,:,:,:,:), allocatable :: mixInput
        real(rp), dimension(:,:,:,:), allocatable :: gJ
        type(Map_T), allocatable, dimension(:) :: Jmaps
        integer :: JStart = JpSt, JShells = JpSh !Coming from cmidefs
        type(mixGrid_T) :: mixGfpd
    end type mhd2mix_T

    ! data for chimp -> gamera conversion
    type chmp2Mhd_T
        !Ni,Nj,Nk,2 array
        !Holds mapping from cell-centered xyz => x1,x1 (projection coordinates)
        !Projection coordinates can be R,phi (cylindrical) or lat,lon (ionospheric)
        real(rp), dimension(:,:,:,:), allocatable :: xyzSquish
        logical , dimension(:,:,:)  , allocatable :: isGood !Good projection or not
        integer :: iMax !Possibly changing i-boundary of squish mapping
    end type chmp2Mhd_T

    ! data for gamera -> chimp conversion
    type mhd2Chmp_T
        logical :: isLonely = .true.
        real(rp) :: Rin
        real(rp) :: lowlatBC

    end type mhd2Chmp_T

    type innerMagBase_T
        contains

        procedure baseInit
        procedure baseAdvance
        procedure baseEval
        procedure baseIO
        procedure baseRestart

        ! functions to be over-written by specific inner magnetosphere implementations
        procedure :: doInit => baseinit
        procedure :: doAdvance => baseAdvance
        procedure :: doEval => baseEval
        procedure :: doIO => baseIO
        procedure :: doRestart => baseRestart

    end type innerMagBase_T

    type voltApp_T

        !Voltron state information
        type(TimeSeries_T) :: tilt
        real(rp) :: time, MJD,tFin
        integer :: ts
        logical :: isSeparate = .false. ! whether Voltron is running in a separate application from gamera

        !Voltron output/restart info
        type (IOClock_T) :: IO
        logical :: isLoud = .true. !Console output
        !Apps
        type(mixApp_T) :: remixApp
        type(mhd2Mix_T) :: mhd2mix
        type(mix2Mhd_T) :: mix2mhd

        type(ebTrcApp_T)  :: ebTrcApp
        type(mhd2Chmp_T)  :: mhd2chmp
        type(chmp2Mhd_T)  :: chmp2mhd
        type(imag2Mix_T)  :: imag2mix

        class(innerMagBase_T), allocatable :: imagApp

        !Shallow coupling information
        real(rp) :: ShallowT
        real(rp) :: ShallowDT

        !Deep coupling information
        real(rp) :: DeepT
        real(rp) :: DeepDT
        logical  :: doDeep = .false. !Whether to do deep coupling
        real(rp) :: rDeep !Radius (in code units) to do deep coupling
        real(rp) :: rTrc  !Radius to do tracing (ebSquish) inside of
        integer  :: iDeep  = 0 !Index of max i shell containing deep coupling radius
        integer  :: imType = 0 !Type of inner magnetosphere model (0 = None)
        integer  :: prType = 0 !Type of projection for coupling   (0 = None)
        logical  :: doQkSquish = .false. !Whether or not to do fast squishing
    end type voltApp_T

    contains

    ! null default subroutines for inner mag base type
    subroutine baseInit(imag,iXML,isRestart,vApp)
        class(innerMagBase_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        type(voltApp_T), intent(inout) :: vApp
    end subroutine

    subroutine baseAdvance(imag,vApp,tAdv)
        class(innerMagBase_T), intent(inout) :: imag
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv
    end subroutine

    subroutine baseEval(imag,x1,x2,x12c,t,imW)
        class(innerMagBase_T), intent(inout) :: imag
        real(rp), intent(in) :: x1,x2,t
        real(rp), intent(in) :: x12C(2,2,2,2)
        real(rp), intent(out) :: imW(NVARIMAG)
        imW = 0.0_rp
    end subroutine

    subroutine baseIO(imag,nOut,MJD,time)
        class(innerMagBase_T), intent(inout) :: imag
        integer, intent(in) :: nOut
        real(rp), intent(in) :: MJD,time
    end subroutine

    subroutine baseRestart(imag,nRes,MJD,time)
        class(innerMagBase_T), intent(inout) :: imag
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time
    end subroutine

end module volttypes

