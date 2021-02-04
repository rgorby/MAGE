module calcdbutils
	use kdefs
	use chmpdefs
	use ebtypes
	implicit none

    real(rp), parameter :: RIon = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880

!Remix holders
    type rmState_T
        real(rp) :: time !CHIMP units
        real(rp), dimension(:,:,:), allocatable :: XY
        integer :: Np,Nth !Remix cap sizes
        integer :: i1=-1,i2=-1 !Bracketing step numbers
        !Data arrays are size Np,Nth (lon,lat)
        real(rp), dimension(:,:), allocatable :: nFac,nSigP,nSigH,nPot
        real(rp), dimension(:,:), allocatable :: sFac,sSigP,sSigH,sPot
    
    end type rmState_T

	INTEGER, parameter :: NORTH=1,SOUTH=2
	INTEGER, parameter :: rSegs=10

!Grid for FAC calculation
    !Np,Nth are phi/theta grid cells
    !Nh = 1/2 for hemisphere
    !Nr = radial segments
    !Arrays are size Np,Nth,Nh=2,Nr x NDIM (for vectors)
    type facGrid_T
    	real(rp), dimension(:,:,:,:,:), allocatable :: XYZcc !XYZ of segment centers
        real(rp), dimension(:,:,:,:,:), allocatable :: Jxyz !Jxyz at segment centers
    	real(rp), dimension(:,:,:,:)  , allocatable :: dV !Differential volume of each segment
    end type facGrid_T

!Grid for ION calculation (Hall+Pederson currents)
    !Np,Nth are phi/theta grid cells
    !Nh = 1/2 for hemisphere
    !Arrays are size Np,Nth,Nh=2 (for vectors)
    type ionGrid_T
        real(rp), dimension(:,:,:,:), allocatable :: XYZcc !XYZ of segment centers
        real(rp), dimension(:,:,:,:), allocatable :: Jxyz !Jxyz at segment centers
        real(rp), dimension(:,:,:)  , allocatable :: dS !Differential surface area of patch

    end type ionGrid_T

!Destination (ground) grid, thin shell (lat/lon/height)
	real(rp) :: dzGG = 60.0 !Default height spacing [km]

	type sphGrid_T
		integer :: NLat,NLon,Nz
		real(rp), dimension(:,:,:,:), allocatable :: xyzI,xyzC !Corner/Center points

	end type sphGrid_T

	type(sphGrid_T) :: gGr !Ground grid
	type(facGrid_T) :: facGrid !FAC grid
    type(ionGrid_T) :: ionGrid !Ionospheric grid

	contains

	subroutine facGridInit(Model,ebState,rmState)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState

        integer :: Np,Nth

        Np = rmState%Np
        Nth = rmState%Nth

        !Allocate and zero out facGrid

        allocate(facGrid%XYZcc(Np,Nth,2,rSegs,NDIM))
        allocate(facGrid%Jxyz (Np,Nth,2,rSegs,NDIM))
        allocate(facGrid%dV   (Np,Nth,2,rSegs))
        facGrid%XYZcc = 0.0
        facGrid%Jxyz  = 0.0
        facGrid%dV    = 0.0

        !TODO: Need to add code here to define facGrid%XYZcc and dV
        !NOTE: XYZ from remix is in units of Rion, not Re

	end subroutine facGridInit

	!Using a rmState (remix data), fill facGrid Jxyz
	subroutine facGridUpdate(Model,ebState,rmState)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState

        !TODO: Write this
        facGrid%Jxyz = 0.0

    end subroutine facGridUpdate

end module calcdbutils
