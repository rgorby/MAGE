module calcdbutils
	use kdefs
	use chmpdefs
	use ebtypes
	implicit none

    real(rp), parameter :: RIon = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880

!Remix holders
    !Holds one hemisphere at one time slice
    type rmHemi_T
        integer :: nStp
        integer :: Np,Nth !Remix cap sizes
        real(rp) :: time
        !Data arrays are size Np,Nth (lon,lat)
        real(rp), dimension(:,:), allocatable :: xFac,xSigP,xSigH,xPot
    end type rmHemi_T

    !Remix state interpolated to specific time using bracketing steps rmHemi_T
    type rmState_T
        real(rp) :: time !CHIMP units
        real(rp), dimension(:,:,:), allocatable :: XY
        integer :: Np,Nth !Remix cap sizes
        integer :: i1=-1,i2=-1 !Bracketing step numbers
        !Data arrays are size Np,Nth (lon,lat)
        real(rp), dimension(:,:), allocatable :: nFac,nSigP,nSigH,nPot
        real(rp), dimension(:,:), allocatable :: sFac,sSigP,sSigH,sPot
        
        !Holds N/S hemisphere at i1/i2
        type(rmHemi_T) :: rmN1,rmN2,rmS1,rmS2

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

        !Individual ground perturbations
        real(rp), dimension(:,:,:,:), allocatable :: dbMAG_xyz,dbION_xyz,dbFAC_xyz
        real(rp), dimension(:,:,:,:), allocatable :: dbMAG_rtp,dbION_rtp,dbFAC_rtp

        !Some parameters for calculation
        integer :: i0=1 !Shell to start at
        real(rp) :: rMax=25.0 !Radius of magnetospheric ball to integrate over [Re]

	end type sphGrid_T

	contains

	subroutine facGridInit(Model,ebState,rmState,facGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(facGrid_T)  , intent(inout) :: facGrid

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

        !CALCDB-TODO: Need to add code here to define facGrid%XYZcc and dV
        !NOTE: XYZ from remix is in units of Rion, not Re

	end subroutine facGridInit

    subroutine ionGridInit(Model,ebState,rmState,ionGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(ionGrid_T)  , intent(inout) :: ionGrid

        integer :: Np,Nth

        Np = rmState%Np
        Nth = rmState%Nth

        allocate(ionGrid%XYZcc(Np,Nth,2,NDIM))
        allocate(ionGrid%Jxyz (Np,Nth,2,NDIM))
        allocate(ionGrid%dS   (Np,Nth,2)) !Surface area per patch

        ionGrid%XYZcc = 0.0
        ionGrid%Jxyz  = 0.0
        ionGrid%dS    = 0.0

        !CALCDB-TODO: Need to add code here to define grid (XYZcc and dS)
        !NOTE: XYZ from remix is in units of Rion, not Re

    end subroutine ionGridInit

	!Using a rmState (remix data), fill facGrid Jxyz
	subroutine facGridUpdate(Model,ebState,rmState,facGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(facGrid_T)  , intent(inout) :: facGrid

        !CALCDB-TODO: Write this
        facGrid%Jxyz = 0.0

    end subroutine facGridUpdate

    !Using a rmState (remix data), fill ionGrid Jxyz
    subroutine ionGridUpdate(Model,ebState,rmState,ionGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(ionGrid_T)  , intent(inout) :: ionGrid

        !CALCDB-TODO: Write this
        ionGrid%Jxyz = 0.0

    end subroutine ionGridUpdate


    !Set rmState given properly set 4 hemispheres and temporal weights
    subroutine hemi2rm(rmState,w1,w2)
        type(rmState_T)  , intent(inout) :: rmState
        real(rp), intent(in) :: w1,w2

        rmState%nFac  = w1*rmState%rmN1%xFac  + w2*rmState%rmN2%xFac 
        rmState%nPot  = w1*rmState%rmN1%xPot  + w2*rmState%rmN2%xPot 
        rmState%nSigP = w1*rmState%rmN1%xSigP + w2*rmState%rmN2%xSigP
        rmState%nSigH = w1*rmState%rmN1%xSigH + w2*rmState%rmN2%xSigH
        rmState%sFac  = w1*rmState%rmS1%xFac  + w2*rmState%rmS2%xFac 
        rmState%sPot  = w1*rmState%rmS1%xPot  + w2*rmState%rmS2%xPot 
        rmState%sSigP = w1*rmState%rmS1%xSigP + w2*rmState%rmS2%xSigP
        rmState%sSigH = w1*rmState%rmS1%xSigH + w2*rmState%rmS2%xSigH


    end subroutine hemi2rm
end module calcdbutils
