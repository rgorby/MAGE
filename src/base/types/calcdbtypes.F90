!Main types for calculating ground db

module calcdbtypes
	use kdefs
	use chmpdefs
	use ebtypes
	implicit none

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
        logical :: doPed,doHall
    end type rmState_T

	INTEGER, parameter :: rSegs=30

!--------
!Native source grids (SM coords)
!Data for SM-XYZ locations/J vectors to contribute to ground delta-B
!facGrid_T, ionGrid_T, and ebGr_T

	!Grid for FAC calculation
	    !Np,Nth are phi/theta grid cells
	    !Nh = 1/2 for hemisphere
	    !Nr = radial segments
	    !Arrays are size Np,Nth,Nh=2,Nr x NDIM (for vectors)
    type facGrid_T
        integer :: Np,Nth,rSegs
    	real(rp), dimension(:,:,:,:,:), allocatable :: XYZcc !XYZ of segment centers
        real(rp), dimension(:,:,:,:,:), allocatable :: Jxyz !Jxyz at segment centers
    	real(rp), dimension(:,:,:,:)  , allocatable :: dV !Differential volume of each segment
        real(rp) :: dt,dp
        real(rp), dimension(:,:,:), allocatable :: pcc,tcc !phi/theta cell-centered

    end type facGrid_T

	!Grid for ION calculation (Hall+Pederson currents)
	    !Np,Nth are phi/theta grid cells
	    !Nh = 1/2 for hemisphere
	    !Arrays are size Np,Nth,Nh=2 (for vectors)
    type ionGrid_T
        integer :: Np,Nth
        real(rp), dimension(:,:,:,:), allocatable :: XYZcc !XYZ of patch centers
        real(rp), dimension(:,:,:,:), allocatable :: Jxyz !Jxyz at patch centers
        real(rp), dimension(:,:,:)  , allocatable :: dS !Differential surface area of patch
        real(rp) :: dt,dp
        real(rp), dimension(:,:,:), allocatable :: pcc,tcc !phi/theta cell-centered
        real(rp), dimension(:,:,:,:), allocatable :: Etp !E field, theta/phi components
        real(rp), dimension(:,:,:,:), allocatable :: hJ  !Hall current, theta/phi components
        real(rp), dimension(:,:,:,:), allocatable :: pJ  !Pede current, theta/phi components
    end type ionGrid_T

!--------
!Biot-Savart grid
!1D collection of points, currents, geometric factors to do discrete BS integral
!Coordinate system (points and vectors) matches that of the ground grid, likely GEO
	type BSGrid_T
		integer :: NumP !Number of points
		!shape XYZcc,Jxyz = NumP,NDIM
		!shape dV = NumP
		real(rp), dimension(:,:), allocatable :: XYZcc !XYZ of centers
		real(rp), dimension(:,:), allocatable :: Jxyz !Vector currents
		real(rp), dimension(:)  , allocatable :: dV !Volume element
        real(rp) :: jScl = 1.0 !Scaling term for current contribution to get nT from BS

        logical :: isActive = .true. !Whether this grid should be used
	end type BSGrid_T

!--------
!Destination (ground) grid, thin shell (lat/lon/height)

	type grGrid_T
		integer :: NLat,NLon,Nz
		real(rp), dimension(:,:,:,:), allocatable :: GxyzI,GxyzC !Corner/Center points, ground coordinates
        real(rp), dimension(:,:,:,:), allocatable :: SMxyzC !Center points, SM
        logical :: doGEO = .true. !Do GEO coordinates on ground
        !Ground geomagnetic information (possibly redundant, meh)
        real(rp), dimension(:,:,:), allocatable :: smlat,smlon,dBn

        !Individual ground perturbations
        real(rp), dimension(:,:,:,:), allocatable :: dbMAG_xyz,dbION_xyz,dbFAC_xyz
        real(rp), dimension(:,:,:,:), allocatable :: dbMAG_rtp,dbION_rtp,dbFAC_rtp

        !Some parameters for calculation
        real(rp) :: rMax=30.0 !Radius of magnetospheric ball to integrate over [Re]

        !Indices
        real(rp) :: SML,SML_MLat,SML_MLon
        real(rp) :: SMU,SMU_MLat,SMU_MLon
        real(rp) :: SMR_00,SMR_06,SMR_12,SMR_18
        real(rp) :: SML_00,SML_06,SML_12,SML_18
        real(rp) :: SMU_00,SMU_06,SMU_12,SMU_18

        !Derived indices
        real(rp) :: SMR,SME,SMO

        

	end type grGrid_T

end module calcdbtypes
