module kaimagtypes
	use kaimagdefs
	use ioclock

	implicit none

	type kaimagGrid_T
		! Convention i is the theta/latitude direction and j is phi/longitude direction
		! Nip = # of physical cells
		! Ni =  Total number of cells namely Nip + 2*Ng
		! Ng = Number of ghost cells
		! Sizes of data structures
		! NB - Even though computation will be done on the r,theta,phi grid we compute xyz for display
		! Coords (x,y,z, r, theta, phi): (Ni+1,Nj+1)
		! Cell centered (xc,yc,zc,rc,thetac,phic) : (Ni,Nj)
		! Face centered (xi,yi,zi,ri,thetai,phii) : (Ni+1,Nj)
		! Volume : (Ni,Nj)
		! Edge Lengths : (Ni+1,Nj)
		integer :: Ni, Nj, Nip1, Njp1
		integer :: Nip, Njp 
		integer :: Ng

		!Local indices of active region
		integer :: is, ie, js, je
		!Local indices with ghost cells 
		integer :: isg, ieg, jsg, jeg 

		!Corner-centered coordiates 
		real(rp), dimension(:,:), allocatable :: x,y,z
		real(rp), dimension(:,:), allocatable :: r, theta, phi 

		!Cell-centered coordinates
		real(rp), dimension(:,:), allocatable :: xc,yc,zc 
		real(rp), dimension(:,:), allocatable :: rc, thetac, phic 

		!i (theta) face centered coordinates
		real(rp), dimension(:,:), allocatable :: xi,yi,zi  
		real(rp), dimension(:,:), allocatable :: ri, thetai, phii 

		!j (phi) face centered coordinates
		real(rp), dimension(:,:), allocatable :: xj,yj,zj 
		real(rp), dimension(:,:), allocatable :: rj, thetaj, phij 

		!Volume
		real(rp), dimension(:,:), allocatable :: volume

		!Edge Lengths
		real(rp), dimension(:,:), allocatable :: di,dj 

		!Interface volumes for reconstruction 
		real(rp), dimension(:,:), allocatable :: vol_iface, vol_jface

		!Uniform
		logical :: uniform
		real(rp) :: dtheta, dphi
		real(rp) :: Rh


	end type kaimagGrid_T

	type kaimagModel_T
		! Run parameters
		character(len=strLen)  :: RunID, StencilID
		real(rp) :: t,tFin,dt,dtold 
		integer :: ts

		!Numerical paramters
		real(rp) :: CFL, PDMB

		!Output info
		type (IOClock_T) :: IO

		!Various problem type variables
		! Cosine Bell
		real(rp) :: t1 = 0.
		real(rp) :: t2 = 0.
		real(rp) :: p1 = 5.0*PI/6.0
		real(rp) :: p2 = 7.0*PI/6.0
		real(rp) :: h_max = 1.0

		! Cosine Bell Velocity 
		real(rp) :: Amp = 10.0
		real(rp) :: Peroid = 5.0 

		!Constant Phi Velocity
		real(rp) ::  Vphi0 = 1.

		!Constant Theta Velocity
		real(rp) ::  Vth0 =  -1.
 

	end type kaimagModel_T

	type kaimagState_T
		real(rp) :: time
		real(rp), dimension(:,:), allocatable :: rho,uc,vc
	end type kaimagState_T

	!StateIC_T
	!Generic initialization routine: ICs, Grid, Model
	abstract interface
		subroutine kaimagStateIC_T(Model,Grid,State,inpXML)
			Import :: kaimagModel_T, kaimagGrid_T, kaimagState_T, XML_Input_T
			type(kaimagModel_T), intent(inout) :: Model
			type(kaimagGrid_T), intent(in) :: Grid
			type(kaimagState_T), intent(inout) :: State 
			type(XML_Input_T), intent(in) :: inpXML
		end subroutine kaimagStateIC_T
	end interface

	!updateV_T
	!Generic Function for updating Velocities
	abstract interface
		subroutine kaimagupdateV_T(Model,Grid,State)
			Import :: kaimagModel_T, kaimagGrid_T, kaimagState_T
			type(kaimagModel_T), intent(in) :: Model
			type(kaimagGrid_T), intent(in) :: Grid
			type(kaimagState_T), intent(inout) :: State 
		end subroutine kaimagupdateV_T
	end interface
	
end module kaimagtypes
