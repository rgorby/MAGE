

module sifgrids
	use xml_input
	use sifdefs
	use siftypes
	use shellgrid


	implicit none

	contains



	subroutine finalizeLLGrid(Grid)
		!! Use a fully-created shell grid to allocate and populate the rest of the grid parameters
		type(sifGrid_T), intent(inout) :: Grid


        associate(shGr=>Grid%shGrid)

        end associate

	end subroutine finalizeLLGrid

	subroutine sifGenUniSphGrid(Grid, iXML)
		type(sifGrid_T)  , intent(inout) :: Grid
		type(XML_Input_T), intent(in)    :: iXML

		real(rp), dimension(:), allocatable :: Theta
	    real(rp), dimension(:), allocatable :: Phi
	    real(rp) :: dTheta, dPhi, thetaL, thetaU
	    integer :: Nt,Np,Ng,Ngn,Ngs,Nge,Ngw
	    integer :: i


	    call iXML%Set_Val(thetaL , "grid/ThetaL", 50.0)  
	    	!! Lower lat boundary [deg]
	    call iXML%Set_Val(thetaU , "grid/ThetaU", 80.0)
	    	!! Upper lat boundary [deg]
		call iXML%Set_Val(Nt, "grid/Nt", 30 )  ! 1 deg resolution
		call iXML%Set_Val(Np, "grid/Np", 360)  ! 1 deg resolution
		call iXML%Set_Val(Ng, "grid/Ng", 4  )  ! Number of ghosts, in every direction for now

		! Turn degrees into radians
		thetaL = thetaL*deg2rad
		thetaU = thetaU*deg2rad
		
		! Probably need to change
		Ngn = Ng
	    Ngs = Ng
	    Nge = Ng
	    Ngw = Ng

		! Allocate arrays
		allocate(Theta(Nt))
		allocate(Phi(Np))

		! Create uniform grids
		dTheta = (thetaU-thetaL)/(Nt-1)
	    dPhi = 2*PI/(Np-1)

	    do i=1,Nt
	    	write(*,*)i
	        Theta(i) = thetaL + (i-1)*dTheta
	    enddo

	    do i=1,Np
	        Phi(i) = (i-1)*dPhi
	    enddo

	    write(*,*)"Orig theta:",Theta
	    write(*,*)"Orig phi:",Phi

	    call GenShellGrid(Grid%shGrid,Theta,Phi,Ngn,Ngs,Nge,Ngw)

	end subroutine

	subroutine sifGenGridFromShGrid(Grid, shGrid)
		type(sifGrid_T)  , intent(inout) :: Grid
		type(ShellGrid_T), intent(inout) :: shGrid

	end subroutine sifGenGridFromShGrid


end module sifgrids