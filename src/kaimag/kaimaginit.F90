module kaimaginit
	use kdefs
	use kaimagtypes
	use kaimagsolver
	use kaimagprob
	use xml_input

	implicit none


	contains

	subroutine createUniformGrid(Model,Grid,xmlInp)
		type(kaimagModel_T), intent(in) :: Model
		type(KaimagGrid_T), intent(inout) :: Grid
		type(XML_Input_T), intent(in) :: xmlInp

		!Defaults that can be overidden by XML inputs
		integer, parameter :: Ni = 128
    	integer, parameter :: Nj = 248
    	integer, parameter :: NO = 8
    	integer, parameter :: NO2 = NO/2

    	real(rp) :: dtheta, dphi, Rh
    	real(rp), dimension(:), allocatable :: theta0, phi0

    	!Center weights for volume reconstruction
    	real(rp), dimension(8), parameter :: Cent8C = [-3,29,-139,533,533,-139,29,-3]/840.0_rp
    	! loop indices 
    	integer :: i,j

    	write(*,*) 'Especially those giving birth!'

    	call xmlInp%Set_Val(Grid%ni,'grid/nt',Ni)
    	call xmlInp%Set_Val(Grid%nj,'grid/np',Nj)
    	call xmlInp%Set_Val(Grid%ng,'grid/ng',NO2)
    	Grid%nip1 = Grid%ni + 1
    	Grid%njp1 = Grid%nj + 1

    	dtheta = 180./(Grid%Ni-Grid%ng*2)*deg2rad
    	dphi = 360./(Grid%Nj-Grid%ng*2)*deg2rad
    	Rh = 1. 

    	Grid%uniform = .true.
		Grid%dtheta = dtheta
		Grid%dphi = dphi
		Grid%Rh = Rh

    	allocate(theta0(Grid%nip1))
    	allocate(phi0(Grid%njp1))
    	

    	do i=1,Grid%nip1
    		theta0(i)=-Grid%ng*dtheta+(i-1)*dtheta
    	enddo 


    	do j=1,Grid%njp1
    		phi0(j)=-Grid%ng*dphi+(j-1)*dphi
    	enddo 


    	!allocate corner gird

    	allocate(Grid%x(Grid%nip1,Grid%njp1))
    	allocate(Grid%y(Grid%nip1,Grid%njp1))
    	allocate(Grid%z(Grid%nip1,Grid%njp1))
    	allocate(Grid%theta(Grid%nip1,Grid%njp1))
    	allocate(Grid%phi(Grid%nip1,Grid%njp1))
    	allocate(Grid%r(Grid%nip1,Grid%njp1))

    	!allocate center grid
    	allocate(Grid%xc(Grid%ni,Grid%nj))
    	allocate(Grid%yc(Grid%ni,Grid%nj))
    	allocate(Grid%zc(Grid%ni,Grid%nj))
    	allocate(Grid%thetac(Grid%ni,Grid%nj))
    	allocate(Grid%phic(Grid%ni,Grid%nj))
    	allocate(Grid%rc(Grid%ni,Grid%nj))

    	!allocate i face grid
    	allocate(Grid%xi(Grid%nip1,Grid%nj))
    	allocate(Grid%yi(Grid%nip1,Grid%nj))
    	allocate(Grid%zi(Grid%nip1,Grid%nj))
    	allocate(Grid%thetai(Grid%nip1,Grid%nj))
    	allocate(Grid%phii(Grid%nip1,Grid%nj))
    	allocate(Grid%ri(Grid%nip1,Grid%nj))

    	!allocate j face grid
    	allocate(Grid%xj(Grid%ni,Grid%njp1))
    	allocate(Grid%yj(Grid%ni,Grid%njp1))
    	allocate(Grid%zj(Grid%ni,Grid%njp1))
    	allocate(Grid%thetaj(Grid%ni,Grid%njp1))
    	allocate(Grid%phij(Grid%ni,Grid%njp1))
    	allocate(Grid%rj(Grid%ni,Grid%njp1))

    	!allocate volume
    	allocate(Grid%volume(Grid%ni,Grid%nj))
    	allocate(Grid%di(Grid%ni,Grid%nj))
    	allocate(Grid%dj(Grid%ni,Grid%nj))

    	!allocate face areas
    	allocate(Grid%vol_iface(Grid%ni,Grid%nj))
    	allocate(Grid%vol_jface(Grid%ni,Grid%nj))

    	! calculate locations of cell corners
    	do i=1,Grid%nip1
    		do j=1,Grid%njp1
    			Grid%x(i,j) = Rh*sin(theta0(i))*cos(phi0(j))
    			Grid%y(i,j) = Rh*sin(theta0(i))*sin(phi0(j))
    			Grid%z(i,j) = Rh*cos(theta0(i))
    			Grid%theta(i,j) = theta0(i)
    			Grid%phi(i,j) = phi0(j)
    			Grid%R(i,j) = Rh
    		enddo 
    	enddo 

    	! calcluate locations of cell centers
    	do i=1,Grid%ni
    		do j=1,Grid%nj
    			Grid%xc(i,j) = 0.25*(Grid%x(i,j)+ &
    										  Grid%x(i+1,j)+ &
    										  Grid%x(i,j+1)+ &
    										  Grid%x(i+1,j+1))
    			Grid%yc(i,j) = 0.25*(Grid%y(i,j)+ &
    										  Grid%y(i+1,j)+ &
    										  Grid%y(i,j+1)+ &
    										  Grid%y(i+1,j+1))
    			Grid%zc(i,j) = 0.25*(Grid%z(i,j)+ &
    										  Grid%z(i+1,j)+ &
    										  Grid%z(i,j+1)+ &
    										  Grid%z(i+1,j+1))
    			Grid%thetac(i,j) = 0.25*(Grid%theta(i,j)+ &
    										  Grid%theta(i+1,j)+ &
    										  Grid%theta(i,j+1)+ &
    										  Grid%theta(i+1,j+1))
    			Grid%phic(i,j) = 0.25*(Grid%phi(i,j)+ &
    										  Grid%phi(i+1,j)+ &
    										  Grid%phi(i,j+1)+ &
    										  Grid%phi(i+1,j+1))
    			Grid%Rc(i,j) = 0.25*(Grid%R(i,j)+ &
    										  Grid%R(i+1,j)+ &
    										  Grid%R(i,j+1)+ &
    										  Grid%R(i+1,j+1))
    		enddo 
    	enddo 

    	! calcluate locations of i interfaces
    	do i=1,Grid%nip1
    		do j=1,Grid%nj
    			Grid%xi(i,j) = 0.5*(Grid%x(i,j)+ &
    										  Grid%x(i,j+1))
    			Grid%yi(i,j) = 0.5*(Grid%y(i,j)+ &
    										  Grid%y(i,j+1))
    			Grid%zi(i,j) = 0.5*(Grid%z(i,j)+ &
    										  Grid%z(i,j+1))
    			Grid%thetai(i,j) = 0.5*(Grid%theta(i,j)+ &
    										  Grid%theta(i,j+1))
    			Grid%phii(i,j) = 0.5*(Grid%phi(i,j)+ &
    										  Grid%phi(i,j+1))
    			Grid%ri(i,j) = 0.5*(Grid%r(i,j)+ &
    										  Grid%r(i,j+1))

    		enddo 
    	enddo

    	! calcluate locations of j interfaces
    	do i=1,Grid%ni
    		do j=1,Grid%njp1
    			Grid%xj(i,j) = 0.5*(Grid%x(i,j)+ &
    										  Grid%x(i+1,j))
    			Grid%yj(i,j) = 0.5*(Grid%y(i,j)+ &
    										  Grid%y(i+1,j))
    			Grid%zj(i,j) = 0.5*(Grid%z(i,j)+ &
    										  Grid%z(i+1,j))
    			Grid%thetaj(i,j) = 0.5*(Grid%theta(i,j)+ &
    										  Grid%theta(i+1,j))
    			Grid%phij(i,j) = 0.5*(Grid%phi(i,j)+ &
    										  Grid%phi(i+1,j))
    			Grid%rj(i,j) = 0.5*(Grid%r(i,j)+ &
    										  Grid%r(i+1,j))

    		enddo 
    	enddo

    	! calculate cell volume and edges
    	do i = 1,Grid%ni
    		do j = 1,Grid%nj 
    			Grid%volume(i,j) = Grid%R(i,j)*Grid%R(i,j)* &
    										 sin(Grid%thetac(i,j))* &
    										 Grid%dtheta*Grid%dPhi
			 	Grid%di(i,j) = Grid%Rc(i,j)*Grid%dtheta
			 	Grid%dj(i,j) = abs(Grid%Rc(i,j)*sin(Grid%thetac(i,j))*Grid%dphi)
			 enddo 
		enddo

		! calculate volumes for high order recontruction
		! right now this is hard coded for 8th order interpolation - will need to fix
		do i = Grid%ng+1,Grid%ni-Grid%ng+1
			do j = Grid%ng+1,Grid%nj-Grid%ng 
				! Grid%vol_iface(i,j) = (-3.*Grid%volume(i-4,j)+29*Grid%volume(i-3,j)- &
				! 			 		   139*Grid%volume(i-2,j)+533*Grid%volume(i-1,j)+ &
				! 			 		   533*Grid%volume(i,j)-139*Grid%volume(i+1,j)+ &
				! 			 		   29*Grid%volume(i+2,j)-3*Grid%volume(i+3,j))/840.0
				Grid%vol_iface(i,j) = dot_product(Cent8C,Grid%volume(i-4:i+3,j))
			enddo
		enddo
		do i = Grid%ng+1,Grid%ni-Grid%ng
			do j = Grid%ng+1,Grid%nj-Grid%ng+1 
				! Grid%vol_jface(i,j) = (-3.*Grid%volume(i,j-4)+29*Grid%volume(i,j-3)- &
				! 			 		   139*Grid%volume(i,j-2)+533*Grid%volume(i,j-1)+ &
				! 			 		   533*Grid%volume(i,j)-139*Grid%volume(i,j+1)+ &
				! 			 		   29*Grid%volume(i,j+2)-3*Grid%volume(i,j+3))/840.0
				Grid%vol_jface(i,j) = dot_product(Cent8C,Grid%volume(i,j-4:j+3))
			enddo
		enddo

		end subroutine createUniformGrid

		subroutine initKaimagModel(Model,xmlInp)
			type(kaimagModel_T), intent(inout) :: Model
			type(XML_Input_T), intent(in) :: xmlInp

			real(rp) :: t,tFin,dt,CFL,PDMB 
			integer :: ts
			character(len=strLen) :: RunID, StencilID


			!Set defaults that can be overidden by XML input
			t = 0.0
			tFin = 5.0
			dt = 1.0

			ts = 1

			!Numerical parameter defaults
			CFL = 0.3 ! 0.3 is the default value when using PDM solver
			PDMB = 1.0 ! numerical diffusion with 1 litle and 0 alot
			RunID = 'DEFAULT'
			StencilID = 'Up7C' ! 7th order upwind is default

			call xmlInp%Set_Val(Model%RunID,'sim/runid',trim(RunID))
			call xmlInp%Set_Val(Model%t,'sim/t',t)
			call xmlInp%Set_Val(Model%tFin,'sim/tFin',tFin)
			call xmlInp%Set_Val(Model%dt,'sim/dt',dt)
			call xmlInp%Set_Val(Model%ts,'sim/ts',ts)
			call xmlInp%Set_Val(Model%CFL,'sim/cfl',CFL)
			call xmlInp%Set_Val(Model%PDMB,'sim/PDMB',PDMB)
			call xmlInp%Set_Val(Model%StencilID,'sim/StencilID',StencilID)

			Model%CFL = min(0.5/(pdmb+0.5) ,CFL) !Set CFL based on PDM
			!Time options
			call Model%IO%init(xmlInp,Model%t,Model%ts)

		end subroutine initKaimagModel

		subroutine setInitialConditions(Model,Grid,State,xmlInp,updateV)
			type(kaimagModel_T), intent(inout) :: Model
			type(kaimagGrid_T), intent(in) :: Grid
			type(kaimagState_T), intent(inout) :: State
			type(XML_Input_T), intent(in) :: xmlInp
			procedure(kaimagupdateV_T), pointer, intent(inout) :: updateV

			procedure(kaimagStateIC_T), pointer :: initState => NULL()
			
			character(len=strLen) :: icStr

			!First we need to allocate the space for the state variables.
			allocate(State%rho(Grid%ni,Grid%nj))
			allocate(State%uc(Grid%ni,Grid%nj))
			allocate(State%vc(Grid%ni,Grid%nj))
			State%rho = 0 
			State%uc = 0
			State%vc = 0 

			call xmlInp%Set_Val(icStr,"sim/icType","BELL")
			call kaimagSetIC_T(initState,icStr)

			call initState(Model,Grid,State,xmlInp)

			Model%t = 0.0
			call updateV(Model,Grid,State)

			Model%dt = CalcDT(Model,Grid,State)
			model%dtold = Model%dt

			call enforceBCs(Model,Grid,State)
			


		end subroutine setInitialConditions






end module kaimaginit  