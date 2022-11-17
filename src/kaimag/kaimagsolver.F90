module kaimagsolver
	use kdefs
	use kaimagdefs
	use kaimagtypes

	implicit none

	contains



	function calcDT(Model,Grid,State)

		type(kaimagModel_T), intent(in) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(in) :: State
		real(rp) :: calcDT 

		
		real(rp) :: dtMin, dtOld, dtijk,u
		integer :: i,j,k
		integer :: is,ie,js,je

		dtMin = HUGE 
		do j=1,Grid%nj 
			do i=1,Grid%ni
				u = sqrt(State%uc(i,j)**2+State%vc(i,j)**2) 
				dtijk = Model%CFL/(abs(u)/Grid%di(i,j)+abs(u)/Grid%dj(i,j))
				dtMin = min(dtijk,dtMin)
			enddo 
		enddo 

		calcDT = dtMin
		!write(*,*) 'calcDT', calcDT, Model%t, Model%tFin
		!Make sure we don't overstep the end of the simulation
		if (Model%t > Model%tFin) then 
			calcDT = max(Model%tFin -Model%t,TINY)
		endif

	end function calcDT

	subroutine EnforceBCs(Model,Grid,State)
		type(kaimagModel_T), intent(in) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State

		integer :: i,is,ie 
		integer :: j,js,je
		integer :: ny,ny2


		! Peroidic in the PHI Direction
		!write(*,*) 'From EnforceBC Left', (Grid%nj-2*Grid%ng+1),(Grid%nj-Grid%ng)
		!write(*,*) State%rho(5,(Grid%nj-2*Grid%ng+1):(Grid%nj-Grid%ng))
		is = Grid%ng+1
		ie = Grid%ni-Grid%ng
		js = 1 
		je = Grid%ng 
		do j=js,je
			do i=is,ie
				State%rho(i,j) = State%rho(i,(Grid%nj-2*Grid%ng+j))
			enddo
		enddo
		!write(*,*) 'From EnforceBC Right ', (Grid%nj-Grid%ng+1),Grid%nj
		!write(*,*) State%rho(5,(Grid%ng+1):(2*Grid%ng))
		do j = js,je
			do i=is,ie
				State%rho(i,(Grid%nj-Grid%ng+j)) = State%rho(i,(Grid%ng+j))
			enddo
		enddo
		! Pole boundary - no need if no pole in grid
		ny = (Grid%nj-2*Grid%ng)
		ny2 = ny/2
		is = 1
		ie = Grid%ng 
		js = Grid%ng+1
		je = Grid%ng+ny2
		do i=is,ie
			do j=js,je
				State%rho(i,j) = State%rho(2*Grid%ng+1-i,j+ny2)
			enddo
		enddo
		js = Grid%ng+ny2
		je = ny + Grid%ng
		do i=is,ie
			do j=js,je
				State%rho(i,j) = State%rho(2*Grid%ng+1-i,j-ny2)
			enddo
		enddo
		!Southern Hemisphere
		is = Grid%ni-Grid%ng+1
		ie = Grid%ni
		js = Grid%ng+1
		je = Grid%ng+ny2
		do i=is,ie
			do j=js,je
				State%rho(i,j) = State%rho(Grid%ni-Grid%ng-(i-is),j+ny2)
			enddo
		enddo
		js = Grid%ng+ny2
		je = ny + Grid%ng
		do i=is,ie
			do j=js,je
				State%rho(i,j) = State%rho(Grid%ni-Grid%ng-(i-is),j-ny2)
			enddo
		enddo

	end subroutine EnforceBCs

	subroutine Advance2DSph(Model,Grid,State,oState,updateV)
		type(kaimagModel_T), intent(inout) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State,oState
		procedure(kaimagupdateV_T), pointer, intent(inout) :: updateV

		integer :: i,j,k
		integer :: is,ie,js,je,ks,ke

		real(rp), dimension(:,:), allocatable :: rho_h,rho_left,rho_right
		real(rp), dimension(:,:), allocatable :: u_left,u_right,v_left,v_right
		real(rp), dimension(:,:), allocatable :: rho_flux_x,rho_flux_y

		!Allocate space for work arrays

		allocate(rho_h     (Grid%ni,Grid%nj))
		allocate(rho_left  (Grid%ni,Grid%nj))
		allocate(rho_right (Grid%ni,Grid%nj))
		allocate(rho_flux_x(Grid%ni,Grid%nj))
		allocate(rho_flux_y(Grid%ni,Grid%nj))
		allocate(u_left    (Grid%ni,Grid%nj))
		allocate(u_right   (Grid%ni,Grid%nj))
		allocate(v_left    (Grid%ni,Grid%nj))
		allocate(v_right   (Grid%ni,Grid%nj))

		!write(*,*) 'Time before updateV',Model%t
		call updateV(Model,Grid,State)

		js = 1 
		je = Grid%nj
		is = 1 
		ie = Grid%ni
		do j=js,je
			do i=is,ie
				rho_h(i,j) = State%rho(i,j) + Model%dt/Model%dtOld/2* & 
				(State%rho(i,j)-oState%rho(i,j))
			enddo 
		enddo
		!Save the old state for the next predictor step
		oState = State
		Model%dtOld = Model%dt

		! First do the Theta direction
		call reconstruct(Model,Grid,rho_h,rho_left,rho_right,THETADIR)
		call reconstruct(Model,Grid,State%uc,u_left,u_right,THETADIR)

		js = Grid%ng+1
		je = Grid%nj-Grid%ng
		is = Grid%ng+1
		ie = Grid%ni-Grid%ng+1
		do j=js,je
			do i=is,ie 
				rho_flux_x(i,j) = (1+sign(1.0_rp,u_left(i,j)))/2*rho_left(i,j)*u_left(i,j) + &
								  (1-sign(1.0_rp,u_right(i,j)))/2*rho_right(i,j)*u_right(i,j)
			enddo
		enddo

		!Now do PHI Direction
		call reconstruct(Model,Grid,rho_h,rho_left,rho_right,PHIDIR)
		call reconstruct(Model,Grid,State%vc,v_left,v_right,PHIDIR)

		js = Grid%ng+1
		je = Grid%nj-Grid%ng+1
		is = Grid%ng+1
		ie = Grid%ni-Grid%ng
		do j=js,je
			do i=is,ie 
				rho_flux_y(i,j) = (1+sign(1.0_rp,v_left(i,j)))/2*rho_left(i,j)*v_left(i,j) + &
								  (1-sign(1.0_rp,v_right(i,j)))/2*rho_right(i,j)*v_right(i,j)
			enddo
		enddo

		js = Grid%ng+1
		je = Grid%nj-Grid%ng
		is = Grid%ng+1
		ie = Grid%ni-Grid%ng
		do j=js,je
			do i=is,ie 
				State%rho(i,j) = oState%rho(i,j) - Model%dt*( &
					(rho_flux_x(i+1,j)*sin(Grid%thetai(i+1,j))- &
					rho_flux_x(i,j)  *sin(Grid%thetai(i,j)))/ &
					Grid%Rc(i,j)/sin(Grid%thetac(i,j))/Grid%dtheta+ &
					(rho_flux_y(i,j+1)-rho_flux_y(i,j))/ &
					Grid%Rc(i,j)/sin(Grid%thetac(i,j))/Grid%dphi)
			enddo
		enddo

	end subroutine Advance2DSph

	subroutine reconstruct(Model,Grid,var,var_left,var_right,DIR)
		type(kaimagModel_T), intent(in) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		real(rp), dimension(:,:), intent(in) :: var 
		real(rp), dimension(:,:), intent(inout) :: var_left, var_right
		integer, intent(in) :: DIR 

		real(rp), dimension(:,:), allocatable :: var_h0, var_h, var_interp
		real(rp), dimension(8), target :: Cent8C = [-3,29,-139,533,533,-139,29,-3]/840.0_rp
    	real(rp), dimension(6), target :: Cent6C = [1,-8,37,37,-8,1]/60.0_rp
    	real(rp), dimension(7), target :: Up7C = [-3,25,-101,319,214,-38,4]/420.0_rp
    	real(rp), dimension(7), target :: High5C = [0,2,-13,47,27,-3,0]/60.0_rp
    	
    	real(rp), pointer :: pStencil(:)
    	real(rp), dimension(8) :: QdV,Q
    	real(rp) :: dVi,Qi
    	
		integer i,j
		integer is,ie, js,je

		allocate(var_h0(Grid%ni,Grid%nj))
		allocate(var_h(Grid%ni,Grid%nj))
		allocate(var_interp(Grid%ni,Grid%nj))

		var_left = 0.0*var
		var_right = 0.0*var

		var_h0 = var
		var_h = var*Grid%volume

		if (Model%StencilID == 'High5C') then
		 	pStencil => High5C
		else 
		 	pStencil => Up7C
		end if

		!Get loop bounds
		if (dir .eq. THETADIR) then
			is = Grid%ng+1
			ie = Grid%ni-Grid%ng+1 
			js = Grid%ng+1 
			je = Grid%nj-Grid%ng 
		else
			is = Grid%ng+1
			ie = Grid%ni-Grid%ng
			js = Grid%ng+1
			je = Grid%nj-Grid%ng+1 
		endif

		!Now loop over interfaces
		do j=js,je
			do i=is,ie
				!Pull stencil
				if (dir .eq. THETADIR) then
					QdV = var_h (i-4:i+3,j)
					Q   = var_h0(i-4:i+3,j)
					dVi = Grid%vol_iface(i,j)
				else
					QdV = var_h (i,j-4:j+3)
					Q   = var_h0(i,j-4:j+3)
					dVi = Grid%vol_jface(i,j)
				endif

				!Left state
				Qi = dot_product(pStencil,QdV(1:7:+1))/dVi
				var_left(i,j) = PDM(Q(3),Q(4),Q(5),Qi,Model%PDMB)

				!Right state
				Qi = dot_product(pStencil,QdV(8:2:-1))/dVi
				var_right(i,j) = PDM(Q(6),Q(5),Q(4),Qi,Model%PDMB)
			enddo
		enddo

		! if (dir .eq. IDIR) then
		! 	is = Grid%ng+1
		! 	ie = Grid%ni-Grid%ng+1 
		! 	js = Grid%ng+1 
		! 	je = Grid%nj-Grid%ng 
		! 	do j=js,je
		! 		do i=is,ie 
		! 			var_interp(i,j) = dot_product(pStencil,var_h(i-4:i+2,j))
		! 			var_interp(i,j) = var_interp(i,j)/Grid%vol_iface(i,j)
		! 			var_left(i,j) = PDM(var_h0(i-2,j),var_h0(i-1,j),&
		! 								var_h0(i,j),var_interp(i,j),Model%PDMB)
		! 		enddo
		! 	enddo
		! 	do j=js,je
		! 		do i=is,ie 
		! 			var_interp(i,j) = dot_product(pStencil,var_h(i+3:i-3:-1,j))
		! 			var_interp(i,j) = var_interp(i,j)/Grid%vol_iface(i,j)
		! 			var_right(i,j) = PDM(var_h0(i+1,j),var_h0(i,j),&
		! 								 var_h0(i-1,j),var_interp(i,j),Model%PDMB)
		! 		enddo
		! 	enddo
		! endif

		! if (dir .eq. JDIR) then
		! 	is = Grid%ng+1
		! 	ie = Grid%ni-Grid%ng
		! 	js = Grid%ng+1
		! 	je = Grid%nj-Grid%ng+1 
		! 	do j=js,je
		! 		do i=is,ie 
		! 			var_interp(i,j) = dot_product(pStencil,var_h(i,j-4:j+2))

		! 			var_interp(i,j) = var_interp(i,j)/Grid%vol_jface(i,j)
		! 			var_left(i,j) = PDM(var_h0(i,j-2),var_h0(i,j-1),&
		! 								var_h0(i,j),var_interp(i,j),Model%PDMB)
		! 		enddo
		! 	enddo
		! 	do j=js,je
		! 		do i=is,ie 
		! 			var_interp(i,j) = dot_product(pStencil,var_h(i,j+3:j-3:-1))
		! 			var_interp(i,j) = var_interp(i,j)/Grid%vol_jface(i,j)
		! 			var_right(i,j) = PDM(var_h0(i,j+1),var_h0(i,j),&
		! 								 var_h0(i,j-1),var_interp(i,j),Model%PDMB)
		! 		enddo
		! 	enddo
		! endif
	end subroutine reconstruct 

    function PDM(q0,q1,q2,qI,pdmb)
        real(rp), intent(in) :: q0,q1,q2,qI,pdmb
        real(rp) :: PDM
        real(rp) :: maxQ,minQ, qN, dq0,dq1,dqL
        real(rp) :: s0,s1,s01

        !Max/min of nearest neighbors
        maxQ = max(q1,q2)
        minQ = min(q1,q2)
        qN = max(minQ,min(qI,maxQ))

        !Local differences/flips
        dq0 = pdmb*(q1-q0)
        dq1 = pdmb*(q2-q1)

        s0 = sign(1.0_rp,dq0)
        s1 = sign(1.0_rp,dq1)

        s01 = abs(s0+s1)

        !Local slopes
        dqL = qN-q1

        !Replace w/ limited value
        PDM = qN - s1*max(0.0,abs(dqL)-s01*abs(dq0))
        
    end function PDM

    subroutine calcError(Model,Grid,State,StateIC)
		type(kaimagModel_T), intent(in) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
    	type(kaimagState_T), intent(in) :: State
    	type(kaimagState_T), intent(in) :: StateIC
    	integer i,j,is,ie,js,je

    	real(rp) :: totalError, initialTotalDensity, finalTotalDensity
    	real(rp) :: massError, initialMass, finalMass
    	real(rp) :: l2norm

    	! compare final density to initial density

		is = Grid%ng+1
		ie = Grid%ni-Grid%ng 
		js = Grid%ng+1 
		je = Grid%nj-Grid%ng
		initialMass = SUM(StateIC%rho(is:ie,js:je)*Grid%volume(is:ie,js:je))
		finalMass =  SUM(State%rho(is:ie,js:je)*Grid%volume(is:ie,js:je))
	    initialTotalDensity = SUM(StateIC%rho(is:ie,js:je))
	    finalTotalDensity = SUM(State%rho(is:ie,js:je))
	    totalError = SUM(ABS(StateIC%rho(is:ie,js:je) - State%rho(is:ie,js:je)))

	 	l2norm = norm2((StateIC%rho(is:ie,js:je) - State%rho(is:ie,js:je))*Grid%volume(is:ie,js:je))
	    write(*,*) 'initial mass', initialMass
	    write(*,*) 'final mass', finalMass
 	    write(*,*) 'Precent loss',ABS((finalMass-initialMass))/initialMass*100.
	    write(*,*) 'Total error', totalError
	    write(*,*) 'Percent Dissipation',totalError/initialTotalDensity*100.
	    write(*,*) 'L2 norm', l2norm

	    return 
	end subroutine calcError
end module kaimagsolver
