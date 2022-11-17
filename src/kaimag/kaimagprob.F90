!contains initial conditions for various problems
module kaimagprob
	use kaimagtypes
	use math
	use xml_input

	implicit none

	contains

	subroutine kaimagSetIC_T(initState,icStr)
		procedure(kaimagStateIC_T), pointer, intent(out) :: initState
		character(len=*), intent(in) :: icStr 
		!procedure(kaimagStateIC_T), pointer, intent(in) :: userInitFunc

		initState => NULL()

		select case (icStr)

		case ("BELL")
			initState => DenBell
		case ("RHO")
			initState => constantDen
		case ("BCCHECK")
			initState => BcCheck

		case default
			write(*,*) 'Unknown problem id, exiting...'
			write(*,*) icStr

		end select

	end subroutine kaimagSetIC_T

	subroutine kaimagSetUpdateV_T(Model,Grid,State,inpXML,updateV,vStr)
		type(kaimagModel_T), intent(inout) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State
		type(XML_Input_T), intent(in) :: inpXML
		procedure(kaimagUpdateV_T), pointer, intent(out) :: updateV
		character(len=*), intent(in) :: vStr 
		!procedure(kaimagStateIC_T), pointer, intent(in) :: userInitFunc

		updateV => NULL()

		select case (vStr)

		case ("BELL")
			call inpXML%Set_Val(Model%Amp,"prob/Amp",Model%Amp)
			call inpXML%Set_Val(Model%Peroid,"prob/Peroid",Model%Peroid)
			updateV => VBell
		case ("VPHI")
			call inpXML%Set_Val(Model%Vphi0,"prob/vphi0",Model%Vphi0)
			updateV => VPhi
		case ("VTHETA")
			call inpXML%Set_Val(Model%Vth0,"prob/vth0",Model%Vth0)
			updateV => VTheta

		case default
			write(*,*) 'Unknown problem id, exiting...'
			write(*,*) vStr

		end select

	end subroutine kaimagSetUpdateV_T


	subroutine DenBell(Model,Grid,State, inpXML)
		type(kaimagModel_T), intent(inout) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State
		type(XML_Input_T), intent(in) :: inpXML

		real(rp) :: ri

		integer :: i,j
		integer :: is,ie,js,je

		write(*,*) 'Initializing cosine bell density'

		call inpXML%Set_Val(Model%t1,"prob/t1",Model%t1)
		call inpXML%Set_Val(Model%t2,"prob/t2",Model%t2)
		call inpXML%Set_Val(Model%p1,"prob/p1",Model%p1)
		call inpXML%Set_Val(Model%p2,"prob/p2",Model%p2)
		call inpXML%Set_Val(Model%h_max,"prob/h_max",Model%h_max)

		is = 1
		ie = Grid%ni 
		js = 1 
		je = Grid%nj 
		do i=is,ie
			do j=js,je
				ri = Grid%rc(i,j)*acos(sin(Model%t1)*sin(PI/2.-Grid%thetac(i,j))+ &
										   cos(Model%t1)*cos(PI/2.-Grid%thetac(i,j))* &
										   cos(Grid%phic(i,j)-Model%p1))
				if (ri < Grid%rc(i,j)/2.) then 
					State%rho(i,j) = State%rho(i,j)+ &
									 Model%h_max/2.0*(1+cos(PI*ri/(Grid%rc(i,j)/2)))
				end if 
			end do 
		enddo 

		do i=is,ie
			do j=js,je
				ri = Grid%rc(i,j)*acos(sin(Model%t1)*sin(PI/2-Grid%thetac(i,j))+ &
										   cos(Model%t1)*cos(PI/2-Grid%thetac(i,j))* &
										   cos(Grid%phic(i,j)-Model%p2))
				if (ri < Grid%rc(i,j)/2.) then 
					State%rho(i,j) = State%rho(i,j) + & 
									 Model%h_max/2.0*(1+cos(PI*ri/(Grid%rc(i,j)/2)))
				end if 
			enddo 
		enddo
	end subroutine DenBell

	subroutine VBell(Model,Grid,State)
		type(kaimagModel_T), intent(in) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State


		real(rp) :: time
		integer :: i,j 

		time = Model%t

		do i=1,Grid%ni
			do j=1,Grid%nj
				State%uc(i,j) = Model%Amp*Grid%Rc(i,j)/Model%Peroid* &
								sin(2.*(Grid%phic(i,j)-2*PI*time/Model%Peroid))* &
								cos(PI/2-Grid%thetac(i,j))*cos(PI*time/Model%Peroid)
				State%vc(i,j) = -Model%Amp*Grid%Rc(i,j)/Model%Peroid* &
								(sin(Grid%phic(i,j)-2*PI*time/Model%Peroid))**2* &
								sin(2*(PI/2-Grid%thetac(i,j)))*cos(PI*time/Model%Peroid) + &
								2*PI*Grid%Rc(i,j)/Model%Peroid*cos(PI/2-Grid%thetac(i,j))
			enddo 
		enddo 
	end subroutine VBell

	subroutine BcCheck(Model,Grid,State,inpXML)
		type(kaimagModel_T), intent(inout) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State
		type(XML_Input_T), intent(in) :: inpXML

		integer :: i,j
		integer :: is,ie,js,je

		write(*,*) 'Initializing boundary condition check values'

		js = Grid%ng+1
		je = Grid%nj-Grid%ng 
		is = Grid%ng+1 
		ie = Grid%ni-Grid%ng 
		do j=js,je
			do i=is,ie 
				State%rho(i,j) = i-Grid%ng + 10000*(j-Grid%ng)
			enddo
		enddo
	end subroutine BcCheck

	subroutine constantDen(Model,Grid,State,inpXML)
		type(kaimagModel_T), intent(inout) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State
		type(XML_Input_T), intent(in) :: inpXML

		real(rp) :: rhoval=1.0

		integer :: i,j
		integer :: is,ie,js,je

		write(*,*) 'Initializing constant density'

		call inpXML%Set_Val(rhoval,'prob/rhoval',rhoval)

		js = Grid%ng+1
		je = Grid%nj-Grid%ng 
		is = Grid%ng+1 
		ie = Grid%ni-Grid%ng 
		do j=js,je
			do i=is,ie 
				State%rho(i,j) = rhoval
			enddo
		enddo

	end subroutine constantDen

	subroutine VPhi(Model,Grid,State)
		type(kaimagModel_T), intent(in) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State


		real(rp) :: time
		integer :: i,j 

		time = Model%t

		do i=1,Grid%ni
			do j=1,Grid%nj
				State%uc(i,j) = 0.0
				State%vc(i,j) = Model%Vphi0*sin(Grid%thetac(i,j))
			enddo 
		enddo 
	end subroutine VPhi

	subroutine VTheta(Model,Grid,State)
		type(kaimagModel_T), intent(in) :: Model
		type(kaimagGrid_T), intent(in) :: Grid
		type(kaimagState_T), intent(inout) :: State

		real(rp) :: time
		integer :: i,j 

		time = 0.0

		! do i=1,Grid%ni
		! 	do j=1,Grid%nj
		! 		State%uc(i,j) = -Model%Amp*Grid%Rc(i,j)/Model%Peroid* &
		! 						sin(2.*(Grid%phic(i,j)-2*PI*time/Model%Peroid))* &
		! 						cos(PI/2-Grid%thetac(i,j))*cos(PI*time/Model%Peroid)
		! 		State%vc(i,j) = 0.0
		! 	enddo 
		! enddo 

		do i=1,Grid%ni
			do j=1,Grid%nj/4
				State%uc(i,j) = Model%Vth0
				State%vc(i,j) = 0
			enddo 
		enddo 
		do i=1,Grid%ni
			do j=Grid%nj/4+1,Grid%nj/2
				State%uc(i,j) = -1.0*Model%Vth0
				State%vc(i,j) = 0
			enddo 
		enddo
		do i=1,Grid%ni
			do j=Grid%nj/2+1,Grid%nj/2+Grid%nj/4
				State%uc(i,j) = Model%Vth0
				State%vc(i,j) = 0
			enddo 
		enddo 
		do i=1,Grid%ni
			do j=Grid%nj/2+Grid%nj/4+1,Grid%nj
				State%uc(i,j) = 1.0*Model%Vth0
				State%vc(i,j) = 0
			enddo 
		enddo  
	end subroutine VTheta

end module kaimagprob

