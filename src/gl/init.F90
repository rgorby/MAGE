module init
    use gltypes
    use math
    use xml_input
    use glutils
    use glio5
    use files    
    !use output
    !use glioH5

    implicit none
    !Initialization defaults (2D Square)
    integer, private :: Nc = 64
    real(rp) :: xMin = -1.0, xMax = 1.0

    !----------------
    ! Notes: 
    !
    !
    !----------------
    contains
      
        !> Assume reasonable defaults if not overidden by user 
        !> for Gibson-Low model


        subroutine setModelDefaults(Model)
            type(glModel_T), intent(inout) :: Model

            ! ; 
            ! ; set up defaults if necessary
            ! ; as in Dove et al 2011
            ! ;
            Model%frontheight = 1.35    ! height of the CME front (in Rsun) *if* apar was 0
            Model%legsang = 25.         ! solid angle of the CME
            Model%apar = 0.05           ! transformation r -> r+a, which, e.g., converts a sphere to a teardrop
            Model%ao = 1.75             ! "the constant of proportionality between the stream function S and the pressure Pi
                                        ! It is later divided by muse (which is negative).
                                        ! so that the amplitude of Ao ends up being the same as the maximum amplitude of Br (bubble coords)."
            Model%alnotrbub = 1.        ! "the eigenvalue parameter of the twisted flux rope model"; "RESTRICTED TO EIGENVALUES OF J5/2
                                        ! BESSEL FUNCTION"; DEFAULT ALNOTRBUB=1 corresponds to smallest eigenvalue -- 5.763854

            Model%sigma = PI/2.         ! the rotation about the [?] surface normal, as in, like a tilt of an AR, default is pi/2
            Model%outScale = 0.1        ! "a parameter which scales the outer field solution since there is a discontinuity around bubble for field"
            Model%pio = 0.              ! constant for the background pressure
                                        ! for 1d8 see last lines, within outarray in giblowprams.pro
                                        ! this is a block of atmosphere parameters outside of the bubble (we probably won't ever use these)
                                        ! as in, Pout=(aa/(bb+1.))*r^(-bb-1.)+(cc/(dd+1.))*r^(-dd-1.)+(ee/(ff+1.))*r^(-ff-1.)
                                        ! Schmit & Gibson 2011
            Model%aa = 1.27*1d8         
            Model%bb = 22.8
            Model%cc = 5.97*1d8
            Model%dd = 16.0
            Model%ee = 1.57*1d8
            Model%ff = 3.87

            Model%bonly = 0.            ! a flag to not calculate plasma parameters (dens, pres, etc), only the magnetic field
            Model%c1bonly = 0.5         ! c[1|2]bonly "only used if BONLY set; then creates a cavity which is C1BONLY fractional density
            Model%c2bonly = 0.9         ! (relative to background) and C2BONLY fractional height (relative to bubble)"
            Model%isothermal = 1.55e6   ! if bonly is set, the model temperature will have this value
            Model%bubbleonly = 0.       ! a flag to only calculate the model for interior of the bubble (and don't do anything for points outside of bubble)

            Model%alpha = 0.            ! acceleration (we probably won't ever make it non-zero)
            Model%velmult = 0.1         ! multiplier for 'velocity'
            Model%phiss = 1.0           ! center of the bubble - included for completeness: this is not an *input*, it's calculated from other parameters
            Model%time = 0.*3600        ! Current Model time in seconds!
            Model%velimpose = 0.        ! if not zero, output model velocity is overriden with this value (we probably won't use it)
            Model%dt = 1.0              ! Model timestep - default to 1.0 second,
            Model%ts = 0
            Model%tfin = 100.           ! defaults to 100 seconds
            Model%alnotrbubsign = 1.0   ! Chirality

            ! Alternate defauult GLpars, keeping here for reference AJM        
            ! GLpars = (/1.35, 25., 0.05, 1.75, 1.0, 1.57080, 0.0, 0.0, 0.0, 22.8, &
            !         5.97, 16.0, 0.0, 3.87, 0.0, 0.5, 0.9, 1.55e6, 0.0, 0.0, &
            !         1.0, 1.0, 0.0, 1.0/)
            !real(rp) :: alnotsign  !introduced to change direction of winding TODO: not sure if this is diff from chiral?
    
        end subroutine setModelDefaults

        !>
        !>
        subroutine initModel(Model, xmlInp)
            type(glModel_T), intent(inout) :: Model
            type(XML_Input_T), intent(inout) :: xmlInp

            !Get RunID
            call xmlInp%Set_Val(Model%RunID,'sim/runid',"Sim")
            call setModelDefaults(Model)

            ! Check State Type "Sphere" or 
            call xmlInp%Set_Val(Model%CoordName, "prob/StateCoord", "Sphere")
            !read in GL parameters to override defaults
            call xmlInp%Set_Val(Model%tfin, "prob/tfin", 100.)
            call xmlInp%Set_Val(Model%dt, "prob/dt", 1.)
            call xmlInp%Set_Val(Model%latitude, "prob/lat", 0.)
            call xmlInp%Set_Val(Model%longitude, "prob/lon", 0.)
            call xmlInp%Set_Val(Model%frontheight, "prob/frontheight", 1.35)
            call xmlInp%Set_Val(Model%legsang, "prob/legsang", 25.)
            call xmlInp%Set_Val(Model%topmorph, "prob/topmorph", 2.5)
            call xmlInp%Set_Val(Model%Bpar, "prob/Bmax", 1.75)
            call xmlInp%Set_Val(Model%alpha, "prob/alpha", 0.)
            call xmlInp%Set_Val(Model%sigma, "prob/orientation", 1.5708)
            call xmlInp%Set_Val(Model%alnotrbubsign, "prob/chiral", 1.0) 
            call xmlInp%Set_Val(Model%vel_fh, "prob/vel_fh", 1.)  !V_CME_statdraw
            call xmlInp%Set_Val(Model%IO%doTimerOut,'output/timer',.false.)
            call xmlInp%Set_Val(Model%isLoud,'sim/isLoud',.true.)
            call xmlInp%Set_Val(Model%isDebug,'sim/isDebug',.true.)
            call xmlInp%Set_Val(Model%isPrecheck,'sim/isPrecheck',.false.)

            ! From precheck
            if (Model%isPrecheck) then
                Model%k = tan(0.5*Model%legsang*mdtor)
                Model%apar = Model%frontheight*(1 - 0.5*Model%k*(Model%topmorph - 2))/(0.5*Model%k*Model%topmorph)
                Model%r0 = 2*Model%frontheight/Model%topmorph
                Model%x0 = Model%r0/Model%k
            end if
            Model%s_eta0 = sqrt(eta0)
            ! from v_fh = frontheight*sqrt(eta)*Rsun_km*velmult; sqrt(eta)*Rsun_km=258.55
            ! can get velmult=Model%vel_fh/(Model%frontheight*s_eta0*Rsun*1d-5)
            Model%s_eta = Model%vel_fh/(Model%frontheight*Rsun*1d-5) ! sqrt(eta0)*velmult
            Model%velmult = Model%s_eta/Model%s_eta0
            ! (this is simply to rewrite the sequence
            !  s_eta0=sqrt(eta0); velmult=Model%vel_fh/(Model%frontheight*Rsun*1d-5*s_eta0); s_eta=s_eta0*velmult)

            Model%alnotrbub_abs = abs(Model%alnotrbub)
    
            if (abs(Model%alnotrbub_abs - alnotrbub0) .gt. 0.01) then
                Model%alnotrbub_abs = 1
            end if
    
            Model%alnotrbub_abs = int(Model%alnotrbub_abs)
            if ((Model%alnotrbub_abs .lt. 1) .or. (Model%alnotrbub_abs .gt. 6)) then
                if (Model%isLoud) write(*,*) 'warning: use eigenvalues between 0 and 6; resetting to 1'
                Model%alnotrbub_abs = 1
            end if
    
            Model%xo = (Model%frontheight + Model%apar)/(1.+tan(Model%legsang/2.*mdtor))
            Model%rbub = Model%frontheight - Model%xo + Model%apar
            Model%alnotrbubuse = alnotrbubuse1(int(Model%alnotrbub_abs))*Model%alnotrbubsign
            Model%alnot = Model%alnotrbubuse/Model%rbub
            Model%gfunot1 = sin(Model%alnotrbubuse)/Model%alnotrbubuse - cos(Model%alnotrbubuse)
            Model%Muse = 8.*pi*(Model%rbub*Model%rbub/3./Model%gfunot1 - 1./Model%alnot/Model%alnot)
            Model%ao = -Model%ao/Model%Muse
            Model%alpha = Model%alpha*Model%velmult*Model%velmult
            Model%eta = eta0*Model%velmult*Model%velmult

            if (Model%isDebug) then
                write(*,"(1X,A14,2X,F)") "frontheight: ", Model%frontheight
                write(*,"(1X,A14,2X,F)") "xo: ", Model%xo
                write(*,"(1X,A14,2X,F)") "apar: ", Model%apar
                write(*,"(1X,A14,2X,F)") "r0: ", Model%r0
                write(*,"(1X,A14,2X,F)") "x0: ", Model%x0
                write(*,"(1X,A14,2X,F)") "k: ", Model%k
                write(*,"(1X,A14,2X,F)") "rbub: ", Model%rbub
                write(*,"(1X,A14,2X,F)") "alnotrbubuse: ", Model%alnotrbubuse
                write(*,"(1X,A14,2X,F)") "alnot: ", Model%alnot
                write(*,"(1X,A14,2X,F)") "gfunot1: ", Model%gfunot1
                write(*,"(1X,A14,2X,F)") "ao: ", -Model%ao*Model%Muse
                write(*,"(1X,A14,2X,F)") "ao scaled: ", Model%ao
                write(*,"(1X,A14,2X,F)") "Muse: ", Model%Muse
                write(*,"(1X,A14,2X,F)") "eta: ", Model%eta
                write(*,"(1X,A14,2X,F)") "velmult: ", Model%velmult
                write(*,"(1X,A14,2X,F)") "s_eta: ", Model%s_eta 
            end if

            !Output/Restart (IOCLOCK)
            call Model%IO%init(xmlInp,Model%time,Model%ts)
            ! TODO: Implement solution for acceleration
            if ( (Model%alpha .gt. 0) .or. (Model%alpha .lt. 0)) then
                if (Model%isLoud) write(*,*) 'Currently not supported, defaulting to constant velocity: alpha = 0.0'
                Model%alpha = 0.0
                !   print *, 'have only set up the alpha/eta = 0.49 case, will assume that one'
                !   Model%alpha = 6.762d-8*Model%velmult*Model%velmult
                !
                !   restore,'/MODELS/GIBLOW/PARAMETERS/phiss_vs_t_0.490000_1.38000d-07.dat'
                !   timearr=timearr/velmult
                !   test = where(abs(timearr-time) eq min(abs(timearr-time)))
                !   phissuse = phissarr(test)
                !  phiss = phissuse(0)
            end if
    
            ! if (alpha .lt. 0) then
            !  print *, 'have not set up decelerating case'
            !  stop
            ! endif
        end subroutine initModel

        !> Initialize Solution and State
        !>
        subroutine initSolutionState(Model, State, Solution)
            type(glModel_T), intent(inout) :: Model
            type(glSolution_T), intent(inout) :: Solution
            type(glState_T), intent(inout) :: State

            call allocState(State)
            call allocSolution(State, Solution)

            ! Initalize Solution
            Solution%dens = 0.
            Solution%pres = 0.
            Solution%temp = 0.
            Solution%b = 0.
            Solution%v = 0.
            Solution%j = 0.
            Solution%inside_mask = 0.
            State%ri = 1
 
        end subroutine initSolutionState

        !>
        !>
        !>
        subroutine initStandaloneSolutionState(Model, State, Solution, xmlInp)
            type(glModel_T), intent(inout) :: Model
            type(glSolution_T), intent(inout) :: Solution
            type(glState_T), intent(inout) :: State
            type(XML_Input_T), intent(inout) :: xmlInp

            real(rp) :: dr, dth, dphi
            real(rp) :: x, y, z, x1, x2, x3
            real(rp) :: rtpBds(6)
            integer :: i, j, k
            ! Set up State grid dimensions       
            call xmlInp%Set_Val(State%Ni,"idir/N",Nc)
            call xmlInp%Set_Val(State%Nj,"jdir/N",Nc)
            call xmlInp%Set_Val(State%Nk,"kdir/N",Nc)
            call xmlInp%Set_Val(rtpBds(1),"idir/min",xMin)
            call xmlInp%Set_Val(rtpBds(2),"idir/max",xMax)
            call xmlInp%Set_Val(rtpBds(3),"jdir/min",xMin)
            call xmlInp%Set_Val(rtpBds(4),"jdir/max",xMax)
            call xmlInp%Set_Val(rtpBds(5),"kdir/min",xMin)
            call xmlInp%Set_Val(rtpBds(6),"kdir/max",xMax)

            call initSolutionState(Model, State, Solution)
            ! Set up State
            ! Set dx's (uniform for now)
            dr = (rtpBds(2)-rtpBds(1))/State%Ni
            dth = (rtpBds(4)-rtpBds(3))/State%Nj
            dphi = (rtpBds(6)-rtpBds(5))/State%Nk

            ! Create r(Ni), theta(Nj,Nk), phipb(Nj,Nk) arrays
              
            do k=1, State%Nk
                do j=1, State%Nj
                    ! Theta
                    x2 = rtpBds(3)+(j-1)*pi*dth + Model%latitude
                    ! Phi
                    x3 = rtpBds(5)+(k-1)*2*pi*dphi - Model%longitude
                    State%thpb(j,k) = x2
                    State%phpb(j,k) = x3
                    do i=1, State%Ni
                        x1 = rtpBds(1) + (i-1)*dr
                        State%r(i) = x1
                        
                        x = x1*sin(x2)*cos(x3)
                        y = x1*sin(x2)*sin(x3)
                        z = x1*cos(x2)
                        ! This grid is for external processing to set up cartesian coordinates 
                        ! for analyzing in visualization routines e.g. Paraview
                        State%xyz(i,j,k,:) = [x,y,z]
                    end do 
                end do
            end do
        end subroutine initStandaloneSolutionState

        !> This init routine is to be used for standalone runs of the Gibson-Low Model
        !> on a State. Used mainly for testing purposes and to export h5 file of the 
        !> GL soluiton on a State for additional analysis or reading into an external 
        !> program. 
        subroutine initGLStandalone(Model, State, Solution, xmlInp)
            type(glModel_T), intent(inout) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            type(XML_Input_T), intent(inout) :: xmlInp

            call initModel(Model, xmlInp)
            call initStandaloneSolutionState(Model, State, Solution, xmlInp)
            ! State Type dicates output coordinate system
            select case (Model%CoordName)
                case ("Sphere")
                    Solution%CoordSystem = SPHERICAL
                    write(*,*) 'Coordinate system is Spherical'
                case ("Cart3D")
                    Solution%CoordSystem = CARTESIAN
                    write(*,*) 'Coordinate system is Cartesian'
                case default
                    write(*,*) 'You must specify a valid StateCoord, "Sphere" or "Cart3D"'
                    stop
            end select
            GLH5File   = genName(Model%RunID, State%Ni, State%Nj, State%Nk, State%Ni+1, State%Nj+1, State%Nk+1)
        end subroutine initGLStandalone

        !> Expectation is that this interface is to be used by Gamera etc. 
        !> User must initalize the glApp components -> Model, State, State, Solution
        !> Model is initialized from xmlInp
        !> Parameters:
        !>  xyz: the (i,j,k,xyz) locations over which to solve the Gibson Low Model
        !>  glApp: instance of the Gibson Low model app
        !>  xmlInp: the xmlInp for this model
        !> 
        subroutine initGLInterfaceXYZ(xyz, glApp,  xmlInp)
            type(glApp_T), intent(inout) :: glApp
            type(XML_Input_T), intent(inout) :: xmlInp
            real(rp), dimension(:,:,:,:), intent(in) :: xyz

            real(rp), dimension(:), allocatable :: r
            real(rp), dimension(:,:), allocatable :: theta, phi
            integer, dimension(2) :: adims
            integer :: rdim, i

            rdim = size(xyz(:,1,1,1))
            adims = shape(xyz(1,:,:,ZDIR))

            allocate(r(rdim))
            allocate(theta(adims(1), adims(2)))
            allocate(phi(adims(1), adims(2)))
            do i=1, rdim
                r(i) = norm2(xyz(i,1,1,:))
            end do
            theta = acos(xyz(1,:,:,ZDIR)/norm2(xyz(1,1,1,:)))
            phi = atan2(xyz(1,:,:,YDIR),xyz(1,:,:,XDIR))

            call initGLInterfaceRTP(r, theta, phi, glApp, xmlInp)

        end subroutine initGLInterfaceXYZ

        !> Expectation is that this interface is to be used by Gamera etc. 
        !> User must initalize the glApp components -> Model, State, State, Solution
        !> Model is initialized from xmlInp
        !> State must minimally be defined by r(i), theta(j,k), phi(j,k) 
        !> 
        subroutine initGLInterfaceRTP(r, theta, phi, glApp,  xmlInp)
            type(glApp_T), intent(inout) :: glApp
            type(XML_Input_T), intent(inout) :: xmlInp
            real(rp), dimension(:), intent(in) :: r
            real(rp), dimension(:,:), intent(in) :: theta, phi
            integer :: rdim
            integer, dimension(2) :: adims

            rdim = size(r)
            adims = shape(theta)

            associate(State=>glApp%State, Model=>glApp%Model, Solution=>glApp%Solution)

                call initModel(glApp%Model, xmlInp)

                call initSolutionState(Model, State, Solution)

                State%Ni = rdim
                State%Nj = adims(1)
                State%Nk = adims(2)
                State%r = r
                State%thpb = theta + Model%latitude
                State%phpb = phi - Model%longitude
            end associate
        end subroutine initGLInterfaceRTP
        
end module init