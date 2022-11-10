module glinit
    use gltypes
    use glsolution
    use math
    use xml_input
    use glutils
    use glio5
    use files    

    implicit none
    !Initialization defaults (2D Square)
    integer, private :: Nc = 64
    real(rp) :: xMin = -1.0, xMax = 1.0

    !----------------------------------------------
    ! Notes: 
    !   Standalone model constructs a psuedo-grid in which to calculate the 
    !   Gibson-Low solution
    !
    !   Interface initialization assumes the construction of state
    !   from external source specifying the r(Ni), theta(Nj, Nk), phi(Nj, Nk)
    !   values. Radius and angles are independent of eachother. 
    !
    !----------------------------------------------
    contains

        !> Set model atmosphere parameters to Schmit & Gibson 2011
        subroutine setSchmitAtmosphere(Model)
            type(glModel_T), intent(inout) :: Model
            Model%ao = 1.75
            ! for 1d8 see last lines, within outarray in giblowprams.pro
            ! this is a block of atmosphere parameters outside of the bubble (we probably won't ever use these)
            ! as in, Pout=(aa/(bb+1.))*r^(-bb-1.)+(cc/(dd+1.))*r^(-dd-1.)+(ee/(ff+1.))*r^(-ff-1.)
            ! Schmit & Gibson 2011 - Isothermal Model 1.55MK
            Model%aa = 1.27*1.0d8         
            Model%bb = 22.8
            Model%cc = 5.97*1.0d8
            Model%dd = 16.0
            Model%ee = 1.57*1.0d8
            Model%ff = 3.87
        end subroutine setSchmitAtmosphere

        !> Set model atmosphere paramaters to zero so that there is 
        !> no outer solution pressure
        !> 
        subroutine setNoAtmosphere(Model)
            type(glModel_T), intent(inout) :: Model
            ! for 1d8 see last lines, within outarray in giblowprams.pro
            ! this is a block of atmosphere parameters outside of the bubble (we probably won't ever use these)
            ! as in, Pout=(aa/(bb+1.))*r^(-bb-1.)+(cc/(dd+1.))*r^(-dd-1.)+(ee/(ff+1.))*r^(-ff-1.)
            Model%aa = 0.0      
            Model%bb = 22.8
            Model%cc = 0.0
            Model%dd = 16.0
            Model%ee = 0.0
            Model%ff = 3.87
        end subroutine setNoAtmosphere
     
        !> Assume reasonable defaults if not overidden by user 
        !> for Gibson-Low model
        subroutine setModelDefaults(Model)
            type(glModel_T), intent(inout) :: Model

            ! ; 
            ! ; set up defaults if necessary
            ! ; as in Dove et al 2011, ao = 1.75
            ! ;
            Model%frontheight = 1.35    ! height of the CME front (in Rsun) *if* apar was 0
            Model%legsang = 25.         ! solid angle of the CME
            Model%apar = 0.05           ! transformation r -> r+a, which, e.g., converts a sphere to a teardrop
            Model%ao = 1.0              ! "the constant of proportionality between the stream function S and the pressure Pi
                                        ! It is later divided by muse (which is negative).
                                        ! so that the amplitude of Ao ends up being the same as the maximum amplitude of Br (bubble coords)."
            Model%alnotrbub = 1.        ! "the eigenvalue parameter of the twisted flux rope model"; "RESTRICTED TO EIGENVALUES OF J5/2
                                        ! BESSEL FUNCTION"; DEFAULT ALNOTRBUB=1 corresponds to smallest eigenvalue -- 5.763854

            Model%sigma = PI/2.         ! the rotation about the [?] surface normal, as in, like a tilt of an AR, default is pi/2
            Model%cmer = -PI/2.         ! Central meridian longitude
            Model%outScale = 0.1        ! "a parameter which scales the outer field solution since there is a discontinuity around bubble for field"
            Model%pio = 0.              ! constant for the background pressure
                                        

            Model%bonly = 0.            ! a flag to not calculate plasma parameters (dens, pres, etc), only the magnetic field
            Model%c1bonly = 0.5         ! c[1|2]bonly "only used if BONLY set; then creates a cavity which is C1BONLY fractional density
            Model%c2bonly = 0.9         ! (relative to background) and C2BONLY fractional height (relative to bubble)"
            Model%isothermal = 1.55e6   ! if bonly is set, the model temperature will have this value
            Model%bubbleonly = 0.       ! a flag to only calculate the model for interior of the bubble (and don't do anything for points outside of bubble)

            Model%alpha = 0.            ! acceleration (we probably won't ever make it non-zero)
            Model%velmult = 1.0         ! multiplier for 'velocity'
            Model%phiss = 1.0           ! center of the bubble - included for completeness: this is not an *input*, it's calculated from other parameters
            Model%time = 0.*3600.        ! Current Model time in seconds!
            Model%velimpose = 0.        ! if not zero, output model velocity is overriden with this value (we probably won't use it)
            Model%dt = 1.0              ! Model timestep - default to 1.0 second,
            Model%ts = 0
            Model%tfin = 100.           ! defaults to 100 seconds
            Model%alnotrbubsign = 1.0   ! Chirality
            
            Model%updateModelTime => updateGLTime

            ! Alternate default GLpars, keeping here for reference AJM - This removes background pressure        
            ! GLpars = (/1.35, 25., 0.05, 1.75, 1.0, 1.57080, 0.0, 0.0, 0.0, 22.8, &
            !         5.97, 16.0, 0.0, 3.87, 0.0, 0.5, 0.9, 1.55e6, 0.0, 0.0, &
            !         1.0, 1.0, 0.0, 1.0/)
        end subroutine setModelDefaults

        !> Allocate and Initalize Model parameters
        !> Calulates necessary geometric and scalar parameters
        !> for Gibson Low model from XML input
        subroutine initModel(Model, inpXML)
            type(glModel_T), intent(inout) :: Model
            type(XML_Input_T), intent(in) :: inpXML

            !Get RunID
            call inpXML%Set_Val(Model%RunID,'sim/runid',"Sim")
            call setModelDefaults(Model)
            call setNoAtmosphere(Model)

            ! Check State Type "Sphere" or 
            call inpXML%Set_Val(Model%CoordName, "prob/StateCoord", "Sphere")
            !read in GL parameters to override defaults
            call inpXML%Set_Val(Model%tfin, "prob/tfin", 100.)
            call inpXML%Set_Val(Model%dt, "prob/dt", 1.)
            call inpXML%Set_Val(Model%latitude, "prob/lat", 0.)
            call inpXML%Set_Val(Model%longitude, "prob/lon", 0.)
            call inpXML%Set_Val(Model%frontheight, "prob/frontheight", 1.35)
            call inpXML%Set_Val(Model%legsang, "prob/legsang", 25.)
            call inpXML%Set_Val(Model%topmorph, "prob/topmorph", 2.5)
            call inpXML%Set_Val(Model%bmax, "prob/Bmax", 0.001)
            call inpXML%Set_Val(Model%alpha, "prob/alpha", 0.)
            call inpXML%Set_Val(Model%apar, "prob/apar", 0.05)
            call inpXML%Set_Val(Model%sigma, "prob/orientation", 1.5708)
            call inpXML%Set_Val(Model%cmer, "prob/cmer", -1.5708)
            call inpXML%Set_Val(Model%alnotrbubsign, "prob/chiral", 1.0) 
            call inpXML%Set_Val(Model%vel_fh, "prob/vel_fh", 1.)
            call inpXML%Set_Val(Model%isLoud,'sim/isLoud',.true.)
            call inpXML%Set_Val(Model%isDebug,'sim/isDebug',.true.)
            call inpXML%Set_Val(Model%isTopomorph,'sim/isTopomorph',.false.)
            call inpXML%Set_Val(Model%isAtmosphere,'sim/isAtmosphere',.false.)
            call inpXML%Set_Val(Model%scaleBmax,'sim/scaleBmax',.false.)
            call inpXML%Set_Val(Model%Tstart_transient,'time/Tstart_transient',0.0)

            if (Model%isAtmosphere) call setSchmitAtmosphere(Model)

            Model%s_eta0 = sqrt(eta0)
            Model%s_eta = Model%vel_fh/(Model%frontheight*Rsolar) 
            call inpXML%Set_Val(Model%velmult, "prob/velmult", Model%s_eta/Model%s_eta0)

            Model%Tstart_transient = Model%Tstart_transient*3600. ! change from [hr] to [s]
            Model%k = tan(0.5*Model%legsang*mdtor)

            if (Model%isTopomorph) then ! Tethered CME
                Model%apar = Model%frontheight*(1. - 0.5*Model%k*(Model%topmorph - 2.))/(0.5*Model%k*Model%topmorph)
                Model%xo = (Model%frontheight + Model%apar)/(1.+tan(Model%legsang/2.*mdtor))
                Model%rbub = Model%frontheight - Model%xo + Model%apar
                Model%r0 = 2*Model%frontheight/Model%topmorph
            else ! Spheromak
                Model%xo = (Model%frontheight + Model%apar)/(1.+tan(Model%legsang/2.*mdtor))
                Model%rbub = Model%frontheight - Model%xo + Model%apar
                Model%r0 =  Model%xo - Model%apar
            end if

            Model%x0 = Model%r0/Model%k
            Model%alnotrbub_abs = abs(Model%alnotrbub)

            !Bessel eigenvalues
            if (abs(Model%alnotrbub_abs - alnotrbub0) .gt. 0.01) then
                Model%alnotrbub_abs = 1
            end if
            Model%alnotrbub_abs = int(Model%alnotrbub_abs)
            if ((Model%alnotrbub_abs .lt. 1) .or. (Model%alnotrbub_abs .gt. 6)) then
                if (Model%isLoud) write(*,*) 'warning: use eigenvalues between 0 and 6; resetting to 1'
                Model%alnotrbub_abs = 1
            end if
            
            Model%alnotrbubuse = alnotrbubuse1(int(Model%alnotrbub_abs))*Model%alnotrbubsign
            Model%alnot = Model%alnotrbubuse/Model%rbub
            Model%gfunot1 = sin(Model%alnotrbubuse)/Model%alnotrbubuse - cos(Model%alnotrbubuse)
            Model%Muse = 8.*pi*(Model%rbub*Model%rbub/3./Model%gfunot1 - 1./Model%alnot/Model%alnot)
            Model%ao = -Model%ao/Model%Muse
            Model%alpha = Model%alpha*Model%velmult*Model%velmult
            Model%eta = eta0*Model%velmult*Model%velmult

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
    
            if (Model%isLoud) then
                call printModelParameters(Model,"after Initialization")
            end if
            ! if (alpha .lt. 0) then
            !  print *, 'have not set up decelerating case'
            !  stop
            ! endif
        end subroutine initModel

        !> Reinitialize the CME model with scaling paramater 
        !> ao, from calcBmax
        !> 
        subroutine reInitModelScale(Model, State, Solution) 
            type(glModel_T), intent(inout) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            real(rp) :: ao 
            ao = calcModelBmax(Model, State, Solution, 100)
            if (Model%isLoud) write(*,"(1X,A14,2X,F)") "ao rescaled: ", ao
            if (Model%isLoud) write(*,"(1X,A14,2X,F)") "ao pre-scaling: ", -Model%ao*Model%Muse
            Model%ao = -ao/Model%Muse
        end subroutine

        !>  Allocate and Initialize Solution and State
        !>
        subroutine initSolutionState(Model, State, Solution, bounds)
            type(glModel_T), intent(inout) :: Model
            type(glSolution_T), intent(inout) :: Solution
            type(glState_T), intent(inout) :: State
            integer, optional, dimension(6), intent(in) :: bounds

            if(present(bounds)) then
                State%is = bounds(1)
                State%ie = bounds(2)
                State%js = bounds(3)
                State%je = bounds(4)
                State%ks = bounds(5)
                State%ke = bounds(6)
            else
                State%is = 1
                State%ie = State%Nip + 1
                State%js = 1
                State%je = State%Njp + 1
                State%ks = 1
                State%ke = State%Nkp + 1
            end if
            
            call allocState(State)
            call allocSolution(State, Solution)

            State%ri = 1 
        end subroutine initSolutionState

        !> Initialize Standalone Solution and State 
        !> grid, bounds, etc. 
        !>
        subroutine initStandaloneSolutionState(Model, State, Solution, inpXML)
            type(glModel_T), intent(inout) :: Model
            type(glSolution_T), intent(inout) :: Solution
            type(glState_T), intent(inout) :: State
            type(XML_Input_T), intent(in) :: inpXML

            real(rp) :: dr, dth, dphi
            real(rp) :: x, y, z, x1, x2, x3
            real(rp) :: rtpBds(6)
            integer :: i, j, k
            ! Set up State grid dimensions       
            call inpXML%Set_Val(State%Nip,"idir/N",Nc)
            call inpXML%Set_Val(State%Njp,"jdir/N",Nc)
            call inpXML%Set_Val(State%Nkp,"kdir/N",Nc)
            call inpXML%Set_Val(rtpBds(1),"idir/min",xMin)
            call inpXML%Set_Val(rtpBds(2),"idir/max",xMax)
            call inpXML%Set_Val(rtpBds(3),"jdir/min",xMin)
            call inpXML%Set_Val(rtpBds(4),"jdir/max",xMax)
            call inpXML%Set_Val(rtpBds(5),"kdir/min",xMin)
            call inpXML%Set_Val(rtpBds(6),"kdir/max",xMax)

            call initSolutionState(Model, State, Solution)
            ! Set up State
            ! Set dx's (uniform for now)
            dr = (rtpBds(2)-rtpBds(1))/State%Nip
            dth = (rtpBds(4)-rtpBds(3))/State%Njp
            dphi = (rtpBds(6)-rtpBds(5))/State%Nkp

            ! Create r(Ni), theta(Nj,Nk), phipb(Nj,Nk) arrays
              
            do k=State%ks, State%ke
                do j=State%js, State%je    
                    ! Phi
                    x3 = rtpBds(5)+(k-1)*dphi
                    State%phpb(j,k) = x3*2.*pi + Model%cmer       
                    ! Theta
                    x2 = rtpBds(3)+(j-1)*dth
                    State%thpb(j,k) = x2*pi                        
                    do i=State%is, State%ie    
                        ! Rho
                        x1 = rtpBds(1) + (i-1)*dr                      
                        State%r(i) = x1
                        ! construct xyz cartesian locations
                        x = x1*sin(x2*pi)*cos(x3*2.*pi)
                        y = x1*sin(x2*pi)*sin(x3*2.*pi)
                        z = x1*cos(x2*pi)
                        ! This grid is for external processing to set up cartesian coordinates 
                        ! for analyzing in visualization routines e.g. Paraview
                        State%xyz(i,j,k,XDIR:ZDIR) = [x,y,z]
                    end do 
                end do
            end do
        end subroutine initStandaloneSolutionState

        !> This init routine is to be used for standalone runs of the Gibson-Low Model
        !> on a State. Used mainly for testing purposes and to export h5 file of the 
        !> GL soluiton on a State for additional analysis or reading into an external 
        !> program. 
        subroutine initGLStandalone(Model, State, Solution, inpXML)
            type(glModel_T), intent(inout) :: Model
            type(glState_T), intent(inout) :: State
            type(glSolution_T), intent(inout) :: Solution
            type(XML_Input_T), intent(in) :: inpXML

            generateCMESolution => generateGLSolution

            call initModel(Model, inpXML)
            !Output/Restart (IOCLOCK)
            call Model%IO%init(inpXML,Model%time,Model%ts)

            call initStandaloneSolutionState(Model, State, Solution, inpXML)
            if (Model%scaleBmax) then 
                call reInitModelScale(Model, State, Solution)
            else
                Model%ao = Model%ao*Model%Bmax
            end if
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
            call setH5File(genName(Model%RunID, State%Nip, State%Njp, State%Nkp, 1, 1, 1))
        end subroutine initGLStandalone

        !> Initialize GL Types for use in external code 
        !> e.g. Gamera Helio
        !>
        !>
        !>
        subroutine initGLInterface(xyz, Model, State, Solution, inpXML, bounds)
            real(rp), dimension(:,:,:,:), intent(in) :: xyz
            class(glModel_T), intent(inout) :: Model
            class(glState_T), intent(inout) :: State
            class(glSolution_T), intent(inout) :: Solution
            type(XML_Input_T), intent(in) :: inpXML
            integer, dimension(6), intent(in) :: bounds

            generateCMESolution => generateGLSolution

            call initModel(Model, inpXML)
            call initSolutionState(Model, State, Solution, bounds)
            call setGLStateXYZ(xyz, Model, State)
            if (Model%scaleBmax) then 
                call reInitModelScale(Model, State, Solution)
            else
                Model%ao = Model%ao*Model%Bmax
            end if
            if (Model%isLoud) then
                call printModelParameters(Model, "after interface init finished.")
            end if
        end subroutine
        
        !> Expectation is that this interface is to be used by Gamera etc. 
        !> User must initalize the glApp components -> Model, State, State, Solution
        !> Model is initialized from inpXML
        !> Parameters:
        !>  xyz: the (i,j,k,xyz) locations over which to solve the Gibson Low Model
        !>  glApp: instance of the Gibson Low model app
        !>  inpXML: the inpXML for this model
        !> 
        subroutine setGLStateXYZ(xyz, Model,  State)
            type(glModel_T), intent(inout)  :: Model
            type(glState_T), intent(inout)  :: State
            real(rp), dimension(:,:,:,:), intent(in) :: xyz
            real(rp), dimension(:), allocatable :: r
            real(rp), dimension(:,:), allocatable :: theta, phi
            integer :: i

            allocate(r(State%is:State%ie))
            allocate(theta(State%js, State%je))
            allocate(phi(State%ks, State%ke))
            do i=State%is, State%ie
                r(i) = norm2(xyz(i,1,1,:))
            end do
            theta = acos(xyz(1,:,:,ZDIR)/norm2(xyz(1,1,1,:)))
            phi = atan2(xyz(1,:,:,YDIR),xyz(1,:,:,XDIR))

            call setGLStateRTP(r, theta, phi, Model, State)

        end subroutine setGLStateXYZ

        !> Expectation is that this interface is to be used by Gamera etc. 
        !> User must initalize the glApp components -> Model, State, State, Solution
        !> Model is initialized from inpXML
        !> State must minimally be defined by r(i), theta(j,k), phi(j,k) 
        !> 
        subroutine setGLStateRTP(r, theta, phi, Model, State)
            type(glModel_T), intent(inout) :: Model
            type(glState_T), intent(inout)  :: State
            real(rp), dimension(:), intent(in) :: r
            real(rp), dimension(:,:), intent(in) :: theta, phi
            ! integer :: rdim
            ! integer, dimension(2) :: adims

            ! rdim = size(r)
            ! adims = shape(theta)

            ! State%Ni = rdim
            ! State%Nj = adims(1)
            ! State%Nk = adims(2)
            State%r = r
            State%thpb = theta + Model%latitude
            State%phpb = phi - Model%longitude
        end subroutine setGLStateRTP
        
end module glinit