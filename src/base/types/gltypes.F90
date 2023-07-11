
!> Gibson-Low Model Types
!>
module gltypes
    use gldefs
    use cmetypes

    implicit none

    !Unit information
    !Holds information about GL MHD scaling & how to scale input values
    ! Not currently used, model is in cgs by default
    type glUnits_T
        real(rp) :: glx0 = 1.0 ! [m]
        real(rp) :: glv0 = 1.0 ! [m/s]
        real(rp) :: glD0 = 1.0 ! [kg/m3]
        real(rp) :: glP0 = 1.0 ! [Pa]
        real(rp) :: glB0 = 1.0 ! [T]
        real(rp) :: glT0 = 1.0 ! [s]
    end type glUnits_T

    !Scaling for output data
    type glOut_T
        real(rp) :: tScl=1.0,dScl=1.0,vScl=1.0,pScl=1.0,bScl=1.0
        real(rp) :: eScl=1.0,jScl=1.0
        character(len=strLen) :: uID = "CODE" !Overall units ID
        character(len=strLen) :: tID = "CODE TIME"
        character(len=strLen) :: dID = "CODE DEN"
        character(len=strLen) :: vID = "CODE VEL"
        character(len=strLen) :: pID = "CODE PRESSURE"
        character(len=strLen) :: bID = "CODE BFIELD"
        character(len=strLen) :: eID = "CODE EFIELD"
        character(len=strLen) :: jID = "CODE CURRENT"
    end type glOut_T

    !> Overall model information
    !> Algorithmic/Run options
    type, extends(baseCMEModel_T) :: glModel_T
        real(rp) :: alnotrbubuse
        real(rp) :: gfunot1
        real(rp) :: alnotrbub_abs  !to preserve sign of alnotrbub
        real(rp) :: alnotrbubsign
        real(rp) :: alnotrbub
        ! Morphology/position
        real(rp) :: frontheight
        real(rp) :: legsang
        real(rp) :: topmorph
        real(rp) :: sigma
        real(rp) :: cmer
        real(rp) :: Bmax
        real(rp) :: cmeV
        real(rp) :: vel_fh
        real(rp) :: frontvolume
        real(rp) :: k 
        real(rp) :: apar 
        real(rp) :: r0
        real(rp) :: x0
        real(rp) :: ao
        real(rp) :: bmagmax
        real(rp) :: alnot
        real(rp) :: muse
        real(rp) :: outScale
        !intermediate parameters
        real(rp) :: xo
        real(rp) :: rbub
        real(rp) :: anot
        real(rp) :: eta
        real(rp) :: phiss
        real(rp) :: pio
        ! Physical Properties
        real(rp) :: aa
        real(rp) :: bb
        real(rp) :: cc
        real(rp) :: dd
        real(rp) :: ee
        real(rp) :: ff
        real(rp) :: bonly
        real(rp) :: c1bonly
        real(rp) :: c2bonly
        real(rp) :: isothermal
        real(rp) :: bubbleonly
        !temporal properties
        real(rp) :: alpha
        real(rp) :: velmult
        real(rp) :: tFin
        real(rp) :: dt
        real(rp) :: velimpose
        real(rp) :: s_eta0
        real(rp) :: s_eta
        character(len=strLen)  :: CoordName
        !> Whether precheck logic should run to calculate parameters
        logical :: isTopomorph
        ! Whether an atmosphere exists outside the bubble
        logical :: isAtmosphere
        ! Scale Bmax - whether a preconditioning is run to find Bmax
        logical :: scaleBmax
        !Unit information
        type (glUnits_T) :: Units
        type (glOut_T)   :: glOut
    end type glModel_T

    !> State information
    type, extends(baseCMEState_T) :: glState_T

        ! Used for standalone Gibson-Low Model run
        ! to define user specified 3D grid dimensions to supply r,theta,phi to 
        ! construct state and to define size of solution ranks
        integer :: Nip,Njp,Nkp,is,js,ks,ie,je,ke

        ! Solution State Radii index (r(i))
        integer :: ri
        ! Count of inside/outside bubble locations 
        real(rp) :: outside_count, inside_count
        ! Coordinate values
        real(rp), dimension(:, :, :, :), allocatable :: xyz ! For Standalone Grid Output
        real(rp), dimension(:), allocatable :: r ! radii of spherical shells to pass to GL - copied to rpb for GL solution
        real(rp), dimension(:, :), allocatable :: rpb, thpb, phpb ! input to giblow_coords
        real(rp), dimension(:, :), allocatable :: rout, thout, phout
        real(rp), dimension(:, :), allocatable :: rcap, thcap, phcap
        real(rp), dimension(:, :), allocatable :: rsquig, rlam, xtr, ytr, ztr, F
        real(rp), dimension(:, :), allocatable :: xtilde, ytilde, ztilde, rat
        real(rp), dimension(:, :), allocatable :: rtilde, thtilde, phtilde
        ! Field Values
        real(rp), dimension(:, :), allocatable :: glpi
        ! Inside Solution Values
        real(rp), dimension(:, :), allocatable :: bpresin, presin, presin_0 
        real(rp), dimension(:, :), allocatable :: pbackin
        real(rp), dimension(:, :), allocatable :: dbackin, densin
        ! Background and Total Values
        real(rp), dimension(:, :), allocatable :: densback, DensbackHEonly
        real(rp), dimension(:, :), allocatable :: ptot, presback
        real(rp), dimension(:, :), allocatable :: bmag
        ! Outside Solution Values
        real(rp), dimension(:, :), allocatable :: streamout, densout, presout
        real(rp), dimension(:, :), allocatable :: brlambout, bthlambout, bphlambout
        real(rp), dimension(:, :), allocatable :: tderivout, tderivR, tderivmu, tderiv
        real(rp), dimension(:, :), allocatable :: jrlambout, jthlambout, jphlambout
        logical, dimension(:, :), allocatable :: cavinside
        ! Physical Solution Values
        real(rp), dimension(:, :), allocatable :: blittlerlamb, blittlethlamb, blittlephlamb
        real(rp), dimension(:, :), allocatable :: jlittlerlamb, jlittlethlamb, jlittlephlamb
        real(rp), dimension(:, :), allocatable :: stream
        real(rp), dimension(:, :), allocatable :: bstrengthphys

        contains
            procedure :: updateGrid => setGLStateXYZ
    end type glState_T

    !> Solution information
    type, extends(baseCMESolution_T) :: glSolution_T
        ! solution bubble r(i) < rbub = Model%frontheight - Model%xo + Model%apar
        integer :: CoordSystem
    end type glSolution_T

    contains
        !>
        !>
        subroutine updateGLTime(Model, time, gT0)
            type(glModel_T) :: Model       
            real(rp), intent(in) :: time, gT0
            Model%time = (time - Model%Tstart_transient)*gT0
            if (Model%isLoud .and. (Model%isDebug)) then
                write(*,"(1X,A14,2X,F13.2)") "Sim time: ", Model%time 
                write(*,"(1X,A14,2X,F13.2)") "Updated time: ", Model%time 
            end if
        end subroutine updateGLTime


        !> Expectation is that this interface is to be used by Gamera etc. 
        !> User must initalize the glApp components -> Model, State, State, Solution
        !> Model is initialized from inpXML
        !> State must minimally be defined by r(i), theta(j,k), phi(j,k) 
        !> 
        subroutine setGLStateRTP(r, theta, phi, Model, State)
            real(rp), dimension(:), intent(in) :: r
            real(rp), dimension(:,:), intent(in) :: theta, phi
            class(*), intent(inout) :: Model
            class(glState_T), intent(inout)  :: State

            select type (Model)
                type is (glModel_T)
                    State%r = r
                    State%thpb = theta + State%currentLatitude
                    State%phpb = phi - State%currentLongitude
            end select

        end subroutine setGLStateRTP

        !> Expectation is that this interface is to be used by Gamera etc. 
        !> User must initalize the glApp components -> Model, State, State, Solution
        !> Model is initialized from inpXML
        !> Parameters:
        !>  State: this instance of the CME State
        !>  xyz: the (i,j,k,xyz) locations over which to solve the Gibson Low Model
        !>  Model: the CME model object
        !> 
        subroutine setGLStateXYZ(State, xyz, Model)
            class(glState_T), intent(inout)  :: State
            real(rp), dimension(:,:,:,:), intent(in) :: xyz
            class(*), intent(inout)  :: Model

            real(rp), dimension(:), allocatable :: r
            real(rp), dimension(:,:), allocatable :: theta, phi
            integer :: i

            allocate(r(State%is:State%ie))
            allocate(theta(State%js:State%je, State%ks:State%ke))
            allocate(phi(State%js:State%je, State%ks:State%ke))
            !write(*,*) "Allocated r, theta, phi"
            do i=State%is, State%ie
                r(i) = norm2(xyz(i,State%js,State%ks,:))
            end do
            theta = acos(xyz(State%is,:,:,ZDIR)/norm2(xyz(State%is,State%js,State%ks,:)))
            phi = atan2(xyz(State%is,:,:,YDIR),xyz(State%is,:,:,XDIR))
            call setGLStateRTP(r, theta, phi, Model, State)
            !write(*,*) "Set r, theta, phi"

        end subroutine setGLStateXYZ

end module gltypes