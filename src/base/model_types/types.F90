!Main data types

module types
    use gdefs
    use xml_input

    implicit none

    ! larger then 0 will write out the massage
    integer :: verbose = -1

!Ring information
!Holds info for ring-averaging about singularity
!Note, always assuming IDIR is ring direction
    type Ring_T
        integer :: PLDIR !Direction (IJK) of direction wrapping pole
        integer :: NumR !Number of rings to average over
        integer, dimension(:), allocatable :: NCh !Number of chunks to average over for each ring
        integer :: Np !Number of cells about pole (assuming power of 2)
        integer :: nSi,nSe !Start/End of slicing index
        character(len=strLen) :: GridID = "NONE" !Grid type: lfm,cyl,sph
        logical :: doS,doE !Whether to do +/- sides of singularity
    end type Ring_T


!Overall model information
    !Algorithmic/Run options
    type Model_T
        character(len=strLen) :: RunID
        integer :: verbose = -1
        integer :: nSpc, nDim, nG !Number of species, dimensions, ghosts
        real(rp) :: CFL, gamma
        real(rp) :: Vd0=0.5 !Coefficient for diffusive electric field, generally 0.5
        real(rp) :: t, tFin, dt
        real(rp) :: dtOut,tOut !Time between file outputs, time of next file output
        real(rp) :: dtRes,tRes !Same for restarts
        integer :: nOut=0,nRes=0 !Output numbering for H5/restarts 
        integer :: ts, tsOut !tsOut = steps between screen output
        logical :: doHall=.false., doMultiF=.false., &
                   doBackground=.false. , doMHD=.false., doArmor=.false.
        logical :: do25D =.false., doTimer=.false., doGrav = .false.
        logical :: doResOut = .false. !Output restarts
        logical :: isRestart=.false.
        logical :: doDivB=.false. !Output DivB
        logical :: useResistivity=.false.
        logical :: doForce=.false. !Use external force

        integer :: nTh=1 !Number of threads per node/group

        !Boris correction information
        logical :: doBoris,doHogs, doHeat,doPsphere
        real(rp) :: cHogH,cHogM !Coefficients for HOGS terms
        real(rp) :: Ca !Max speed for Boris correction
        real(rp) :: hRate,hTau  ! heating rate and cadence

        !Ring average information
        logical :: doRing = .false.
        type (Ring_T) :: Ring
        
        !MPI information
        !NumRi,NumRj,NumRk = # of ranks in each dimension
        !Ri,Rj,Rj = Placement of this rank
        logical :: isMPI = .false.
        integer :: NumRi=0,NumRj=0,NumRk=0
        integer :: Ri=0,Rj=0,Rk=0
        
        !Background field function pointer
        procedure(VectorField_T), pointer, nopass :: B0 => NULL()

        !Gravitational potential function pointer
        procedure(ScalarFun_T), pointer, nopass :: Phi => NULL()

        !User hack function pointers
        procedure(HackFlux_T), pointer, nopass :: HackFlux => NULL()
        procedure(HackE_T)   , pointer, nopass :: HackE    => NULL()
        procedure(HackStep_T), pointer, nopass :: HackStep => NULL()

    end type Model_T

    type bcHolder
        class(*), allocatable :: p
    end type bcHolder

!Overall grid information
    !Nip = # of physical cells, Ni = Nip+2*Ng
    !Sizes of data structures
    !Coords (x,y,z): (Ni+1,Nj+1,Nk+1)
    !Total coords (xyz): (Ni+1,Nj+1,Nk+1,NDIM)
    !Face-center coords (xfc): (Ni+1,Nj+1,Nk+1,NDIM,NDIM)
        !xfc(i,j,k,:,IDIR) = xyz coordinates of i,j,k I-Face
    !Cell-center edge lengths (di,dj,dk) = (Ni,Nj,Nk)
    !Cell volumes (volume) = (Ni,Nj,Nk)
    !Face ares (face) = (Ni+1,Nj+1,Nk+1,NDIM)
    !Edge lengths (edge) = (Ni,Nj,Nk,NDIM)
    !Face transforms (Tf) = (Ni+1,Nj+1,Nk+1,9,NDIM)
    !Edge velocity transforms (Te) = (Ni,Nj,Nk,9,NDIM)
    !Edge mag transforms (Teb) = (Ni,Nj,Nk,4,NDIM)
    type Grid_T
        integer :: Nip,Njp,Nkp
        integer :: Ni,Nj,Nk

        !Local indices, active region
        integer :: is,ie,js,je,ks,ke
        !Local indices, including ghosts
        integer :: isg,ieg,jsg,jeg,ksg,keg

        !Local indices to set domain for DT calculation
        integer :: isDT,jsDT,ksDT
        integer :: ieDT,jeDT,keDT

        !Local indices to set domain for Bxyz half-step calculation
        !Can be trimmed to set own Bxyz's
        integer :: isMG,jsMG,ksMG
        integer :: ieMG,jeMG,keMG
        
        !Shift between global and local indexes, always zero for
        !single process job
        integer ::ijkShift(3)

        !Corner-centered xyz-coordinates
        real(rp), dimension(:,:,:), allocatable :: x,y,z
        real(rp), dimension(:,:,:,:), allocatable :: xyz !4D array with all corners
        
        !Volume-centered xyz-coordinates (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), dimension(:,:,:,:), allocatable :: xyzcc
        
        !Face-centered coordinates
        !Ni,Nj,Nk,3,3
        !xfc(i,j,k,XDIR,IDIR) = x-coord of I-Face
        real(rp), dimension(:,:,:,:,:), allocatable :: xfc
        
        !TODO, combine these Courant # edge lengths
        real(rp), dimension(:,:,:), allocatable :: di,dj,dk
        !Cell volumes
        real(rp), dimension(:,:,:), allocatable :: volume

        !Transform matrices for face-centered vectors
        !Transform from xyz vector to face-centered triad
        !Ni+1,Nj+1,Nk+1,3*3,3
        real(rp), dimension(:,:,:,:,:), allocatable :: Tf
        
        !Transform matrices for edge-centered triad
        !Ni+1,Nj+1,Nk+1,3*3,3 (could be x,x,x,6,3)
        real(rp), dimension(:,:,:,:,:), allocatable :: Te
        
        !B field edge transform (xnqi,ynqi,xnqj,ynqj)
        !Ni,Nj,Nk,4,NDIM
        real(rp), dimension(:,:,:,:,:), allocatable :: Teb

        !Face areas
        real(rp), dimension(:,:,:,:), allocatable :: Face
        !Edge lengths in each direction, used for E field
        real(rp), dimension(:,:,:,:), allocatable :: Edge

        !Boundary conditions
        !Inner/outer i,j,k directions
        type(bcHolder) :: externalBCs(MAXBC)

        !Data for background field incorporation
        !fcB0(i,j,k,D,F) is the D component of the Bxyz on face F
        !edgB0(i,j,k,1:2,E) is the 1:2 component of b1/b2 on edge E
        !B0(i,j,k,d) is the d (xyz) component of B0
        !bFlux0(i,j,k,F) is the magnetic flux through face F
        real(rp), dimension(:,:,:,:,:), allocatable :: fcB0
        real(rp), dimension(:,:,:,:,:), allocatable :: edgB0
        real(rp), dimension(:,:,:,:), allocatable :: B0,bFlux0
        real(rp), dimension(:,:,:,:), allocatable :: dpB0 !Stress from background field

        real(rp), dimension(:,:,:,:), allocatable :: gxyz !Gravitational acceleration (cell-centered)
    end type Grid_T

!State information
    type State_T
        !Size (local grid) Ni,Nj,Nk,NVAR,(nSpc+1)
        real(rp), dimension(:,:,:,:,:), allocatable :: Gas
        !Size (local grid) Ni+1,Nj+1,Nk+1,nDim
        real(rp), dimension(:,:,:,:), allocatable :: magFlux 
        real(rp), dimension(:,:,:,:), allocatable :: Efld
        !TODO: Does this need to be in each state variable?  
        real(rp), dimension(:,:,:,:), allocatable :: Deta
        real(rp), dimension(:,:,:,:), allocatable :: Bxyz
        real(rp), dimension(:,:,:,:), allocatable :: eqMap   ! instantaneous equatorial mapping
        real(rp), dimension(:,:,:), allocatable :: eqPres    ! instantaneous equatorial pressure
        real(rp), dimension(:,:,:), allocatable :: eqDen    ! instantaneous equatorial density (for plasmasphere)
        
        real(rp) :: time
    end type State_T


    !StateIC_T
    !Generic initialization routine: ICs, Grid, Model
    abstract interface
        subroutine StateIC_T(Model,Grid,State,inpXML)
            Import :: Model_T, Grid_T, State_T, strLen, XML_Input_T
            type(Model_T), intent(inout) :: Model
            type(Grid_T), intent(inout) :: Grid
            type(State_T), intent(inout) :: State
            type(XML_Input_T), intent(in) :: inpXML
        end subroutine StateIC_T

    end interface

    !ScalarFun_T
    !Generic scalar function
    abstract interface
        subroutine ScalarFun_T(x,y,z,F)
            Import :: rp
            real(rp), intent(in) :: x,y,z
            real(rp), intent(out) :: F
        end subroutine ScalarFun_T
    end interface

    !VectorField_T
    !Generic vector field
    abstract interface
        subroutine VectorField_T(x,y,z,Ax,Ay,Az)
            Import :: rp
            real(rp), intent(in) :: x,y,z
            real(rp), intent(out) :: Ax,Ay,Az
        end subroutine VectorField_T
    end interface

    !GasIC_T
    !Generic (x,y,z) -> primitive hydro variables
    abstract interface
        subroutine GasIC_T(x,y,z,D,Vx,Vy,Vz,P)
            Import :: rp
            real(rp), intent(in) :: x,y,z
            real(rp), intent(out) :: D,Vx,Vy,Vz,P
        end subroutine GasIC_T
    end interface

    !BC_T
    !Boundary condition type
    !Contained by Grid structure, takes Model/Grid/State
    abstract interface
        subroutine BC_T(Model,Grid,State)
            Import :: Model_T, Grid_T, State_T
            type(Model_T), intent(in) :: Model
            type(Grid_T), intent(in) :: Grid
            type(State_T), intent(inout) :: State

        end subroutine BC_T
    end interface

    !HackFlux_T
    !User hack function for fluxes
    !Gives user opportunity to change gas/mag fluxes before application
    !Contained by Model structure, takes Model/Grid/gFlx/mFlx
    abstract interface
        subroutine HackFlux_T(Model,Gr,gFlx,mFlx)
            Import :: rp,NDIM,NVAR,BLK,Model_T, Grid_T
            type(Model_T), intent(in) :: Model
            type(Grid_T), intent(in) :: Gr
            real(rp), intent(inout) :: gFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
            real(rp), intent(inout), optional :: mFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM,1:NDIM)

        end subroutine HackFlux_T
    end interface

    !HackE_T
    !User hack electric field
    !Gives user opportunity to change State%Efld before application
    !Contained by Model structure, takes Model/Grid/State
    abstract interface
        subroutine HackE_T(Model,Grid,State)
            Import :: Model_T, Grid_T, State_T
            type(Model_T), intent(in) :: Model
            type(Grid_T), intent(in) :: Grid
            type(State_T), intent(inout) :: State

        end subroutine HackE_T
    end interface

    !HackStep_T
    !User-defined function to be called after update w/ latest state
    !Contained by model structure, takes Model/Grid/State and can edit Grid/State
    abstract interface
        subroutine HackStep_T(Model,Grid,State)
            Import :: Model_T, Grid_T, State_T
            type(Model_T), intent(in) :: Model
            type(Grid_T), intent(inout) :: Grid
            type(State_T), intent(inout) :: State

        end subroutine HackStep_T
    end interface

end module types
