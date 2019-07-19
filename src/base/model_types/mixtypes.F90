module mixtypes
  use mixdefs

  implicit none

  type mixParams_T
     ! constants to be read from XML input file

     ! conductance model
     integer  :: euv_model_type
     integer  :: et_model_type
     real(rp) :: alpha
     real(rp) :: beta
     real(rp) :: r
     real(rp) :: f107
     real(rp) :: pedmin 
     real(rp) :: hallmin
     real(rp) :: sigma_ratio
     real(rp) :: ped0
     logical  :: const_sigma
     logical  :: do_ramp
     logical  :: apply_cap 

     ! solver
     integer :: maxitr
     integer :: mr
     real(rp) :: tol_abs
     real(rp) :: tol_rel
     real(rp) :: llbc_value

     ! state
     real(rp) :: tilt
     integer :: hemisphere

     ! grid
     real(rp) :: LowLatBoundary
     integer :: Np, Nt

     ! init
     logical :: init_from_file

     ! IO
     real(rp) :: dtOut
  end type mixParams_T

  type mixState_T
     real(rp), dimension(:,:,:),allocatable :: Vars
     integer :: hemisphere=-1
     real(rp) :: tilt
  end type mixState_T

  type mixGrid_T
     integer :: Np, Nt
     
     real(rp), dimension(:,:), allocatable :: x,y,t,p,r
     real(rp), dimension(:,:), allocatable :: dt,dp
     real(rp), dimension(:,:), allocatable :: ft,fp
     real(rp), dimension(:,:), allocatable :: dtdt,dpdp
     real(rp), dimension(:,:), allocatable :: cosd      
     real(rp), dimension(:,:,:,:), allocatable :: Interpolant
     integer, dimension(:,:), allocatable :: mask
  end type mixGrid_T

  type Map_T
     real(rp), dimension(:,:,:), allocatable :: M  ! map matrix
     integer, dimension(:,:), allocatable :: I1, J1  ! indices of grid
                                                     ! cells on G1 to
                                                     ! which cells on
                                                     ! G2 belong where
                                                     ! G1(G2) are
                                                     ! grids to(from)
                                                     ! which we
                                                     ! interpolate
  end type Map_T

  type Solver_T
     real(rp), dimension(:), allocatable :: RHS,data
     integer,dimension(:), allocatable :: II,JJ,rowI
     real(rp), dimension(:,:), allocatable :: F11,F22,F12 ! (assuming F21=-F12)
     integer :: nnz  ! number of non-zeros
     ! low lat boundary array
     real(rp),dimension(:),allocatable :: LLBC 
     ! solution storage
     real(rp),dimension(:),allocatable :: solution
  end type Solver_T

  ! used to store an entire instance of MIX (e.g., one per hemisphere)
  type mixIon_T
     type(mixState_T)       :: St
     type(mixGrid_T)        :: G, Gfpd    ! G - primary MIX grid used for the solver; Mfpd -- flipped into Gamera space for interpolation
     type(mixParams_T)      :: P
     type(Solver_T)         :: S
     type(Map_T)            :: M, Mfpd    ! M is the map to map in MIX space; Mfpd is flipped to map in Gamera space
  end type mixIon_T

end module mixtypes
