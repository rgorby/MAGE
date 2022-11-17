module siftypes

    use helpertypes

    implicit none

    type kmUnits_T
        real(rp) :: V0   = 1.0  ! [m/s]
        real(rp) :: N0   = 1.0  ! [kg/m3]
        real(rp) :: T0   = 1.0  ! [s]
        real(rp) :: Eta0 = 1.0  ! [N/bVol]
        real(rp) :: Pot0 = 1.0  ! [Volts]
    end type kmUnits_T
    
    type sifModel_T
        character(len=strLen) :: RunID
        integer :: nSpc, nG, nB  ! Number of species, # ghosts, # coupling boundary cells
        real(rp) :: t0, tFin, dt
        logical  :: fixedTimestep  ! fixed or dynamic timestep

        ! Misc flags
        logical :: isMPI
        logical :: isRestart
        logical :: isLoud       ! For debug
        logical :: writeGhosts  ! For debug

        ! Plasmasphere settings
        logical :: doPlasmasphere
        ! Extra paramrs for efilling rate, determining initial profile, etc.

        ! Lambda bounds controls
        logical :: doDynamicLambdaRanges  ! Dynamic lambda ranges so we only evolve certain energies at certain L shells

        ! Grid settings
        integer  :: gType  ! Enum of grid type
        integer  :: Np, Nt, Nl ! # grid cells in phi, theta, and lambda channels
        real(rp) :: latBndL, latBndU  ! Upper and lower latitue boundaries

        ! Detailed information
        type(planet_T) :: planet  ! Planet info like radius, mag. moment, etc.
        !type(precip_T) :: precip  ! Precipitation model info 
        !type(waveModel_T) :: wModel  ! Wave model info (Shanshan)

    end type sifModel_T


    type sifGrid_T
        integer :: gType  ! Enum of grid type

        integer :: Nip, Njp, Nkp ! Number of physical cells, Ni = Nip+2*Ng
        integer :: Ni,Nj,Nk  ! Number of total cells
        integer :: Ng  ! Number of ghost cells

        integer :: is,ie,js,je,ks,ke  ! Active region indices
        integer :: isg,ieg,jsg,jeg,ksg,keg  ! Region including ghosts

        ! Bounds
        real(rp) :: latBndL, latBndU  ! Lower/upper latitude bound

        ! MPI things
        integer :: NumRk  ! Number of ranks in energy space
        integer :: Rk=0  ! Number of this rank

        ! Corner-centered lat/lon-coordinates (isg:ieg+1,jsg:jeg+1, 2)
        real(rp), dimension(:,:,:)  , allocatable :: latlon
        ! Volume-centered lat/lon-coordinates (isg:ieg,jsg:jeg, 2)
        real(rp), dimension(:,:,:)  , allocatable :: latloncc 
        ! Face-centered coordinates
        ! (Ni,Nj,2,2)  i, j, lat/lon, upper/lower
        ! llfc(i,j,k,XDIR,IDIR) = x-coord of I-Face
        real(rp), dimension(:,:,:,:), allocatable :: llfc

        ! (Nj) I of last active cell for each j slice 
        real(rp), dimension(:), allocatable :: iBnd

        ! Corner-centered lambda channel values
        real(rp), dimension(:), allocatable :: lami
        ! Cell-centered cell-centered lamba channel values
        real(rp), dimension(:), allocatable :: lamc

    end type sifGrid_T


    type sifState_T
        real(rp) :: time, MJD, dt

        ! (Ni, Nj, Nk) etas
        real(rp), dimension(:,:,:), allocatable :: eta
        ! (Ni+1, Nj+1, 2) Edge-centered normal velocities
        real(rp), dimension(:,:,:), allocatable :: iVel
        ! (Ni, Nj, 2) Cell-centered velocities
        real(rp), dimension(:,:,:), allocatable :: cVel
        ! (Ni, Nj, Ns, 2) Lower/upper lambda boudnaries for each species when doDynamicLambdaRanges=True
        real(rp), dimension(:,:,:,:), allocatable :: kBnds


        ! Variables coming from MHD flux tube tracing, size (Ni, Nj, Ns)
        real(rp), dimension(:,:,:), allocatable :: Pavg
        real(rp), dimension(:,:,:), allocatable :: Davg
        ! (Ni, Nk, NDIM)
        real(rp), dimension(:,:,:), allocatable :: Bmin
        real(rp), dimension(:,:,:), allocatable :: xyzMin
        ! (Ni, Nj)
        integer , dimension(:,:), allocatable :: topo    ! Topology (open, closed)
        integer , dimension(:,:), allocatable :: active  ! (active, buffer, inactive)
        real(rp), dimension(:,:), allocatable :: espot  ! electro-static potential from REMIX
        real(rp), dimension(:,:), allocatable :: latc  ! Latitude  of conjugate points
        real(rp), dimension(:,:), allocatable :: lonc  ! Longitude of conjugate points
        real(rp), dimension(:,:), allocatable :: bvol  ! Flux-tube volume [!!units]

        ! Varibles coming from RCM, size (Ni, Nj, Ns)
        real(rp), dimension(:,:), allocatable :: precipFlux  ! Precipitation fluxes [!!units]
        real(rp), dimension(:,:), allocatable :: precipAvgE  ! Precipitation avg energy [!!units]
        ! (Ni, Nj, Ns+?)
        ! Last dimension will be D/P of different populations (not necessarily same as species)
        ! Example: Total, hot protons, cold protons, electrons, other species
        real(rp), dimension(:,:,:), allocatable :: Den    ! Density  [!!units]
        real(rp), dimension(:,:,:), allocatable :: Press  ! Pressure [!!units]


    end type sifState_T


end module siftypes








