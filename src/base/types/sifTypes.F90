module siftypes

    use helpertypes
    use shellgrid
    use xml_input

    use sifdefs

    implicit none

    !  var*Var0 = CODE units --> In/Out units
    type kmUnits_T
        real(rp) :: V0   = 1.0  ! [m/s]
        real(rp) :: N0   = 1.0  ! [kg/m3]
        real(rp) :: T0   = 1.0  ! [s]
        real(rp) :: Eta0 = 1.0  ! [#/(Bs*Rx^2)] --> [#/Wb]
        real(rp) :: Pot0 = 1.0  ! [Volts]
    end type kmUnits_T


    type SIFSpecies_T
        !! Container for all info related to a specific species
        
        character(len=strLen) :: name
        integer :: N
        integer :: flav
        integer :: kStart, kEnd
        real(rp) :: fudge
        real(rp), dimension(:), allocatable :: alami
            !! Lambda channel cell interfaces/edges

    end type SIFSpecies_T

    type SIFIO_T
        character(len=strLen) :: SIFH5
        integer :: nOut = 0  ! Next output time
        real(rp) :: intScl = 1.0  ! Scaling for intensisty [TODO]
    end type SIFIO_T


    type sifModel_T

        ! Misc. bookkeeping stuff
        character(len=strLen) :: RunID = ""    
        character(len=strLen) :: configFName
            !! Filename of the .h5 config that holds lambda grid, wavemodel values, etc.
        type(SIFIO_T) :: SIFIO
        
        integer :: nSpc, nG, nB
            !! Number of species, # ghosts, # coupling boundary cells
        real(rp) :: t0, tFin, dt
            !! Start and end time, delta coupling time (may be replaced/ignored later by voltron setting a dynamic coupling time)
        logical  :: fixedTimestep
            !! Fixed or dynamic timestep

        ! https://media0.giphy.com/media/XfDPdSRhYFUhIU7EPw/giphy.gif
        logical :: isMPI
        logical :: isRestart
        logical :: isLoud       ! For debug
        logical :: writeGhosts  ! For debug

        ! Plasmasphere settings
        logical :: doPlasmasphere
            !! Use for now to determine if we should be doing plasmasphere stuff
            !! Likely, in the future, we will determine automatically by the presence of a flavor 0
        ! TODO: Extra params for efilling rate, determining initial profile, etc.

        ! Lambda bounds controls
        logical :: doDynamicLambdaRanges
            !! Dynamic lambda ranges so we only evolve certain energies at certain L shells

        ! Detailed information
        type(planet_T) :: planet  
            !! Planet info like radius, mag. moment, etc.
        !type(precip_T) :: precip  ! Precipitation model info (Shanshan and Dong)
        !type(waveModel_T) :: wModel  ! Wave model info (Shanshan)


        procedure(sifStateIC_T), pointer, nopass :: initState => NULL()
        procedure(sifUpdateV_T), pointer, nopass :: updateV => NULL()
        
        ! Hack functions
        ! e.g. hackDP2Eta

    end type sifModel_T


    type sifGrid_T
        integer :: gType  ! Enum of grid type

        type(ShellGrid_T) :: shGrid

        ! Flags
        logical :: ignoreConfigMismatch
            !! In the case that the config file has more species than Model%nSpc,
            !! SIF will complain and die if this is false, but will carry on if true

        ! MPI things
        integer :: NumRk  ! Number of ranks in energy space
        integer :: Rk=0   ! Number of this rank

        ! Face-centered coordinates
        !!! Maybe don't need because shGrid has centers in lat and lon directions already
        ! (Np/Nt,2,2)  i/j, lon/lat dir, upper/lower
        ! llfc(i,j,XDIR,IDIR) = x-coord of I-Face
        !real(rp), dimension(:,:,:,:), allocatable :: latfc
        !real(rp), dimension(:,:,:,:), allocatable :: lonfc

        ! (Nj) I of last active cell for each j slice
        integer, dimension(:), allocatable :: iBnd

        ! Species / lambda stuff
        integer :: Nk  ! Total number of channels for all species
        integer :: nSpc  ! Model has the main copy of this, but helpful to keep here too
        type(SIFSpecies_T), dimension(:), allocatable :: spc
            !! Collection of SIFSpecies that contain all relevant species info, including alami
        real(rp), dimension(:), allocatable :: alamc
            !! Cell-centered lamba channel values
        real(rp), dimension(:,:,:,:), allocatable :: kBnds
            !! (Ni, Nj, Ns, 2) Lower/upper lambda boudnaries for each species when doDynamicLambdaRanges=True

    end type sifGrid_T


    type sifState_T
        real(rp) :: time, MJD, dt

        real(rp), dimension(:,:,:), allocatable :: eta
            !! (Ni, Nj, Nk) etas
        real(rp), dimension(:,:,:,:), allocatable :: iVel
            !! (Ni+1, Nj+1, Nk, 2) Edge-centered normal velocities
        real(rp), dimension(:,:,:,:), allocatable :: cVel
            !! (Ni, Nj, Nk, 2) Cell-centered velocities


        !> Variables coming from MHD flux tube tracing, size (Ni, Nj, Ns)
        real(rp), dimension(:,:,:), allocatable :: Pavg
        real(rp), dimension(:,:,:), allocatable :: Davg
        !> (Ni, Nk, NDIM)
        real(rp), dimension(:,:,:), allocatable :: Bmin  ! [nT]
        real(rp), dimension(:,:,:), allocatable :: xyzMin ! [Rx]
        !> (Ni, Nj)
        integer , dimension(:,:), allocatable :: topo    ! Topology (0=open, 1=closed)
        integer , dimension(:,:), allocatable :: active  ! (-1=inactive, 0=buffer, 1=active)
        real(rp), dimension(:,:), allocatable :: espot  ! electro-static potential from REMIX [kV]
        real(rp), dimension(:,:), allocatable :: latc  ! Latitude  of conjugate points
        real(rp), dimension(:,:), allocatable :: lonc  ! Longitude of conjugate points
        real(rp), dimension(:,:), allocatable :: bvol  ! Flux-tube volume [Rx/nT]

        !> Varibles coming from SIF, size (Ni, Nj, Nk)
        real(rp), dimension(:,:,:), allocatable :: precipFlux  ! Precipitation fluxes [!!units]
        real(rp), dimension(:,:,:), allocatable :: precipAvgE  ! Precipitation avg energy [!!units]
        ! (Ni, Nj, Ns)
        ! Last dimension will be D/P of different populations (not necessarily same as species)
        ! Example: Total, hot protons, cold protons, electrons, other species
        real(rp), dimension(:,:,:), allocatable :: Den    ! Density  [#/cc]
        real(rp), dimension(:,:,:), allocatable :: Press  ! Pressure [nPa]
        real(rp), dimension(:,:,:), allocatable :: vAvg   ! Average cell velocity [km/s]


    end type sifState_T


!------
! Higher-level types, using above types
!------

    type sifApp_T
        type(sifModel_T) :: Model
        type(sifGrid_T ) :: Grid
        type(sifState_T) :: State
    end type sifApp_T

    abstract interface
        subroutine sifStateIC_T(Model,Grid,State,inpXML)
            Import :: sifModel_T, sifGrid_T, sifState_T, strLen, XML_Input_T
            type(sifModel_T) , intent(in) :: Model
            type(sifGrid_T)  , intent(in) :: Grid
            type(sifState_T) , intent(inout) :: State
            type(XML_Input_T), intent(in) :: inpXML
        end subroutine sifStateIC_T

    end interface

    abstract interface
        subroutine sifUpdateV_T(Model,Grid,State)
            Import :: sifModel_T, sifGrid_T, sifState_T
            type(sifModel_T) , intent(in) :: Model
            type(sifGrid_T)  , intent(in) :: Grid
            type(sifState_T) , intent(inout) :: State
        end subroutine sifUpdateV_T

    end interface


end module siftypes








