module siftypes

    use helpertypes
    use shellgrid
    use xml_input
    use ioclock


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
        
        !> These are all specified by sifconfig.h5
        character(len=strLen) :: name
        integer :: N
        integer :: flav
            !! Species flavor
        integer :: numNuc_p
            !! Number of protons in nucleus
        integer :: numNuc_n
            !! Number of neutrons in nucleus
        integer :: q
            !! Net charge of species
        real(rp) :: fudge
            !! Strong-scattering limit
        real(rp), dimension(:), allocatable :: alami
            !! Lambda channel cell interfaces/edges

        !> These are calculated after read-in
        integer :: kStart, kEnd
            !! Start and end indices for this species in Nk-size arrays
        real(rp) :: amu
            !! Species mass in amu (even the electrons)
        logical :: isElectron
            !! Is this species actually an electron? Determine on initialization and store it here

    end type SIFSpecies_T


    type sifModel_T

        ! Misc. bookkeeping stuff
        character(len=strLen) :: RunID = ""    
        character(len=strLen) :: configFName
            !! Filename of the .h5 config that holds lambda grid, wavemodel values, etc.
        character(len=strLen) :: SIFH5
            !! Filename of the h5 file we output to
        
        integer :: nSpc, nSpcMHD, nG
            !! Number of species in sif
            !! Number of species/fluids in MHD
            !! # ghosts
        real(rp) :: t0, tFin, dt
            !! Start and end time, delta coupling time (may be replaced/ignored later by voltron setting a dynamic coupling time)
        

        ! https://media0.giphy.com/media/XfDPdSRhYFUhIU7EPw/giphy.gif
        logical  :: fixedTimestep
            !! Fixed or dynamic timestep
        logical :: isMPI
        logical :: isRestart
        logical :: isLoud       ! For debug
        logical :: writeGhosts  ! For debug
        logical :: doClockConsoleOut
            !! If we are driving, output clock info
        logical :: doFatOutput
            !! Output extra 3D arrays

        ! Plasmasphere settings
        logical :: doPlasmasphere
            !! Use for now to determine if we should be doing plasmasphere stuff
            !! Likely, in the future, we will determine automatically by the presence of a flavor 0
        ! TODO: Extra params for refilling rate, determining initial profile, etc.

        ! Some constants
        real(rp) :: tiote  ! Ion temp over electron temp. In the future, should be fancier
        real(rp) :: worthyFrac  ! Fracton that a channel must contribute to pressure or density for its i shell to be evolved

        ! Lambda controls
        real(rp) :: kappa
            !! Kappa value, used in case of Kappa mapping from moments to eta channels
        logical :: doDynamicLambdaRanges
            !! Dynamic lambda ranges so we only evolve certain energies at certain L shells

        ! Detailed information
        type(planet_T) :: planet  
            !! Planet info like radius, mag. moment, etc.
        !type(precip_T) :: precip  ! Precipitation model info (Shanshan and Dong)
        !type(waveModel_T) :: wModel  ! Wave model info (Shanshan)

        character(len=strLen) :: icStr
        procedure(sifStateIC_T     ), pointer, nopass :: initState => NULL()
        procedure(sifUpdateV_T     ), pointer, nopass :: updateV   => NULL()
        !> TODO: Retire this and just have a single function with certain options like maxwellian or kappa
        procedure(sifDP2EtaMap_T   ), pointer, nopass :: dp2etaMap => NULL()

    end type sifModel_T


    type sifGrid_T
        integer :: gType  ! Enum of grid type

        type(ShellGrid_T) :: shGrid
        real(rp), dimension(:), allocatable :: delTh
            !! (Ni+1) Delta theta between cell centers. For cell i, delTh(i) = lower theta del, delTh(i+1) = higher theta del
        real(rp), dimension(:), allocatable :: delPh
            !! (Nj+1) Delta phi between cell centers. For cell j, delPh(j) = lower phi del, delPh(j+1) = higher phi del
        real(rp), dimension(:,:), allocatable :: Bmag
            !! Magnitude of B field [nT] at ionosphere

        integer :: nB ! Number of buffer cells between open region and active domain

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
        integer, dimension(:), allocatable :: k2spc
            !! Nk length mapping of k value to corresponding species index

    end type sifGrid_T


    type sifState_T
        real(rp) :: t, dt
            !! Current time and last coupling dt made
        real(rp), dimension(:), allocatable :: dtk
            !! Time step for every lambda channel
        real(rp) :: mjd
            !! Current mjd
        integer :: ts, tss
            !! Current coupling timestep and sub-stepping timestep
        type(IOClock_T) :: IO
            !! Timers for IO operations

        real(rp), dimension(:,:,:), allocatable :: eta
            !! (Ni, Nj, Nk) etas
        real(rp), dimension(:,:,:), allocatable :: pEff
            !! (Ni, Nj, Nk) Effective potential [V] (ExB + corot + gradient-curvature)
        logical, dimension(:,:), allocatable :: activeShells
            !! (Ni, Nk) I shells that should be evolved for a given lambda
        real(rp), dimension(:,:,:,:), allocatable :: iVel
            !! (Ni+1, Nj+1, Nk, 2) Edge-centered normal velocities
        real(rp), dimension(:,:,:,:), allocatable :: cVel
            !! (Ni, Nj, Nk, 2) Cell-centered velocities, [Rp/s] in ionosphere


        !> Variables coming from MHD flux tube tracing, size (Ni, Nj, Ns)
        real(rp), dimension(:,:,:), allocatable :: Pavg
        real(rp), dimension(:,:,:), allocatable :: Davg
        !> (Ni, Nk, NDIM), Bmin vector, [nT]
        real(rp), dimension(:,:,:), allocatable :: Bmin
        !> (Ni+1, Nk+1, NDIM) bMin xyz coordinates [Rx]
        real(rp), dimension(:,:,:), allocatable :: xyzMin
        !> (Ni, Nk, NDIM) cell-centered bMin xyz coordinates [Rx]
        real(rp), dimension(:,:,:), allocatable :: xyzMincc
        !> (Ni+1, Nk+1) corner values
        integer , dimension(:,:), allocatable :: topo    ! Topology (0=open, 1=closed)
        real(rp), dimension(:,:), allocatable :: thcon  ! Co-latitude  of conjugate points
        real(rp), dimension(:,:), allocatable :: phcon  ! Longitude of conjugate points
        !> (Ni, Nj)
        integer , dimension(:,:), allocatable :: active  ! (-1=inactive, 0=buffer, 1=active)
        integer , dimension(:,:), allocatable :: OCBDist  ! Cell distance from open-closed boundary
        real(rp), dimension(:,:), allocatable :: espot  ! electro-static potential from REMIX [kV]
        real(rp), dimension(:,:), allocatable :: bvol  ! Flux-tube volume [Rx/nT]

        !> Varibles coming from SIF, size (Ni, Nj, Nk)
        real(rp), dimension(:,:,:), allocatable :: precipFlux  ! Precipitation fluxes [!!units]
        real(rp), dimension(:,:,:), allocatable :: precipAvgE  ! Precipitation avg energy [!!units]
        ! (Ni, Nj, Nspc+1) (First index is bulk)
        ! Last dimension will be D/P of different populations (not necessarily same as species)
        ! Example: Total, hot protons, cold protons, electrons, other species
        real(rp), dimension(:,:,:), allocatable :: Den    ! Density  [amu/cc]
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

        subroutine sifUpdateV_T(Model,Grid,State)
            Import :: sifModel_T, sifGrid_T, sifState_T
            type(sifModel_T) , intent(in) :: Model
            type(sifGrid_T)  , intent(in) :: Grid
            type(sifState_T) , intent(inout) :: State
        end subroutine sifUpdateV_T

        function sifDP2EtaMap_T(Model,D,kT,vm,amin,amax) result (etaK)
            Import :: rp, sifModel_T
            type(sifModel_T), intent(in) :: Model
            real(rp), intent(in) :: D,kT,vm,amin,amax
            real(rp) :: etaK
        end function sifDP2EtaMap_T
    end interface




end module siftypes








