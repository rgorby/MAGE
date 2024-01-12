module raijutypes

    use helpertypes
    use shellgrid
    use xml_input
    use ioclock
    use kronos

    use raijudefs

    use raijudefs

    implicit none

    !  var*Var0 = CODE units --> In/Out units
    type kmUnits_T
        real(rp) :: V0   = 1.0  ! [m/s]
        real(rp) :: N0   = 1.0  ! [kg/m3]
        real(rp) :: T0   = 1.0  ! [s]
        real(rp) :: Eta0 = 1.0  ! [#/(Bs*Rx^2)] --> [#/Wb]
        real(rp) :: Pot0 = 1.0  ! [Volts]
    end type kmUnits_T

!------
! Species
!------
    type raijuSpecies_T
        !! Container for all info related to a specific species
        
        ! These are all specified by raijuconfig.h5
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
        logical :: mapExtraToPsph
            !! Whether any eta under species' lowest bound gets added to plasmasphere

        ! These are calculated after read-in
        integer :: kStart, kEnd
            !! Start and end indices for this species in Nk-size arrays
        real(rp) :: amu
            !! Species mass in amu (even the electrons)
        integer :: spcType
            !! Enum of species type
        
    end type raijuSpecies_T

!------
! Precipitation models
!------
    type eLossWM_T
        !! Parameters used in electron wave model from Dedong Wang and Shanshan Bao
        
        ! -- Model -- !
        real(rp) :: NpsphHigh = def_NpsphHigh
        real(rp) :: NpsphLow = def_NpsphLow
        real(rp) :: ChorusLMax = def_ChorusLmax
        real(rp) :: PsheetLMin = def_PsheetLmin
        real(rp) :: ChorusEMin = def_ChorusEMin

        logical :: doOutput = .false.
            !! Whether or not we will be asked to provide detailed info in the output file
        
        
        type(TimeSeries_T) :: KpTS
            !! Kp data from wind file

        ! -- Grid -- !
        ! Chorus info
        integer :: Nkp, Nmlt, Nl, Ne
            !! Number of bins for Kp, MLT, L shell, and Energy
        real(rp), allocatable, dimension(:) :: Kp1D
            !! 1D array of Kp dimension for Tau4D
        real(rp), allocatable, dimension(:) :: MLT1D
            !! 1D array of MLT dimension for Tau4D
        real(rp), allocatable, dimension(:) :: L1D
            !! 1D array of L shell dimension for Tau4D [Re]
        real(rp), allocatable, dimension(:) :: Energy1D
            !! 1D array of energy dimension for Tau4D [MeV]
        real(rp), allocatable, dimension(:,:,:,:) :: Tau4D
            !! Tau(Kp, MLT, L, E) table electron scattering table [seconds]

        ! -- State -- !
        real(rp), allocatable, dimension(:,:) :: wPS
        real(rp), allocatable, dimension(:,:) :: wHISS
        real(rp), allocatable, dimension(:,:) :: wCHORUS

    end type eLossWM_T

!------
! Main Model, Grid, State
!------
    type raijuModel_T

        ! Misc. bookkeeping stuff
        character(len=strLen) :: RunID = ""    
        character(len=strLen) :: configFName
            !! Filename of the .h5 config that holds lambda grid, wavemodel values, etc.
        character(len=strLen) :: raijuH5
            !! Filename of the h5 file we output to
        
        integer :: nSpc, nSpcMHD, nG
            !! Number of species in raiju
            !! Number of species/fluids in MHD
            !! # ghosts
        real(rp) :: t0, tFin, dt
            !! Start and end time, delta coupling time (may be replaced/ignored later by voltron setting a dynamic coupling time)
        
        ! Solver params
        integer :: maxItersPerSec
        real(rp) :: CFL
            !! CFL condition used in deciding time step
        integer :: maxOrder
            !! Max reconstruction order we will use when calculating etas at interfaces
        real(rp) :: PDMB
            !! 0 < PDMB < 1 used in PDM limiter. 1 = very diffusive

        ! https://media0.giphy.com/media/XfDPdSRhYFUhIU7EPw/giphy.gif
        logical :: isSA
            !! Is RAIJU running in statnd-alone mode. If so, we shouldn't expect any coupling information to exist
        logical :: fixedTimestep
            !! Fixed or dynamic timestep
        logical :: isMPI
        logical :: isRestart
        logical :: doLosses  ! Whether or not to calculate eta losses during advance
        logical :: isLoud       ! For debug
        logical :: writeGhosts  ! For debug
        logical :: doClockConsoleOut
            !! If we are driving, output clock info
        logical :: doFatOutput
            !! Output extra 3D arrays
        logical :: doDebugOutput
            !! Dump lots of otherwise unnecessary stuff
        logical :: doGeoCorot
            !! If true, calc corotation potential from Geopack
            !! If false, calc corotation potential assuming dipole and rotational axes are aligned
        logical :: doExcesstoPsph
            !! Allow mapping of excess H+ to plasmasphere channel

        ! Plasmasphere settings
        logical :: doPlasmasphere
            !! Use for now to determine if we should be doing plasmasphere stuff
            !! Likely, in the future, we will determine automatically by the presence of a flavor 0
        ! TODO: Extra params for refilling rate, determining initial profile, etc.

        ! Some constants
        real(rp) :: tiote  ! Ion temp over electron temp. In the future, should be fancier

        ! Active shell settings
        logical :: doActiveShell
            !! Use activeShell logic to try to boost dt
        real(rp) :: worthyFrac  
            !! Fracton that a channel must contribute to pressure or density for its i shell to be evolved

        ! Lambda controls
        real(rp) :: kappa
            !! Kappa value, used in case of Kappa mapping from moments to eta channels
        logical :: doDynamicLambdaRanges
            !! Dynamic lambda ranges so we only evolve certain energies at certain L shells

        ! Detailed information
        type(planet_T) :: planet  
            !! Planet info like radius, mag. moment, etc.

        ! Losses
        logical :: doSS, doCC, doCX, doFLC
            !! (Ions) Do strong scattering / coulomb collisions / charge exchange / field-line curvature
        integer :: eLossModel
            !! Enumerator indicating active electon loss model
        procedure(raijuELossRate_T), pointer, nopass :: eLossRateFn => NULL()
            !! Pointer to electron loss function
        type(eLossWM_T) :: eLossWM
            !! Container for electron Wave Model data

        character(len=strLen) :: icStr
        procedure(raijuStateIC_T     ), pointer, nopass :: initState => NULL()
        procedure(raijuUpdateV_T     ), pointer, nopass :: updateV   => NULL()
        !> TODO: Retire this and just have a single function with certain options like maxwellian or kappa
        procedure(raijuDP2EtaMap_T   ), pointer, nopass :: dp2etaMap => NULL()

    end type raijuModel_T


    type raijuGrid_T
        integer :: gType  ! Enum of grid type

        type(ShellGrid_T) :: shGrid
        real(rp), dimension(:), allocatable :: delTh
            !! (Ngi+1) [radians] Delta theta between cell centers. For cell i, delTh(i) = lower theta del, delTh(i+1) = higher theta del
        real(rp), dimension(:), allocatable :: delPh
            !! (Ngj+1) [radians] Delta phi between cell centers. For cell j, delPh(j) = lower phi del, delPh(j+1) = higher phi del
        real(rp), dimension(:,:), allocatable :: areaCC
            !! (Ngi, Ngj) [Rp^2] Area of each cell
        real(rp), dimension(:,:,:), allocatable :: areaFace
            !! (Ngi+1, Ngj+1, 2) [Rp^2] Estimated area at each interface. (i,j,1) = theta dir, (i,j,2) = phi dir
        real(rp), dimension(:,:,:), allocatable :: lenFace
            !! (Ngi+1, Ngj+1, 2) [Rp] arc length of each face
        real(rp), dimension(:,:), allocatable :: Bmag
            !! (Ngi, Ngj) [nT] Magnitude of B field at ionosphere (cell-centered)
        real(rp), dimension(:,:), allocatable :: cosdip
            !! (Ngi) Cosine of the dip angle at ionosphere (cell-centered)
        real(rp), dimension(:,:,:), allocatable :: Br
            !! (Ngi+1, Ngj+1, 2) [nT] Radial/normal component of B field [nT] at ionosphere (cell faces)

        integer :: nB ! Number of buffer cells between open region and active domain

        ! Flags
        logical :: ignoreConfigMismatch
            !! In the case that the config file has more species than Model%nSpc,
            !! raiju will complain and die if this is false, but will carry on if true

        ! MPI things
        integer :: NumRk  ! Number of ranks in energy space
        integer :: Rk=0   ! Number of this rank

        ! (Nj) I of last active cell for each j slice
        integer, dimension(:), allocatable :: iBnd

        ! Species / lambda stuff
        integer :: Nk  ! Total number of channels for all species
        integer :: nSpc  ! Model has the main copy of this, but helpful to keep here too
        type(raijuSpecies_T), dimension(:), allocatable :: spc
            !! Collection of raijuSpecies that contain all relevant species info, including alami
        real(rp), dimension(:), allocatable :: alamc
            !! Cell-centered lamba channel values
        integer, dimension(:), allocatable :: k2spc
            !! Nk length mapping of k value to corresponding species index

    end type raijuGrid_T


    type raijuState_T
        real(rp) :: t, dt
            !! Current time and last coupling dt made
        real(rp), dimension(:), allocatable :: dtk
            !! Time step for every lambda channel
        integer, dimension(:), allocatable :: nStepk
            !! Number of steps each channel has been evolved
        real(rp) :: mjd
            !! Current mjd
        integer :: ts, tss
            !! Current coupling timestep and sub-stepping timestep
        type(IOClock_T) :: IO
            !! Timers for IO operations

        ! -- Solver values -- !
        real(rp), dimension(:,:,:), allocatable :: eta
            !! (Ngi, Ngj, Nk) [#/cc * Rp/nT] etas
        real(rp), dimension(:,:,:), allocatable :: eta_half
            !! (Ngi, Ngj, Nk) [#/cc * Rp/nT] etas  0.5*dt forward from t
        real(rp), dimension(:,:,:), allocatable :: eta_last
            !! (Ngi, Ngj, Nk) [#/cc * Rp/nT] etas from previous time step, used for halt-time calculation
        logical, dimension(:,:), allocatable :: activeShells
            !! (Ngi, Nk) i shells that should be evolved for a given lambda
        real(rp), dimension(:,:,:), allocatable :: pEff
            !! (Ngi+1, Ngj+1, Nk) [V] Effective potential (ExB + corot + gradient-curvature)
            !! Not actually used in calculations, but helpful for output
        real(rp), dimension(:,:,:), allocatable :: gradPotE, gradPotCorot, gradVM
            !! (Ngi+1, Ngj+1,2) [V/m] Th/phi gradient of the ionospheric potential, corotation potential, and flux tube volume raised to -2/3
            !! units of gradVM are [FTV^(2/3)/m]. Multiplying by lambda gives units of [V/m]
        real(rp), dimension(:,:,:,:), allocatable :: iVel
            !! (Ngi+1, Ngj+1, Nk, 2) [m/s] Edge-centered normal velocities
        real(rp), dimension(:,:,:,:), allocatable :: cVel
            !! (Ngi, Ngj, Nk, 2) [m/s] Cell-centered velocities


        ! -- Variables coming from MHD flux tube tracing, size (Ni, Nj, Ns) -- !
        real(rp), dimension(:,:,:), allocatable :: Pavg
            !! (Ngi, Ngj, Ns) [nPa] Average pressure along flux tube
        real(rp), dimension(:,:,:), allocatable :: Davg
            !! (Ngi, Ngj, Ns) [#/cc] Average density along flux tube
        real(rp), dimension(:,:,:), allocatable :: Bmin
            !! (Ngi, Ngj, NDIM) [nT] Bmin vector
        real(rp), dimension(:,:,:), allocatable :: xyzMin
            !! (Ngi+1, Ngj+1, 3) [Rp] bMin xyz coordinates
        real(rp), dimension(:,:,:), allocatable :: xyzMincc
            !! (Ngi, Ngj, 3) [Rp] cell-centered bMin xyz coordinates
        
        ! (Ngi+1, Ngj+1) corner values
        integer , dimension(:,:), allocatable :: topo   
            !! (Ngi+1, Ngj+1) Topology (0=open, 1=closed)
        real(rp), dimension(:,:), allocatable :: thcon
            !! (Ngi+1, Ngj+1) Co-latitude  of conjugate points
        real(rp), dimension(:,:), allocatable :: phcon
            !! (Ngi+1, Ngj+1) Longitude of conjugate points
        real(rp), dimension(:,:), allocatable :: espot
            !! (Ngi+1, Ngj+1) [kV] electro-static potential
        real(rp), dimension(:,:), allocatable :: bvol
            !! (Ngi+1, Ngj+1) [Rp/nT] Flux-tube volume
        
        ! (Ngi, Ngj) cell-centered values
        integer , dimension(:,:), allocatable :: active
            !! (Ngi, Ngj) (-1=inactive, 0=buffer, 1=active)
        integer , dimension(:,:), allocatable :: active_last
            !! (Ngi, Ngj) Active domain from the previous coupling time step, used for half-step eta calculation
        integer , dimension(:,:), allocatable :: OCBDist
            !! (Ngi, Ngj) Cell distance from open-closed boundary

        ! (Ngi, Ngj, Nk) Varibles coming from RAIJU
        real(rp), dimension(:,:,:), allocatable :: lossRates
            !! (Ngi, Ngj, Nk) [1/s] Loss rates for each grid and lambda point. Generally stays the same over coupling time so we store them all here
        real(rp), dimension(:,:,:), allocatable :: lossRatesPrecip
            !! (Ngi, Ngj, Nk) [1/s] Loss rates that result in precipitation. Should be <= lossRates
        real(rp), dimension(:,:,:), allocatable :: precipType_ele
            !! (Ngi, Ngj, Nk) Prepication type used for electrons
        real(rp), dimension(:,:,:), allocatable :: precipNFlux
            !! (Ngi, Ngj, Nk) [#/cm^2/s] Precipitation number fluxes
        real(rp), dimension(:,:,:), allocatable :: precipEFlux
            !! (Ngi, Ngj, Nk) [erg/cm^2/s] Precipitation energy fluxes
        
        ! (Ngi, Ngj, Nspc+1) (First Nspc index is bulk) Moments
        ! Last dimension will be D/P of different populations (not necessarily same as species)
        ! Example: Total, hot protons, cold protons, electrons, other species
        real(rp), dimension(:,:,:), allocatable :: Den    
            !! (Ngi, Ngj, Nspc+1) Density  [amu/cc]
        real(rp), dimension(:,:,:), allocatable :: Press
            !! (Ngi, Ngj, Nspc+1) Pressure [nPa]
        real(rp), dimension(:,:,:), allocatable :: vAvg
            !! (Ngi, Ngj, Nspc+1) Average cell velocity [km/s]


        !> Only used when debugging
        real(rp), dimension(:,:,:,:), allocatable :: etaFaceReconL
        real(rp), dimension(:,:,:,:), allocatable :: etaFaceReconR   
        real(rp), dimension(:,:,:,:), allocatable :: etaFacePDML
        real(rp), dimension(:,:,:,:), allocatable :: etaFacePDMR
        real(rp), dimension(:,:,:,:), allocatable :: etaFlux   

    end type raijuState_T


!------
! Higher-level types, using above types
!------

    type raijuApp_T
        type(raijuModel_T) :: Model
        type(raijuGrid_T ) :: Grid
        type(raijuState_T) :: State
    end type raijuApp_T

!------
! Interfaces
!------

    abstract interface
        subroutine raijuStateIC_T(Model,Grid,State,inpXML)
            Import :: raijuModel_T, raijuGrid_T, raijuState_T, strLen, XML_Input_T
            type(raijuModel_T) , intent(in) :: Model
            type(raijuGrid_T)  , intent(in) :: Grid
            type(raijuState_T) , intent(inout) :: State
            type(XML_Input_T), intent(in) :: inpXML
        end subroutine raijuStateIC_T

        subroutine raijuUpdateV_T(Model,Grid,State)
            Import :: raijuModel_T, raijuGrid_T, raijuState_T
            type(raijuModel_T) , intent(in) :: Model
            type(raijuGrid_T)  , intent(in) :: Grid
            type(raijuState_T) , intent(inout) :: State
        end subroutine raijuUpdateV_T

        function raijuDP2EtaMap_T(Model,D,kT,vm,amin,amax) result (etaK)
            Import :: rp, raijuModel_T
            type(raijuModel_T), intent(in) :: Model
            real(rp), intent(in) :: D,kT,vm,amin,amax
            real(rp) :: etaK
        end function raijuDP2EtaMap_T

        function raijuELossRate_T(Model,Grid,State,i,j,k) result (lossRate2)
            Import :: rp, raijuModel_T, raijuGrid_T, raijuState_T
            type(raijuModel_T) , intent(inout) :: Model
            type(raijuGrid_T)  , intent(in) :: Grid
            type(raijuState_T) , intent(in) :: State
            integer, intent(in) :: i,j,k
            real(rp), dimension(2) :: lossRate2
        end function raijuELossRate_T
    end interface






end module raijutypes








