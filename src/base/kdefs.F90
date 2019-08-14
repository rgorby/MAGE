!Main kaiju definitions/globals
module kdefs
    use, intrinsic :: iso_fortran_env
    implicit none

    !Define variable precisions
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: ip = INT64
    
    !Compute precision, IO precision
    integer, parameter :: rp = dp
    integer, parameter :: iop = sp !Precision for IO

!Globals
    real(rp), parameter :: pi = 4.0D0*atan(1.0D0)
    real(rp), parameter :: TINY=1.0e-12, HUGE=1.0e+12
    real(rp), parameter :: rad2deg = 180.0/PI
    real(rp), parameter :: deg2rad = PI/180.0
    
!Defaults
    integer, parameter :: strLen = 1000 !Default size for strings
    integer, parameter :: ALIGN = 64 !Compiler-directive alignment
    integer, parameter :: NDIM=3 !Number of dimensions, size of vector
    integer, parameter :: rseed0=31337 !Default random seed
    
!Physical constants
    real(rp), parameter :: Mu0 = 4*PI*1.0e-7 ! Tm/A
    real(rp), parameter :: mp   = 1.6726D-24     ![g]    proton mass
    real(rp), parameter :: me   = 9.1094D-28     ![g]    electron mass
    real(rp), parameter :: Kbltz = 1.38e-16      ![cm^2 g /s^2/K=erg/K] Boltzmann constant
    
    real(rp), parameter :: heFrac= 1.16D0       ! Accounts for 4% helium
    real(rp), parameter :: kev2erg= 1.602D-9   ! 1 keV = 1.602e-9 ergs
    real(rp), parameter :: erg2kev= 1.0D0/kev2erg 
    real(rp), parameter :: eCharge= 1.602D-19  ! Charge of electron
    !CGS Constants
    real(rp), parameter :: vc_cgs = 2.9979E10  ![cm/s], speed of light
    real(rp), parameter :: qe_cgs = 4.8032E-10 ![CGS], electron charge
    real(rp), parameter :: Re_cgs = 6.38E8     ![cm]    Earth's radius
    real(rp), parameter :: Me_cgs = 9.1094E-28 ![g] Electron mass

!Planetary constants
    !Earth
    real(rp), parameter :: EarthM0g = 0.31 !Gauss
    real(rp), parameter :: REarth = 6.38e6 ! m
    real(rp), parameter :: RionE  = 6.5    ! Earth Ionosphere radius in 1000 km
    !Saturn
    real(rp), parameter :: SaturnM0g = 0.21 !Gauss
    real(rp), parameter :: RSaturnXE = 9.5  !Rx = X*Re
    !Jupiter
    real(rp), parameter :: JupiterM0g = -4.8 !Gauss
    real(rp), parameter :: RJupiterXE = 11.0 !Rx = X*Re
    !Mercury
    real(rp), parameter :: MercuryM0g = 0.00345  !Gauss
    real(rp), parameter :: RMercuryXE = 0.31397  !Rx = X*Re
    !Neptune
    real(rp), parameter :: NeptuneM0g = 0.142  !Gauss
    real(rp), parameter :: RNeptuneXE = 3.860  !Rx = X*Re

!Numbered accessors
    !Directions
    enum, bind(C) 
        enumerator :: IDIR=1, JDIR, KDIR
    endenum 
    enum, bind(C) 
        enumerator :: XDIR=1, YDIR, ZDIR
    endenum 

    !Conserved variables
    enum, bind(C)
        enumerator :: DEN=1,MOMX,MOMY,MOMZ,ENERGY
    endenum
    !Primitive variables
    enum, bind(C)
        enumerator :: VELX=MOMX,VELY,VELZ,PRESSURE
    endenum
    
!Color options for funsies
character, parameter :: ANSIESCAPE = char(27) !Escape character
integer, parameter :: ANSILEN = 5
character(ANSILEN), parameter :: &
    ANSIRED    = ANSIESCAPE // '[31m', &
    ANSIGREEN  = ANSIESCAPE // '[32m', &
    ANSIYELLOW = ANSIESCAPE // '[33m', &
    ANSIBLUE   = ANSIESCAPE // '[34m', &
    ANSIPURPLE = ANSIESCAPE // '[35m', &
    ANSICYAN   = ANSIESCAPE // '[36m', &
    ANSIWHITE  = ANSIESCAPE // '[37m', &
    ANSIRESET  = ANSIESCAPE // '[0m'

    contains

    ! !Print out basic configuration info
    ! subroutine printConfigStamp()
    !     write(*,*) 'Kaiju configuration'
    !     write(*,'(2a)') 'Compiler = ', compiler_version()
    !     write(*,'(2a)') 'Compiler flags = ', compiler_options()
    !     !write(*,*) 'Git hash = ', GITCOMMITHASH
    !     !write(*,*) 'Git branch = ', GITBRANCH

    ! end subroutine printConfigStamp

end module kdefs
