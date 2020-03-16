!Main kaiju definitions/globals
module kdefs
    use, intrinsic :: iso_fortran_env
    implicit none

    !Define variable precisions
    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    !integer, parameter :: ip = INT64
    integer, parameter :: ip = INT32
    
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
    real(rp), parameter :: Kbltz = 1.38e-16      ![cm^2 g /s^2/K=erg/K] Boltzmann constant
    real(rp), parameter :: mec2 = 0.511 ! [MeV] electron rest mass

    real(rp), parameter :: heFrac= 1.16D0       ! Accounts for 4% helium
    real(rp), parameter :: kev2erg= 1.602D-9   ! 1 keV = 1.602e-9 ergs
    real(rp), parameter :: erg2kev= 1.0D0/kev2erg 
    real(rp), parameter :: eCharge= 1.602D-19  ! Charge of electron
    
    !CGS Constants
    real(rp), parameter :: vc_cgs = 2.9979E+10 ![cm/s], speed of light
    real(rp), parameter :: qe_cgs = 4.8032E-10 ![CGS], electron charge
    real(rp), parameter :: Re_cgs = 6.38E8     ![cm]    Earth's radius
    real(rp), parameter :: Me_cgs = 9.1094E-28 ![g] Electron mass
    real(rp), parameter :: Mp_cgs = 1.6726D-24 ![g] Proton mass

    !MKS Constants
    real(rp), parameter :: vc_mks = vc_cgs*(1.0e-2) ![m/s], Speed of light
    !Helper conversions
    real(rp), parameter :: G2nT = 1.0E+5 !Gauss->nT
    real(rp), parameter :: G2T = 1.0E-4 !Gauss->T
    real(rp), parameter :: kev2J = 1.60218E-16 !keV->J
    real(rp), parameter :: Re_km = Re_cgs*(1.0e-2)*(1.0e-3) !km

!Planetary constants
    !Earth
    real(rp), parameter :: EarthM0g = 0.31 !Gauss
    real(rp), parameter :: REarth = Re_cgs*1.0e-2 !m

    real(rp), parameter :: RionE  = 6.5    ! Earth Ionosphere radius in 1000 km
    real(rp), parameter :: EarthPsi0 = 92.4 ! Corotation potential [kV]
    
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

    !Generate name of output file based on tiling
    function genName(caseName,Ri,Rj,Rk,i,j,k,useOldStyle) result(fName)
        character(len=*), intent(in) :: caseName
        integer, intent(in) :: Ri,Rj,Rk,i,j,k
        logical, optional, intent(in) :: useOldStyle
        character(len=strLen) :: fName

        character(len=strLen) :: fHd,fRn,fijk

        if(present(useOldStyle) .and. useOldStyle) then
            fName = genName_old(caseName,Ri,Rj,Rk,i,j,k)
            return
        endif

        if (Ri > 1 .or. Rj > 1 .or. Rk > 1) then
            write(fHd ,'(a,a)') trim(caseName), '_'
            write(fRn ,'(I0.4,a,I0.4,a,I0.4,a)') Ri,'_',Rj,'_',Rk,'_'
            write(fijk,'(I0.4,a,I0.4,a,I0.4,a)') i-1,'_',j-1,'_',k-1,'.gam.h5'

            fName = trim(fHd) // trim(fRn) // trim(fijk)
        else
            if(index(caseName,'.h5') == 0) then
                ! assume does not have extension
                write(fHd ,'(a,a)') trim(caseName), '.gam.h5'
                fName = trim(fHd)
            else
                ! already has extension
                fName = trim(caseName)
            endif
        endif
        !write(*,*) 'ijk / file = ',i,j,k,trim(fName)
    end function genName

    ! Generate old omega format name based on tiling
    function genName_old(caseName,Ri,Rj,Rk,i,j,k) result(fName)
        character(len=*), intent(in) :: caseName
        integer, intent(in) :: Ri,Rj,Rk,i,j,k
        character(len=strLen) :: fName

        character(len=strLen) :: fHd,fRn,fijk,fTl
        integer :: n

        if (Ri > 1 .or. Rj > 1 .or. Rk > 1) then
            n = (j-1) + (i-1)*Rj + (k-1)*Ri*Rj
            write(fHd ,'(a,a)') trim(caseName), '_'
            write(fRn ,'(I0.4,a,I0.4,a,I0.4,a)') Ri,'_',Rj,'_',Rk,'_'
            write(fijk,'(I0.4,a,I0.4,a,I0.4,a)') i-1,'_',j-1,'_',k-1,'_'
            write(fTl ,'(I0.12,a)') n,'.h5'

            fName = trim(fHd) // trim(fRn) // trim(fijk) // trim(fTl)
        else
            fName = caseName
        endif
        !write(*,*) 'ijk / file = ',i,j,k,trim(fName)
    end function genName_old

end module kdefs
