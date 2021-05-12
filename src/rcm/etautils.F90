!Utilities for D/P <=> eta mapping

MODULE etautils
  USE kdefs, ONLY : TINY,Me_cgs,Mp_cgs
  USE rcmdefs
  USE rcm_precision
  USE rice_housekeeping_module
  USE constants, ONLY : mass_proton,mass_electron,nt,ev,tiote,boltz
  USE Rcm_mod_subs, ONLY : kcsize,alamc,ikflavc
  USE rcm_mhd_interfaces, ONLY : rcmPScl
  USE conversion_module, ONLY : erfexpdiff
  implicit none

  real(rp), private :: density_factor = 0.0 !module private density_factor using planet radius
  real(rp), private :: pressure_factor = 0.0
  logical , private :: doRescaleDef = .true. !Whether to rescale D,P => eta
  logical , private :: doKapDef     = .false. !Whether to do kappa by default

  real(rprec), private :: sclmass(RCMNUMFLAV) !xmass prescaled to proton
  !Kind of hacky limits to Ti/Te ratio
  real(rp), private, parameter :: TioTeMax = 20.0
  real(rp), private, parameter :: TioTeMin = 0.25
  real(rp), private, parameter :: kapDefault = 6.0

  contains

  !Set density/pressure factors using planet radius
  subroutine SetFactors(Rx)
  	real(rp), intent(in) :: Rx

  	pressure_factor = 2./3.*ev/Rx*nt
  	density_factor = nt/Rx

    !Set scaled mass by hand here to avoid precision issues
    sclmass(RCMELECTRON) = Me_cgs/Mp_cgs
    sclmass(RCMPROTON) = 1.0

  end subroutine SetFactors

  !Simple functions to access XXX factors
  function GetDensityFactor() result(df)
  	real(rp) :: df
  	df = density_factor
  end function GetDensityFactor

  function GetPressureFactor() result(pf)
  	real(rp) :: pf
  	pf = pressure_factor
  end function GetPressureFactor

  !Convert single eta to density (RC/plasmasphere) and pressure
  subroutine eta2DP(eta,vm,Drc,Dpp,Prc,doCharge0)
    REAL(rprec), intent(in)  :: eta(kcsize)
    REAL(rprec), intent(in)  :: vm
    REAL(rprec), intent(out) :: Drc,Dpp,Prc
    logical    , intent(in), optional :: doCharge0
    
    integer :: klow
    logical :: doC0

    if (present(doCharge0)) then
      doC0 = doCharge0
    else
      doC0 = .false.
    endif

    
    
    !Set lowest RC channel
    if (use_plasmasphere) then
      klow = 2
    else
      klow = 1
    endif

    Drc = 0.0
    Dpp = 0.0
    Prc = 0.0

    if (vm <= 0) return

    !Do RC channels
    Prc = IntegratePressure(eta,vm,klow,kcsize)
    Drc = IntegrateDensity (eta,vm,klow,kcsize,doC0)

    !Handle plasmasphere
    if (use_plasmasphere) then
    	Dpp = density_factor*sclmass(RCMPROTON)*eta(1)*vm**1.5
    else
    	Dpp = 0.0
    endif

  end subroutine eta2DP

  !Integrate pressure from eta between channels k1,k2
  function IntegratePressure(eta,vm,k1,k2) result(P)
    REAL(rprec), intent(in)  :: eta(kcsize)
    REAL(rprec), intent(in)  :: vm
    integer    , intent(in)  :: k1,k2
    REAL(rprec) :: P
    integer :: k

    P = 0.0
    if (vm <= 0) return

    do k=k1,k2
      !Pressure calc in pascals
      P = P + pressure_factor*ABS(alamc(k))*eta(k)*vm**2.5
    enddo
  end function IntegratePressure

  !Integrate density from eta between channels k1,k2 (neglect cold species)
  !doCharge0: optional argument, whether to attempt to mock up charge neutrality mass in electron regions
  function IntegrateDensity(eta,vm,k1,k2,doCharge0) result(D)
    REAL(rprec), intent(in)  :: eta(kcsize)
    REAL(rprec), intent(in)  :: vm
    integer    , intent(in)  :: k1,k2
    logical    , intent(in), optional :: doCharge0

    real(rp) :: D
    logical :: doC0
    REAL(rprec) :: Di,De
    integer :: k

    if (present(doCharge0)) then
      doC0 = doCharge0
    else
      doC0 = .false.
    endif

    Di = 0.0
    De = 0.0
    if (vm <= 0) return
    do k=k1,k2
      !Density calc 
      if (alamc(k) > TINY) then 
        !Hot ion contribution
        Di = Di + density_factor*sclmass(ikflavc(k))*eta(k)*vm**1.5
      else if (alamc(k) < TINY) then
        !Cold ion counterparts to hot electrons
        De = De +  density_factor*sclmass(RCMPROTON)*eta(k)*vm**1.5
      endif

    enddo !k loop

    if (doC0) then
      D = max(Di,De) !Include mass from cold counterparts to hot electrons if needed
    else
      D = Di
    endif

  end function IntegrateDensity

  !Get both ion and electron pressures
  subroutine IntegratePressureIE(eta,vm,iP,eP)
    REAL(rprec), intent(in)  :: eta(kcsize)
    REAL(rprec), intent(in)  :: vm
    REAL(rprec), intent(out) :: iP,eP

    INTEGER(iprec) :: k
    iP = 0.0
    eP = 0.0
    if (vm <= 0) return
    do k=1,kcsize !Include psphere b/c it won't contribute
      if (abs(alamc(k))<TINY) cycle

      if (alamc(k)>TINY) then
        !Ion pressure
        iP = iP + pressure_factor*ABS(alamc(k))*eta(k)*vm**2.5
      else
        !Elec pressure
        eP = eP + pressure_factor*ABS(alamc(k))*eta(k)*vm**2.5
      endif
    enddo

  end subroutine IntegratePressureIE
  
  !Get Ti/Te for a given eta
  ! Ti/Te = Pi/Pe b/c Ni=Ne
  function GetTioTe(eta,vm) result(TiovTe)
    REAL(rprec), intent(in)  :: eta(kcsize)
    REAL(rprec), intent(in)  :: vm
    REAL(rprec) :: TiovTe
    REAL(rprec) :: eP,iP
    INTEGER(iprec) :: k

    TiovTe = 0.0
    if (vm <= 0) return
    iP = 0.0
    eP = 0.0
    call IntegratePressureIE(eta,vm,iP,eP)
    if (eP>TINY) then
      TiovTe = iP/eP
    else
      TiovTe = 0.0
    endif
  end function GetTioTe

  subroutine MaxVsKap(Drc,Prc,vm)
    REAL(rprec), intent(in)  :: Drc,Prc,vm

    REAL(rprec), dimension(kcsize) :: etaMax,etaKap
    REAL(rprec) :: Dm,Dk,Pm,Pk,Dpp

    call DP2eta(Drc,Prc,vm,etaMax,doRescaleO=.false.,doKapO=.false.)
    call DP2eta(Drc,Prc,vm,etaKap,doRescaleO=.false.,doKapO=.true. )

    call eta2DP(etaMax,vm,Dm,Dpp,Pm)
    call eta2DP(etaKap,vm,Dk,Dpp,Pk)

    write(*,*) 'Max/Kap: D/P = ', Drc*1.0e-6,Dm*1.0e-6,Dk*1.0e-6,Prc*1.0e+9,Pm*1.0e+9,Pk*1.0e+9

  end subroutine MaxVsKap
  !Convert given single density/pressure to eeta
  !Optional flag to rescale moments or provide different Ti/Te
  SUBROUTINE DP2eta(Drc,Prc,vm,eta,doRescaleO,tioteO,doKapO,kapO)
    USE conversion_module, ONLY : almmax,almmin,erfexpdiff
    REAL(rprec), intent(in)  :: Drc,Prc,vm
    REAL(rprec), intent(out) :: eta(kcsize)
    logical    , intent(in), optional :: doRescaleO,doKapO
    REAL(rprec), intent(in), optional :: tioteO,kapO

    REAL(rprec) :: fac,TiovTe,Pion,Pele,kap
    logical :: doRescale,doKap

    eta = 0.0
    kap = 0.0

    if ( (vm<0) .or. (Drc<TINY) ) return
    if (present(doRescaleO)) then
      doRescale = doRescaleO
    else
      doRescale = .true.
    endif

    if (present(doKapO)) then
      doKap = doKapO
    else
      doKap = doKapDef
    endif

    if (doKap) then
      if (present(kapO)) then
        kap = kapO
      else
        kap = kapDefault
      endif
    endif

    if (present(tioteO)) then !Use specified Ti/Te
      TiovTe = tioteO
      call ClampTioTe(TiovTe) !Ensure reasonable number
    else !Use default from defs
      TiovTe = tiote
    endif
    fac = TiovTe/(1.0+TiovTe)

    Pion = Prc*TiovTe/(1.0+TiovTe) !Desired ion pressure
    Pele = Prc*   1.0/(1.0+TiovTe) !Desired elec pressure

    call DPP2eta(Drc,Pion,Pele,vm,eta,doRescale,doKap,kap)

  END SUBROUTINE DP2eta

  !Like DP2eta but take both desired ion and electron pressure
  !Optional flag to rescale pressure
  SUBROUTINE DPP2eta(Drc,Pion,Pele,vm,eta,doRescaleO,doKapO,kapO)
    USE conversion_module, ONLY : almmax,almmin
    REAL(rprec), intent(in)  :: Drc,Pion,Pele,vm
    REAL(rprec), intent(out) :: eta(kcsize)
    logical    , intent(in), optional :: doRescaleO,doKapO
    real(rp)   , intent(in), optional :: kapO

    REAL(rprec) :: Tk,ti,te,prcmI,prcmE,kap
    REAL(rprec) :: pcon,psclI,psclE
    INTEGER(iprec) :: k,klow
    logical :: isIon,doRescale,doKap

    eta = 0.0
    if ( (vm<0) .or. (Drc<TINY) ) return
    if (present(doRescaleO)) then
      doRescale = doRescaleO
    else
      doRescale = .true.
    endif

    if (present(doKapO)) then
      doKap = doKapO
    else
      doKap = doKapDef
    endif

    if (doKap) then
      if (present(kapO)) then
        kap = kapO
      else
        kap = kapDefault
      endif
    endif

    !Set lowest RC channel
    if (use_plasmasphere) then
      klow = 2
    else
      klow = 1
    endif

    !Get ion/electron temperature
    ti = Pion/Drc/boltz
    te = Pele/Drc/boltz

    prcmI = 0.0 !Cumulative ion pressure
    prcmE = 0.0 !Cumulative electron pressure

    !Loop over non-zero channels
    do k=klow,kcsize
      !Get right temperature
      IF (ikflavc(k) == RCMELECTRON) THEN  ! electrons
        Tk = te
        isIon = .false.
      ELSE IF (ikflavc(k) == RCMPROTON) THEN ! ions (protons)
        Tk = ti
        isIon = .true.
      ENDIF

      if (Tk<TINY) then
        eta(k) = 0.0
        cycle
      endif

      if (doKap) then
        eta(k) = Kappa2Eta  (Drc,vm,Tk,almmin(k),almmax(k),alamc(k),kap)
      else
        eta(k) = Maxwell2Eta(Drc,vm,Tk,almmin(k),almmax(k),alamc(k))
      endif

      !Pressure contribution from this channel
      pcon = pressure_factor*ABS(alamc(k))*eta(k)*vm**2.5

      if (isIon) then
        prcmI = prcmI + pcon
      else
        prcmE = prcmE + pcon
      endif

    enddo !k loop

    if (.not. doRescale) return !We're done here

    !Now rescale to get desired Pi and Pe
    !NOTE: This will affect density
    !Check if pressures are above TINY nPa
    if (prcmI*rcmPScl > TINY) then
      psclI = Pion/prcmI
    else
      psclI = 0.0
    endif

    if (prcmE*rcmPScl > TINY) then
      psclE = Pele/prcmE
    else
      psclE = 0.0
    endif

    !Loop over channels and rescale      
    do k=klow,kcsize
      IF (ikflavc(k) == RCMELECTRON) THEN  ! electrons
        eta(k) = psclE*eta(k)
      ELSE IF (ikflavc(k) == RCMPROTON) THEN ! ions (protons)
        eta(k) = psclI*eta(k)
      ELSE
        write(*,*) 'Unknown species!'
        eta(k) = 0.0
      ENDIF
    enddo

  END SUBROUTINE DPP2eta

  SUBROUTINE ClampTioTe(TiovTe)
    REAL(rprec), intent(inout)  :: TiovTe
    if (TiovTe<TioTeMin) TiovTe = TioTeMin
    if (TiovTe>TioTeMax) TiovTe = TioTeMax
  END SUBROUTINE ClampTioTe

!======
  !Specific PSD types
  !General form: Maxwell2Eta(Drc,vm,Tk,almin,almax,almc)
  !almin/max/c are min/max and center of lambda bin
  !Drc = Density [#/m3]
  !vm  = (nt/re)^0.667
  !Tk  = Temperature [K]

  function Maxwell2Eta(Drc,vm,Tk,almin,almax,almc) result(etak)
    real(rp), intent(in) :: Drc,vm,Tk,almin,almax,almc
    real(rp) :: etak

    real(rp) :: A0,xp,xm

    A0 = (Drc/density_factor)/(vm**1.5)
    xp = SQRT(ev*ABS(almax)*vm/boltz/Tk)
    xm = SQRT(ev*ABS(almin)*vm/boltz/Tk)

    !Use quad prec calc of erf/exp differences, Pembroke+ Eqn B5
    etak = erfexpdiff(A0,xp,xm)

  end function Maxwell2Eta

  !NOTE: This is just copied from other RCM code, seems to be eqn 3.12 from 10.1007/s11214-013-9982-9

  function Kappa2Eta(Drc,vm,Tk,almin,almax,almc,kapO) result(etak)
    real(rp), intent(in) :: Drc,vm,Tk,almin,almax,almc
    real(rp), intent(in), optional :: kapO
    real(rp) :: etak

    real(rp) :: kap,ftv,Tkev,Pnpa,tV,pV,kap15,kapgam
    real(rp) :: trans,tran23,A0,Px,Tx,Kx,dLamx,kArg

    if (present(kapO)) then
      kap = kapO
    else
      kap = kapDefault
    endif
    !Start by converting to same units as RCM code
    !tV = ti x ftv^(2/3) = keV [Re/nT]^2/3
    !pV = pi x ftv^(5/3) = nPa [Re/nT]^5/3

    ftv = vm**(-3.0/2) ! Re/nT
    Tkev = Tk*Kbltz*erg2kev
    Pnpa = (Tk*boltz*Drc)*(1.0e+9) !P => nPa

    tV = (Tkev)*(ftv**(2.0/3))
    pV = (Pnpa)*(ftv**(5.0/3))

    kap15 = kap-1.5
    kapgam = gamma(kap+1.0)/gamma(kap-0.5)

    trans = 1/density_factor
    tran23 = trans**(2.0/3)
    A0 = (2.0/sqrt(PI))*(kap15**-1.5)*kapgam
    Px = (pV*(1.0e-9*trans**(5.0/3)))
    Tx = ((tV*kev2J*tran23)**-2.5)
    Kx = sqrt(abs(almc))*tran23*(kev2J*1.0e-3)
    dLamx = (almax-almin)*tran23*(kev2J*1.0e-3)

    kArg = 1.0 + abs(almc)/(kap15*tV*1.0e+3)

    etak = A0*Px*Tx*Kx*dLamx * ((kArg)**(-kap-1.0))
  end function Kappa2Eta

END MODULE etautils
