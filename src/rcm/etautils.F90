!Utilities for D/P <=> eta mapping

MODULE etautils
  USE kdefs, ONLY : TINY,Me_cgs,Mp_cgs
  USE rcmdefs
  USE rcm_precision
  USE rice_housekeeping_module
  USE constants, ONLY : mass_proton,mass_electron,nt,ev,tiote,boltz
  USE Rcm_mod_subs, ONLY : kcsize,alamc,ikflavc
  USE rcm_mhd_interfaces, ONLY : rcmPScl

  implicit none

  real(rp), private :: density_factor = 0.0 !module private density_factor using planet radius
  real(rp), private :: pressure_factor = 0.0
  logical , private :: doRescaleDef = .true. !Whether to rescale D,P => eta
  real(rprec), private :: sclmass(RCMNUMFLAV) !xmass prescaled to proton
  !Kind of hacky limits to Ti/Te ratio
  real(rp), private, parameter :: TioTeMax = 20.0
  real(rp), private, parameter :: TioTeMin = 0.25
  
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
  subroutine eta2DP(eta,vm,Drc,Dpp,Prc)
    REAL(rprec), intent(in)  :: eta(kcsize)
    REAL(rprec), intent(in)  :: vm
    REAL(rprec), intent(out) :: Drc,Dpp,Prc

    integer :: klow
    
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
    Drc = IntegrateDensity (eta,vm,klow,kcsize)

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

  !Integrate density from eta between channels k1,k2
  function IntegrateDensity(eta,vm,k1,k2) result(D)
    REAL(rprec), intent(in)  :: eta(kcsize)
    REAL(rprec), intent(in)  :: vm
    integer    , intent(in)  :: k1,k2
    REAL(rprec) :: D
    integer :: k

    D = 0.0
    if (vm <= 0) return
    do k=k1,k2
      !Density calc 
      if (alamc(k) > 0.0) then ! only add the ion contribution
        D = D + density_factor*sclmass(ikflavc(k))*eta(k)*vm**1.5
      endif
    enddo !k loop

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

  !Convert given single density/pressure to eeta
  !Optional flag to rescale moments or provide different Ti/Te
  SUBROUTINE DP2eta(Drc,Prc,vm,eta,doRescaleO,tioteO)
    USE conversion_module, ONLY : almmax,almmin,erfexpdiff
    REAL(rprec), intent(in)  :: Drc,Prc,vm
    REAL(rprec), intent(out) :: eta(kcsize)
    logical    , intent(in), optional :: doRescaleO
    REAL(rprec), intent(in), optional :: tioteO 

    REAL(rprec) :: fac,TiovTe,Pion,Pele
    logical :: doRescale

    eta = 0.0
    if ( (vm<0) .or. (Drc<TINY) ) return
    if (present(doRescaleO)) then
      doRescale = doRescaleO
    else
      doRescale = .true.
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

    call DPP2eta(Drc,Pion,Pele,vm,eta,doRescale)

  END SUBROUTINE DP2eta

  !Like DP2eta but take both desired ion and electron pressure
  !Optional flag to rescale pressure
  SUBROUTINE DPP2eta(Drc,Pion,Pele,vm,eta,doRescaleO)
    USE conversion_module, ONLY : almmax,almmin,erfexpdiff
    REAL(rprec), intent(in)  :: Drc,Pion,Pele,vm
    REAL(rprec), intent(out) :: eta(kcsize)
    logical    , intent(in), optional :: doRescaleO

    REAL(rprec) :: Tk,ti,te,A0,prcmI,prcmE
    REAL(rprec) :: xp,xm,pcon,psclI,psclE
    INTEGER(iprec) :: k,klow
    logical :: isIon,doRescale

    eta = 0.0
    if ( (vm<0) .or. (Drc<TINY) ) return
    if (present(doRescaleO)) then
      doRescale = doRescaleO
    else
      doRescale = .true.
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

    A0 = (Drc/density_factor)/(vm**1.5)
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

      xp = SQRT(ev*ABS(almmax(k))*vm/boltz/Tk)
      xm = SQRT(ev*ABS(almmin(k))*vm/boltz/Tk)
      !Use quad prec calc of erf/exp differences
      eta(k) = erfexpdiff(A0,xp,xm)
      
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

END MODULE etautils
