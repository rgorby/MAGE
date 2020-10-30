!Utilities for D/P <=> eta mapping

MODULE etautils
  USE kdefs, ONLY : TINY
  USE rcm_precision
  USE rice_housekeeping_module
  USE constants, ONLY : mass_proton,mass_electron,nt,ev,tiote,boltz
  USE Rcm_mod_subs, ONLY : kcsize,alamc,ikflavc

  implicit none

  real(rp), private :: density_factor = 0.0 !module private density_factor using planet radius
  real(rp), private :: pressure_factor = 0.0


  contains

  !Set density/pressure factors using planet radius
  subroutine SetFactors(Rx)
  	real(rp), intent(in) :: Rx

  	pressure_factor = 2./3.*ev/Rx*nt
  	density_factor = nt/Rx
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

    integer :: k,klow
    real(rprec) :: sclmass(RCMNUMFLAV) !xmass prescaled to proton

    !Set lowest RC channel
    if (use_plasmasphere) then
      klow = 2
    else
      klow = 1
    endif

    !Set scaled mass by hand here to avoid precision issues
    sclmass(RCMELECTRON) = mass_electron/mass_proton
    sclmass(RCMPROTON) = 1.0

    Drc = 0.0
    Dpp = 0.0
    Prc = 0.0

    if (vm <= 0) return

    !Do RC channels
    do k=klow,kcsize
    	!Pressure calc in pascals
    	Prc = Prc + pressure_factor*ABS(alamc(k))*eta(k)*vm**2.5

    	!Density calc (ring current)
      if (alamc(k) > 0.0) then ! only add the ion contribution
      	Drc = Drc + density_factor*sclmass(ikflavc(k))*eta(k)*vm**1.5
      endif
    enddo !k loop

    !Handle plasmasphere
    if (use_plasmasphere) then
    	Dpp = density_factor*sclmass(RCMPROTON)*eta(1)*vm**1.5
    else
    	Dpp = 0.0
    endif

  end subroutine eta2DP

  !Convert given single density/pressure to eeta
  SUBROUTINE DP2eta(Drc,Prc,vm,eta)
    USE conversion_module, ONLY : almmax,almmin,erfexpdiff
    REAL(rprec), intent(in)  :: Drc,Prc,vm
    REAL(rprec), intent(out) :: eta(kcsize)

    REAL(rprec), PARAMETER :: fac = tiote/(1.+tiote)

    REAL(rprec) :: Tk,ti,te,A0,prcmI,prcmE,pmhdI,pmhdE
    REAL(rprec) :: xp,xm,pcon,psclI,psclE
    INTEGER(iprec) :: k,klow
    logical :: isIon

    eta = 0.0
    if ( (vm<0) .or. (Drc<TINY) ) return

    !Set lowest RC channel
    if (use_plasmasphere) then
      klow = 2
    else
      klow = 1
    endif

    !Get ion/electron temperature
    ti = fac*Prc/Drc/boltz
    te = ti/tiote

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

    !Now rescale eeta channels to conserve pressure integral between MHD/RCM
    !In particular, we separately conserve ion/electron contribution to total pressure
    pmhdI = Prc*tiote/(1.0+tiote) !Desired ion pressure
    pmhdE = Prc*  1.0/(1.0+tiote) !Desired elec pressure

    psclI = pmhdI/prcmI
    psclE = pmhdE/prcmE

    !Loop over channels and rescale      
    do k=klow,kcsize
      IF (ikflavc(k) == RCMELECTRON) THEN  ! electrons
        eta(k) = psclE*eta(k)
      ELSE IF (ikflavc(k) == RCMPROTON) THEN ! ions (protons)
        eta(k) = psclI*eta(k)
      ENDIF
    enddo

  END SUBROUTINE DP2eta

END MODULE etautils
