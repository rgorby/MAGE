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

    real(rprec), private :: sclmass(RCMNUMFLAV) !xmass prescaled to proton
    integer    , private, dimension(RCMNUMFLAV,2) :: flavorBnds

    real(rp), private, parameter :: kapDefault = 6.0

    real(rp), private, dimension(RCMNUMFLAV) :: kapDefs

    contains

    !Set density/pressure factors using planet radius
    subroutine SetFactors(Rx)
        real(rp), intent(in) :: Rx
        integer :: n,k

        pressure_factor = 2./3.*ev/Rx*nt
        density_factor = nt/Rx

        !Set scaled mass by hand here to avoid precision issues
        sclmass(RCMELECTRON) = Me_cgs/Mp_cgs
        sclmass(RCMPROTON) = 1.0

        !Get flavor bounds
        do n=1,RCMNUMFLAV
        !NOTE: Doing stupid code to avoid findloc

            !Find first value
            do k=1,kcsize
                if (ikflavc(k) == n) exit
            enddo
            flavorBnds(n,1) = k
            !Find last value
            do k=flavorBnds(n,1),kcsize
                if (ikflavc(k) /= n) exit
            enddo
            flavorBnds(n,2) = k-1
        enddo

        !Do fixes
        if (use_plasmasphere) then
            !Bump up electron min to avoid plasmasphere channel
            flavorBnds(RCMELECTRON,1) = flavorBnds(RCMELECTRON,1) + 1
        endif

        !Set default kappa values if using kappa
        kapDefs(RCMELECTRON) = 4.0
        kapDefs(RCMPROTON  ) = kapDefault

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

    !Get Pk - Pressure contribution from each channel
    subroutine eta2Pk(eta,vm,Pk)
        REAL(rprec), intent(in)  :: eta(kcsize)
        REAL(rprec), intent(in)  :: vm
        REAL(rprec), intent(out) :: Pk(kcsize)

        integer :: k
        REAL(rprec) :: dP

        Pk = 0.0
        if (vm <= 0) return

        do k=1,kcsize
            dP = pressure_factor*ABS(alamc(k))*eta(k)*vm**2.5
            Pk(k) = dP
        enddo
    end subroutine eta2Pk

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
            doRescale = doRescaleDef
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

        REAL(rprec) :: Tk,ti,te,kap
        INTEGER(iprec) :: k,klow,n,k1,k2
        logical :: isIon,doRescale,doKap
        REAL(rprec), dimension(RCMNUMFLAV) :: Ds,Ps

        !DEBUG:
        !REAL(rprec) :: Dt,Pit,Pet,Db,Pib,Peb,Da,Pia,Pea

        eta = 0.0
        if ( (vm<0) .or. (Drc<TINY) ) return
        if (present(doRescaleO)) then
            doRescale = doRescaleO
        else
            doRescale = doRescaleDef
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
                !Replace kap-default w/ species specific value
                kap = kapDefs(ikflavc(k))
                eta(k) = Kappa2Eta  (Drc,vm,Tk,almmin(k),almmax(k),alamc(k),kap)
            else
                eta(k) = Maxwell2Eta(Drc,vm,Tk,almmin(k),almmax(k),alamc(k))
            endif

        enddo !k loop

        if (.not. doRescale) return !We're done here

    !Now rescale to get desired Pi and Pe, it's a whole thing
        !Set target density/pressure
        Ds(RCMELECTRON:RCMPROTON) = Drc !Both use Drc
        Ps(RCMELECTRON) = Pele
        Ps(RCMPROTON  ) = Pion

        ! !DEBUG:
        ! Db =  IntegrateDensity(eta,vm,klow,kcsize)
        ! call IntegratePressureIE(eta,vm,Pib,Peb)
        ! Dt = Drc
        ! Pet = Pele
        ! Pit = Pion

        do n=1,RCMNUMFLAV
            !Rescale each flavor
            k1 = flavorBnds(n,1)
            k2 = flavorBnds(n,2)
            call RescaleEta(eta,k1,k2,Ds(n),Ps(n),vm)
        enddo

        ! Da =  IntegrateDensity(eta,vm,klow,kcsize)
        ! call IntegratePressureIE(eta,vm,Pia,Pea)

        ! !$OMP CRITICAL
        ! write(*,*) '---'
        ! write(*,*) 'Target: ', Dt,Pet,Pit
        ! write(*,*) "Before: ", Db/Dt,Peb/Pet,Pib/Pit
        ! write(*,*) "After : ", Da/Dt,Pea/Pet,Pia/Pit
        ! write(*,*) '---'
        ! !$OMP END CRITICAL

    END SUBROUTINE DPP2eta

    !Rescale eta(k1:k2) to have moments D,P
    !Adapted from code in different RCM
    SUBROUTINE RescaleEta(eta,k1,k2,D,P,vm)
        REAL(rprec), intent(inout) :: eta(kcsize)
        integer, intent(in) :: k1,k2
        real(rprec), intent(in) :: D,P,vm

        real(rprec) :: A,B
        integer :: k,kmax

        !Calculate coefficients to rescale
        !eta' = eta x (A + |lam| B)
        !between k1,kmax and 0 between kmax,k2

        call GetRescaleAB(eta,k1,k2,D,P,vm,A,B,kmax)

        do k=k1,kmax
            eta(k) = eta(k) * (A + B*abs(alamc(k)))
        enddo
        do k=kmax+1,k2
            eta(k) = 0.0
        enddo

    END SUBROUTINE RescaleEta

    !Get rescaling coefficients
    SUBROUTINE GetRescaleAB(eta,k1,k2,D,P,vm,A,B,kmax)
        REAL(rprec), intent(in) :: eta(kcsize)
        real(rprec), intent(in) :: D,P,vm
        integer, intent(in) :: k1,k2
        real(rprec), intent(out) :: A,B
        integer, intent(out) :: kmax

        real(rprec) :: etaP,Sig0,Sig1,Sig2,siggy,minScl
        real(rprec) :: AlphaD,AlphaP,DoA,PoA,Mu
        integer :: iflv

        iflv = ikflavc(k1)
        if (iflv == RCMELECTRON) then
            Mu = 1.0
        else
            Mu = sclmass(iflv)
        endif

        kmax = k2
        AlphaD = density_factor*Mu*(vm**1.5)
        AlphaP = pressure_factor  *(vm**2.5)

        DoA = D/AlphaD
        PoA = P/AlphaP

        do kmax=k2,k1+1,-1
            !Get trial A,B coefficients
            Sig0 = sum(eta(k1:kmax))
            Sig1 = sum( eta(k1:kmax)*abs(alamc(k1:kmax)) )
            Sig2 = sum( eta(k1:kmax)*abs(alamc(k1:kmax))*abs(alamc(k1:kmax)) )

            siggy = Sig2*Sig0 - Sig1**2.0
            if (siggy < TINY) cycle

            A = (1/siggy)*( Sig2*DoA - Sig1*PoA )
            B = PoA/Sig2 - A*Sig1/Sig2

            !Check A,B: Require A + |lam|B >= 0 forall lam
            minScl = minval( A + B*abs(alamc(k1:kmax)) )
            if (minScl <= 0) then
                !Bad A,B keep trying
                cycle
            else
                !This is good, let's get out of here
                return
            endif
        enddo

    !If still here, nothing worked. Use single moment rescaling
        B = 0.0
        kmax = k2
        etaP = IntegratePressure(eta,vm,k1,k2)

        if (etaP*rcmPScl > TINY) then
            A = P/etaP
        else
            A = 0.0
        endif

    END SUBROUTINE GetRescaleAB

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

    !NOTE: This is just adaped from other RCM code, seems to be eqn 3.12 from 10.1007/s11214-013-9982-9
    function Kappa2Eta(Drc,vm,Tk,almin,almax,almc,kapO) result(etak)
        real(rp), intent(in) :: Drc,vm,Tk,almin,almax,almc
        real(rp), intent(in), optional :: kapO
        real(rp) :: etak

        real(rp) :: kap,kap15,Tev,E0_ev,E_ev
        real(rp) :: A0,kapgam,kapbar,kArg,delscl

        if (present(kapO)) then
            kap = kapO
        else
            kap = kapDefault
        endif
        
        kap15 = kap-1.5
        Tev = Tk*Kbltz*erg2kev*(1.0e+3) !temperature in eV
        E0_ev = Tev*kap15/kap

        E_ev = abs(almc)*vm !eV

        A0 = (2.0/sqrt(PI))*(Drc/density_factor)/(vm**1.5)
        !TODO: Double check this extra factor of 2
        !A0 = 2.0*A0

        kapgam = gamma(kap+1.0)/gamma(kap-0.5)
        !TODO: Figure out k_0 vs. kappa where kappa - 1.5 = kappa_0
        kapbar = kap15 !Should be kap-3/2 or kap
        !kapbar = kap !Using kap consistent w/ 10.1002/2015JA021166

        kArg = 1.0 + (E_ev/E0_ev)/kapbar
        delscl = vm*(almax-almin)/E0_ev
        etak = A0*kapgam/(kapbar**1.5) * sqrt(E_ev/E0_ev)*delscl*((kArg)**(-kap-1.0))

    end function Kappa2Eta

END MODULE etautils
