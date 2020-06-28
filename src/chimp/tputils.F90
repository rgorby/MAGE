!Simple helper routines for particle calculations

module tputils
    use chmpdefs
    use chmpunits
    use tptypes
    use ebtypes
    use ebinterp
    use math
    use gciter
    use streamline
    
    implicit none

    !Overloader for p->gamma
    interface p2Gam
        module procedure prt2Gam, pv2Gam
    end interface

    !Overloader for adiabaticity
    interface eGC
        module procedure eGC_Mu,eGC_p
    end interface

    logical, private :: doGyroDT = .true. !Use fraction of gyrofrequency for FO timesteps
    logical, private :: doLoudKill = .false. !Output when killing TP
    
    contains


 
!---------------------------------
!Timestep calculation

    function dtFO(Model,pFO,E,B) result(dt)
        type(chmpModel_T), intent(in) :: Model
        real(rp), dimension(NDIM), intent(in) :: pFO,E,B
        real(rp) :: dt
        real(rp) :: dtO,dtE
        real(rp) :: MagF, gamma, Omega, F(NDIM)

        gamma = p2Gam(pFO,Model%m0)

        if (doGyroDT) then
            !dt ~ 2pi/Omega
            Omega = abs(Model%q0)*norm2(B)/Model%m0
            Omega = Omega/gamma !Relativistic correction
            dtO = 2*PI/Omega
            dtE = norm2(pFO)/norm2(Model%q0*E)
            dt = min(dtO,dtE)

        else
            !dt ~ p/F
            F = Model%q0*(E + cross(pFO,B)/(Model%m0*gamma) )
            MagF = max(norm2(F),TINY)
            dt = norm2(pFO)/MagF
        endif

        dt = max(Model%epsht*dt,TINY)
    end function dtFO

!---------------------------------
!Adiabaticity parameter (Model,B,JacB) & either P or Mu

    function eGC_Mu(Model,B,JacB,Mu) result(eps)
        type(chmpModel_T), intent(in) :: Model
        real(rp), dimension(NDIM), intent(in) :: B
        real(rp), intent(in) :: Mu
        real(rp), dimension(NDIM,NDIM) :: JacB
        real(rp) :: eps
        real(rp) :: MagB,rgScl,bScl
        MagB = norm2(B)
        rgScl = Model%m0*sqrt(2*Mu*MagB/Model%m0)/( abs(Model%q0)*MagB )
        bScl = MagB/sqrt(sum(JacB**2.0))
        eps = rgScl/bScl
        
    end function eGC_Mu

    function eGC_p(Model,B,JacB,p) result(eps)
        type(chmpModel_T), intent(in) :: Model
        real(rp), dimension(NDIM), intent(in) :: B,p
        real(rp), dimension(NDIM,NDIM) :: JacB
        real(rp) :: eps
        real(rp) :: MagB,rgScl,bScl

        MagB = norm2(B)
        rgScl = norm2(p)/( abs(Model%q0)*MagB )
        bScl = MagB/sqrt(sum(JacB**2.0))
        eps = rgScl/bScl
    end function eGC_p
!---------------------------------
!Energy calculations, p->K [code] or K [keV]
    !K = (gamma-1)*m*c2 = m*(gamma*v)^2/(gamma+1)
    !Old: K = (Model%m0*mec2*1.0e+3)*(gamma -1)
    !Returns kinetic energy [keV] of a particle
    function prt2kev(Model,prt) result(K)
        type(prt_T), intent(in) :: prt
        type(chmpModel_T), intent(in) :: Model
        
        real(rp) :: K
        real(rp) :: gamma
        if (prt%isGC) then
            gamma = prt%Q(GAMGC)
        else
            gamma = p2Gam(prt%Q(PXFO:PZFO),Model%m0)
        endif
        K = (Model%m0*mec2*1.0e+3)*(gamma-1.0)

    end function prt2kev

    !Convert momentum vector to kinetic energy [kev]
    function p2kev(Model,p) result(K)
        type(chmpModel_T), intent(in) :: Model
        real(rp), intent(in) :: p(NDIM)
        real(rp) :: K

        K = (mec2*1.0e+3)*p2K(p,Model%m0)

    end function
    !Turn momentum into kinetic energy (code units)
    function p2K(p,m0) result(K)
        real(rp), dimension(NDIM), intent(in) :: p
        real(rp), intent(in) :: m0
        real(rp) :: K
        real(rp) :: gamma,u2
        gamma = p2Gam(p,m0)
        u2 = dot_product(p/m0,p/m0)
        K = m0*u2/(gamma+1.0)

    end function p2K


!---------------------------------
!Simple conversions from p to Gamma
    function prt2Gam(prt,m0)
        type(prt_T), intent(in) :: prt
        real(rp), intent(in) :: m0
        real(rp) :: prt2Gam
        real(rp) :: p(NDIM)
        if (prt%isGC) then
            prt2Gam = prt%Q(GAMGC)
        else
            prt2Gam = p2Gam(prt%Q(PXFO:PZFO),m0)
        endif
    end function prt2Gam

    function pv2Gam(p,m0) result(Gam)
        real(rp), dimension(NDIM), intent(in) :: p
        real(rp), intent(in) :: m0
        real(rp) :: Gam

        Gam = sqrt(1 + dot_product(p/m0,p/m0))
    end function pv2Gam

!---------------------------------
!Gets invariant Mu for particle or r/p

    !Returns Mu for particle (calculate for FO, just use for GC)
    function prtMu(prt,t,Model,ebState)
        type(prt_T), intent(in) :: prt
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in) :: t
        real(rp) :: prtMu
        if (.not. prt%isIn) then
            prtMu = -1.0
        endif
        if (prt%isGC) then
            prtMu = prt%Q(MUGC)
        else
            prtMu = CalcMu(prt%Q(XPOS:ZPOS),prt%Q(PXFO:PZFO),t,Model,ebState)
        endif        
    end function
    !Calculate Mu for a particle @ (r,p,t)
    !Return -1 if unconverged
    function CalcMu(r,p,t,Model,ebState,isConvO,RgcO,bO) result(Mu)
        real(rp), dimension(NDIM), intent(in) :: r,p
        real(rp), intent(in) :: t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out), optional :: isConvO
        real(rp), dimension(NDIM), intent(out), optional :: RgcO,bO

        real(rp) :: Mu

        logical :: isConv
        real(rp) :: MagB, p11
        real(rp), dimension(NDIM) :: Rgc, pPerp,E,B,vExB

        Rgc = CalcGC(r,p,t,Model,ebState,isConv)
        if (isConv) then
            !Calculate fields at this position
            call ebFields(r,t,Model,ebState,E,B,vExB=vExB)
            !Split momentum into // and perp
            call SplitP(Model,p,E,B,vExB,p11,pPerp)
            !Calculate mu using perp momentum (in ExB frame)
            MagB = max(norm2(B),TINY)
            Mu = dot_product(pPerp,pPerp)/(2*Model%m0*MagB)
        else
            !Unconverged
            Mu = -1
        endif
        if (present(isConvO)) isConvO = isConv
        if (present(RgcO)) RgcO = Rgc
        if (present(bO)) bO = B
    end function CalcMu
!---------------------------------
!Particle formulation flipping (GC<->FO)
    subroutine gc2fo(Model,prt,t,ebState)
        type(chmpModel_T), intent(in) :: Model
        type(prt_T), intent(inout) :: prt
        real(rp), intent(in) :: t
        type(ebState_T), intent(in)   :: ebState

        real(rp), dimension(NDIM) :: Rgc,E,B,vExB,Ro,pFO,pPerp
        real(rp), dimension(NDIM) :: xhat,yhat,bhat,ebhat,phat
        real(rp) :: bAbs,ebx,eby,psi,absPp,pAbs

        Rgc = prt%Q(XPOS:ZPOS) !Save GC position
        !Get field values @ GC
        call ebFields(Rgc,t,Model,ebState,E,B,vExB=vExB,ijkO=prt%ijk0)

        !Get coordinate system
        call MagTriad(Rgc,B,xhat,yhat,bhat)
        bAbs = norm2(B)
        ebhat = normVec(vExB)

        !Calculate phase angle such that phat is orthogonal to ebhat
        !phat = cos(psi)*xhat + sin(psi)*yhat
        ebx = dot_product(ebhat,xhat)
        eby = dot_product(ebhat,yhat)
        psi = atan2(-ebx,eby)

        !Get perp. p (only gyro)
        !Can use mu, because mu is already corrected for ExB
        phat = cos(psi)*xhat + sin(psi)*yhat
        absPp = sqrt(2*Model%m0*prt%Q(MUGC)*bAbs)
        pPerp = absPp*phat

        !Calculate FO properties
        Ro = cross(bhat,pPerp)/(Model%q0*bAbs) !Gyroradius
        !Calculate provisional FO momentum
        pFO = bhat*prt%Q(P11GC) + pPerp + Model%m0*vExB
        !Constrain magnitude of momentum to conserve energy
        pAbs = Model%m0*sqrt(prt%Q(GAMGC)**2.0 - 1.0) !Original |p|
        pFO = normVec(pFO)*pAbs

        !Set FO properties
        prt%isGC = .false.
        prt%Q(XPOS:ZPOS) = Rgc + Ro !FO position
        prt%Q(PXFO:PZFO) = pFO

        !Get new timestep
        prt%ddt = dtFO(Model,pFO,E,B)
    end subroutine gc2fo

    subroutine fo2gc(Model,prt,t,ebState)
        type(chmpModel_T), intent(in) :: Model
        type(prt_T), intent(inout) :: prt
        real(rp), intent(in) :: t
        type(ebState_T), intent(in)   :: ebState

        real(rp), dimension(NDIM) :: rFO,pFO,Rgc,B,E,bhat
        real(rp) :: Mu
        logical :: isConv

        !Get GC and Mu
        rFO = prt%Q(XPOS:ZPOS)
        pFO = prt%Q(PXFO:PZFO)
        
        Mu = CalcMu(rFO,pFO,t,Model,ebState,isConv,Rgc,B)
        if (.not. isConv) then
            write(*,*) 'Killing particle ', prt%id
            call KillParticle(Model,prt)
            return
        endif
        bhat = normVec(B)
        
        !Store values
        prt%isGC = .true.
        prt%Q(XPOS:ZPOS) = Rgc
        prt%Q(P11GC) = dot_product(bhat,pFO)
        prt%Q(GAMGC) = pv2Gam(pFO,Model%m0)
        prt%Q(MUGC) = Mu
        !Use FO first timestep (lazy)
        E = 0.0
        prt%ddt = dtFO(Model,pFO,E,B)
        
    end subroutine fo2gc
!---------------------------------
!Projection routines

    !Project to equator and calculate K/alpha
    function tpProject(prt,t,Model,ebState,isGood,tOpt) result(tpEq)
        type(prt_T), intent(in) :: prt
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in) :: t
        logical, intent(out) :: isGood
        integer, intent(out), optional :: tOpt

        real(rp) :: tpEq(NVAREQ) !New EQ values
        real(rp) :: B0,Beq,a0,aArg,dE
        real(rp), dimension(NDIM) :: x0,xeq,E,B,vExB,p0
        integer :: OCb,ijk0(NDIM)

        OCb = -1
        tpEq = 0.0
        isGood = .false.
        if (.not. prt%isIn) return

        !Trace local field line and save topology info
        x0 = prt%Q(XPOS:ZPOS)
        ijk0 = prt%ijk0
        call getMagEQ(Model,ebState,x0,t,xeq,Beq,OCb)
        if (present(tOpt)) then
            tOpt = OCb
        endif

        !Check topology and only project for *CLOSED* field lines
        if (OCb < 2) then
            !TP is on open/IMF line, return failure
            return
        endif

        !Find field at TP current location
        call ebFields(x0,t,Model,ebState,E,B,ijkO=ijk0,vExB=vExB)
        B0 = norm2(B) !Field strength at current location
        
        !Check for null fields
        if ( (B0<=TINY) .or. (Beq<=TINY) ) return

        !Check for forbidden region
        a0 = prt%alpha !Current pitch angle
        aArg = sqrt(Beq/B0)*abs(sin(a0))
        if (aArg > 1) return

        !If we're still here, this must be okay
        isGood = .true.

        !Now set values
        !Assuming magnetic equator close enough to Z=0
        tpEq(EQX:EQY) = xeq(XPOS:YPOS) 
        tpEq(EQTIME) = t
        !Get momentum in ExB and lab frame
        if (prt%isGC) then
            !Create momentum using random gyrophase
            p0 = ConP(prt%Q(GAMGC),prt%alpha,x0,B) - Model%m0*vExB
        else
            !FO particle is just momentum
            p0 = prt%Q(PXFO:PZFO) - Model%m0*vExB
        endif

        !Calculate kinetic energy in ExB frame and lab
        tpEq(EQKEV) = prt2kev(Model,prt)
        tpEq(EQKEB) = p2kev(Model,p0)

        !Last, calculate equatorial pitch angle
        tpEq(EQALP) = asin(aArg)
        if (tpEq(EQALP) < 0) tpEq(EQALP) = tpEq(EQALP) + PI

        contains
        
        !Create full momentum vector from gamma/alpha
        function ConP(gamma,alpha,r,B)
            real(rp), intent(in) :: gamma,alpha
            real(rp), dimension(NDIM), intent(in) :: r,B
            real(rp), dimension(NDIM) :: ConP

            real(rp), dimension(NDIM) :: bhat,xhat,yhat,pxy
            real(rp) :: pMag,p11,psi
            integer :: pSgn

            call MagTriad(r,B,xhat,yhat,bhat)
            pMag = Model%m0*sqrt(gamma**2.0 - 1.0)
            if (alpha <= PI/2) then
                pSgn = +1
            else
                pSgn = -1
            endif
            psi = genRand(0.0_rp,2*PI)

            !Get new momentum
            p11 = pSgn*pMag*sqrt( 1 - sin(alpha)**2.0 )
            pxy = pMag*sin(alpha)*( cos(psi)*xhat + sin(psi)*yhat )
            ConP = p11*bhat + pxy

        end function ConP
    end function tpProject

!---------------------------------
!Various particle book keeping routines
    subroutine KillParticle(Model,prt)
        type(chmpModel_T), intent(in) :: Model
        type(prt_T), intent(inout) :: prt

        if (doLoudKill) then
            !$OMP CRITICAL
            write(*,*) '-------'
            write(*,*) 'Killing particle ', prt%id
            write(*,*) 'Position = ', prt%Q(XPOS:ZPOS)
            write(*,*) 'Energy (keV) = ', prt2kev(Model,prt)
            write(*,*) 'Mu = ', prt%Q(MUGC)
            write(*,*) 'ddt = ', prt%ddt
            write(*,*) '-------'
            !$OMP END CRITICAL
        endif

        prt%isIn = .false.


    end subroutine KillParticle

    !Rotate 3D momentum (or any vector) to xy plane
    !Not projection, magnitude stays same
    subroutine p2xy(p)
        real(rp), intent(inout) :: p(NDIM)

        real(rp) :: MagP,pxy(NDIM)

        MagP = norm2(p)
        p(ZDIR) = 0.0
        pxy = normVec(p)*MagP
        p = pxy
    end subroutine p2xy
    
end module tputils