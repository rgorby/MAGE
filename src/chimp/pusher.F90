
module pusher
    use chmpdefs
    use tptypes
    use ebtypes
    use gridloc
    use pxing
    use gcutils
    use wpicalc
    implicit none

    real(rp), parameter :: dtX = 2.0 !Max increase in timestep, dtNew <= dtX*dtOld
    integer, parameter :: Nrk = 4 !Number of RK steps
    logical :: doGCKill = .FALSE. !Kill particles in GC that violate adiabaticity
    logical :: doRK4 = .true. !Do RK4 or RK2
    !logical, private, parameter :: doBoris = .true. !Do Boris vs. Higuera-Cary pusher
    
    contains
    !Advance particle prt using ebState fields
    !Advance from t->t+dt
    subroutine PushTP(prt,t,dt,Model,ebState)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(prt_t), intent(inout) :: prt
        real(rp), intent(in) :: t,dt
        
        type(prt_t) :: oprt
        real(rp) :: dtCum,dtRem,ddt
        real(rp) :: t1,t2
        logical :: isGood,doKill

        dtCum = 0.0 !How far we've advanced
        do while ( (dtCum<dt) .and. prt%isIn )

        !---------------------------
        !Setup substep
            !Save old particle state
            oprt = prt
            !Use precalculated particle timestep
            ddt = prt%ddt
            dtRem = dt-dtCum
            if (dtRem < ddt) ddt = dtRem !Don't overshoot

        !---------------------------
        !Take step
            !Split based on GC/FO
            if (prt%isGC) then
                isGood = StepGC(prt,t+dtCum,ddt,Model,ebState)
            else
                isGood = StepFO(prt,t+dtCum,ddt,Model,ebState)
            endif
            if (isGood) dtCum = dtCum + ddt

        !---------------------------
        !Do per substep operations here

            !Check for lost
            if (prt%isIn) prt%isIn = inDomain(prt%Q(XPOS:ZPOS),Model,ebState%ebGr)

            !Check for suicide conditions
            doKill = ( prt%ddt <= TINY ) .or. ( prt2kev(Model,prt)<=Model%MinK )
            if (doKill) then
                call KillParticle(Model,prt)
            endif

            !Check for wave particle interaction
            if ( Model%doWPI .and. isGood .and. prt%isIn) then
                call AdvanceWPI(oprt,prt,t+dtCum,Model,ebState)
            endif

            !Check for equatorial crossing
            if (isGood .and. prt%isIn) then     
                !OldP/t1 and prt/t2
                t2 = t+dtCum
                t1 = t2-ddt
                call ChkEQX(oprt,prt,t1,t2,Model,ebState)
            endif
            !Try to upgrade FO->GC
            if ( Model%isDynamic .and. (.not. prt%isGC) .and. prt%isIn) then
                !Test adiabaticity (and update alpha)
                call Upgrade2GC(prt,t+dtCum,Model,ebState)
            endif

        enddo
    !---------------------------
    !Do per step operations here

    end subroutine PushTP

    !Both substep routines
    !isGood = SubstepXX(prt,t,ddt,m0,q0,Model,ebState)

    function StepGC(prt,t,ht,Model,ebState) result(isGood)
        type(prt_T), intent(inout) :: prt
        real(rp), intent(in) :: t,ht
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState

        logical :: isGood,isOuts(Nrk)
        integer :: n, ijk(NDIM)
        real(rp) :: htOld,htNew, p11,pMag, eGCs(Nrk)
        real(rp), dimension(NVARTP) :: xGC,dxGC1,dxGC2,dxGC3,dxGC4,dQ,oQ

        xGC = prt%Q !GC configuration
        ijk = prt%ijk0

        dxGC1 = ht*DerivGC(xGC          ,t       ,ijk,Model,ebState,eGCs(1),isOuts(1))
        if (doRK4) then
            !Do RK4
            dxGC2 = ht*DerivGC(xGC+0.5*dxGC1,t+0.5*ht,ijk,Model,ebState,eGCs(2),isOuts(2))
            dxGC3 = ht*DerivGC(xGC+0.5*dxGC2,t+0.5*ht,ijk,Model,ebState,eGCs(3),isOuts(3))
            dxGC4 = ht*DerivGC(xGC+    dxGC3,t+    ht,ijk,Model,ebState,eGCs(4),isOuts(4),htGC=htNew)
        else
            !Do RK2
            dxGC2 = ht*DerivGC(xGC+0.5*dxGC1,t+0.5*ht,ijk,Model,ebState,eGCs(2),isOuts(2),htGC=htNew)
            eGCs(3:4) = 0.0
            isOuts(3:4) = .false.
        endif

        !Check if particle left during RK substeps
        !Mark particle as out, don't use any steps
        if (any(isOuts)) then
            isGood = .false.
            prt%isIn = .false.
            return
        endif

        !Test for good step
        if (maxval(eGCs)>=Model%epsgc) then
            !Field variation is too large for GC, throw away everything
            isGood = .false.
            if (Model%isDynamic) then
                !Flip to FP
                !write(*,*) 'Downgrading particle ', prt%id
                call gc2fo(Model,prt,t,ebState)
                return
            else if (doGCKill) then
                !GC only integrator & kill particle
                call KillParticle(Model,prt)
                return
            endif

        endif


        !If still here then the step was good (or ignoring it)
        isGood = .true.

        !If still here, perform update
        associate(Q=>prt%Q)
        oQ = Q

        if (doRK4) then
            dQ = dxGC1/6 + dxGC2/3 + dxGC3/3 + dxGC4/6
            
        else
            dQ = dxGC2
        endif
        Q = Q + dQ

    !Test update: gamma
        if ( Q(GAMGC) <= (1.0 + TINY) ) then
            !Bad update, redo with smaller timestep
            Q = oQ
            prt%ddt = 0.5*prt%ddt            
            isGood = .false.
            return
            !call FixGC(prt,t+ht,Model,ebState)
        endif

    !Test update: p11
        p11 = Q(P11GC)
        pMag = Model%m0*sqrt( Q(GAMGC)**2.0 - 1.0 )

        if ( abs(p11) > pMag ) then
            !Constrain p11
            if (p11 >= 0.0) then
                p11 = +pMag
            else
                p11 = -pMag
            endif
            Q(P11GC) = p11

        endif

        !Otherwise finish update
        prt%Ngc = prt%Ngc+1

        !Update timestep
        htOld = prt%ddt
        prt%ddt = min(Model%epsht*htNew,dtX*htOld)

        !Update pitch angle
        prt%alpha = acos( p11/max(pMag,abs(p11)) )

        end associate

        if (Model%do2D) then
            write(*,*) '2D Integration not implemented for GC yet'
            stop
        endif
    end function StepGC
    
    !Full orbit integrator, synchronized integrator w/ matching x/v time states (ie not leapfrog)
    !See excellent description by Ripperda++ 2018, 10.3847/1538-4365/aab114

    !General form
    !1: x^n+1/2 = x^n + u^n/gamma^n x dt/2
    !2: Evaluate E,B at x^n+1/2,t^n+1/2
    !3: u^- = u^n + q dt/2m E^n+1/2
    !4: Rotation to get u^+ and u^n+1
    !5: x^n+1 = x^n+1/2 + u^n+1 / gamma^n+1 x dt/2

    function StepFO(prt,t,ht,Model,ebState) result(isGood)
        type(prt_T), intent(inout) :: prt
        real(rp), intent(in) :: t,ht
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState

        logical :: isGood
        real(rp), dimension(NDIM) :: X_n,X_n12,X_np1
        real(rp), dimension(NDIM) :: P_n,P_n12,P_np1
        real(rp), dimension(NDIM) :: U_n,U_m,U_p,U_np1
        real(rp), dimension(NDIM) :: E_n12,B_n12,tvec,svec,tauvec
        real(rp) :: Gam_n,Gam_np1,Gam_m,Gam_p,s,tau2,ht2
        real(rp) :: Ust,Sig,htNew,htOld

        associate(Q=>prt%Q,q0=>Model%q0,m0=>Model%m0)
    !0: Setup
        X_n = Q(XPOS:ZPOS)
        P_n = Q(PXFO:PZFO)
        Gam_n = p2Gam(P_n,m0)
        U_n = P_n/m0

        ht2 = 0.5*ht

    !1: 
        !Advance position to half time
        X_n12 = X_n + (U_n/Gam_n)*ht2
    !2:
        call ebFields(X_n12,t+ht2,Model,ebState,E_n12,B_n12,ijkO=prt%ijk0)
    !3:
        U_m = U_n + (q0/m0)*E_n12*ht2
    !4:
        Gam_m = sqrt(1 + dot_product(U_m,U_m))

        if (doBoris) then
            tvec = B_n12 * (q0/m0) * ht2 / Gam_m
            svec = 2*tvec/(1+dot_product(tvec,tvec))
            U_p = U_m + cross( U_m + cross(U_m,tvec), svec )
            U_np1 = U_p + (q0/m0)*ht2*E_n12
            !Gam_p = sqrt(1 + U_p**2.0) !Don't actually need this
        else
            !Do Higuera-Cary
            tauvec = (q0/m0)*B_n12*ht2
            Ust = dot_product(U_m,tauvec)
            tau2 = dot_product(tauvec,tauvec)
            Sig = Gam_m**2.0 - tau2
            Gam_p = sqrt( Sig + sqrt(Sig**2.0 + 4*(tau2 + Ust**2.0)) )/sqrt(2.0)
            tvec = tauvec/Gam_p
            s = 1.0/(1.0+dot_product(tvec,tvec))

            U_p = s*( U_m + dot_product(U_m,tvec)*tvec + cross(U_m,tvec) )
            U_np1 = U_p + (q0/m0)*ht2*E_n12 + cross(U_m,tvec) !Note extra term relative to Boris

        endif

    !5:
        P_np1 = m0*U_np1
        Gam_np1 = p2Gam(P_np1,m0)
        X_np1 = X_n12 + (U_np1/Gam_np1)*ht2
        

    !Finish up
        !Store update
        isGood = .true.
        Q(XPOS:ZPOS) = X_np1
        Q(PXFO:PZFO) = P_np1
        prt%Nfo = prt%Nfo+1

        !Get new timestep/alpha
        !NOTE: Using fields from half-step to avoid extra field evaluation
        !Use half substep estimation of momentum in Lorentz force
        !Also using half-substep lag for calculating alpha (doesn't feed back)
        P_n12 = 0.5*( P_n + P_np1 )
        htNew = dtFO(Model,P_n12,E_n12,B_n12)
        prt%alpha = angVec(B_n12,P_n12) !Half substep lag

        !Store timestep (dtFO multiplies by epsdt)
        htOld = prt%ddt
        prt%ddt = min(htNew,dtX*htOld)

        !We're done here, will test for FO->GC upgrade upstairs
        end associate

    end function StepFO

    ! !Full orbit integrator, Boris push w/ matching x/v time states (ie not leapfrog)
    ! !x^* = x^n + 0.5*h v^n
    ! !v^- = v^n + 0.5*h E^*
    ! !w/ Q^* = Q(x=x*,t=t^n+1/2), interpolant at half-step position/time
    ! !v^+ = Boris[v^-,B^*], ie Boris rotation
    ! !v^n+1 = v^+ + 0.5*h E^*
    ! !x^n+1 = x^* + 0.5*h*v^n+1

    ! function StepFO(prt,t,ht,Model,ebState) result(isGood)
    !     type(prt_T), intent(inout) :: prt
    !     real(rp), intent(in) :: t,ht
    !     type(chmpModel_T), intent(in) :: Model
    !     type(ebState_T), intent(in)   :: ebState

    !     logical :: isGood
    !     real(rp), dimension(NDIM) :: r0,rS,p0,E,B
    !     real(rp), dimension(NDIM) :: pMin,pPos,pF,rF,pHf
    !     real(rp) :: gamma,ht2,htNew,htOld

    !     associate(Q=>prt%Q,q0=>Model%q0,m0=>Model%m0)
    !     r0 = Q(XPOS:ZPOS)
    !     p0 = Q(PXFO:PZFO)

    !     ht2 = 0.5*ht
    !     !Advance position to half time
    !     rS = r0 + ht2*p0/(m0*p2Gam(p0,m0))

    !     !Get fields at rS,t+h/2
    !     call ebFields(rS,t+ht2,Model,ebState,E,B,ijkO=prt%ijk0)

    !     !Do half E
    !     pMin = p0 + q0*ht2*E

    !     !Do Boris rotation for B field
    !     pPos = Boris(Model,pMin,B,ht)

    !     !Finish E push
    !     pF = pPos + q0*ht2*E

    !     if (Model%do2D) then
    !         !Force momentum to xy plane
    !         call p2xy(pF)
    !     endif

    !     !Finish position push
    !     rF = rS + ht2*pF/(m0*p2Gam(pF,m0))

    ! !Store update
    !     isGood = .true.
    !     Q(XPOS:ZPOS) = rF
    !     Q(PXFO:PZFO) = pF
    !     prt%Nfo = prt%Nfo+1

    ! !Get new timestep/alpha
    !     !NOTE: Using fields from half-step to avoid extra field evaluation
    !     !Use half substep estimation of momentum in Lorentz force
    !     !Also using half-substep lag for calculating alpha (doesn't feed back)
    !     pHf = 0.5*(pF+p0)
    !     htNew = dtFO(Model,pHf,E,B)
    !     prt%alpha = angVec(B,pHf) !Half substep lag

    !     !Store timestep (dtFO multiplies by epsdt)
    !     htOld = prt%ddt
    !     prt%ddt = min(htNew,dtX*htOld)

    !     !We're done here, will test for FO->GC upgrade upstairs
        
    !     end associate
    ! end function StepFO

    ! !Boris rotation, ie p^- -> p^+
    ! !Note, using relativistic formulation where rotation angle depends on gamma
    ! !See chapter 15 of Birdsall & Langdon
    ! function Boris(Model,pMin,B,ht) result(pPos)
    !     type(chmpModel_T), intent(in) :: Model
    !     real(rp), dimension(NDIM), intent(in) :: pMin,B
    !     real(rp), intent(in) :: ht
    !     real(rp), dimension(NDIM) :: pPos
    !     real(rp) :: gamma
    !     real(rp), dimension(NDIM) :: Tv,Sv,pbi

    !     !Use updated gamma to account for first half of E push
    !     gamma = p2Gam(pMin,Model%m0)

    !     Tv = 0.5*ht*Model%q0*B/(Model%m0*gamma)
    !     Sv = 2*Tv/(1+dot_product(Tv,Tv))

    !     pbi = pMin + cross(pMin,Tv) !Bisector
    !     pPos = pMin + cross(pbi,Sv)

    ! end function Boris

end module pusher

