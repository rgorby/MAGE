!Various routines necessary for guiding center formulation
!Main work is calculating time derivatives of GC variables

module gcutils
    use chmpdefs
    use math
    use tptypes
    use ebtypes
    use ebinterp
    use tputils
    implicit none

    contains

    !Calculates time derivatives ("velocity") of GC variables
    function DerivGC(xGC,t,ijk,Model,ebState,epsGC,isOut,htGC) result(vGC)
        type(chmpModel_T), intent(in) :: Model
        real(rp), intent(in) :: t,xGC(NVARTP)
        integer, intent(inout) :: ijk(NDIM)
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(out) :: epsGC
        logical , intent(out) :: isOut
        real(rp), intent(out), optional :: htGC

        real(rp) :: vGC(NVARTP)

        type(gcFields_T) :: gcFields
        real(rp), dimension(NDIM) :: r,E,B,vExB,vxb,bhatdot,cbhat
        real(rp), dimension(NDIM) :: bhat,gMagB,DotVeb,Bstar,Estar,gKinEB
        real(rp), dimension(NDIM,NDIM) :: JacVeb
        real(rp) :: MagB,Bst11,m0,q,Mu,p11,gamma,ebGam
        real(rp) :: htgcT,htgcB

        vGC = 0.0

    !Get atomic field information
        !Start by getting all necessary fields
        r = xGC(XPOS:ZPOS) !GC position
        isOut = .not. inDomain(r,Model,ebState%ebGr)
        if (isOut) then
            return
        endif
        
        !If still here, then point is in domain
        call ebFields(r,t,Model,ebState,E,B,ijk,vExB,gcFields)

    !Do main calculation for star fields
        associate( JacB=>gcFields%JacB, JacE=>gcFields%JacE,DotB=>gcFields%DotB,DotE=>gcFields%DotE )

        !Main TP parameters
        m0 = Model%m0
        q  = Model%q0
        Mu  = xGC(MUGC)
        p11 = xGC(P11GC)
        gamma = max(1.0,xGC(GAMGC))

        !Main derived field parameters
        MagB = norm2(B)
        bhat = normVec(B)

        if (MagB < TINY) then
            !If field is too weak then GC can't work
            isOut = .true.
            return
        endif
        !TP/EB params
        ebGam = gamma - 0.5*dot_product(vExB,vExB) !c=1
        ebGam = max(1.0,ebGam)

        !Start getting derivative terms
        !gMagB = gradient(|B|), vector
        !      = bhat \cdot \grad \vec{B}
        gMagB = VdT(bhat,JacB)

        !Curl of bhat
        !Curl(bhat) = (1/MagB)*( Curl(B) + bhat \cross gMagB)
        cbhat = (1/MagB)*( Jac2Curl(JacB) + cross(bhat,gMagB) )

        !Now get Jacobian of ExB velocity, JacVeb
        JacVeb = ( (MagB**2.0)*( matmul(xMat(E),JacB) - matmul(xMat(B),JacE) ) &
               - 2*Dyad(VdT(B,JacB),cross(E,B)) )/(MagB**4.0)

        !Prep for E*
        !Gradient of ExB kinetic energy
        gKinEB = m0*matmul(transpose(JacVeb),vExB)

        !Time deriv of ExB velocity
        DotVeb = ( (MagB**2.0)*(cross(DotE,B) + cross(E,DotB)) &
               - 2*dot_product(B,DotB)*cross(E,B) )/(MagB**4.0)

        bhatdot = (1/MagB)*(DotB - dot_product(DotB,bhat)*bhat)

        !Calculate B*
        !B* = B + (c/q)*p11*curl(bhat) + (c/q)*curl(Peb)
        Bstar = B + (p11/q)*cbhat + (m0/q)*Jac2Curl(JacVeb)

        !Calculate E*
        Estar = E - Mu*gMagB/(q*ebGam) &
              - gKinEB/q               &
              - p11*bhatdot/q          &
              - m0*DotVeb/q

    !Turn star fields into time derivatives
        Bst11 = dot_product(bhat,Bstar) !Parallel comp. of B*

        vGC(XPOS:ZPOS) = cross(Estar,bhat)/Bst11    &
                       + p11*Bstar/(m0*ebGam*Bst11)
        vGC(P11GC)     = q*dot_product(Bstar,Estar)/Bst11
        vGC(GAMGC)     = p11*vGC(P11GC)/(m0*m0*ebGam)         &
                       + Mu*dot_product(DotB,bhat)/(m0*ebGam) &
                       + dot_product(DotVeb,vExB)             &
                       + (1/m0)*dot_product(vGC(XPOS:ZPOS), (Mu*gMagB/ebGam) + gKinEB)

    !Finalize
        !Calculate adiabaticity parameter
        !bScl = sum(JacB**2.0)/(MagB**3.0)
        !epsGC = sqrt(2*m0*Mu*bScl)/abs(q)

        epsGC = eGC(Model,B,JacB,Mu)

        !Calculate new timestep if variable is present
        if ( present(htGC) ) then
            !Timestep, dt ~ |p|/F, F = Lorentz(vGC) + dP11/dt
            htgcT = m0*sqrt(gamma**2.0-1) ! |P|

            vxb = q*cross(vGC(XDIR:ZDIR),B) !Lorentz force on GC velocity
            htgcB = sqrt( norm2(vxb)**2.0 + norm2(q*E)**2.0 + vGC(P11GC)**2.0 )
            htGC = htgcT/max(htgcB,TINY)

        endif

        end associate

    end function DerivGC

    !For an FO particle, test adiabaticity, update alpha and upgrade to GC if possible
    subroutine Upgrade2GC(prt,t,Model,ebState)
        type(prt_T), intent(inout) :: prt
        real(rp), intent(in) :: t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState

        type(gcFields_T) :: gcFields
        real(rp), dimension(NDIM) :: r,p,E,B,vExB
        real(rp) :: eAdb

        if (prt%isGC) then
            write(*,*) 'This should not happen, Upgrade2GC w/ GC'
            stop
        endif

        r = prt%Q(XPOS:ZPOS)
        p = prt%Q(PXFO:PZFO)

        !Get fields and calculate adiabaticity
        call ebFields(r,t,Model,ebState,E,B,prt%ijk0,vExB,gcFields)
        eAdb = eGC(Model,B,gcFields%JacB,p)

        !Update alpha no matter what
        prt%alpha = angVec(B,p)

        !Test adiabaticity and upgrade if possible
        if (eAdb<Model%epsgc) then
            !write(*,*) 'Upgrading particle ', prt%id
            !Flip to GC
            call fo2gc(Model,prt,t,ebState)
        endif
        
    end subroutine Upgrade2GC

    !Fix degenerate gamma in GC
    !Reconstruct gamma from momenta and ExB drift
    subroutine FixGC(prt,t,Model,ebState)
        type(prt_T), intent(inout) :: prt
        real(rp), intent(in) :: t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState

        real(rp) :: ebGam2,pMag,vxb,MagF,dt
        real(rp), dimension(NDIM) :: r,E,B,vExB

        associate(Q=>prt%Q)

        r = Q(XPOS:ZPOS)
        call ebFields(r,t,Model,ebState,E,B,prt%ijk0,vExB)

        ebGam2 = 1.0 + (Q(P11GC)/Model%m0)**2.0 + 2*Q(MUGC)*norm2(B)/Model%m0
        Q(GAMGC) = sqrt(ebGam2) + 0.5*dot_product(vExB,vExB)

        !Set new timestep
        pMag = Model%m0*sqrt( Q(GAMGC)**2.0 - 1.0 )
        vxb = pMag*norm2(B)/(Model%m0*Q(GAMGC))
        MagF = Model%q0*( norm2(E) + vxb)
        dt = pMag/MagF

        prt%ddt = max(Model%epsht*dt,TINY)

        end associate
    end subroutine FixGC
end module gcutils
