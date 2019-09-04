!Various routines to handle metrics at particle crossings
module wpi
    use chmpdefs
    use tptypes
    use ebtypes
    use tputils
    
    implicit none

    contains
    !Main subroutine for wave particle inteactions
    subroutine PerformWPI(prt,t,tau,Model,ebState)
        type(prt_t), intent(inout) :: prt
        real(rp), intent(in) :: t,tau
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        
        real(rp), dimension(NDIM) :: r,p,E,B,xhat,yhat,bhat
        real(rp), dimension(NDIM) :: p11,pxy,vExB
        real(rp) :: gamma,ebGam,aNew,pNew,pMag,MagB,Mu,p11Mag
        integer :: pSgn=1

        !Change in pitch angle and momentum due to different waves 
        real(rp) :: dAwhistler=0.0,dPwhistler=0.0

        if (Model%do2D) then
            !Trap here and quickly grab values
            prt%Qeq(EQX:EQY) = prt%Q(XPOS:YPOS)
            prt%Qeq(EQTIME)  = t
            prt%Qeq(EQKEV)   = prt2kev(Model,prt)
            prt%Qeq(EQALP)   = prt%alpha
            prt%Qeq(EQKEB)   = 0.0
            return
        endif

        !Check to see if waves are present
        call ChkWhistler(Model,ebState,prt,tau,dAwhistler,dPwhistler)

        !Get local coordinate system
        r = prt%Q(XPOS:ZPOS)
        call ebFields(r,t,Model,ebState,E,B,ijkO=prt%ijk0,vExB=vExB)
        call MagTriad(r,B,xhat,yhat,bhat)
        MagB = max(norm2(B),TINY)

        !Update the pitch angle 
        aNew = prt%alpha + dAwhistler
        prt%alpha = aNew

        !Updating the momentum vairables
        !FIXME: dp should be 3 dimensional, need to update momentum/energy components
        !This is calculating new momentum from alpha, not necessary when have dp
        !Look into pSgn to see if it is necessary in this application
        prt%Q(PSITP) = pNew
        if (prt%isGC) then
            !Scatter GC particle
            !Energy is unchanged, p11 and Mu change
            gamma = prt%Q(GAMGC)
            ebGam = max(1.0,gamma - 0.5*dot_product(vExB,vExB))
            pMag = Model%m0*sqrt(ebGam**2.0 - 1.0)
            p11Mag = pSgn*pMag*sqrt( 1 - sin(aNew)**2.0 )
            Mu = ( (Model%m0**2)*(ebGam**2.0 - 1.0) - p11Mag*p11Mag) / (2*Model%m0*MagB)

            prt%Q(P11GC) = p11Mag
            prt%Q(MUGC ) = Mu 
            prt%Q(GAMGC) = gamma
        else
            p = prt%Q(PXFO:PZFO)
            gamma = pv2Gam(p,Model%m0)
            pMag = Model%m0*sqrt(gamma**2.0 - 1.0)
            
            !Get new momentum
            p11 = pSgn*pMag*sqrt( 1 - sin(aNew)**2.0 )
            pxy = pMag*sin(aNew)*( cos(pNew)*xhat + sin(pNew)*yhat )

            p = p11*bhat + pxy
            prt%Q(PXFO:PZFO) = p
            
        endif

    end subroutine PerformWPI

    !Check to see if particle interacts with whistler waves
    subroutine ChkWhistler(Model,ebState,prt,tau,dalpha,dp)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(prt_t), intent(inout) :: prt
        real(rp), intent(in) :: tau
        real(rp), intent(inout) :: dalpha,dp

        logical :: doWave = .false. !Particle is in location where waves are occurring 

        !See if particle meets criteria for whistler waves
        !FIXME: for now have waves everywhere
        doWave = .true. 

        !perform wave particle interaction
        if (doWave) then 
            !Determine the wave that partilce is in resonance with
            !FIXME: need to add this when include a form for Daa

            !Calculate changes to pitch angle and momentum
            call WhistlerWPI(Model,ebState,prt,tau,dalpha,dp)
        endif

    end subroutine ChkWhistler

    subroutine WhistlerWPI(Model,ebState,prt,tau,dalpha,dp)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(prt_t), intent(inout) :: prt
        real(rp), intent(in) :: tau
        real(rp), intent(inout) :: dalpha,dp

        real(rp) :: Daa,rand
        integer :: aSgn

        !Calculate Daa from the local plasma
        Daa = 1.0

        !Calculate dAlpha from Daa
        dalpha = sqrt(2.0*tau*Daa)

        !Randomly determine sign of dAlpha
        rand=genRand(0.0_rp,1.0_rp)
        if (rand > 0.5) then
            aSgn = +1
        else
            aSgn = -1
        endif 

        dalpha = aSgn*dalpha

        !Get Change in momentum
        dp = 0.0

    end subroutine WhistlerWPI

end module wpi
