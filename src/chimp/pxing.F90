!Various routines to handle metrics at particle crossings
module pxing
    use chmpdefs
    use tptypes
    use ebtypes
    use tputils
    
    implicit none

    real(rp) :: zEq = 0.0 !Defined Z for equator

    contains
    !Checks for equatorial crossing and saves EQX data into nPrt
    !oPrt,oT particle data/time of previous state
    subroutine ChkEQX(oPrt,nPrt,oT,nT,Model,ebState)
        type(prt_t), intent(in)    :: oPrt
        type(prt_t), intent(inout) :: nPrt
        real(rp), intent(in) :: oT,nT
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState

        real(rp) :: oZp,nZp,dZ,w1,w2,req
        real(rp), dimension(NDIM) :: xeq,E,B,vExB,pExB

        if (Model%do2D) then
            !Trap here and quickly grab values
            nPrt%Qeq(EQX:EQY) = nPrt%Q(XPOS:YPOS)
            nPrt%Qeq(EQTIME)  = nT
            nPrt%Qeq(EQKEV)   = prt2kev(Model,nPrt)
            nPrt%Qeq(EQALP)   = nPrt%alpha
            nPrt%Qeq(EQKEB)   = 0.0
            return
        endif
        
        oZp = oPrt%Q(ZPOS)-zEq
        nZp = nPrt%Q(ZPOS)-zEq
        if (oZp*nZp <=0) then
        !Have a crossing
            req = 0.5*( norm2(oPrt%Q(XPOS:YPOS)) + norm2(nPrt%Q(XPOS:YPOS)) )
            !Do PA scattering first if necessary
            if (Model%doEQScat .and. (req <= Model%reqScat) ) then
                call PAScat(Model,ebState,nPrt,nT)
            endif

            !Get weights from Z distance for interpolation
            if ( max(abs(oZp),abs(nZp)) > TINY ) then
                dZ = abs(oZp-nZp)
                w1 = abs(nZp)/dZ
                w2 = abs(oZp)/dZ
            else
                dZ = 1.0
                w1 = 0.0
                w2 = 1.0
            endif

            !Interpolate quantities at crossing
            nPrt%Qeq(EQX:EQY) = w1*oPrt%Q(XPOS:YPOS) + w2*nPrt%Q(XPOS:YPOS)
            nPrt%Qeq(EQTIME)  = w1*oT + w2*nT
            nPrt%Qeq(EQKEV)   = w1*prt2kev(Model,oPrt) + w2*prt2kev(Model,nPrt)
            nPrt%Qeq(EQALP)   = w1*oPrt%alpha + w2*nPrt%alpha

            !Handle ExB correction to particle energy
            !FIXME: Need to handle cases where nPrt and oPrt have different isGC
            if (Model%imeth == IFO) then
                !Interp exb velocity at crossing
                xeq = [nPrt%Qeq(EQX),nPrt%Qeq(EQX),zEq]
                call ebFields(xeq,nPrt%Qeq(EQTIME),Model,ebState,E,B,ijkO=nPrt%ijk0,vExB=vExB)
                pExB = Model%m0*vExB
                nPrt%Qeq(EQKEB) = w1*p2kev(Model,oPrt%Q(PXFO:PZFO)-pExB) + w2*p2kev(Model,nPrt%Q(PXFO:PZFO)-pExB)

            else
                !FIXME: For now not bothering correcting for GC electrons
                nPrt%Qeq(EQKEB)   = nPrt%Qeq(EQKEV)
            endif

        endif !Found EQX
        
    end subroutine ChkEQX

    !Randomly change pitch angle
    subroutine PAScat(Model,ebState,prt,t)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(prt_t), intent(inout) :: prt
        real(rp), intent(in) :: t
        
        real(rp), dimension(NDIM) :: r,p,E,B,xhat,yhat,bhat
        real(rp), dimension(NDIM) :: p11,pxy,vExB
        real(rp) :: gamma,ebGam,aNew,pNew,pMag,MagB,Mu,p11Mag
        integer :: pSgn

        !Get local coordinate system
        r = prt%Q(XPOS:ZPOS)
        call ebFields(r,t,Model,ebState,E,B,ijkO=prt%ijk0,vExB=vExB)
        call MagTriad(r,B,xhat,yhat,bhat)
        MagB = max(norm2(B),TINY)

        !Get new alpha/psi
        aNew = genRand(0.0_rp,1*PI)
        pNew = genRand(0.0_rp,2*PI)

        prt%alpha = aNew

        if (aNew <= PI/2) then
            pSgn = +1
        else
            pSgn = -1
        endif

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

    end subroutine PAScat

end module pxing