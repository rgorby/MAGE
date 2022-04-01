!Various field manipulation utilities
module ebutils
    use chmpdefs
    use ebtypes
    use math
    
    implicit none

    contains

    !Jacobian(ADIR,BDIR) = dF_{A}/dB
    
    !Turn Jacobian matrix into curl vector
    function Jac2Curl(JacA) result(CurlA)
        real(rp), dimension(NDIM,NDIM), intent(in) :: JacA
        real(rp), dimension(NDIM) :: CurlA

        CurlA(XDIR) = JacA(ZDIR,YDIR) - JacA(YDIR,ZDIR)
        CurlA(YDIR) = JacA(XDIR,ZDIR) - JacA(ZDIR,XDIR)
        CurlA(ZDIR) = JacA(YDIR,XDIR) - JacA(XDIR,YDIR)
    end function Jac2Curl

    !Calculate divergence of field from Jacobian
    function Jac2Div(JacA) result(DivA)
        real(rp), dimension(NDIM,NDIM), intent(in) :: JacA
        real(rp) :: DivA

        DivA = JacA(XDIR,XDIR) + JacA(YDIR,YDIR) + JacA(ZDIR,ZDIR)
        
    end function Jac2Div
    !Get local magnetic triad from B
    subroutine MagTriad(r,B,xhat,yhat,bhat)
        real(rp), dimension(NDIM), intent(in) :: r,B
        real(rp), dimension(NDIM), intent(out) :: xhat,yhat,bhat

        bhat = normVec(B)
        if ( abs(norm2(bhat)-1.0) > TINY ) then
            !Something went wrong, force to zhat
            bhat = [0.0_rp,0.0_rp,1.0_rp]
        endif
        yhat = normVec( cross(r,bhat) )
        xhat = cross(yhat,bhat)

    end subroutine MagTriad

    !Return radius of curvature in code length based on B and JacB
    !Assuming B and jB are both in code units
    function getRCurv(B,jB) result(rcurv)
        real(rp), intent(in), dimension(NDIM)      :: B
        real(rp), intent(in), dimension(NDIM,NDIM) :: jB

        real(rp) :: MagB, invrad, rcurv
        real(rp), dimension(NDIM,NDIM) :: Jacbhat
        real(rp), dimension(NDIM) :: bhat, gMagB

        bhat = normVec(B)
        MagB = norm2(B)

        !Start getting derivative terms
        !gMagB = gradient(|B|), vector
        !      = bhat \cdot \grad \vec{B}
        gMagB = VdT(bhat,jB)

        Jacbhat = ( MagB*jB - Dyad(gMagB,B) )/(MagB*MagB)

        invrad = norm2(VdT(bhat,Jacbhat))
        if (invrad>TINY) then
           rcurv = 1.0/invrad
        else
           rcurv = -TINY
        endif
    end function getRCurv


!Split momentum into 3 components: p11,pPerp,pExB
    subroutine SplitP(Model,p,E,B,vExB,p11,pPerp)
        type(chmpModel_T), intent(in) :: Model
        real(rp), dimension(NDIM), intent(in) :: p,E,B,vExB
        real(rp), intent(out) :: p11, pPerp(NDIM)

        real(rp), dimension(NDIM) :: bhat

        bhat = normVec(B)
        p11 = dot_product(p,bhat)
        pPerp = p - bhat*p11 - Model%m0*vExB

    end subroutine SplitP

    
    !Guarantee E.B = 0
    subroutine CleanE(E,B)
        real(rp), intent(inout) :: E(NDIM)
        real(rp), intent(in) :: B(NDIM)

        real(rp) :: EdB,BdB
        BdB = dot_product(B,B)
        EdB = dot_product(E,B)
        if (BdB >= TINY) then
            E = E - EdB*B/BdB
        endif
    end subroutine CleanE

    !Get ExB drift velocity
    function EBDrift(E,B) result(vExB)
        real(rp), dimension(NDIM), intent(in) :: E,B
        real(rp), dimension(NDIM) :: vExB

        real(rp) :: BdB

        BdB = dot_product(B,B)
        if (BdB>=TINY) then
            vExB = cross(E,B)/BdB
        else
            vExB = 0.0
        endif

    end function EBDrift
    
    !Takes i,j,k cell index and returns active cell ip,jp,kp of mirror
    !Map in i,k,j order
    subroutine ijk2Active(Model,Grid,i,j,k,ip,jp,kp)
        type(chmpModel_T), intent(in)    :: Model
        type(   ebGrid_T), intent(in)    :: Grid
        integer, intent(in) :: i,j,k
        integer, intent(out) :: ip,jp,kp

        integer :: Np,Np2

        Np  = Grid%Nkp
        Np2 = Grid%Nkp/2

        !Start w/ i index, do mirror back into active
        if (i < Grid%is) then
            ip = Grid%is + (Grid%is-i) - 1
        elseif (i > Grid%ie) then
            ip = Grid%ie - (i-Grid%ie) + 1
        else
            ip = i
        endif

        !Next do k, map via periodicity
        if (k < Grid%ks) then
            kp = Grid%ke - (Grid%ks-k) + 1
        elseif (k > Grid%ke) then
            kp = Grid%ks + (k-Grid%ke) - 1
        else
            kp = k
        endif

        !Finally do j
        if (j < Grid%js) then
            jp = Grid%js + (Grid%js-j) - 1
            kp = k+Np2
            if (kp>Np) kp = kp-Np
        elseif (j > Grid%je) then
            jp = Grid%je - (j-Grid%je) + 1
            kp = k+Np2
            if (kp>Np) kp = kp-Np
        else
            jp = j
        endif
    end subroutine ijk2Active
end module ebutils