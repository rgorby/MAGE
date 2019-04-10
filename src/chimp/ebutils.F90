!Various field manipulation utilities
module ebutils
    use chmpdefs
    use ebtypes
    use math
    
    implicit none

    contains

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
        if (norm2(bhat)<1) then
            !Something went wrong, force to zhat
            bhat = [0.0_rp,0.0_rp,1.0_rp]
        endif
        yhat = normVec( cross(r,bhat) )
        xhat = cross(yhat,bhat)

    end subroutine MagTriad

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

end module ebutils