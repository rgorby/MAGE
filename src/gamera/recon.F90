module recon
    use gamtypes
    use gamutils
    use math
    implicit none

    !GetLR_T
    !Performs stencil->LR calculation
    !Q = n+1/2 state to reconstruct on
    !nQ = n state to limit on
    abstract interface
        subroutine GetLR_T(dV,nQ,Q,Vi,Ql,Qr)
            import :: rp,vecLen,recLen,limLen
            real(rp), intent(in), dimension(vecLen,limLen) :: nQ
            real(rp), intent(in), dimension(vecLen,recLen) :: dV,Q
            real(rp), intent(in), dimension(vecLen)  :: Vi
            real(rp), intent(inout), dimension(vecLen) :: Ql,Qr
        end subroutine GetLR_T
    end interface

    !Note, be careful how interpWgt is initialized to ensure sum(:) = 1 exactly
    real(rp), dimension(recLen), parameter :: interpWgt  = [-3,29,-139,533,533,-139,29,-3]/840.0_rp
    real(rp), dimension(recLen), parameter :: interpWgt6 = [ 0, 1,  -8, 37, 37,  -8, 1, 0]/60.0_rp

    !Set choice of LR method in init
    procedure(GetLR_T), pointer :: GetLRs
    
    real(rp) :: pdmb !Set via initModel
    real(rp), dimension(8), parameter :: Cent8C = [-3,29,-139,533,533,-139,29,-3]/840.0_rp
    real(rp), dimension(7), parameter :: Up7C   = [-3,25,-101,319,214,-38,4]/420.0_rp
    real(rp), dimension(5), parameter :: High5C = [2,-13,47,27,-3]/60.0_rp
    real(rp), dimension(6), parameter :: Cent6C = [1,-8,37,37,-8,1]/60.0_rp
    
    contains


    !Given brickette of volumes,conserved variables (n+1/2) and limit conserved variables (n)
    !Return L/Rs of primitive variables
    subroutine BlockStateLRs(Model,VolB,nConB,ConB,Wl,Wr)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in ), dimension(vecLen,recLen)      :: VolB
        real(rp), intent(in ), dimension(vecLen,limLen,NVAR) :: nConB
        real(rp), intent(in ), dimension(vecLen,recLen,NVAR) :: ConB
        real(rp), intent(out), dimension(vecLen,NVAR) :: Wl,Wr

        !Hold primitives
        real(rp), dimension(vecLen,recLen,NVAR) :: PrimB
        !DIR$ attributes align : ALIGN :: PrimB
        real(rp), dimension(vecLen,limLen,NVAR) :: nPrimB
        !DIR$ attributes align : ALIGN :: nPrimB

        !DIR$ ASSUME_ALIGNED nConB: ALIGN
        !DIR$ ASSUME_ALIGNED ConB: ALIGN
        !DIR$ ASSUME_ALIGNED VolB: ALIGN
        !DIR$ ASSUME_ALIGNED Wl: ALIGN
        !DIR$ ASSUME_ALIGNED Wr: ALIGN

        !Convert to primitives
        call BlockCon2Prim(nConB,nPrimB,limLen)
        call BlockCon2Prim( ConB, PrimB,recLen)

        !Reconstruct and limit
        call BlockLRs(VolB,nPrimB,PrimB,Wl,Wr,NVAR)

        contains
            subroutine BlockCon2Prim(inCon,oPrim,xLen)
                integer , intent(in ) :: xLen
                real(rp), intent(in ) :: inCon(vecLen,xLen,NVAR)
                real(rp), intent(out) :: oPrim(vecLen,xLen,NVAR)
                integer :: i,n
                real(rp) :: D,Mx,My,Mz,E,KinE,P

                do n=1,xLen
                    do i=1,vecLen
                        D = max(inCon(i,n,DEN),dFloor)
                        Mx = inCon(i,n,MOMX)
                        My = inCon(i,n,MOMY)
                        Mz = inCon(i,n,MOMZ)
                        E  = inCon(i,n,ENERGY)
                        KinE = 0.5*(Mx**2.0 + My**2.0 + Mz**2.0)/D

                        P = max((Model%gamma-1)*(E-KinE),pFloor)
                        
                        oPrim(i,n,DEN)  = D
                        oPrim(i,n,VELX) = Mx/D
                        oPrim(i,n,VELY) = My/D
                        oPrim(i,n,VELZ) = Mz/D
                        oPrim(i,n,PRESSURE) = P
                    enddo
                enddo

            end subroutine BlockCon2Prim
    end subroutine BlockStateLRs

    !Computes LR's using volume-weighted interpolation and splitting
    subroutine BlockLRs(VolB,nQb,Qb,Ql,Qr,NumV)
        integer , intent(in) :: NumV
        real(rp), intent(in), dimension(vecLen,limLen,NumV) :: nQb
        real(rp), intent(in), dimension(vecLen,recLen,NumV) :: Qb
        real(rp), intent(in), dimension(vecLen,recLen) :: VolB
        real(rp), intent(out),dimension(vecLen,NumV) :: Ql,Qr

        integer :: i,nv
        real(rp), dimension(vecLen) :: Vi !Volume interpolant
        !DIR$ attributes align : ALIGN :: Vi
        !DIR$ ASSUME_ALIGNED Qb: ALIGN
        !DIR$ ASSUME_ALIGNED VolB: ALIGN
        !DIR$ ASSUME_ALIGNED Ql: ALIGN
        !DIR$ ASSUME_ALIGNED Qr: ALIGN
        

        !Calculate volume interpolants
        !Use 8th order central for volume always
        call Central8(VolB,Vi)

        !Get LRs for each variable
        do nv=1,NumV
            !For each variable pass Qb,VolB,Vi
            call GetLRs(VolB,nQb(:,:,nv),Qb(:,:,nv),Vi,Ql(:,nv),Qr(:,nv))
        enddo

    end subroutine BlockLRs

    !Central 8/PDM LRs
    subroutine Cen8LRs(dV,nQ,Q,Vi,Ql,Qr)
        real(rp), intent(in), dimension(vecLen,limLen) :: nQ
        real(rp), intent(in), dimension(vecLen,recLen) :: dV,Q
        real(rp), intent(in), dimension(vecLen)  :: Vi
        real(rp), intent(inout), dimension(vecLen) :: Ql,Qr

        integer :: i,n
        real(rp), dimension(vecLen,recLen) :: QdV !Volume-weighted quantity
        real(rp), dimension(vecLen) :: Qi !Interpolated quantity

        !DIR$ ASSUME_ALIGNED dV: ALIGN
        !DIR$ ASSUME_ALIGNED Q : ALIGN
        !DIR$ ASSUME_ALIGNED Vi: ALIGN
        !DIR$ ASSUME_ALIGNED Ql: ALIGN
        !DIR$ ASSUME_ALIGNED Qr: ALIGN
        !DIR$ ASSUME_ALIGNED nQ : ALIGN

        !Volume-weight
        do n=1,recLen
            do i=1,vecLen
                QdV(i,n) = dV(i,n)*Q(i,n)
            enddo
        enddo

        !Reconstruct and unweight
        call Central8(QdV,Qi)
        do i=1,vecLen
            Qi(i) = Qi(i)/Vi(i)
        enddo

        !Split into LRs
        call pdmLR(nQ,Qi,Ql,Qr)

    end subroutine Cen8LRs


    !Upwind 7/PDM LRs
    subroutine Up7LRs(dV,nQ,Q,Vi,Ql,Qr)
        real(rp), intent(in), dimension(vecLen,limLen) :: nQ
        real(rp), intent(in), dimension(vecLen,recLen) :: dV,Q
        real(rp), intent(in), dimension(vecLen)  :: Vi
        real(rp), intent(inout), dimension(vecLen) :: Ql,Qr

        integer :: i,n
        real(rp), dimension(vecLen,recLen) :: QdV !Volume-weighted quantity

        !DIR$ attributes align : ALIGN :: QdV
        !DIR$ ASSUME_ALIGNED dV: ALIGN
        !DIR$ ASSUME_ALIGNED Q : ALIGN
        !DIR$ ASSUME_ALIGNED Vi: ALIGN
        !DIR$ ASSUME_ALIGNED Ql: ALIGN
        !DIR$ ASSUME_ALIGNED Qr: ALIGN
        !DIR$ ASSUME_ALIGNED nQ : ALIGN

        !Volume-weight
        do n=1,recLen
            do i=1,vecLen
                QdV(i,n) = dV(i,n)*Q(i,n)
            enddo
        enddo

        !High-order LRs, unweight
        do i=1,vecLen
            !High-order interpolation, unweight
            !Do it ugly to avoid temporary array creation
            Ql(i) = Up7(QdV(i,1),QdV(i,2),QdV(i,3),QdV(i,4),QdV(i,5),QdV(i,6),QdV(i,7))/Vi(i)
            Qr(i) = Up7(QdV(i,8),QdV(i,7),QdV(i,6),QdV(i,5),QdV(i,4),QdV(i,3),QdV(i,2))/Vi(i)

            Ql(i) = PDM(nQ(i,1),nQ(i,2),nQ(i,3),Ql(i))
            Qr(i) = PDM(nQ(i,4),nQ(i,3),nQ(i,2),Qr(i))

        enddo

    end subroutine Up7LRs


    !8th order central interpolation
    subroutine Central8(Qb,Qi)
        real(rp), intent(in) :: Qb(vecLen,recLen)
        real(rp), intent(out) :: Qi(vecLen)

        integer :: i,n
        !DIR$ ASSUME_ALIGNED Qb: ALIGN
        !DIR$ ASSUME_ALIGNED Qi : ALIGN

        do i=1,vecLen
            Qi(i) = dot_product(interpWgt,Qb(i,1:recLen))
        enddo

    end subroutine Central8

    !6th order central interpolation
    subroutine Central6(Qb,Qi)
        real(rp), intent(in) :: Qb(vecLen,recLen)
        real(rp), intent(out) :: Qi(vecLen)

        integer :: i,n
        !DIR$ ASSUME_ALIGNED Qb: ALIGN
        !DIR$ ASSUME_ALIGNED Qi : ALIGN

        do i=1,vecLen
            Qi(i) = dot_product(interpWgt6,Qb(i,1:recLen))
        enddo

    end subroutine Central6

    subroutine pdmLR(Qb,Qi,Ql,Qr)
        real(rp), intent(in) :: Qb(vecLen,limLen), Qi(vecLen)
        real(rp), intent(out) :: Ql(vecLen),Qr(vecLen)

        real(rp) :: v0,v1,v2,v3, maxV,minV,vN
        real(rp) :: dv0,dv1,dv2, s0,s1,s2, q0,q1, dvL,dvR
        integer :: i

        do i=1,vecLen
            !Grab limiter values
            v0 = Qb(i,1)
            v1 = Qb(i,2)
            v2 = Qb(i,3)
            v3 = Qb(i,4)

            !Max/Min of nearest neighbors
            maxV = max(v1,v2)
            minV = min(v1,v2)
            vN = max(minV,min(Qi(i),maxV)) 

            !Local differences (signed) and flips
            dv0 = pdmb*(v1-v0)
            dv1 = pdmb*(v2-v1)
            dv2 = pdmb*(v3-v2)
    
            s0 = sign(1.0_rp,dv0)
            s1 = sign(1.0_rp,dv1)
            s2 = sign(1.0_rp,dv2)
    
            q0 = abs(s0+s1)
            q1 = abs(s1+s2)
    
            !Local slopes
            dvL = vN-v1
            dvR = v2-vN

            !Limited L/R values
            Ql(i) = vN - s1*max(0.0,abs(dvL)-q0*abs(dv0))
            Qr(i) = vN + s1*max(0.0,abs(dvR)-q1*abs(dv2))

        enddo
    end subroutine pdmLR

    !PDM Left
    function PDM(q0,q1,q2,qI)
        real(rp), intent(in) :: q0,q1,q2,qI
        real(rp) :: PDM
        real(rp) :: maxQ,minQ, qN, dq0,dq1,dqL
        real(rp) :: s0,s1,s01

        !Max/min of nearest neighbors
        maxQ = max(q1,q2)
        minQ = min(q1,q2)
        qN = max(minQ,min(qI,maxQ))

        !Local differences/flips
        dq0 = pdmb*(q1-q0)
        dq1 = pdmb*(q2-q1)

        s0 = sign(1.0_rp,dq0)
        s1 = sign(1.0_rp,dq1)

        s01 = abs(s0+s1)

        !Local slopes
        dqL = qN-q1

        !Replace w/ limited value
        PDM = qN - s1*max(0.0,abs(dqL)-s01*abs(dq0))
        
    end function PDM

    function Up7(q0,q1,q2,q3,q4,q5,q6) 
        real(rp), intent(in) :: q0,q1,q2,q3,q4,q5,q6
        real(rp) :: Up7

        Up7 =  q0*Up7C(1) + q1*Up7C(2) + q2*Up7C(3) + q3*Up7C(4) &
             + q4*Up7C(5) + q5*Up7C(6) + q6*Up7C(7)

    end function Up7

    !Smoothness detector for length-5 stencil
    function isSmooth5(Q) result(isSmooth)
        real(rp), intent(in) :: Q(-2:+2)
        logical :: isSmooth

        real(rp) :: D1,D2,D3,D4
        logical :: isMin,isMax,isLim

        D1 = Q(-1)-Q(-2)
        D2 = Q( 0)-Q(-1)
        D3 = Q(+1)-Q( 0)
        D4 = Q(+2)-Q(+1)

        isMax = (D1>0) .and. (D2>0) .and. (D3<0) .and. (D4<0)
        isMin = (D1<0) .and. (D2<0) .and. (D3>0) .and. (D4>0)
        isLim = ( abs(D1) > abs(D2) ) .and. ( abs(D3) < abs(D4) )

        if ( (isMax .or. isMin) .and. isLim ) then
            !If max or min *AND* gradient constraint
            isSmooth = .true.
        else
            isSmooth = .false.
        endif
    end function isSmooth5

    !Smoothness detector for length-7 stencil
    function isSmooth7(Q) result(isSmooth)
        real(rp), intent(in) :: Q(-3:+3)
        logical :: isSmooth

        real(rp) :: D1,D2,D3,D4,D5,D6
        logical :: SmL,SmC,SmR

        D1 = Q(-2)-Q(-3)
        D2 = Q(-1)-Q(-2)
        D3 = Q( 0)-Q(-1)
        D4 = Q(+1)-Q( 0)
        D5 = Q(+2)-Q(+1)
        D6 = Q(+3)-Q(+2)

        isSmooth = .false.

        !Check left (D1-D2 / D3-D4)
        SmL = isExtreme(D1,D2,D3,D4)

        !Check center (D2-D3 / D4-D5)
        SmC = isExtreme(D2,D3,D4,D5)

        !Check right (D3,D4,D5,D6)
        SmR = isExtreme(D3,D4,D5,D6)

        if (SmL .or. SmC .or. SmR) then
            isSmooth = .true.
        else
            isSmooth = .false.
        endif
    end function isSmooth7

    !Checks individual block of differences
    function isExtreme(D1,D2,D3,D4)
        real(rp), intent(in) :: D1,D2,D3,D4
        logical :: isExtreme
        logical :: isMin,isMax,isLim

        isMax = (D1>0) .and. (D2>0) .and. (D3<0) .and. (D4<0)
        isMin = (D1<0) .and. (D2<0) .and. (D3>0) .and. (D4>0)
        isLim = ( abs(D1) > abs(D2) ) .and. ( abs(D3) < abs(D4) )
        if ( (isMax .or. isMin) .and. isLim ) then
            isExtreme = .true.
        else
            isExtreme = .false.
        endif

    end function isExtreme
    
end module recon
