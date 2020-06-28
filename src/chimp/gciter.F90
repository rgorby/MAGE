!FO->GC iteration
module gciter
    use chmpdefs
    use chmpunits
    use tptypes
    use ebtypes
    use ebinterp
    use math

    implicit none

    real(rp), private :: Ak0 = 1.0/64 !Minimum step in line search
    real(rp), private :: nScl = 1.0 !Error scaling for new method
    real(rp), private :: rScl = 1.0e-1 !Random perturbation

    !General GC Iteration routine
    abstract interface
        function GCIter_T(xFO,pFO,t,Model,ebState,isConvO,NitO) result(Rgc)
            Import :: rp,NDIM,chmpModel_T,ebState_T
            real(rp), dimension(NDIM), intent(in) :: xFO,pFO
            real(rp), intent(in) :: t
            type(chmpModel_T), intent(in) :: Model
            type(ebState_T), intent(in)   :: ebState
            logical, intent(out), optional :: isConvO
            integer, intent(out), optional :: NitO
            real(rp) :: Rgc(NDIM)
        end function GCIter_T
    end interface

    !procedure(GCIter_T), pointer :: CalcGC => ChimpGC
    !procedure(GCIter_T), pointer :: CalcGC => TestGC
    procedure(GCIter_T), pointer :: CalcGC


    contains

    !Direct GC calculation (no iteration)
    function DirectGC(xFO,pFO,t,Model,ebState,isConvO,NitO) result(Rgc)
        real(rp), dimension(NDIM), intent(in) :: xFO,pFO
        real(rp), intent(in) :: t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out), optional :: isConvO
        integer, intent(out), optional :: NitO
        real(rp) :: Rgc(NDIM)

        integer, dimension(NDIM) :: ijk
        real(rp), dimension(NDIM) :: vFO,E,B,vExB
        real(rp), dimension(NDIM) :: bhat,vP,r0,rE,r1
        real(rp), dimension(NDIM) :: beta, zeta, gMagB,gOm
        real(rp), dimension(NDIM,NDIM) :: JacB,Jbh,M,vv,zz
        real(rp) :: v11,MagB,gamma,Om,bCof
        logical :: isIn
        type(gcFields_T) :: gcFields


        call locate(xFO,ijk,Model,ebState%ebGr,isIn)
        if (.not. isIn) then
            Rgc = xFO
            return
        endif

        !Get field values/derivatives at position
        call ebFields(xFO,t,Model,ebState,E,B,ijk,vExB,gcFields)
        MagB = norm2(B)
        bhat = normVec(B)

        !Get particle properties
        gamma = sqrt( 1 + dot_product(pFO/Model%m0,pFO/Model%m0) )
        Om = Model%q0*MagB/(gamma*Model%m0) !Relativistic gyrofreq
        vFO = pFO/(gamma*Model%m0)
        v11 = dot_product(vFO,bhat)
        vP = vFO - v11*bhat - vExB

        !Rgc = xFO + r0 + rE + r1

        r0 = cross(vP,bhat)/Om
        rE = cross(vExB,bhat)/Om

        !Now do next order correction
        JacB = gcFields%JacB
        gMagB = VdT(bhat,JacB) !gradient of |B|
        Jbh = ( MagB*JacB - Dyad(gMagB,B) )/(MagB**2.0)

        gOm = Model%q0*gMagB/(gamma*Model%m0) !Gradient of gyrofrequency

        M = ( Om*Jbh - Dyad(bhat,gOm) )/(Om**2.0)

        beta = v11*bhat + vP/4.0
        zeta = cross(vP,bhat)

        r1 = ( dCross(beta,zeta,M) + dCross(zeta,beta,M) )/Om
        r1 = r1 + v11*VdT(vP,Jbh)/(Om**2.0)

        bCof = v11*dot_product(vP,VdT(bhat,Jbh)) &
             + 0.125*(tColon(Dyad(vP,vP)-Dyad(zeta,zeta),Jbh))
        r1 = r1 + bhat*bCof/(Om**2.0)

        Rgc = xFO + r0 + rE + r1
        
        if (present(isConvO)) isConvO=.true.
        if (present(NitO)) NitO = 1

    end function DirectGC

    !Fast direct GC calculation (no iteration, only first term)
    function FastGC(xFO,pFO,t,Model,ebState,isConvO,NitO) result(Rgc)
        real(rp), dimension(NDIM), intent(in) :: xFO,pFO
        real(rp), intent(in) :: t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out), optional :: isConvO
        integer, intent(out), optional :: NitO
        real(rp) :: Rgc(NDIM)

        integer, dimension(NDIM) :: ijk
        real(rp), dimension(NDIM) :: vFO,E,B,vExB
        real(rp), dimension(NDIM) :: bhat,vP,r0,rE,r1
        
        
        real(rp) :: v11,MagB,gamma,Om
        logical :: isIn

        call locate(xFO,ijk,Model,ebState%ebGr,isIn)
        if (.not. isIn) then
            Rgc = xFO
            return
        endif

        !Get field values/derivatives at position
        call ebFields(xFO,t,Model,ebState,E,B,ijk,vExB)
        MagB = norm2(B)
        bhat = normVec(B)

        !Get particle properties
        gamma = sqrt( 1 + dot_product(pFO/Model%m0,pFO/Model%m0) )
        Om = Model%q0*MagB/(gamma*Model%m0) !Relativistic gyrofreq
        vFO = pFO/(gamma*Model%m0)
        v11 = dot_product(vFO,bhat)
        vP = vFO - v11*bhat - vExB

        !Rgc = xFO + r0 + rE + 0*r1

        r0 = cross(vP,bhat)/Om
        rE = cross(vExB,bhat)/Om

        Rgc = xFO + r0 + rE
        
        if (present(isConvO)) isConvO=.true.
        if (present(NitO)) NitO = 1

    end function FastGC

    !Hopefully final GC iteration function
    function ChimpGC(xFO,pFO,t,Model,ebState,isConvO,NitO) result(Rgc)
        real(rp), dimension(NDIM), intent(in) :: xFO,pFO
        real(rp), intent(in) :: t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out), optional :: isConvO
        integer, intent(out), optional :: NitO
        real(rp) :: Rgc(NDIM)

        logical :: isIn,isConv,doIter
        integer :: n,ijk(NDIM)
        real(rp) :: dR0,Err,newErr
        real(rp), dimension(NDIM) :: newRgc,newRg,Rg
        real(rp), dimension(NDIM) :: E,B,vExB
        Rgc = 0.0
    !Initialize
        !Get local scale using total p/B
        call locate(xFO,ijk,Model,ebState%ebGr,isIn)
        if (.not. isIn) return
        call ebFields(xFO,t,Model,ebState,E,B,vExB=vExB,ijkO=ijk)

        dR0 = norm2(pFO)/(Model%q0*norm2(B))
        !Setup iteration
        n = 0
        Rgc = xFO
        Rg = GyroVec(E,B,vExB)

        Err = HUGE !Guarantee 1 iteration
        doIter = .true.
        isConv = .false.

    !Iteration loop    
        do while (doIter)
            !Calculate Rgc^{n+1} via Rg^n
            newRgc = xFO-Rg

            !Get Rg^n+1 for next iteration/testing error
            !Check point
            isIn = inDomain(newRgc,Model,ebState%ebGr)
            if (.not. isIn) then
                doIter = .false.
                cycle
            endif

            call ebFields(newRgc,t,Model,ebState,E,B,vExB=vExB,ijkO=ijk)
            newRg = GyroVec(E,B,vExB)

            !Calculate error in new Rgc
            newErr = norm2(newRgc-xFO+newRg)/dR0

            !Test whether new iterate is better
            if (newErr < Err) then
                !Take value/setup next loop
                Rgc = newRgc
                Rg  = newRg
                Err = newErr
                n = n+1
            else
                !This was worse, keep old value and bail
                doIter = .false.
                isConv = .false.
                cycle
            endif

            !Test exit conditions
            if (Err <= Model%TolGC) then
                doIter = .false.
                isConv = .true.
            else if (n >= Model%MaxIter) then
                doIter = .false.
                isConv = .false.
            endif

        enddo

        if (present(isConvO)) isConvO=isConv
        if (present(NitO)) NitO = n

        contains

        function GyroVec(E,B,vExB) result(Rg)
            real(rp), dimension(NDIM), intent(in) :: E,B,vExB
            real(rp), dimension(NDIM) :: Rg
            real(rp) :: MagB, p11
            real(rp), dimension(NDIM) :: bhat,pPerp
            !Jump to ExB frame and get perp momentum
            call SplitP(Model,pFO,E,B,vExB,p11,pPerp)
            bhat = normVec(B)
            MagB = max(norm2(B),TINY)
            Rg = cross(bhat,pPerp)/(Model%q0*MagB)

        end function GyroVec

    end function ChimpGC


    !Calculates guiding center given current position
    function TestGC(r,p,t,Model,ebState,isConvO,NitO) result(Rgc)
        real(rp), dimension(NDIM), intent(in) :: r,p
        real(rp), intent(in) :: t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out), optional :: isConvO
        integer, intent(out), optional :: NitO

        real(rp) :: Rgc(NDIM)


        !real(rp), dimension(NDIM) :: fpRgc,bRgc,ndRgc,kRgc
        real(rp), dimension(NDIM) :: dRgc,iRgc
        logical :: dCon,iCon
        integer :: dN,iN
        !real(rp), dimension(NDIM) :: B,pP,bhat
        !?integer :: Nfp,Nb,Nnd,Nk
        !ogical :: fpCon,bCon,ndCon,kCon
        real(rp) :: rMag

        rMag = norm2(Rho(r,p,t,Model,ebState))

        ! B = fldInterp(r,t,Model,ebState,BFLD)
        ! bhat = normVec(B)
        ! pP = p - dot_product(p,bhat)*bhat
        ! K  = (mec2*1.0e+3)*p2K(p ,Model%m0)
        ! Kp = (mec2*1.0e+3)*p2K(pP,Model%m0)
        ! !$OMP CRITICAL
        
        ! fpRgc = fixpGC(r,p,t,Model,ebState,fpCon,Nfp)
        ! ndRgc = ndGC(r,p,t,Model,ebState,ndCon,Nnd)
        ! kRgc = ChimpGC(r,p,t,Model,ebState,kCon,Nk)
        dRgc = DirectGC(r,p,t,Model,ebState,dCon,dN)
        iRgc =  ChimpGC(r,p,t,Model,ebState,iCon,iN)

        !$OMP critical
        write(*,*) '------------------'
        write(*,'(a,4f14.8)') 'r / |r| = ', r,rMag
        write(*,*) 'Direct : '
        write(*,'(a,3f14.8)')   '   r     = ', dRgc
        write(*,'(a,I0,es9.2)') '   N/err = ', dN,ObjF(dRgc,r,p,t,Model,ebState)/rMag
        write(*,'(a,f14.8)')    '   dR    = ', norm2(dRgc-r)
        write(*,*) 'isConv = ', dCon

        write(*,*) 'Chimp : '
        write(*,'(a,3f14.8)')   '   r     = ', iRgc
        write(*,'(a,I0,es9.2)') '   N/err = ', iN,ObjF(iRgc,r,p,t,Model,ebState)/rMag
        write(*,'(a,f14.8)')    '   dR    = ', norm2(iRgc-r)
        write(*,*) 'isConv = ', iCon

        write(*,*) '------------------'
        !$OMP end critical

        Rgc = iRgc
        
    end function TestGC

    !Objective function
    !Evaluate error using RGC as guiding center for xFO,pFO
    function ObjF(Rgc,xFO,pFO,t,Model,ebState)
        real(rp), dimension(NDIM), intent(in) :: Rgc,xFO,pFO
        real(rp), intent(in) :: t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState

        real(rp) :: ObjF
        real(rp) :: rhovec(NDIM) !Rho vector for GC
        rhovec = Rho(Rgc,pFO,t,Model,ebState)
        ObjF = norm2(Rgc-xFO+rhovec) !Objective function value

    end function ObjF

    !GC vector, xFO = Rgc + Rho
    function Rho(Rgc,p,t,Model,ebState)
        real(rp), dimension(NDIM), intent(in) :: Rgc,p
        real(rp), intent(in) :: t
        real(rp) :: Rho(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState

        real(rp), dimension(NDIM) :: E,B,vExB,bhat,pPerp
        real(rp) :: p11,MagB

        call ebFields(Rgc,t,Model,ebState,E,B,vExB=vExB)
        call SplitP(Model,p,E,B,vExB,p11,pPerp)
        bhat = normVec(B)
        MagB = max(norm2(B),TINY)
        rho = cross(bhat,pPerp)/(Model%q0*MagB)
    end function Rho

    ! !Calculates guiding center given current position
    ! function bfgsGC(xFO,pFO,t,Model,ebState,isConvO,NitO) result(Rgc)
    !     real(rp), dimension(NDIM), intent(in) :: xFO,pFO
    !     real(rp), intent(in) :: t
    !     type(chmpModel_T), intent(in) :: Model
    !     type(ebState_T), intent(in)   :: ebState
    !     logical, intent(out), optional :: isConvO
    !     integer, intent(out), optional :: NitO
    !     real(rp) :: Rgc(NDIM)

    !     logical :: doIter,isIn
    !     real(rp) :: dR,Ak,err
    !     real(rp), dimension(NDIM) :: xk,gObj,Pk,Sk,Yk,newG
    !     real(rp), dimension(NDIM,NDIM) :: Bk, Bki
    !     integer :: n,i

    !     !Calculate dR for differencing
    !     dR = norm2(Rho(xFO,pFO,t,Model,ebState))

    !     !Set initial values
    !     n = 0
    !     xk = xFO
    !     gObj = GradO(xk,dR,xFO,pFO,t) !Initial value of gradient
    !     Bk = 0.0 !Set Hessian to identity tensor initially
    !     do i=1,NDIM
    !         Bk(i,i) = 1.0*norm2(gObj)/dR
    !     enddo
    !     !err = ObjF(xk,xFO,pFO,t)/dR
    !     !write(*,*) 'n0 err = ', err
        
    !     doIter = .true.

    !     do while (doIter)
    !         !Invert Bk
    !         call matinv3(Bk,Bki)
    !         !Get direction
    !         Pk = -matmul(Bki,gObj)

    !         !Do line search for alpha
    !         !Ak = 0.05
    !         Ak = LineSearch(xk,Pk,xFO,pFO,t)

    !         Sk = Ak*Pk
    !         xk = xk + Sk !New position
    !         !write(*,*) 'Bk = ',Bk
    !         !write(*,*) 'Bki = ',Bki
    !         !write(*,*) 'Pk/Sk = ',Pk,Sk
    !         newG = GradO(xk,dR,xFO,pFO,t) !New gardient
    !         Yk = newG - gObj
    !         gObj = newG
    !         !write(*,*) 'Yk = ',Yk
    !         call StepHessian(Bk,Yk,Sk)

    !         !err = ObjF(xk,xFO,pFO,t)/dR
    !         err = nScl*norm2(Sk)/dR
    !         isIn = inDomain(xk,Model,ebState%ebGr)

    !         if (err < Model%TolGC .or. .not. isIn) then
    !             doIter = .false.
    !         else
    !             n = n+1
    !         endif

    !     enddo

    !     if (present(NitO)) NitO = n
    !     if (present(isConvO)) isConvO = .true.
    !     !write(*,'(a,I0,es9.2)') '   N/Err = ', n,err
    !     Rgc = xk

    !     contains

    !     !Local functions
    !     function LineSearch(xk,Pk,x0,p0,t) result(Ak)
    !         real(rp), dimension(NDIM), intent(in) :: xk,Pk,x0,p0
    !         real(rp), intent(in) :: t
    !         real(rp) :: Ak

    !         real(rp) :: f0,fAk,aScl
    !         logical :: doIter
    !         Ak = 1.0
    !         doIter = .true.
    !         do while(doIter)
    !             !Evaluate optimality condition
    !             f0 = ObjF(xk       ,x0,p0,t,Model,ebState)
    !             fAk = ObjF(xk+Ak*Pk,x0,p0,t,Model,ebState)

    !             aScl = 1.0-0.5*Ak
    !             if (fAk < aScl*f0) then
    !                 !Done
    !                 doIter = .false.
    !             else
    !                 Ak = Ak/2.0
    !             endif
    !             if (Ak<=Ak0) doIter=.false.
    !         enddo
    !     end function LineSearch
    !     !Update Hessian
    !     subroutine StepHessian(Bk,Yk,Sk)
    !         real(rp), dimension(NDIM,NDIM), intent(inout) :: Bk
    !         real(rp), dimension(NDIM), intent(in) :: Yk,Sk

    !         real(rp), dimension(NDIM,NDIM) :: oBk
    !         integer :: i,j
    !         real(rp) :: up1,up2,scl1,scl2 !Denominators
    !         !Save original Bk
    !         oBk = Bk
    !         do i=1,NDIM
    !             do j=1,NDIM
    !                 Bk(i,j) = oBk(i,j)
    !                 scl1 = Yk(j)*Sk(i)
    !                 scl2 = Sk(j)*oBk(i,j)*Sk(i)
    !                 if ( abs(scl1) > TINY .and. abs(scl2) > TINY) then
    !                     up1 = Yk(i)*Yk(j)
    !                     up2 = oBk(i,j)*Sk(i)*Sk(j)*oBk(i,j)
    !                     Bk(i,j) = oBk(i,j) + up1/scl1 - up2/scl2
    !                 endif
    !             enddo
    !         enddo
    !         !Bk = oBk
    !     end subroutine StepHessian

    !     !Gradient of objective function @ r
    !     !dR is distance for finite difference
    !     function GradO(r,rmag,x0,p0,t) result(gObj)
    !         real(rp), dimension(NDIM), intent(in) :: r,x0,p0
    !         real(rp), intent(in) :: rmag,t
    !         real(rp) :: gObj(NDIM)

    !         integer :: i
    !         real(rp) :: dR(NDIM)
    !         real(rp) :: fP,fM
    !         gObj = 0.0
    !         do i=1,NDIM
    !             dR = 0.0
    !             dR(i) = rmag
    !             fP = ObjF(r+dR,x0,p0,t,Model,ebState)
    !             fM = ObjF(r-dR,x0,p0,t,Model,ebState)
    !             gObj(i) = (fP-fM)/(2*rmag)
    !         enddo
            
    !     end function GradO
    ! end function bfgsGC

    ! !Calculates guiding center given current position
    ! function fixpGC(r,p,t,Model,ebState,isConvO,NitO) result(Rgc)
    !     real(rp), dimension(NDIM), intent(in) :: r,p
    !     real(rp), intent(in) :: t
    !     type(chmpModel_T), intent(in) :: Model
    !     type(ebState_T), intent(in)   :: ebState
    !     logical, intent(out), optional :: isConvO
    !     integer, intent(out), optional :: NitO
    !     real(rp) :: Rgc(NDIM)

    !     logical :: isIn,isConv
    !     real(rp), dimension(NDIM) :: E,B,vExB,bhat,obhat,pPerp,dR
    !     real(rp) :: bErr,rErr,MagB,p11,oErr
    !     real(rp) :: oRgc(2,NDIM) !Keep track of 2 previous iterations
    !     integer :: n

    !     !Initialize iteration
    !     Rgc = r !First guess
    !     oRgc(1,:) = 0.0
    !     oRgc(2,:) = Rgc
    !     n = 0
    !     obhat = 0.0
    !     bErr = Model%TolGC+1 !Guarantee 1 iteration
    !     oErr = Model%TolGC+1 !Guarantee 1 iteration
    !     ! write(*,*) '-----------------------'
    !     ! write(*,'(a,3f14.8)') '   r = ', r
    !     ! write(*,'(a,3f14.8)') '   p = ', p

    !     do while (bErr>Model%TolGC .and. n<Model%MaxIter)
    !     !do while (oErr>Model%TolGC .and. n<Model%MaxIter)
    !         !Check point
    !         isIn = inDomain(r,Model,ebState%ebGr)
    !         if (.not. isIn) then
    !             isConv = .false.
    !             exit
    !         endif
    !         !Get local fields and split momentum
    !         call FOFields(Rgc,t,Model,ebState,E,B,vExB)
    !         call SplitP(Model,p,E,B,vExB,p11,pPerp)
    !         bhat = normVec(B)
    !         MagB = max(norm2(B),TINY)
    !         dR = cross(bhat,pPerp)/(Model%q0*MagB)

    !         !Update GC
    !         Rgc = r - dR

    !         !Check error and setup next iteration
    !         bErr = norm2(bhat-obhat)
    !         rErr = norm2(Rgc-oRgc(1,:)) !Check for repeating
    !         obhat = bhat
    !         n = n+1
    !         if (rErr <= TINY) then
    !             !End here and take average of 2 Rgc's
    !             Rgc = 0.5*(Rgc+oRgc(2,:))
    !             n = Model%MaxIter+1 !End convergence and signal failure
    !             bErr = HUGE
    !         endif
    !         ! oErr = ObjF(Rgc,r,p,t,Model,ebState)/norm2(dR)
    !         ! write(*,*) '---------------'
    !         ! write(*,*) 'n / berr / err = ',n,bErr,oErr
    !         ! write(*,*) '---------------'

    !         ! if (n>200) then
    !             ! write(*,'(a,I0,2E10.3)') 'n / bErr / rErr = ',n,bErr,rErr
    !             ! write(*,'(a,E10.3)') '    dRgc = ', norm2(Rgc-oRgc(2,:))
    !             ! write(*,'(a,3f14.8)') '   Rgc = ', Rgc
    !             ! write(*,*) ''
    !             ! write(*,'(a,3f14.8)') '  bhat = ', bhat
    !             ! write(*,*) 'normbhat = ', norm2(bhat)
    !             ! write(*,'(a,3f14.8)') ' |pP| = ', normVec(pPerp)
    !             ! write(*,'(a,3f14.8)') '  B = ', B
    !             ! write(*,'(a,3f14.8)') '    dR = ', dR
    !         ! endif
    !         oRgc(1,:) = oRgc(2,:)
    !         oRgc(2,:) = Rgc
    !     enddo

    !     if (bErr < Model%TolGC) then
    !         isConv = .true.
    !     else
    !         isConv = .false.
    !         ! !$OMP CRITICAL
    !         ! write(*,*) '!!! GC Convergence Failure !!!'

    !         ! !$OMP END CRITICAL
    !     endif
    !     !write(*,'(a,I0,es9.2)') '   N/Err = ', n,bErr
    !     !write(*,*) 'Conv/N = ',isConv,n
    !     !write(*,*) '-----------------------'
    !     if (present(isConvO)) isConvO=isConv
    !     if (present(NitO)) NitO = n
    ! end function fixpGC

    ! !Calculates guiding center given current position
    ! !Uses Nelder Mead iteration
    ! function ndGC(xFO,pFO,t,Model,ebState,isConvO,NitO) result(Rgc)
    !     real(rp), dimension(NDIM), intent(in) :: xFO,pFO
    !     real(rp), intent(in) :: t
    !     type(chmpModel_T), intent(in) :: Model
    !     type(ebState_T), intent(in)   :: ebState
    !     logical, intent(out), optional :: isConvO
    !     integer, intent(out), optional :: NitO
    !     real(rp) :: Rgc(NDIM)

    !     integer :: i
    !     real(rp) :: A,Gam,P,Sig
    !     real(rp) :: rmag,err
    !     logical :: isConv,doIter
    !     type(Simplex_T) :: Smp
    !     !Variables for reflection/expansion/contraction points/values
    !     real(rp), dimension(NDIM) :: X0,X1,Xi,Xr,Xe,Xc
    !     real(rp) :: Fr,Fe,Fc,F1,F2,F3,F4

    !     !Hard-wired ND coefficients
    !     A = 1.0
    !     Gam = 2.0
    !     P = 0.5
    !     Sig = 0.5

    !     doIter = .true.
    !     !Create initial simplex (sets rmag), not sorted
    !     call InitSmp(Smp)

    !     !err = Smp%F(1)/rmag
    !     !do while ( err>Model%TolGC .and. Smp%n<Model%MaxIter)
    !     do while (doIter)
    !     !Initialize/test convergence
    !         !Sort simplex
    !         call SortSmp(Smp)
    !         !Save values for simplicity
    !         F1 = Smp%F(1)
    !         F2 = Smp%F(2)
    !         F3 = Smp%F(3)
    !         F4 = Smp%F(4)
    !         err = F1/rmag
    !         !call DebugSmp(Smp)
    !         Smp%n = Smp%n + 1

    !         !Test for exit conditions
    !         if (err <= Model%TolGC) then
    !             isConv = .true.
    !             doIter = .false.
    !             cycle
    !         endif
    !         if (Smp%n > Model%MaxIter) then
    !             isConv = .false.
    !             doIter = .false.
    !             cycle
    !         endif

    !     !Centroid
    !         !Sum points on triangle opposite worst point
    !         X0 = sum(Smp%x(:,1:NDIM),dim=2)/NDIM
    !     !Reflection
    !         Xr = X0 + A*(X0 - Smp%x(:,NDIM+1))
    !         Fr = ObjF(Xr,xFO,pFO,t,Model,ebState)
            
    !         !Test reflection point
    !         if ( (F1 <= Fr) .and. (Fr < F3) ) then
    !             !Replace worst (NDIM+1)
    !             Smp%F(4)   = Fr
    !             Smp%x(:,4) = Xr
    !             cycle
    !         endif
    !     !Expansion
    !         !Test reflection point
    !         if (Fr < F1) then
    !             !This is the best so far, calculate expanded point
    !             Xe = X0 + Gam*(Xr-X0)
    !             Fe = ObjF(Xe,xFO,pFO,t,Model,ebState)

    !             !Test expanded vs reflected
    !             if (Fe < Fr) then
    !                 !Replace worst with expanded
    !                 Smp%F(4)   = Fe
    !                 Smp%x(:,4) = Xe
    !             else 
    !                 !Replace worst with reflected
    !                 Smp%F(4)   = Fr
    !                 Smp%x(:,4) = Xr

    !             endif
    !             !Cycle either way
    !             cycle
    !         endif !Fr<F1
    !     !Contraction
    !         !If still here then Fr>=F3
    !         !Calculate contraction point using better of Xr and X4
    !         if (Fr < F4) then
    !             !Outside contraction
    !             Xc = X0 + P*(Xr - X0)
    !         else !Inside contraction
    !             Xc = X0 + P*(Smp%x(:,NDIM+1) - X0)
    !         endif
    !         Fc = ObjF(Xc,xFO,pFO,t,Model,ebState)
    !         !Test contraction point
    !         if (Fc < F4) then
    !             !Better than worst, replace with contracted
    !             Smp%F(4)   = Fc
    !             Smp%x(:,4) = Xc

    !             cycle
    !         endif
    !     !Shrink
    !         !Everything is terrible, replace all but best point with shrink
    !         !Replace point and re-evaluate function
    !         X1 = Smp%x(:,1)
    !         !write(*,*) 'Shrink!'
    !         do i=1,NDIM
    !             Xi = Smp%x(:,i+1) !Old point
    !             Xi = X1 + Sig*(Xi-X1)
    !             Smp%x(:,i+1) = Xi
    !             Smp%F(  i+1) = ObjF(Xi,xFO,pFO,t,Model,ebState)
    !         enddo
    !     enddo

    !     !Finalize
    !     if (isConv) then
    !         Rgc = Smp%x(:,1) !Take best value
    !     else
    !         Rgc = 0.0
    !         Rgc = Smp%x(:,1) !Take best value
    !     endif
        
    !     if (present(isConvO)) isConvO=.true.
    !     if (present(NitO)) NitO = Smp%n

    !     contains
    !     !Local functions
    !     subroutine DebugSmp(Smp)
    !         type(Simplex_T), intent(in) :: Smp

    !         write(*,'(a,I0,es9.2)') 'n/Err = ', Smp%n,Smp%F(1)/rmag
    !         write(*,'(a,3f14.8)') '   X1 = ', Smp%x(:,1)
    !         write(*,'(a,3f14.8)') '   X2 = ', Smp%x(:,2)
    !         write(*,'(a,3f14.8)') '   X3 = ', Smp%x(:,3)
    !         write(*,'(a,3f14.8)') '   X4 = ', Smp%x(:,4)

    !         write(*,'(a,4es9.2)') '  E14 = ', Smp%F(:)/rmag

    !     end subroutine DebugSmp
    !     subroutine InitSmp(Smp)
    !         type(Simplex_T), intent(inout) :: Smp

    !         real(rp),dimension(NDIM) :: RhoV,B,bhat
    !         integer :: i

    !         B = fldInterp(xFO,t,Model,ebState,BFLD)
    !         bhat = normVec(B)
    !         RhoV = Rho(xFO,pFO,t,Model,ebState)
    !         rmag = norm2(RhoV)

    !         ! Smp%x(XDIR:ZDIR,1) = xFO !Initial point
    !         ! do i=1,NDIM
    !         !     Smp%x(XDIR:ZDIR,i+1) = xFO
    !         !     !Add to d component
    !         !     Smp%x(i,i+1) = Smp%x(i,i+1) + rmag
    !         ! enddo

    !         Smp%x(:,1) = xFO + rmag*bhat
    !         Smp%x(:,2) = xFO - rmag*bhat
    !         Smp%x(:,3) = xFO + RhoV
    !         Smp%x(:,4) = xFO - RhoV
    !         !Evaluate objective function at points
    !         do i=1,NDIM+1
    !             Smp%F(i) = ObjF(Smp%x(:,i),xFO,pFO,t,Model,ebState)
    !         enddo
    !         Smp%n = 0
    !     end subroutine InitSmp

    !     !Reorder points in simplex to ensure strict monotonicity
    !     subroutine SortSmp(Smp)
    !         type(Simplex_T), intent(inout) :: Smp

    !         integer :: i,i0
    !         real(rp) :: F0,x0(NDIM)

    !         do i=1,NDIM
    !             !Find best point between i:NDIM+1
    !             i0 = minloc(Smp%F(i:NDIM+1),dim=1) + (i-1)
                
    !             !Swap positions i/i0
    !             F0 = Smp%F(i)
    !             x0 = Smp%x(:,i)
    !             Smp%F(  i) = Smp%F(i0)
    !             Smp%x(:,i) = Smp%x(:,i0)
    !             Smp%F(  i0) = F0
    !             Smp%x(:,i0) = x0

    !         enddo
    !     end subroutine SortSmp
    ! end function ndGC
end module gciter
