!Various low-level math routines
module quadrature
    use types
    use math
    
    implicit none
    !Constants for Gaussian integration
    integer, parameter :: gNp = 12
    

    !Zeros of 12th order Legendre polynomial
    real(rp), dimension(gNp), parameter :: gA   = &
              [ 0.1252334085114689154724_dp, 0.3678314989981801937527_dp, &
                0.5873179542866174472967_dp, 0.7699026741943046870369_dp, &
                0.9041172563704748566785_dp, 0.9815606342467192506906_dp, &
               -0.1252334085114689154724_dp,-0.3678314989981801937527_dp, &
               -0.5873179542866174472967_dp,-0.7699026741943046870369_dp, &
               -0.9041172563704748566785_dp,-0.9815606342467192506906_dp ]

    !Gaussian integration coefficients
    real(rp), dimension(gNp), parameter :: gWgt = &
              [ 0.2491470458134027850006_dp, 0.233492536538354808761_dp , &
                0.203167426723065921749_dp , 0.160078328543346226335_dp , &
                0.1069393259953184309603_dp, 0.0471753363865118271946_dp, &
                0.2491470458134027850006_dp, 0.233492536538354808761_dp , &
                0.203167426723065921749_dp , 0.160078328543346226335_dp , &
                0.1069393259953184309603_dp, 0.0471753363865118271946_dp]

    !Alternative Gaussian order for speed
    ! integer, parameter :: gNp = 6
    ! real(rp), dimensin(gNp), parameter ::    gA = [0.6612093864662645, 0.2386191860831969, 0.9324695142031521, &
    !                                                  -0.6612093864662645,-0.2386191860831969,-0.9324695142031521]
    ! real(rp), dimension(gNp), parameter :: gWgt = [0.3607615730481386, 0.4679139345726910, 0.1713244923791704, &
    !                                                   0.3607615730481386, 0.4679139345726910, 0.1713244923791704]                              

    contains

    !Gaussian Volume integral of function W: (x,y,z) -> (D,Vx,Vy,Vz,P) (GasIC_T)
    !Returns vector of integrated primitive quantities
    !xyzC are cell corners 8xNDIM
    !Assuming specific corner order
    !(-1,-1,-1) -> (+1,-1,-1) -> (+1,+1,-1) -> (-1,+1,-1) -> (-1,-1,+1) -> (+1,-1,+1) -> (+1,+1,+1) -> (-1,+1,+1)
    function GaussianVolumeIntegral(xyzC,W) result(wInt)
        real(rp), dimension(8,NDIM), intent(in) :: xyzC
        procedure(GasIC_T), pointer, intent(in) :: W
        real(rp), dimension(NVAR) :: wInt

        real(rp), dimension(NVAR) :: Wp
        real(rp), dimension(NDIM) :: xp
        real(rp), dimension(NDIM,NDIM) :: Jacob
        integer, dimension(8) :: iS,jS,kS
        real(rp), dimension(8) :: iScl,jScl,kScl
        integer :: i,j,k
        real(rp) :: dJ,wijk
        
        
        wInt = 0.0

        iS = [-1,+1,+1,-1,-1,+1,+1,-1]
        jS = [-1,-1,+1,+1,-1,-1,+1,+1]
        kS = [-1,-1,-1,-1,+1,+1,+1,+1]

        do i=1,gNp
            do j=1,gNp
                do k=1,gNp
                    wijk = gWgt(i)*gWgt(j)*gWgt(k)!Weight of this point
                    
                    iScl = [1-gA(i),1+gA(i),1+gA(i),1-gA(i),1-gA(i),1+gA(i),1+gA(i),1-gA(i)]
                    jScl = [1-gA(j),1-gA(j),1+gA(j),1+gA(j),1-gA(j),1-gA(j),1+gA(j),1+gA(j)]
                    kScl = [1-gA(k),1-gA(k),1-gA(k),1-gA(k),1+gA(k),1+gA(k),1+gA(k),1+gA(k)]

                    !Calculate Gaussian evaluation point
                    xp(XDIR) = 0.125*dot_product(xyzC(:,XDIR),iScl*jScl*kScl)
                    xp(YDIR) = 0.125*dot_product(xyzC(:,YDIR),iScl*jScl*kScl)
                    xp(ZDIR) = 0.125*dot_product(xyzC(:,ZDIR),iScl*jScl*kScl)

                    !Calculate Jacobian at this point
                    Jacob(XDIR,1) = 0.125*dot_product(xyzC(:,XDIR),iS*jScl*kScl)
                    Jacob(YDIR,1) = 0.125*dot_product(xyzC(:,YDIR),iS*jScl*kScl)
                    Jacob(ZDIR,1) = 0.125*dot_product(xyzC(:,ZDIR),iS*jScl*kScl)

                    Jacob(XDIR,2) = 0.125*dot_product(xyzC(:,XDIR),iScl*jS*kScl)
                    Jacob(YDIR,2) = 0.125*dot_product(xyzC(:,YDIR),iScl*jS*kScl)
                    Jacob(ZDIR,2) = 0.125*dot_product(xyzC(:,ZDIR),iScl*jS*kScl)

                    Jacob(XDIR,3) = 0.125*dot_product(xyzC(:,XDIR),iScl*jScl*kS)
                    Jacob(YDIR,3) = 0.125*dot_product(xyzC(:,YDIR),iScl*jScl*kS)
                    Jacob(ZDIR,3) = 0.125*dot_product(xyzC(:,ZDIR),iScl*jScl*kS)

                    !Get determinant of Jacobian
                    dJ = abs(DetJ(Jacob))

                    !Evaluate for variables at this point
                    call W(xp(XDIR),xp(YDIR),xp(ZDIR),Wp(1),Wp(2),Wp(3),Wp(4),Wp(5))
                    wInt = wInt + dJ*wijk*Wp

                enddo
            enddo
        enddo

    end function GaussianVolumeIntegral

    !TODO: Clean up these routines to remove redundant code

    !Calculates face flux of given vector field
    subroutine GaussianFaceFlux(f0,f1,f2,f3,Axyz,flx)
        real(rp), dimension(NDIM), intent(in) :: f0,f1,f2,f3
        procedure(VectorField_T), pointer, intent(in) :: Axyz
        real(rp), intent(out) :: flx

        real(rp), dimension(NDIM) :: dx,dy,dz,etaV,psiV, epV, N
        real(rp) :: area, eta, psi, w, xp,yp,zp, Axp,Ayp,Azp
        integer :: i,j

        flx = 0.0

        !Calculate face displacements from f0
        call faceDeltas(f0,f1,f2,f3,dx,dy,dz)

        !Loop over 12x12 grid of points on face
        do i=1,gNp
            do j=1,gNp
                !Displacement of the i,j point
                eta = 0.5*(1+gA(i))
                psi = 0.5*(1+gA(j))
                w = gWgt(i)*gWgt(j)

                !Calculate evaluation point and evaluate
                xp = f0(XDIR) + dx(1)*eta + dx(2)*psi + dx(3)*eta*psi
                yp = f0(YDIR) + dy(1)*eta + dy(2)*psi + dy(3)*eta*psi
                zp = f0(ZDIR) + dz(1)*eta + dz(2)*psi + dz(3)*eta*psi

                call Axyz(xp,yp,zp,Axp,Ayp,Azp)

                !Calculate displacement vectors and cross for area
                etaV(1) = dx(1) + dx(3)*psi
                etaV(2) = dy(1) + dy(3)*psi
                etaV(3) = dz(1) + dz(3)*psi

                psiV(1) = dx(2) + dx(3)*eta
                psiV(2) = dy(2) + dy(3)*eta
                psiV(3) = dz(2) + dz(3)*eta

                area = 0.25*norm2(cross(etaV,psiV))*w
                N = area*cross(etaV,psiV)/norm2(cross(etaV,psiV))
                flx = flx + dot_product(N,[Axp,Ayp,Azp])
            enddo
        enddo

    end subroutine GaussianFaceFlux

    !Gaussian face integral of Axyx on face defined by f0,f1,f2,f3
    !Note: Order of corner points important, f0=SW,f1=SE,f2=NW,f3=NE
    !Returns F,F2,Fij
    !F = face integral of A
    !F2 = face integral of Ai*Ai
    !Fij = face integral of cross: Axy,Ayz,Azx
    subroutine GaussianFaceIntegral(f0,f1,f2,f3,Axyz,fInt,fInt2,fIntX)
        real(rp), dimension(NDIM), intent(in) :: f0,f1,f2,f3
        procedure(VectorField_T), pointer, intent(in) :: Axyz
        real(rp), dimension(NDIM), intent(out) :: fInt,fInt2,fIntX
        real(rp), dimension(NDIM) :: dx,dy,dz,etaV,psiV, epV
        real(rp) :: area, eta, psi, w, xp,yp,zp, Axp,Ayp,Azp
        integer :: i,j

        fInt = 0.0
        fInt2 = 0.0
        fIntX = 0.0

        !Calculate face displacements from f0
        call faceDeltas(f0,f1,f2,f3,dx,dy,dz)

        !Loop over 12x12 grid of points on face
        do i=1,gNp
            do j=1,gNp
                !Displacement of the i,j point
                eta = 0.5*(1+gA(i))
                psi = 0.5*(1+gA(j))
                w = gWgt(i)*gWgt(j)

                !Calculate evaluation point and evaluate
                xp = f0(XDIR) + dx(1)*eta + dx(2)*psi + dx(3)*eta*psi
                yp = f0(YDIR) + dy(1)*eta + dy(2)*psi + dy(3)*eta*psi
                zp = f0(ZDIR) + dz(1)*eta + dz(2)*psi + dz(3)*eta*psi

                call Axyz(xp,yp,zp,Axp,Ayp,Azp)

                !Calculate displacement vectors and cross for area
                etaV(1) = dx(1) + dx(3)*psi
                etaV(2) = dy(1) + dy(3)*psi
                etaV(3) = dz(1) + dz(3)*psi

                psiV(1) = dx(2) + dx(3)*eta
                psiV(2) = dy(2) + dy(3)*eta
                psiV(3) = dz(2) + dz(3)*eta

                area = 0.25*norm2(cross(etaV,psiV))*w
                !Accumulate integrals into total

                !Calculate face integrals of vector field
                fInt(XDIR) = fInt(XDIR) + area*Axp
                fInt(YDIR) = fInt(YDIR) + area*Ayp
                fInt(ZDIR) = fInt(ZDIR) + area*Azp

                !Calculate face integral of squares
                fInt2(XDIR) = fInt2(XDIR) + area*Axp*Axp
                fInt2(YDIR) = fInt2(YDIR) + area*Ayp*Ayp
                fInt2(ZDIR) = fInt2(ZDIR) + area*Azp*Azp

                !Calculate face integral of cross terms
                !Fij(1) = xy, Fij(2) = yz, Fij(3) = zx
                fIntX(1) = fIntX(1) + area*Axp*Ayp
                fIntX(2) = fIntX(2) + area*Ayp*Azp
                fIntX(3) = fIntX(3) + area*Azp*Axp
                
            enddo
        enddo

    end subroutine GaussianFaceIntegral

    subroutine GaussianFaceStress(f0,f1,f2,f3,Axyz,Mbb)
        real(rp), dimension(NDIM), intent(in) :: f0,f1,f2,f3
        procedure(VectorField_T), pointer, intent(in) :: Axyz
        real(rp), dimension(NDIM), intent(out) :: Mbb

        real(rp), dimension(NDIM) :: dx,dy,dz,etaV,psiV, epV
        real(rp), dimension(NDIM) :: dS,B0
        real(rp) :: area, eta, psi, w, xp,yp,zp, Axp,Ayp,Azp
        integer :: i,j

        Mbb = 0.0

        !Calculate face displacements from f0
        call faceDeltas(f0,f1,f2,f3,dx,dy,dz)

        !Loop over 12x12 grid of points on face
        do i=1,gNp
            do j=1,gNp
                !Displacement of the i,j point
                eta = 0.5*(1+gA(i))
                psi = 0.5*(1+gA(j))
                w = gWgt(i)*gWgt(j)

                !Calculate evaluation point and evaluate
                xp = f0(XDIR) + dx(1)*eta + dx(2)*psi + dx(3)*eta*psi
                yp = f0(YDIR) + dy(1)*eta + dy(2)*psi + dy(3)*eta*psi
                zp = f0(ZDIR) + dz(1)*eta + dz(2)*psi + dz(3)*eta*psi

                call Axyz(xp,yp,zp,Axp,Ayp,Azp)

                !Calculate displacement vectors and cross for area
                etaV(1) = dx(1) + dx(3)*psi
                etaV(2) = dy(1) + dy(3)*psi
                etaV(3) = dz(1) + dz(3)*psi

                psiV(1) = dx(2) + dx(3)*eta
                psiV(2) = dy(2) + dy(3)*eta
                psiV(3) = dz(2) + dz(3)*eta

                area = 0.25*norm2(cross(etaV,psiV))*w

                !Get local dS
                !dS  = area*cross(etaV,psiV)/norm2(cross(etaV,psiV))
                dS = 0.25*w*cross(etaV,psiV)

                !Accumulate integral into total
                ! 0.5(B0^2)*I - B0*B0
                B0 = [Axp,Ayp,Azp]

                
                Mbb = Mbb + 0.5*dot_product(B0,B0)*dS &
                          - B0*dot_product(B0,dS)
                
            enddo
        enddo

    end subroutine GaussianFaceStress
    !Gaussian face system (N,T1,T2) on face defined by f0,f1,f2,f3
    !Note: Order of corner points important, f0=SW,f1=SE,f2=NW,f3=NE
    !Returns N,T1,T2

    subroutine GaussianFaceSystem(f0,f1,f2,f3,N,T1,T2)
        real(rp), dimension(NDIM), intent(in) :: f0,f1,f2,f3
        real(rp), dimension(NDIM), intent(inout) :: N,T1,T2

        real(rp), dimension(NDIM) :: dx,dy,dz,etaV,psiV, epV
        real(rp) :: area, eta, psi, w, xp,yp,zp, Axp,Ayp,Azp
        integer :: i,j


        !Initialize
        N = 0
        T1 = 0
        T2 = 0

        !Calculate face displacements from f0
        call faceDeltas(f0,f1,f2,f3,dx,dy,dz)

        !Loop over 12x12 grid of points on face
        do i=1,gNp
            do j=1,gNp
                !Displacement of the i,j point
                eta = 0.5*(1+gA(i))
                psi = 0.5*(1+gA(j))
                w = gWgt(i)*gWgt(j)

                !Calculate evaluation point and evaluate
                xp = f0(XDIR) + dx(1)*eta + dx(2)*psi + dx(3)*eta*psi
                yp = f0(YDIR) + dy(1)*eta + dy(2)*psi + dy(3)*eta*psi
                zp = f0(ZDIR) + dz(1)*eta + dz(2)*psi + dz(3)*eta*psi


                !Calculate displacement vectors and cross for area
                etaV(1) = dx(1) + dx(3)*psi
                etaV(2) = dy(1) + dy(3)*psi
                etaV(3) = dz(1) + dz(3)*psi

                psiV(1) = dx(2) + dx(3)*eta
                psiV(2) = dy(2) + dy(3)*eta
                psiV(3) = dz(2) + dz(3)*eta

                area = 0.25*norm2(cross(etaV,psiV))*w
                !Accumulate integrals into total
                N  = N  + area*cross(etaV,psiV)/norm2(cross(etaV,psiV))
                T1 = T1 + area*etaV/norm2(etaV)
                T2 = T2 + area*psiV/norm2(psiV)
            enddo
        enddo

        !Use accumulated T1/T2 to calculate normal
        T1 = T1/norm2(T1)
        T2 = T2/norm2(T2)
        N = cross(T1,T2)/norm2(cross(T1,T2))

        !Recalculate T1 to enforce orthogonality of triad
        T1 = cross(T2,N)/norm2(cross(T2,N))
        
    end subroutine GaussianFaceSystem

    !Gaussian edge integral of Axyz between e1/e2
    function GaussianEdgeIntegral(e1,e2,Axyz) result(eInt)
        real(rp), dimension(NDIM), intent(in) :: e1,e2
        procedure(VectorField_T), pointer, intent(in) :: Axyz
        real(rp), dimension(NDIM) :: eInt

        integer :: i
        real(rp), dimension(NDIM) :: dR, dR2, xBar, xp,Ap

        dR = e2-e1 !Displacement vector
        dR2 = 0.5*(e2-e1)
        xBar = 0.5*(e2+e1)

        eInt = 0.0
        do i=1,gNp
            !Evaluation point
            xp = xBar + gA(i)*dR2

            !Evaluate at each point and add to total integral
            call Axyz(xp(XDIR),xp(YDIR),xp(ZDIR),Ap(XDIR),Ap(YDIR),Ap(ZDIR))
            eInt = eInt+0.5*gWgt(i)*Ap

        enddo
        

    end function GaussianEdgeIntegral

    !Face displacements, assuming standard face corner ordering
    subroutine faceDeltas(f0,f1,f2,f3,dx,dy,dz)
        real(rp), dimension(NDIM), intent(in) :: f0,f1,f2,f3
        real(rp), dimension(NDIM), intent(out) :: dx,dy,dz
        !Calculate face displacements from f0
        dx(IDIR) = f1(XDIR) - f0(XDIR)
        dx(JDIR) = f2(XDIR) - f0(XDIR)
        dx(KDIR) = f3(XDIR) + f0(XDIR) - f2(XDIR) - f1(XDIR)

        dy(IDIR) = f1(YDIR) - f0(YDIR)
        dy(JDIR) = f2(YDIR) - f0(YDIR)
        dy(KDIR) = f3(YDIR) + f0(YDIR) - f2(YDIR) - f1(YDIR)

        dz(IDIR) = f1(ZDIR) - f0(ZDIR)
        dz(JDIR) = f2(ZDIR) - f0(ZDIR)
        dz(KDIR) = f3(ZDIR) + f0(ZDIR) - f2(ZDIR) - f1(ZDIR)

    end subroutine faceDeltas
end module quadrature
