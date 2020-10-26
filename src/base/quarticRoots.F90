module quarticRoots
use kdefs 

implicit none 

    integer, parameter :: NROOTS = 4   
   
contains   
    ! Reference: 
    ! Peter Strobach (2010), Journal of Computational and Applied Mathematics 234
    !    http://www.sciencedirect.com/science/article/pii/S0377042710002128 
    ! Solves for the x1-x4 roots of the quartic equation y(x)=ax^4+bx^3+cx^2+dx+e. 
    subroutine fastQuarticSolver(coef,x)  
        complex(rp), dimension(NROOTS+1), intent(in) :: coef
        complex(rp), dimension(NROOTS), intent(inout) :: x
        complex(rp), dimension(NROOTS) :: depCoef
        complex(rp) :: a,b,c,d ! depressed coefficients of quartic equation
        complex(rp), dimension(2,2) :: quadCoef1,quadCoef2,quadCoef ! in form [alpha, beta][gamma,delta]
        complex(rp), dimension(2) :: rootTest
        complex(rp), dimension(NROOTS) :: e1,e2 !error vector of quadratic coefficients
        complex(rp), dimension(NROOTS) :: xcf !holds roots from closed form solution of quartic roots
        
        real(rp) :: eps,epsTest1,epsTest2,epsilon1,epsilon2
        integer(rp) :: maxIter,iter,i,j,ei
        real(rp), dimension(4) :: epsArr1,epsArr2 ! holds last four values of error in coefficients
        complex(rp) :: alpha,beta,x1,x2
        logical :: isGood

        !setting iteration values
        maxIter = 16
        iter=1
        isGood=.false.

        ! Changing coeff. into depressed form used by Strobach (2010) -- x^4+ax^3+bx^2+cx+d
        a = coef(2)/coef(1); b = coef(3)/coef(1); c = coef(4)/coef(1); d = coef(5)/coef(1)
        depCoef = [a,b,c,d]
        eps = epsilon(1.0_rp) ! machine error used to check if error is close to zero

        ! Check if there are two double roots or one quartic root
        alpha = a/2.0
        beta = (b-alpha**2.0)/2.0
        epsTest1 = c-2.0*alpha*beta
        epsTest2 = d-beta**2.0
      
        if (abs(epsTest1) < eps .and. abs(epsTest2) < eps) then
            x(1:2) = quadraticSolve([alpha,beta])
            x(3:4) = x(1:2)
            return
        endif

        ! Check if there is one cubic root and one simple root (two possible cubic roots to check)
        rootTest = quadraticSolve([a/2.0,b/6.0])
        do i = 1, size(rootTest)
            x1 = rootTest(i)
            x2 = -a-3.0*x1
            epsTest1 = c +x1**2.0*(x1+3.0*x2)
            epsTest2 = d-x2*x1**3.0
            if (abs(epsTest1) < eps .and. abs(epsTest2) < eps) then
                x(1:3) = x1
                x(4) = x2
                return
            endif
        end do

        ! generalized for 4 simple roots or a double root and two simple roots
        ! Initializing arrays holding error in solution, setting to be large
        epsArr1(:) = HUGE
        epsArr2(:) = HUGE

        !will solve for roots by approximating quartic as two quadratics (will solve two solutions simultaneously, only one will converge)
        !determining initial quadratic coefficients 
        xcf = closedFormRoots(depCoef)
        call isort(xcf)  ! put into descending order according to magnitude - abs(xcf)
        
        quadCoef1(1,1) = -real(xcf(1)+xcf(2))  ! alpha01
        quadCoef1(1,2) = real(xcf(1)*xcf(2))   ! beta01
        call fastGammaDelta(depCoef,quadCoef1) ! gamma01, delta01

        quadCoef2(1,1) = -real(xcf(2)+xcf(3))  ! alpha02
        quadCoef2(1,2) = real(xcf(2)*xcf(3))   ! beta02
        call fastGammaDelta(depCoef,quadCoef2) ! gamma02, delta02

        ! Calculating initial error in quadratic coefficients
        e1(1) = a-quadCoef1(1,1)-quadCoef1(2,1)
        e1(2) = b-quadCoef1(1,2)-quadCoef1(1,1)*quadCoef1(2,1)-quadCoef1(2,2)
        e1(3) = c-quadCoef1(1,2)*quadCoef1(2,1)-quadCoef1(1,1)*quadCoef1(2,2)
        e1(4) = d-quadCoef1(1,2)*quadCoef1(2,2)

        e2(1) = a-quadCoef2(1,1)-quadCoef2(2,1)
        e2(2) = b-quadCoef2(1,2)-quadCoef2(1,1)*quadCoef2(2,1)-quadCoef2(2,2)
        e2(3) = c-quadCoef2(1,2)*quadCoef2(2,1)-quadCoef2(1,1)*quadCoef2(2,2)
        e2(4) = d-quadCoef2(1,2)*quadCoef2(2,2)

        !fixed point type iterative refinement or backward optimizer loop for quad coefficients
        do while (iter <= maxIter .and. .not.isGood)
            
            call backwardOptimizer(depCoef,quadCoef1,e1,epsilon1)
            call backwardOptimizer(depCoef,quadCoef2,e2,epsilon2)

            ! Check to see if either branch has converged any(eArr == ee)
            if(epsilon1 < eps .or. any(epsArr1 == epsilon1)) then
                isGood=.true.
                quadCoef = quadCoef1
            end if

            if(epsilon2 < eps .or. any(epsArr2 == epsilon2)) then
                isGood=.true.
                quadCoef = quadCoef2
            end if

            !putting epsilon into error array to compare for next iteration
            ei = modulo(iter,4); if(ei == 0) ei=4
            epsArr1(ei) = epsilon1; epsArr2(ei) = epsilon2
            iter = iter+1
        end do

        ! If hit max iterations than use coefficientsthat had smallest error
        if (iter == maxIter+1) then
            if (epsilon1 < epsilon2) then
                quadCoef = quadCoef1
            else
                quadCoef = quadCoef2
            end if
        end if
      
        x(1:2) = quadraticSolve(quadCoef(1,:))
        x(3:4) = quadraticSolve(quadCoef(2,:))

        !removing error in imaginary part for real roots
        do i = 1, NROOTS
            !if(abs(aimag(x(i))) < TINY*abs(x(i))) x(i) = real(x(i))
            if(abs(aimag(x(i))) < TINY) x(i) = real(x(i))
        end do

    end subroutine fastQuarticSolver 

    ! Ferrari's method from https://en.wikipedia.org/wiki/Quartic_function
    ! note change in notation from Strobach (2010), using similar notaion as wiki to make it easier to compare
    function closedFormRoots(coef) result(roots)
        complex(rp), dimension(NROOTS), intent(in) :: coef
        complex(rp), dimension(NROOTS) :: roots
        complex(rp) :: b,c,d,e ! w/ a=1, coefficients of quartic equation
        complex(rp) :: p,q,r   !depressed quartic coefficients
        complex(rp) :: cubB,cubC,cubD,cubP,cubQ,t,m,S  

        ! in form x^4+bx^3+cx^2+dx+e=0 
        b = coef(1); c = coef(2); d = coef(3); e = coef(4)

        p = (8.0*c-3.0*b**2.0)/8.0
        q = (b**3.0-4.0*b*c+8.0*d)/8.0
        r = (-3.0*b**4/256.0)+e-(b*d/4.0)+(c*b**2.0)/16.0

        ! solving for m (from eq 1a) using Cardano's formula 
        ! rewritten so cubic coeff: a=1, b=p, c=(p^2-4r)/4, d=-q^2/8
        cubB = p
        cubC = (p**2.0-4.0*r)/4.0
        cubD = -q**2.0/8.0
        ! Depressed cubic coefficients
        cubP = cubC-cubB**2.0/3.0 
        cubQ = (2.0*cubB**3.0-9.0*p*cubC+27.0*cubD)/27.0
        t = (-cubQ/2.0+sqrt(cubQ**2.0/4.0+cubP**3.0/27.0))**(1.0/3.0) &
            + (-cubQ/2.0-sqrt(cubQ**2.0/4.0+cubP**3.0/27.0))**(1.0/3.0)
        m = t-p/3.0 ! m is largest source of error in the close form solution of roots

        S = sqrt(2.0*m)/2.0
        roots(1) = -b/4.0-S+0.5*sqrt(-4*S**2-2*p+q/S)
        roots(2) = -b/4.0-S-0.5*sqrt(-4*S**2-2*p+q/S)
        roots(3) = -b/4.0+S+0.5*sqrt(-4*S**2-2*p-q/S)
        roots(4) = -b/4.0+S-0.5*sqrt(-4*S**2-2*p-q/S)

    end function closedFormRoots

    ! Table 3 of Strobach(2010) to calculate gamma and delta from alpha0, beta0
    subroutine fastGammaDelta(coef,quadCoefs)
        complex(rp), dimension(4), intent(in) :: coef
        complex(rp), dimension(2,2), intent(inout) :: quadCoefs
        complex(rp) :: a,b,c,d
        complex(rp) :: alpha0,beta0,gamma0,delta0
        complex(rp) :: phi1,phi2,c1,c2,L1,L2,L3,y1,y2

        a = coef(1); b = coef(2); c = coef(3); d = coef(4)
        alpha0 = quadCoefs(1,1)
        beta0  = quadCoefs(1,2)

        phi1=1+alpha0**2+beta0**2
        phi2=alpha0*(1+beta0)
        c1=a-alpha0+alpha0*(b-beta0)+beta0*c
        c2=b-beta0+alpha0*c+beta0*d
        L1=sqrt(phi1)
        L3=phi2/L1
        L2=sqrt(phi1-phi2**2/phi1)
        y1=c1/L1
        y2=(c2-y1*L3)/L2

        delta0=y2/L2
        gamma0=(y1-delta0*L3)/L1
        quadCoefs(2,1) = gamma0; quadCoefs(2,2) = delta0

    end subroutine fastGammaDelta

    subroutine backwardOptimizer(coef,quadCoefs,e,epsilon)
        complex(rp), dimension(NROOTS), intent(in) :: coef
        complex(rp), dimension(2,2), intent(inout) :: quadCoefs
        complex(rp), dimension(NROOTS), intent(inout) :: e ! Coefficient Error
        real(rp), intent(inout) :: epsilon

        complex(rp) :: a,b,c,d ! coefficients of quartic equation
        complex(rp) :: alpha,beta,gamma,delta ! coefficients of 2 quadratics approximation
        complex(rp) :: U23,U33,L43,U44
        complex(rp) :: x1,x2,x3,x4,y1,y2,y3,y4

        a = coef(1); b = coef(2); c = coef(3); d = coef(4)
        
        alpha = quadCoefs(1,1); beta  = quadCoefs(1,2)
        gamma = quadCoefs(2,1); delta = quadCoefs(2,2)

        U23=alpha-gamma
        U33=beta-delta-gamma*U23
        L43=-delta*U23/U33
        U44=beta-delta-L43*U23

        x1=e(1)
        x2=e(2)-gamma*x1
        x3=e(3)-delta*x1-gamma*x2
        x4=e(4)-delta*x2-L43*x3
        y4=x4/U44
        y3=(x3-U23*y4)/U33
        y2=x2-U23*y3-y4
        y1=x1-y3

        alpha=alpha+y1
        beta=beta+y2
        gamma=gamma+y3
        delta=delta+y4

        !udating coefficients
        quadCoefs(1,1) = alpha; quadCoefs(1,2) = beta
        quadCoefs(2,1) = gamma; quadCoefs(2,2) = delta

        !calculating error in coefficients
        e(1)=a-alpha-gamma
        e(2)=b-beta-alpha*gamma-delta
        e(3)=c-beta*gamma-alpha*delta
        e(4)=d-beta*delta   
        epsilon=sum(abs(e))

    end subroutine backwardOptimizer

    ! Chapter 5.6 from Numerical Recipes
    !Assume in form x^2+b*x+c=0
    function quadraticSolve(coef) result(x)
        complex(rp), dimension(2), intent(in) :: coef
        complex(rp), dimension(2) :: x
        complex(rp) :: b,c
        complex(rp) :: q,sqt
        real(rp)  :: a,sgn=1.0

        a = 1.0; b = coef(1); c = coef(2)

        sqt = sqrt(b**2.0-4.0*a*c)

        if(real(conjg(b)*sqt) < 0.0) sgn = -1.0 

        q = -0.5*(b+sgn*sqt) 

        x = [q/a,c/q]

    end function quadraticSolve

    ! insertion sort in descending order (only use for array size < 20)
    subroutine isort(x) 
        complex(rp), dimension(:), intent(inout) :: x
        complex(rp) :: temp
        integer :: i,j

        do j = 1, size(x) 
            temp = x(j)
            i=j
            do while (i>1 .and. abs(x(i-1)) < abs(temp))
                x(i)=x(i-1)
                i = i-1
            enddo
            x(i)=temp
        end do
    end subroutine isort
   
end module quarticRoots 