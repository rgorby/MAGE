!contains wave mode and wave model specific functions
module wpifuns
    use chmpdefs
    use kdefs
    use tputils
    use wpitypes
    use xml_input
    use strings

    implicit none

    real(rp) :: memp = Me_cgs/Mp_cgs

    !dispersion related functions 
    procedure(VgFun_T) , pointer :: Vg => NULL()
    procedure(ResFun_T), pointer :: ResRoots => NULL()

    !function provides wave frequency spectrum 
    procedure(wsFun_T) , pointer :: waveSpec => NULL()

       integer, parameter :: NROOTS = 4 !Number of possible resonant roots

    !Group Velocity function type
    !Returns Vg for a given resonant wave in normalized units
    abstract interface
        function VgFun_T(wave,astar,x,y) result(Vg)
            import :: rp,wave_T
            type(wave_T), intent(in) :: wave
            real(rp), intent(in) :: astar,x,y
            real(rp) :: Vg
        end function VgFun_T
    end interface

    !Resonant root function type
    !Returns resonant root in unitless variables (divided by gyrofrequency)
    abstract interface
        subroutine ResFun_T(Model,wave,prt,astar,x,y) 
            import :: rp,prt_T,chmpModel_T,wave_T
            type(chmpModel_T), intent(in) :: Model
            type(prt_T), intent(in) :: prt
            type(wave_T), intent(in) :: wave
            real(rp), intent(in) :: astar
            real(rp), dimension(:), allocatable, intent(out) :: x,y
        end subroutine ResFun_T
    end interface

    !wave frequency function type
    !Returns wave spectral density for a given wave number
    abstract interface
        function wsFun_T(wModel,x) result(Ws)
            import :: rp,wModel_T
            type(wModel_T), intent(in) :: wModel
            real(rp), intent(in) :: x
            real(rp) :: Ws
        end function wsFun_T
    end interface

    contains

    !Initialize the wave
    subroutine initWPI(Model,wModel,wave,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(wModel_T), intent(inout) :: wModel
        type(wave_T), intent(inout) :: wave
        type(XML_Input_T), intent(inout) :: inpXML

        !initializing the wave model
        call inpXML%Set_Val(wModel%model,'wpi/wmodel',"Gauss")

        select case (trim(toUpper(wModel%model)))
        !-------
        case("GAUSSIAN","GAUSS")
            call inpXML%Set_Val(wModel%xm,'wpi/xm',0.25)
            call inpXML%Set_Val(wModel%Dx,'wpi/dx',0.2)
            call inpXML%Set_Val(wModel%B1,'wpi/b1',0.01)

            !normalizing
            wModel%B1 = wModel%B1*inBScl

            waveSpec => gaussWS 
        case default
            write(*,*) '<Unknown wave mode. Exiting....>'
            stop
        end select

        !initializing the wave
        call inpXML%Set_Val(wave%mode,'wpi/mode',"eWhistler")

        select case (trim(toUpper(wave%mode)))
        !-------
        case("EWHISTLER","E_WHISTLER")
            wave%mode = "eRW"
            wave%s = 1.0
            wave%lam = -1.0
            wave%emin = 0.0
            Vg => epVg
            ResRoots => epResRoots
        case default
            write(*,*) '<Unknown wave mode. Exiting....>'
            stop
        end select

    end subroutine initWPI

    !!!!!!!!!!!!!!!!! Wave Model Related Functions !!!!!!!!!!!!!!!!!!!!
    !assuming gaussian wave frequency spectrum
    function gaussWS(wModel,x) result(Ws)
        type(wModel_T), intent(in) :: wModel
        real(rp), intent(in) :: x
        real(rp) :: xm,dx,sigma,Ws

        !Setting up wave spectrum (assumed to be Gaussian) Ws^~(w) 
        !Eq 31/32 from Summers 2005 without normalization constant from B1 (included in Daa)
        sigma = 2.0   ! 95% of wave spectrum will be xm-dx < x < xm+dx (value used in Summers 2005 and Horne 2003)
        dx = wModel%Dx/sigma
        xm = wModel%xm 
        Ws = EXP(-(x-xm)**2/dx**2)/(SQRT(PI)*ERF(sigma)*dx)

    end function gaussWS

    !!!!!!!!!!!!!!!!! wave related functions !!!!!!!!!!!!!!!!!!!!
    !variable notation similar to Summers et al. 2005

    !Calculates the unitless wave number of the resonant root from the resonance criteria
    subroutine resCrit(Model,wave,prt,astar,xj,yj)
        type(chmpModel_T), intent(in) :: Model
        type(wave_T), intent(in) :: wave
        type(prt_t), intent(in) :: prt
        real(rp), intent(in) :: astar
        real(rp), dimension(:), allocatable, intent(inout) :: xj, yj
        real(rp) :: a,mu,K,beta

        a = wave%s*wave%lam/prt2Gam(prt,Model%m0)
        mu = cos(prt%alpha)
        K = prt2kev(Model,prt)/(Model%m0*mec2*1.0e+3) ! normalized energy 
        beta = sqrt(K*(K+2.0))/(K+1.0) !beta = v/c
        yj = (xj+a)/(beta*mu) !resonance criteria [see Eq 24 of Summers 2005 for notation]

    end subroutine resCrit

    !Solving the generalized resonance condition for particles with 90 deg pitch angles (A2/3 of Summers 2005)
    subroutine res90deg(Model,wave,prt,astar,xjs,yjs)
        type(chmpModel_T), intent(in) :: Model
        type(wave_T), intent(in) :: wave
        type(prt_t), intent(in) :: prt
        real(rp), intent(in) :: astar
        real(rp), dimension(2), intent(out) :: xjs,yjs
        real(rp) :: a,b,s,y0

        b = (1.0+memp)/astar
        s = wave%s
        a = s*wave%lam/prt2Gam(prt,Model%m0)

        y0 = abs(a)*sqrt(1.+b/((a+s)*(s*memp-a)))
        
        xjs = [-a,-a]
        yjs = [y0,-y0]

    end subroutine res90deg

    ! Calculating the value of the critical root to remove the singularity if necessary (Appendix B of Summers 2005)
    subroutine criticalRoot(Model,wave,prt,astar,xjs,yjs)
        type(chmpModel_T), intent(in) :: Model
        type(wave_T), intent(in) :: wave
        type(prt_t), intent(in) :: prt
        real(rp), intent(in) :: astar
        real(rp), dimension(2), intent(out) :: xjs, yjs
        complex(rp), dimension(NROOTS) :: roots
        real(rp), dimension(NROOTS+1) :: coef
        real(rp), allocatable :: xc(:),yc(:)
        real(rp) :: a,b,s,K,beta
        real(rp) :: b0,b1,b2,b3,b4

        b = (1.0+memp)/astar
        s = wave%s
        a = s*wave%lam/prt2Gam(prt,Model%m0)

        K = prt2kev(Model,prt)/(Model%m0*mec2*1.0e+3) ! normalized energy 
        beta = sqrt(K*(K+2.0))/(K+1.0) !v/c

        b0 = 1.0
        b1 = 2.0*s*(-1.0+memp)+b/a
        b2 = 1.0-4.0*memp+memp**2.+b*s*(-1.0+memp)/(2.*a)
        b3 = -s*(-1.0+memp)*(b+4.*memp)/2.0
        b4 = memp*(b+memp)

        coef = [b0,b1,b2,b3,b4]
        !roots = np.roots(coef)

        ! Keeping roots that are positive, below the gyrofrequency (xj<1), and real (others are non-physical)
        xc = pack(roots, (real(roots)>0 .and. real(roots)<1 .and. aimag(roots) == 0)) 
        !Should only be one root
        if (size(xc) > 1) then
            ! Should only be one resonant wave
            write(*,*) 'wpiCalc:criticalRoot:: Too many resonant roots, w/|Ome|: ', xc
            stop
        end if
        yc = xc*sqrt(1.0-b/((xc-s)*(xc+s*memp))) !get two waves, one in each direction, field & anti-field alligned 

        xjs = [xc(1),xc(1)]
        yjs = [-yc(1),yc(1)]
        
    end subroutine criticalRoot 
    
    !Generalized dw/dk for plasma of electrons and protons (Appendix C of Summers 2005)
    function epVg(wave,astar,x,y) result(Fxy)
        type(wave_T), intent(in) :: wave
        real(rp), intent(in) :: astar,x,y
        real(rp) :: b,s,gx,Vg,Fxy
        real(rp) :: c1,c2,c3,c4

        b = (1.0+memp)/astar
        s = wave%s

        c1 = 2.0*s*(-1.0+memp)
        c2 = 1.0-4.0*memp+memp**2.0
        c3 = -s*(memp-1.0)*(b+4.0*memp)/2.0
        c4 = memp*(b+memp)

        gx = x**4.+c1*x**3.+c2*x**2.+c3*x+c4
        Fxy = y*((x-s)**2.0)*((x+s*memp)**2.0)/(x*gx)

    end function epVg

    !Solving resonance for generalized plasma of protons and electrons (Appendix A of Summers 2005)
    subroutine epResRoots(Model,wave,prt,astar,xjs,yjs)
        type(chmpModel_T), intent(in) :: Model
        type(prt_T), intent(in) :: prt
        type(wave_T), intent(in) :: wave
        real(rp), intent(in) :: astar
        complex(rp), dimension(NROOTS) :: roots
        real(rp), dimension(:), allocatable, intent(out) :: xjs, yjs
        real(rp) :: a,b,s,mu,beta,denom,alpha,gamma,K
        real(rp) :: a0,a1,a2,a3,a4 !polunomial coefficients

        alpha = prt%alpha
        gamma = prt2Gam(prt,Model%m0)
        K = prt2kev(Model,prt)/(Model%m0*mec2*1.0e+3) !in code units

        b = (1.0+memp)/astar
        s = wave%s
        a = s*wave%lam/gamma

        mu = cos(alpha)
        beta = sqrt(K*(K+2.0))/(K+1.0)

        denom = 1.0-((beta*mu)**2.0)

        a0 = 1.0
        a1 = (2.0*a+s*(-1.0+memp)-((beta*mu)**2.0)*s*(-1+memp))/denom
        a2 = (a**2.0+2.0*a*s*(-1.0+memp)-memp+(b+memp)*(beta*mu)**2.0)/denom
        a3 = ((a**2.0)*s*(-1.0+memp)-2.0*a*memp)/denom
        a4 = (-memp*a**2)/denom

        !!!!FIXME: root solver and picking out correct roots!!!
        !roots = np.roots(coef)

        ! Keeping roots that are positive, below the gyrofrequency (xj<1), and real (others are non-physical)
        xjs = pack(roots, (real(roots)>0 .and. real(roots)<1 .and. aimag(roots) == 0))
        call resCrit(Model,wave,prt,astar,xjs,yjs)

    end subroutine epResRoots

    !Returns kinetic energy of particle needed to resonate with given wave (takes in dimensionless w and k)
    !Assumes all energy is in the 11 direction
    function Kres_whistleR(x,astar) result(Kres)
        real(rp), intent(in) :: x,astar
        real(rp) :: y,gamma,Kres

        !calculating minimum wavenumber (ymin) of given wave from dispersion relation
        y = sqrt(x/(astar*(1.-x))) ! Take wave moving in same direction as B
        gamma = (x - y*sqrt(y**2.-x**2.+1.))/(x**2.-y**2.) 
        Kres = (gamma-1.) !returns min kinetic energy (K) in chimp units

    end function Kres_whistleR

end module wpifuns
