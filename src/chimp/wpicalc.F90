!Various routines to handle metrics at particle crossings
module wpicalc
    use chmpdefs
    use tptypes
    use ebtypes
    use tputils
    use wpitypes
    use wpifuns
    use xml_input
    use chmpunits
    
    implicit none

    contains

    !FIXME: include capability of particle interaction with multiple waves
    !Main subroutine for wave particle inteactions
    subroutine PerformWPI(prt,t,dt,Model,ebState)
        type(prt_t), intent(inout) :: prt
        real(rp), intent(in) :: t,dt ! Note: dt is substep size for particle (ddt in pusher.F90)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(wave_T) :: wave
        type(wModel_T) :: wModel
        
        real(rp), dimension(NDIM) :: r,p,E,B,xhat,yhat,bhat
        real(rp), dimension(NDIM) :: p11,pxy
        real(rp), dimension(NVARMHD) :: Qmhd
        real(rp) :: gamma,aNew,pNew,pMag,MagB,Mu,p11Mag,K,rho
        real(rp) :: Ome,wp,astar,xj,yj,Daa
        real(rp) :: dtCum,dtRem,ddt ! substep used in wpi calculations, not same as ddt in pusher.F90 
        integer :: pSgn=1

        real(rp) :: dAlim = 0.05  !Limiting change in pith-angle to be below this value each wpi 
        real(rp) :: da=0.0,dp=0.0 !Change in pitch angle and momentum due to wpi

        if (Model%do2D) then
            !Trap here and quickly grab values
            prt%Qeq(EQX:EQY) = prt%Q(XPOS:YPOS)
            prt%Qeq(EQTIME)  = t
            prt%Qeq(EQKEV)   = prt2kev(Model,prt)
            prt%Qeq(EQALP)   = prt%alpha
            prt%Qeq(EQKEB)   = 0.0
            return
        endif

        ! pulling wave information
        wave = ebState%ebWave
        wModel = ebState%ebWmodel

        !Get local coordinate system
        r = prt%Q(XPOS:ZPOS)
        call ebFields(r,t,Model,ebState,E,B,ijkO=prt%ijk0,vExB=vExB)
        call MagTriad(r,B,xhat,yhat,bhat)
        MagB = max(norm2(B),TINY)

        !Get MHD density
        if (Model%doMHD) then
            Qmhd = mhdInterp(r,t,Model,ebState)
            rho = Qmhd(DEN)
        else 
            write(*,*) 'PerformWPI: Exiting.... no access to full MHD variables, set doMHD to true'
            stop
        endif

        !Check to see if waves are present
        call ChkWave(wModel,r,doWave)
        if (.not. doWave) return

        !Calulating the ratio of nonrelativistic gyrofrequency to plasma frequency at prt's location
        Ome = magB 
        wpe = sqrt(4*PI*rho)*qe_cgs*L0/(vc_cgs*sqrt(Me_cgs)) ! constants are to normalize plasma frequency
        astar = Ome**2/wpe**2

        !calculate minimum energy needed to resonate
        wave%emin = K_whistleR(w0+Dw,B0,n0)

        !check if particle above minimum required energy
        K = prt2kev(prt,Model%m0)/(Model%m0*mec2*1.0e+3)
        if (K .lt. wave%emin) return

        !Setting adaptive time-step/particle sub-stepping to limit da and reduce error
        dtCum = 0.0 !How far we've advanced
        do while (dtCum<dt)
            ! Determine wj,kj that particle resonants with
            call Resonance(Model,wave,wModel,prt,astar,xj,yj)
            if (xj == 999) exit ! no resonant waves present

            Daa = DiffCoef(Model,wave,wModel,prt,astar,MagB,xj,yj) ! Calculate the diffusion coeffcient

            dtRem = dt - dtCum
            if (sqrt(2.0*dt*Daa) > abs(dAlim)) then
                ddt = dAlim**2/(2.*Daa)
                if (dtRem < ddt) ddt = dtRem ! So don't overshoot global step
            else
                ddt = dt
            endif

            dtCum += ddt

            ! Calculate the resulting change in pitch angle and energy of the particle
            call DiffCurve(Model,wave,prt,Daa,ddt,xj,yj,da,dp)

            !Update the pitch angle
            aNew = prt%alpha + da
            prt%alpha = aNew

            if (aNew <= PI/2) then
                pSgn = +1
            else
                pSgn = -1
            endif

            !Updating the momentum components
            gamma = prt2Gam(prt,Model%m0)
            pMag = Model%m0*sqrt(gamma**2.0 - 1.0)+dp !dp is scalar and change in total momentum
            if (prt%isGC) then
                !Scatter GC particle
                p11Mag = pSgn*pMag*sqrt( 1 - sin(aNew)**2.0 )
                Mu = ( (Model%m0**2)*(gamma**2.0 - 1.0) - p11Mag*p11Mag) / (2*Model%m0*MagB)

                prt%Q(P11GC) = p11Mag
                prt%Q(MUGC ) = Mu 
                prt%Q(GAMGC) = gamma
            else                
                !Update momentum for FO particle
                p11 = pSgn*pMag*sqrt( 1 - sin(aNew)**2.0 )
                ! Assume gyrophase remians the same after scattering
                psi = prt%Q(PSITP) 
                pxy = pMag*sin(aNew)*( cos(psi)*xhat + sin(psi)*yhat )

                p = p11*bhat + pxy
                prt%Q(PXFO:PZFO) = p
                
            endif
        enddo

    end subroutine PerformWPI

    !Check to see if waves are present at location 
    subroutine ChkWave(wModel,r,doWave)
        type(wModel_T), intent(in) :: wModel
        real(rp), dimension(NDIM), intent(in) :: r !Particle location
        logical, intent(inout) :: doWave = .false. !Holds if waves are present 

        !FIXME: for now have waves everywhere
        doWave = .true. 

    end subroutine ChkWave

    ! Determines the wave in resosonance with the particle from all possible roots
    subroutine Resonance(Model,wave,wModel,prt,astar,xj,yj)
        type(chmpModel_T), intent(in) :: Model
        type(wave_T), intent(in) :: wave
        type(wModel_T), intent(in) :: wModel
        type(prt_t), intent(in) :: prt
        real(rp), intent(in) :: astar
        real(rp), intent(inout) :: xj,yj
        real(rp), dimension(:), allocatable :: xjs, yjs
        real(rp), dimension(2) :: xcs, ycs
        real(rp) :: pa,mu,K,beta

        pa = prt%alpha
        if (wave%mode .eq. "eRW".and. pa .eq. PI/2.) then
            call res90deg(Model,wave,prt,astar,xjs,yjs)
        else   
            call ResRoots(Model,wave,prt,astar,xjs,yjs)
        end if

        ! taking desired root from all physical roots
        call selectRoot(wModel,prt,xjs,yjs,xj,yj)

        if (xj == 999) return ! no resonant roots

        ! Check to see if wave is a critical root (causes singularity in Daa if not fixed)
        ! can occur for electrons and R-mode waves and protons with L-mode waves
        mu = cos(pa)
        K = prt2kev(prt,Model%m0)/(Model%m0*mec2*1.0e+3) 
        beta = sqrt(K*(K+2.0))/(K+1.0) ! beta = v/c

        if (beta*mu == Vg(wave,astar,xj,yj)) then
            call criticalRoot(Model,wave,prt,astar,xcs,ycs)
            call selectRoot(wModel,prt,xcs,ycs,xj,yj)
        end if

    end subroutine Resonance

    ! Determining if resonant waves are present in the wave model at location
    ! and considering only waves propagating in opposite direction of particle
    subroutine selectRoot(wModel,prt,xjs,yjs,xj,yj)
        type(wModel_T), intent(in) :: wModel
        type(prt_t), intent(in) :: prt
        real(rp), intent(inout) :: xj,yj
        real(rp), dimension(:), allocatable, intent(in) :: xjs, yjs
        real(rp), dimension(:), allocatable :: xpres, ypres
        real(rp) :: pa,xMin,xMax
        
        pa = prt%alpha
        xMin = wModel%xm-wModel%Dx
        xMax = wModel%xm+wModel%Dx
        if (pa >= PI/2.) then
            xpres = pack(xjs, (xjs >= xMin .and. xjs <= xMax .and. yjs < 0))
            ypres = pack(yjs, (xjs >= xMin .and. xjs <= xMax .and. yjs < 0))
        else
            xpres = pack(xjs, (xjs >= xMin .and. xjs <= xMax .and. yjs > 0))
            ypres = pack(yjs, (xjs >= xMin .and. xjs <= xMax .and. yjs > 0))
        end if

        ! Check if no resonant waves are present
        if (size(xpres) == 0) then  
            xj = 999
            yj = 999
        else if (size(xpres) > 1) then
            ! Should only be one resonant wave
            write(*,*) 'Too many resonant roots, w/|Ome|: ', xpres
            stop
        else
            xj = xpres(1)
            yj = ypres(1)
        end if
    end subroutine selectRoot

    !!!!!!!!!!!FIXME: need to add capability to sum over mulitple roots!!!!!!!!!!!!
    !Calculates diffusion coefficient assuming wave spectrum is Gaussian (Summers 2005 eq 33)
    function DiffCoef(Model,wave,wModel,prt,astar,B0,xj,yj) return(Daa)
        type(chmpModel_T), intent(in) :: Model
        type(wave_T), intent(in) :: wave
        type(wModel_T), intent(in) :: wModel
        type(prt_t), intent(in) :: prt
        real(rp), intent(in) :: astar,B0,xj,yj
        real(rp) :: pa,K,beta,R,Ome,DScl,Fxy

        pa= prt%alpha
        K = prt2kev(prt,Model%m0)/(Model%m0*mec2*1.0e+3) 
        beta = sqrt(K*(K+2.0))/(K+1.0) ! beta = v/c
        R = (wModel%B1/B0)**2  !ratio of the wave amplitude to background field strength
        Ome = B0 ! normalized non-relativistic electron gyrofrequency, has Om^2/Ome therefore dont need sign of q

        DScl = (PI/2.0)*abs(Ome)*(K+1)**-2.0

        Fxy = Vg(wave,astar,xj,yj)

        Daa = R*(1.0-xj*cos(pa)/(yj*beta))**2*(abs(Fxy)/abs(beta*cos(pa)-Fxy))*waveSpec(wModel,xj)

        Daa = DScl*Daa

    end function Daa

    !Calculates the change in pitch-angle and corresponding change in momentum along the diffusion curve 
    subroutine DiffCurve(Model,wave,prt,Daa,dt,xj,yj,da,dp) 
        type(chmpModel_T), intent(in) :: Model
        type(wave_T), intent(in) :: wave
        type(prt_t), intent(in) :: prt
        real(rp), intent(in) :: Daa,astar,xj,yj
        real(rp), intent(inout) :: da,dp
        real(rp) :: eta,u,gamu,pa,gamma,pMag,E,A,B

        eta = genRand(-1.0_rp,1.0_rp) !generating random number between -1 and 1 for random walk eq for da
        da = sqrt(2.0*dt*Daa)*eta

        u = xj/yj ! phase velocity of the wave normalized by c
        gamu = sqrt(1-u**2.)**-1.0
        
        pa = prt%alpha
        gamma = prt2Gam(prt,Model%m0)
        pMag = Model%m0*sqrt(gamma**2.0-1.0)
        E = prt2kev(prt,Model%m0)/(Model%m0*mec2*1.0e+3)+1.0 !full energy of the particle including rest mass

        A = gamu**2.0*(pMag*cos(pa)-u*E)*pMag*sin(pa)-pMag**2.0*sin(pa)*cos(pa)
        B = gamu**2.0*(pMag*cos(pa)-u*E)*(cos(pa)-u*pMag/E)+pMag*sin(pa)**2.0

        dp = da*A/B

    end subroutine DiffCurve

    !Calculates the unitless wave number of the resonant root from the resonance criteria
    subroutine resCrit(Model,wave,prt,astar,xj,yj)
        type(chmpModel_T), intent(in) :: Model
        type(wave_T), intent(in) :: wave
        type(prt_t), intent(in) :: prt
        real(rp), intent(in) :: astar,xj
        real(rp), intent(inout) :: yj
        real(rp) :: a,mu,K,beta,K

        a = wave%s*wave%lam/prt2Gam(prt,Model%m0)
        mu = cos(prt%alpha)
        K = prt2kev(prt,Model%m0)/(Model%m0*mec2*1.0e+3) ! normalized energy 
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
        real(rp), allocatable :: xc(:),yc(:)
        real(rp) :: a,b,s,K,beta
        real(rp) :: b0,b1,b2,b3,b4

        b = (1.0+memp)/astar
        s = wave%s
        a = s*wave%lam/prt2Gam(prt,Model%m0)

        K = prt2kev(prt,Model%m0)/(Model%m0*mec2*1.0e+3) ! normalized energy 
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
        yc = xc*np.sqrt(1.0-b/((xc-s)*(xc+s*memp))) !get two waves, one in each direction, field & anti-field alligned 

        xjs = [xc(1),xc(1)]
        yjs = [-yc(1),yc(1)]
        
    end subroutine criticalRoot 

end module wpicalc
