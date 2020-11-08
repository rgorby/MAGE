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
        
        real(rp), dimension(NDIM) :: r,p,E,B,xhat,yhat,bhat,req
        real(rp), dimension(NDIM) :: p11,pxy,vExB
        real(rp), dimension(NVARMHD) :: Qmhd
        real(rp) :: gamma,gamNew,aNew,pMag,MagB,Mu,p11Mag,K,rho,psi
        real(rp) :: Ome,wpe,astar,xHigh,xj,yj,Daa,inWScl
        real(rp) :: dtCum,dtRem,ddt ! substep used in wpi calculations, not same as ddt in pusher.F90 
        real(rp) :: zEq = 0.0 !Defined Z for equator
        integer :: pSgn

        real(rp) :: dAlim  !Limiting change in pith-angle to be below this value each wpi to reduce error in diff. Curve
        real(rp) :: da,dp !Change in pitch angle and momentum due to wpi
        real(rp) :: K0,K1,dK ! energy before and after wpi

        logical :: doWave

        if (.not.prt%isIn) return !shouldn't be here if particle is not in domain

        if (Model%do2D) then
            !Trap here and quickly grab values
            prt%Qeq(EQX:EQY) = prt%Q(XPOS:YPOS)
            prt%Qeq(EQTIME)  = t
            prt%Qeq(EQKEV)   = prt2kev(Model,prt)
            prt%Qeq(EQALP)   = prt%alpha
            prt%Qeq(EQKEB)   = 0.0
            return
        endif

        !initializing some variables
        pSgn = 1; dAlim = 0.05; da=0.0; dp=0.0; dK=0.0

        ! pulling wave information
        wave = ebState%ebWave
        wModel = ebState%ebWmodel

        !Get local coordinate system
        r = prt%Q(XPOS:ZPOS)
        req = [prt%Qeq(EQX),prt%Qeq(EQX),zEq] ! defining magnetic equator at z=0
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
        call ChkWave(wModel,r,req,doWave)
        if (.not. doWave) return ! no waves so exit

        !Calulating the ratio of nonrelativistic gyrofrequency to plasma frequency at prt's location
        !normalization constant to put wpe in code units
        inWScl = 114704.0207604 !(qe_cgs*L0/(vc_cgs*sqrt(Me_cgs)))^2
        Ome = magB 
        wpe = sqrt(4*PI*rho) 
        astar = Ome**2/((wpe**2)*inWScl)

        !calculate minimum energy needed to resonate
        xHigh = wModel%xm+wModel%Dx
        wave%emin = Kres_whistleR(xHigh,astar)

        !check if particle above minimum required energy
        K = prt2kev(Model,prt)/(Model%m0*mec2*1.0e+3)
        if (K .lt. wave%emin) return

        !Setting adaptive time-step/particle sub-stepping to limit da and reduce error
        dtCum = 0.0 !How far we've advanced
        do while (dtCum<dt)
            ! Determine wj,kj that particle resonants with
            call Resonance(Model,wave,wModel,prt,astar,xj,yj)
            if (xj == 999) return ! no resonant waves present

            Daa = DiffCoef(Model,wave,wModel,prt,astar,MagB,xj,yj) ! Calculate the diffusion coeffcient

            dtRem = dt - dtCum
            if (sqrt(2.0*dt*Daa) > abs(dAlim)) then
                ddt = dAlim**2/(2.*Daa)
                if (dtRem < ddt) ddt = dtRem ! So don't overshoot global step
            else
                ddt = dt
            endif

            dtCum = dtCum + ddt

            ! Calculate the resulting change in pitch angle and energy of the particle
            call DiffCurve(Model,wave,prt,Daa,ddt,xj,yj,da,dp)

            !Update the pitch angle
            aNew = prt%alpha + da
            prt%alpha = aNew

            !reflecting pitch angle at 0° and 180° (to keep in this range)
            if (aNew .lt. 0) then 
                aNew = abs(aNew)
            else if (aNew .gt. PI) then
                aNew = 2.0*PI-aNew
            endif

            if (aNew <= PI/2) then
                pSgn = +1
            else
                pSgn = -1
            endif

            !Updating the momentum components
            gamma = prt2Gam(prt,Model%m0)
            pMag = Model%m0*sqrt(gamma**2.0 - 1.0)+dp !dp is scalar and change in total momentum
           
            !updating wpi diagnostic variables
            prt%dAwpi = prt%dAwpi + da
            K0 = prt2kev(Model,prt)

            if (prt%isGC) then
                !Scatter GC particle
                p11Mag = pSgn*pMag*sqrt( 1 - sin(aNew)**2.0 )
                gamNew = sqrt(1+(pMag/Model%m0)**2.0)
                Mu = ( (Model%m0**2)*(gamNew**2.0 - 1.0) - p11Mag*p11Mag) / (2*Model%m0*MagB)

                prt%Q(P11GC) = p11Mag
                prt%Q(MUGC ) = Mu 
                prt%Q(GAMGC) = gamNew

                !updating wpi K diagnostic
                K1 = prt2kev(Model,prt)
                dK = (K1 - K0)
                prt%dKwpi = prt%dKwpi + dK

                ! catch to see if momentum becomes nan and when updated
                if (isnan(pMag) .or. isnan(p11Mag) .or. isnan(gamNew)) then
                    write(*,*) 'PERFORMWPI:: ERROR:: A NAN OCCURRED IN THE MOMENTUM UPDATE'
                    write(*,*) 'Resonant wave: xj,yj:  ',xj,yj 
                    write(*,*) 'Daa,  da,  dp: ', Daa,da,dp
                    write(*,*) 'new ps and p of tp: ', aNew,pMag
                    write(*,*) 'new p11 and gamma of tp: ', p11Mag,gamNew
                    stop
                endif
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
    subroutine ChkWave(wModel,r,req,doWave)
        type(wModel_T), intent(in) :: wModel
        real(rp), dimension(NDIM), intent(in) :: r,req !Particle location and coordinates of last equatorial crossing
        logical, intent(inout) :: doWave !Holds if waves are present 
        real(rp), dimension(2) :: Lbds, LONbds,MLATbds ! Hold the range where waves are present
        real(rp) :: lat,lon,L,rMag
        logical :: inL,inLat,inLon,inWaves

        !FIXME: need more accurate wave model
        ! Setting simple limits in L, MLT and MLAT for now
        Lbds(1)    = 4.5
        Lbds(2)    = 8.0
        LONbds(1)  = 4*PI/5 ! 144 degrees, 21 MLT
        LONbds(1)  = -PI/12 ! -15 degrees, 11 MLT
        MLATbds(1) = -PI/9  ! -20 degrees
        MLATbds(2) = PI/9   ! -20 degrees

        rMag = norm2(r)

        lat = asin(r(3)/rMag)
        lon = atan2(r(2),r(1))
        L = norm2(req)
        !L = rMag/cos(lat)**2 !FIXME: assumes dipolar field

        ! check if particle is in presence of waves
        inL = (L .gt. Lbds(1) .and. L .le. Lbds(2))
        inLon = (lon .gt. LONbds(1) .or. lon .lt. LONbds(2))
        inLat = (lat .gt. MLATbds(1) .and. lat .lt. MLATbds(2))
        inWaves = inL .and. inLat .and. inLon

        if (inWaves) then
            doWave = .true.
        else   
            doWave = .false.
        end if 

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
        real(rp), dimension(:), allocatable :: xcs, ycs
        real(rp) :: pa,mu,K,beta,pa90p,pa90m

        pa90p = 0.5*PI+5.0E-6
        pa90m = 0.5*PI-5.0E-6
        pa = prt%alpha
        if (wave%mode .eq. "eRW".and. pa .ge. pa90m .and. pa .le. pa90p) then
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
        K = prt2kev(Model,prt)/(Model%m0*mec2*1.0e+3) 
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
            xpres = pack(xjs, (xjs >= xMin .and. xjs <= xMax .and. yjs > 0))
            ypres = pack(yjs, (xjs >= xMin .and. xjs <= xMax .and. yjs > 0))
        else
            xpres = pack(xjs, (xjs >= xMin .and. xjs <= xMax .and. yjs < 0))
            ypres = pack(yjs, (xjs >= xMin .and. xjs <= xMax .and. yjs < 0))
        end if
        
        ! Check if no resonant waves are present
        if (size(xpres) == 0) then  
            xj = 999
            yj = 999
        else if (size(xpres) > 1) then
            ! Should only be one resonant wave
            write(*,*) 'Particle gamma and pitch angle: ', prt%Q(GAMGC),pa
            write(*,*) 'All roots (w/|Ome|): ', xjs
            write(*,*) 'All roots (kc/|Ome|): ', yjs
            write(*,*) 'Too many resonant roots (w/|Ome|, kc/|Ome|): ', xpres, ypres
            stop
        else
            xj = xpres(1)
            yj = ypres(1)
        end if

    end subroutine selectRoot

    !!!!!!!!!!!FIXME: need to add capability to sum over mulitple roots!!!!!!!!!!!!
    !Calculates diffusion coefficient assuming wave spectrum is Gaussian (Summers 2005 eq 33)
    function DiffCoef(Model,wave,wModel,prt,astar,B0,xj,yj) result(Daa)
        type(chmpModel_T), intent(in) :: Model
        type(wave_T), intent(in) :: wave
        type(wModel_T), intent(in) :: wModel
        type(prt_t), intent(in) :: prt
        real(rp), intent(in) :: astar,B0,xj,yj
        real(rp) :: pa,K,beta,R,Ome,DScl,Fxy,Daa

        pa= prt%alpha
        K = prt2kev(Model,prt)/(Model%m0*mec2*1.0e+3) 
        beta = sqrt(K*(K+2.0))/(K+1.0) ! beta = v/c
        R = (wModel%B1/B0)**2  !ratio of the wave amplitude to background field strength
        Ome = B0 ! normalized non-relativistic electron gyrofrequency, has Om^2/Ome therefore dont need sign of q

        DScl = (PI/2.0)*abs(Ome)*(K+1)**(-2.0)

        Fxy = Vg(wave,astar,xj,yj)

        Daa = R*(1.0-xj*cos(pa)/(yj*beta))**2*(abs(Fxy)/abs(beta*cos(pa)-Fxy))*waveSpec(wModel,xj)

        Daa = DScl*Daa

    end function DiffCoef

    !Calculates the change in pitch-angle and corresponding change in momentum along the diffusion curve 
    subroutine DiffCurve(Model,wave,prt,Daa,dt,xj,yj,da,dp,constDA) 
        type(chmpModel_T), intent(in)    :: Model
        type(wave_T),      intent(in)    :: wave
        type(prt_t),       intent(in)    :: prt
        real(rp),          intent(in)    :: Daa,dt,xj,yj
        real(rp),          intent(inout) :: da,dp
        real(rp), intent(in), optional   :: constDA ! use a constant change in pitch angle if desired
        real(rp) :: eta,u,gamu,pa,gamma,pMag,E,A,B

        if (present(constDA)) then
            da = constDA
        else
            eta = genRand(-1.0_rp,1.0_rp) !generating random number between -1 and 1 for random walk eq for da
            da = sqrt(2.0*dt*Daa)*eta
        end if

        u = xj/yj ! phase velocity of the wave normalized by c
        gamu = sqrt(1-u**2.)**(-1.0)
        
        pa = prt%alpha
        gamma = prt2Gam(prt,Model%m0)
        pMag = Model%m0*sqrt(gamma**2.0-1.0)
        E = prt2kev(Model,prt)/(Model%m0*mec2*1.0e+3)+1.0 !full energy of the particle including rest mass

        A = gamu**2.0*(pMag*cos(pa)-u*E)*pMag*sin(pa)-pMag**2.0*sin(pa)*cos(pa)
        B = gamu**2.0*(pMag*cos(pa)-u*E)*(cos(pa)-u*pMag/E)+pMag*sin(pa)**2.0

        dp = da*A/B

    end subroutine DiffCurve

end module wpicalc
