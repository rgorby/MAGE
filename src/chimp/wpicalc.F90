!Various routines to handle metrics at particle crossings
module wpicalc
    use chmpdefs
    use tptypes
    use ebtypes
    use tputils
    use wpitypes
    use wpihelper
    use earthhelper
    use xml_input
    use chmpunits
    use quarticRoots
    
    implicit none

    ! sets subcycling limit in change pitch-angle to be below this value 
    ! to reduce error and not allow large resonance times 
    real(rp), private :: dAalim = 0.0087_rp 

    real(rp), parameter :: dtW = 2.0 !Max increase in timestep, dtNew <= dtX*dtOld

    !scales ratio of e- gyro to plasma frequency in code units to real value
    real(rp), parameter :: inWScl = (qe_cgs*Re_cgs/vc_cgs/sqrt(Me_cgs))**2.0 ! Assumes L0 = Re_cgs

    contains

    !FIXME: include capability of particle interaction with multiple waves
    !Main subroutine for wave particle inteactions
    subroutine PerformWPI(prt,t,dt,Model,ebState,constB0,constNe)
        type(prt_t), intent(inout) :: prt
        real(rp), intent(in) :: t,dt ! Note: dt is substep size for particle (ddt in pusher.F90)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), optional, intent(in)   :: constB0,constNe ! use constant values for B0 and ne, usefull for testing
        type(wave_T) :: wave
        type(wModel_T) :: wModel
        
        real(rp), dimension(NDIM) :: r,p,E,B,vxb,xhat,yhat,bhat,req
        real(rp), dimension(NDIM) :: p11,pxy,vExB
        real(rp), dimension(NVARMHD) :: Qmhd
        complex(rp), dimension(2) :: quadCoef,dtLim ! used to set subcycle dt given dAlim
        real(rp) :: ebGam,gamNew,aNew,pNew,MagB,Mu,p11Mag,K,rho,psi,L,phi
        real(rp) :: Ome,wpe,astar,xHigh,xj,yj,dGDda,Daa
        real(rp) :: dtCum,dtRem,ddt,daMax ! substep used in wpi calculations, not same as ddt in pusher.F90 
        real(rp) :: zEq = 0.0 !Defined Z for equator
        real(rp) :: aCoef,bCoef !coefficients to LangevinEq
        real(rp) :: p0,G,gamma,ebGam0,p0eb,pMagEB
        integer  :: pSgn

        real(rp) :: pa,Kprt ! pitch-angle and energy of particle in subcycle loop
        real(rp) :: da,dp,dK,dgam !Changes due to wpi
        real(rp) :: K0,K1,mu0,p110,Pperp20,Pperp21 ! energy before and after wpi
        real(rp) :: htwpiB,q,B11,vWPI,htNew,htOld
        real(rp) :: vGC(NVARTP)

        ! helps sets subcycling limit in change pitch-angle to be below this value 
        ! to reduce error in change in energy not to be to large
        real(rp) :: dAlim,dAplim,deltaPdeltaA,u,epsp,deltaP

        logical :: doWave

        if (.not.prt%isIn) return !shouldn't be here if particle is not in domain

        !initializing some variables
        pSgn = 1; da=0.0; dK=0.0; epsp = 0.005_rp ! change in momentum is below this limit
        vGC = 0.0

        ! pulling wave information
        wave = ebState%ebWave
        wModel = ebState%ebWmodel

        !Get local coordinate system
        r = prt%Q(XPOS:ZPOS)
        req = [prt%Qeq(EQX),prt%Qeq(EQY),zEq] ! defining magnetic equator at z=0

        !Check to see if waves are present
        call ChkWave(wModel,r,req,prt%alpha,doWave)
        if (.not. doWave) return ! no waves so exit

        !Get MHD density
        if (present(constNe)) then
            rho = constNe
        else 
            if (Model%doMHD) then
                Qmhd = mhdInterp(r,t,Model,ebState)
                rho = Qmhd(DEN)
            else 
                ! Pull density from Gallagher plasmasphere
                L = norm2(req)
                phi = atan2(req(2),req(1))
                rho = GallagherRP(L,phi)  ! use Gallagher from earthhelper
            endif
        endif

        if (rho < TINY) return ! no plasma to interact with so exit

        !pull background magnetic field strength
        if (present(constB0)) then
            MagB = constB0
            vExB = 0.0_rp
            E = 0.0_rp
        else 
            call ebFields(r,t,Model,ebState,E,B,ijkO=prt%ijk0,vExB=vExB)
            call MagTriad(r,B,xhat,yhat,bhat)
            MagB = max(norm2(B),TINY)
        endif
        !Calulating the ratio of nonrelativistic gyrofrequency to plasma frequency at prt's location
        Ome = MagB 
        wpe = sqrt(4*PI*rho) 
        astar = Ome**2/((wpe**2)*inWScl)

        !calculate minimum energy needed to resonate
        xHigh = wModel%xm+wModel%Dx
        wave%emin = Kres_whistleR(xHigh,astar,Model%m0)

        !check if particle above minimum required energy
        K = pGC2K(Model,prt)
        if (K < wave%emin) return

        !Setting adaptive time-step/particle sub-stepping to limit da and reduce error
        dtCum = 0.0 !How far we've advanced
        ! ddt = dt
        do while (dtCum<dt)
            ! Get pitch-angle and energy of the particle
            pa = prt%alpha
            Kprt = pGC2K(Model,prt)
            ! Determine wj,kj that particle resonants with
            call Resonance(wave,wModel,Model%m0,Kprt,pa,astar,xj,yj)
            if (xj == 999) return ! no resonant waves present
            prt%xj = xj
            prt%yj = yj

            Daa = DiffCoef(wave,wModel,Model%m0,K,pa,astar,MagB,xj,yj) ! Calculate the diffusion coeffcient
            dGDda = derivGD(wave,wModel,Model%m0,K,pa,astar,MagB)

            gamma = prt2Gam(prt,Model%m0)
            p0 = Model%m0*sqrt(gamma**2.0-1.0) 
            G = sin(pa)*p0*p0

            !coefficients for Langevin Equation
            aCoef = dGDda/G
            bCoef = sqrt(2.0*Daa)/p0

            daMax = abs(aCoef)*dt+bCoef*sqrt(dt)
            ! daMax = abs(aCoef)*ddt+bCoef*sqrt(ddt)

            !limit da to change p by epsp
            u = xj/yj ! phase velocity of the wave normalized by c
            deltaPdeltaA = dpda(pa,p0,Model%m0,u)
            dAplim = epsp*p0/deltaPdeltaA !also have dAlim in case dpda is very small

            dtRem = dt - dtCum
            if (abs(daMax) > abs(dAalim) .or. abs(daMax) > abs(dAplim)) then
                dAlim = min(dAalim,dAplim)!use value that is the smallest
                if (aCoef < TINY) then
                    ddt = (dAlim/bCoef)**2
                else 
                    quadCoef = [-(2.*aCoef*abs(dAlim)+bCoef**2)/aCoef**2,(dAlim/aCoef)**2.0]
                    dtLim = quadraticSolve(quadCoef)
                    ddt = minval(real(dtLim)) ! only real roots but only taking real part
                end if
            else
                ddt = dt
            endif
            if (dtRem < ddt) ddt = dtRem ! So don't overshoot global step
            dtCum = dtCum + ddt

            ! Calculate the resulting updated pitch angle and energy of the particle
            call LangevinEq(wave,wModel,Model,prt,aCoef,bCoef,ddt,astar,aNew,pNew)

            !reflecting pitch angle at 0° and 180° (to keep in this range)
            if (aNew < 0) then 
                aNew = abs(aNew)
            else if (aNew > PI) then
                aNew = 2.0*PI-aNew
            endif

            if (aNew <= PI/2) then
                pSgn = +1
            else
                pSgn = -1
            endif
           
            !Update the pitch angle
            prt%alpha = aNew

            !updating wpi diagnostic variables
            da = aNew - pa
            prt%dAwpi = prt%dAwpi + da
            K0 = prt2kev(Model,prt)

            if (prt%isGC) then
                !Scatter GC particle
                K1 = sqrt(pNew**2.+Model%m0**2.) - Model%m0
                dK = K1-Kprt
                dgam = dK/Model%m0
                gamNew = gamma+dgam
                ! gamNew = sqrt(1+(pNew/Model%m0)**2.0)

                p11Mag = pSgn*pNew*sqrt( 1 - sin(aNew)**2.0 )
                
                ! Mu update occurs in ExB frame 
                ebGam = max(1.0,gamNew - 0.5*dot_product(vExB,vExB))
                pMagEB = Model%m0*sqrt(ebGam**2.0 - 1.0)
                ! p11eb = pSgn*pMagEB*sqrt( 1 - sin(aNew)**2.0 )

                mu0 = prt%Q(MUGC )
                ebGam0 = max(1.0,gamma - 0.5*dot_product(vExB,vExB))
                p0eb = Model%m0*sqrt(ebGam0**2.0 - 1.0)
                p110 = prt%Q(P11GC)
                Pperp20 = p0eb*p0eb-p110*p110
                Pperp21 = pMagEB*pMagEB-p11Mag*p11Mag

                !Old 
                ! Mu = ( (Model%m0**2)*(ebGam**2.0 - 1.0) - p11Mag*p11Mag) / (2*Model%m0*MagB)
                if (Pperp20 < TINY) then
                    !particle will most likely precipitate
                    Mu = ( Pperp21) / (2*Model%m0*MagB)
                else
                    Mu = (Pperp21/Pperp20)*mu0
                endif

                prt%Q(P11GC) = p11Mag
                prt%Q(MUGC ) = Mu 
                prt%Q(GAMGC) = gamNew

                !updating wpi K diagnostic
                K1 = prt2kev(Model,prt)
                dK = (K1 - K0)
                prt%dKwpi = prt%dKwpi + dK

                ! catch to see if momentum becomes nan and when updated
                if (isnan(pNew) .or. isnan(p11Mag) .or. isnan(gamNew)) then
                    write(*,*) 'PERFORMWPI:: ERROR:: A NAN OCCURRED IN THE MOMENTUM UPDATE'
                    write(*,*) 'Resonant wave: xj,yj:  ',xj,yj 
                    write(*,*) 'Daa,  da,  dK: ', Daa,da,dK
                    write(*,*) 'new pa and p of tp: ', aNew,pNew
                    write(*,*) 'new p11 and gamma of tp: ', p11Mag,gamNew
                    stop
                endif

                !Calculate new Timestep, dt ~ |p|/F, F = Lorentz(vGC) + dP11/dt
                ! htgcT = m0*sqrt(gamma**2.0-1) ! |P|
                ! q  = Model%q0
                !!!!!!! Need to use effective B* and E* as in gcutils:DerivGC?!!!!!!!!!!!!!
                ! B11 = dot_product(bhat,B) !Parallel comp. of B
                ! ebGam = max(1.0,gamNew - 0.5*dot_product(vExB,vExB))
                ! vGC(XPOS:ZPOS) = cross(E,bhat)/B11 + p11Mag*B/(Model%m0*ebGam*B11)
                ! vGC(P11GC)     = q*dot_product(B,E)/B11
                ! vxb = q*cross(vGC(XDIR:ZDIR),B) !Lorentz force on GC velocity
                ! htwpiB = sqrt( norm2(vxb)**2.0 + norm2(q*E)**2.0 + vGC(P11GC)**2.0 )
                ! htNew = Model%epsht*(pNew/max(htwpiB,TINY))
            else                
                !Update momentum for FO particle
                p11 = pSgn*pNew*sqrt( 1 - sin(aNew)**2.0 )
                ! Assume gyrophase remians the same after scattering
                psi = prt%Q(PSITP) 
                pxy = pNew*sin(aNew)*( cos(psi)*xhat + sin(psi)*yhat )

                p = p11*bhat + pxy
                prt%Q(PXFO:PZFO) = p
                
            endif
            prt%Nwpi = prt%Nwpi + 1
        enddo
        !update particle time-step
        ! if (pNew > p0) then
        !     ddt = prt%ddt*(vMag_old/vMag_new)
        ! end if

        ! htOld = prt%ddt
        ! prt%ddt = min(htNew,dtW*htOld)
    end subroutine PerformWPI

    !Check to see if waves are present at location 
    subroutine ChkWave(wModel,r,req,pa,doWave)
        type(wModel_T), intent(in) :: wModel
        real(rp), dimension(NDIM), intent(in) :: r,req !Particle location and coordinates of last equatorial crossing
        real(rp), intent(in)   :: pa ! pitch-angle of particle 
        logical, intent(inout) :: doWave !Holds if waves are present 
        real(rp), dimension(2) :: Lbds, LONbds,MLATbds ! Hold the range where waves are present
        real(rp) :: lat,lon,L,rMag
        logical :: inL,inLat,inLon,inWaves,offEq,inNorth,inSouth

        !FIXME: need more accurate wave model
        ! Setting simple limits in L, MLT and MLAT for now
        Lbds(1)    = 4.5_rp
        Lbds(2)    = 8.0_rp
        ! LONbds(1)  = 4*PI/5 ! 144 degrees, 21 MLT
        ! LONbds(2)  = -PI/12 ! -15 degrees, 11 MLT
        ! MLATbds(1) = -PI/9  ! -20 degrees
        ! MLATbds(2) = PI/9   ! -20 degrees

        !!!!!!!!!!!!!!! for dipole testing !!!!!!!!!!!!!!!!!
        MLATbds(1) = -PI/18  ! 10 degrees
        MLATbds(2) = PI/18   ! 10 degrees
        inLon = .true. ! uncomment line below when remove!!!!!!!!!!!!!

        rMag = norm2(r)

        lat = asin(r(3)/rMag)
        lon = atan2(r(2),r(1))
        L = norm2(req)

        ! check if particle is in presence of waves
        inL = (L >= Lbds(1) .and. L <= Lbds(2))
        ! inLon = (lon >= LONbds(1) .or. lon <= LONbds(2))
        inLat = (lat >= MLATbds(1) .and. lat <= MLATbds(2))

        ! only include waves propagating off equator, resonance is with waves propagating 
        ! in opposite direction of prt, therefore particle needs to be moving toward equator
        inNorth = (pa >= PI/2 .and. lat >= 0)
        inSouth = (pa <= PI/2 .and. lat < 0)
        offEq = (inNorth .or. inSouth)

        inWaves = inL .and. inLat .and. inLon .and. offEq

        if (inWaves) then
            doWave = .true.
        else   
            doWave = .false.
        end if 

    end subroutine ChkWave

    ! Determines the wave in resosonance with the particle from all possible roots
    subroutine Resonance(wave,wModel,m0,K,pa,astar,xj,yj,includeAllWaves)
        type(wave_T), intent(in) :: wave
        type(wModel_T), intent(in) :: wModel
        real(rp), intent(in) :: m0,K,pa,astar
        real(rp), intent(inout) :: xj,yj
        logical, optional, intent(in) :: includeAllWaves
        real(rp), dimension(:), allocatable :: xjs, yjs
        real(rp), dimension(:), allocatable :: xcs, ycs
        real(rp) :: mu,beta,pa90p,pa90m,Km
        logical :: doAllWaves

        ! pa90p = 0.5*PI+5.0E-6
        ! pa90m = 0.5*PI-5.0E-6
        pa90p = 0.5*PI+5.0E-6
        pa90m = 0.5*PI-5.0E-6

        if (wave%mode .eq. "eRW".and. pa .ge. pa90m .and. pa .le. pa90p) then
            call res90deg(wave,m0,K,astar,xjs,yjs)
        else   
            call ResRoots(wave,m0,K,pa,astar,xjs,yjs)
        end if

        ! taking desired root from all physical roots
        if (present(includeAllWaves)) then
            doAllWaves = includeAllWaves
        else
            doAllWaves = .false.
        endif
        call selectRoot(wModel,pa,xjs,yjs,xj,yj,doAllWaves)

        if (xj == 999) return ! no resonant roots

        ! Check to see if wave is a critical root (causes singularity in Daa if not fixed)
        ! can occur for electrons and R-mode waves and protons with L-mode waves
        mu = cos(pa) 
        Km = K/m0
        beta = sqrt(Km*(Km+2.0))/(Km+1.0) ! beta = v/c

        if (beta*mu == Vg(wave,astar,xj,yj)) then
            call criticalRoot(wave,m0,K,astar,xcs,ycs)
            call selectRoot(wModel,pa,xcs,ycs,xj,yj,doAllWaves)
        end if

    end subroutine Resonance

    ! Determining if resonant waves are present in the wave model at location
    ! and considering only waves propagating in opposite direction of particle
    subroutine selectRoot(wModel,pa,xjs,yjs,xj,yj,doAllWaves)
        type(wModel_T), intent(in) :: wModel
        real(rp), intent(in)       :: pa
        real(rp), intent(inout)    :: xj,yj
        logical,  intent(in)       :: doAllWaves
        real(rp), dimension(:), allocatable, intent(in) :: xjs, yjs
        real(rp), dimension(:), allocatable :: xpres, ypres
        real(rp) :: xMin,xMax
        
        if (doAllWaves) then
            xMin = 0
            xMax = 1
        else
            xMin = wModel%xm-wModel%Dx
            xMax = wModel%xm+wModel%Dx
        end if
 
        ! particle resonantes with wave propogating in opposite direction
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
            write(*,*) 'Particle pitch angle: ', pa
            write(*,*) 'All roots (w/|Ome|): ', xjs
            write(*,*) 'All roots (kc/|Ome|): ', yjs
            write(*,*) 'Too many resonant roots (w/|Ome|, kc/|Ome|): ', xpres, ypres
            stop
        else
            xj = xpres(1)
            yj = ypres(1)
        end if

    end subroutine selectRoot

    !Calculates the derivative of the diffusion coefficient
    !Equation 5.7.7 from Numerical Recipes 
    function derivGD(wave,wModel,m0,K,pa,astar,B0) result(dGDdx)
        type(wave_T), intent(in) :: wave
        type(wModel_T), intent(in) :: wModel
        real(rp), intent(in) :: m0,K,pa,astar,B0
        real(rp) :: eps,h,xStp,yStp,ap,am,pMag
        real(rp) :: Dp,Dm,dGDdx
        logical :: doAllWaves

        doAllWaves = .true.

        pMag = sqrt(K*(K+2.0*m0))

        eps = epsilon(1.0_rp) !machine error
        if (pa < dAalim) then
            h = eps*dAalim !to avoid setting h to zero
        else
            h = eps*pa
        end if

        ap = pa + h
        call Resonance(wave,wModel,m0,K,ap,astar,xStp,yStp,doAllWaves)
        Dp = DiffCoef(wave,wModel,m0,K,ap,astar,B0,xStp,yStp)

        am = pa - h
        call Resonance(wave,wModel,m0,K,am,astar,xStp,yStp,doAllWaves)
        Dm = DiffCoef(wave,wModel,m0,K,am,astar,B0,xStp,yStp)

        dGDdx = pMag**2.*(sin(ap)*Dp - sin(am)*Dm)/(2.0*h)

    end function derivGD

    !Calculates the change in pitch-angle and corresponding change in momentum along the diffusion curve 
    subroutine LangevinEq(wave,wModel,Model,prt,aCoef,bCoef,dt,astar,a1,p1,constDA) 
        type(wave_T),   intent(in) :: wave
        type(wModel_T), intent(in) :: wModel
        type(chmpModel_T), intent(in)    :: Model
        type(prt_t),       intent(in)    :: prt
        real(rp),          intent(in)    :: aCoef,bCoef,dt,astar
        real(rp),          intent(inout) :: a1,p1   ! new pitch-angle and momentum of particle after wpi
        real(rp), intent(in), optional   :: constDA ! use a constant change in pitch angle if desired
        logical, dimension(6)            :: isRoot ! array to check if any partial step of integrator has no resonant roots
        real(rp) :: eta,da
        real(rp) :: gamma

        if (present(constDA)) then
            da = constDA
        else
            eta = genRand(-1.0_rp,1.0_rp) !generating random number between -1 and 1 for random walk eq for da
            da = aCoef*dt + bCoef*eta*sqrt(dt)
        end if

        !catch to not go past alpha = 90 degrees
        ! if ((prt%alpha < PI/2 .and. prt%alpha+da > PI/2) .or. (prt%alpha > PI/2 .and. prt%alpha+da < PI/2)) then
        !     da = aSgn*abs(a1-PI/2)
        !     write(*,*) 'RKF45: da over 90, fixing:  ',a1,da
        ! end if

        !push back when get to alpha = 90 degrees so dont cross 
        ! if ((prt%alpha < PI/2 .and. prt%alpha+da > PI/2) .or. (prt%alpha > PI/2 .and. prt%alpha+da < PI/2)) then
        !push back when get to alpha = 70 degrees so dont cross 
        ! if ((prt%alpha < 1.22173 .and. prt%alpha+da > 1.22173) .or. (prt%alpha > 2.0944 .and. prt%alpha+da < 2.0944)) then
        ! if ((prt%alpha < 0.698132 .and. prt%alpha+da > 0.698132) .or. (prt%alpha > 2.44346 .and. prt%alpha+da < 2.44346)) then

        !keep between 85 and 70
        ! if ((prt%alpha > 1.22173 .and. prt%alpha+da < 1.22173) .or. (prt%alpha < 1.48353 .and. prt%alpha+da > 1.48353)) then
        ! if ((prt%alpha > 1.6581 .and. prt%alpha+da < 1.6581) .or. (prt%alpha < 1.9199 .and. prt%alpha+da > 1.9199)) then
        !     da = -1.0*da
        ! end if

        isRoot = .false.
        do while(.not.all(isRoot))
            call rkf45(wave,wModel,prt,Model%m0,astar,da,isRoot,p1)
            if (.not.all(isRoot)) then
                da = da/2.0_rp ! if da too large, reduce 
                write(*,*) "lowering da", isRoot
            end if
        enddo

        !update pitch angle with final da
        a1 = prt%alpha + da

    end subroutine LangevinEq

    !resonant diffusion surface equation for dp/da for fixed wave phase speed
    function dpda(pa,pMag,m0,u) result(f)
        real(rp), intent(in)    :: pa,pMag,m0,u
        real(rp) :: gamu,E,A,B,f

        gamu = sqrt(1-u**2.)**(-1.0)
        
        E =  sqrt(pMag**2.+m0**2.) !full energy of the particle including rest mass (code units)

        A = gamu**2.0*(pMag*cos(pa)-u*E)*pMag*sin(pa)-pMag**2.0*sin(pa)*cos(pa)
        B = gamu**2.0*(pMag*cos(pa)-u*E)*(cos(pa)-u*pMag/E)+pMag*sin(pa)**2.0

        f = A/B

    end function dpda

    ! adaptive step RK following the Runge-Kutta-Fehlberg Method (RKF45)
    ! http://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
    subroutine rkf45(wave,wModel,prt,m0,astar,da,isRoot,p1) 
        type(wave_T),   intent(in) :: wave
        type(wModel_T), intent(in) :: wModel
        type(prt_t),    intent(in) :: prt
        real(rp),       intent(in) :: m0,astar,da
        logical, dimension(6), intent(inout) :: isRoot
        real(rp),       intent(out):: p1

        real(rp) :: u,a1,aStp,pStp,kStp,xStp,yStp,uStp,y1
        real(rp) :: tol,h,s,dah,daRem,unity,aSgn,gamma
        real(rp) :: k1,k2,k3,k4,k5,k6 
        integer  :: Nmax,i
        logical :: doAllWaves

        tol = 1.0e-8
        h   = da/100.0 !da & h holds +/- direction 
        unity = 1.0
        doAllWaves = .true.

        aSgn = sign(unity,da)
        a1 = prt%alpha

        gamma = prt2Gam(prt,m0)
        p1 = m0*sqrt(gamma**2.0-1.0)

        ! u = prt%xj/prt%yj ! phase velocity of the wave normalized by c

        i  = 0
        dah = 0
        do while(abs(dah) < abs(da)) 
            daRem = aSgn*(abs(da) - abs(dah))
            if (abs(daRem) < abs(h))then 
                h = daRem ! So don't overshoot 
            end if
            !Apply Runge-Kutta-Fehlberg formulas to find next value of p
            aStp = a1
            pStp = p1
            kStp = sqrt(pStp**2.+m0**2.) - m0
            call Resonance(wave,wModel,m0,kStp,aStp,astar,xStp,yStp,doAllWaves)
            if (xStp == 999) then 
                isRoot(1) = .false.
            else
                isRoot(1) = .true.
            endif
            uStp = xStp/yStp
            k1 = h * dpda(a1, p1, m0, uStp) 

            aStp = a1+0.25*h
            pStp = p1+0.25*k1
            kStp = sqrt(pStp**2.+m0**2.) - m0
            call Resonance(wave,wModel,m0,kStp,aStp,astar,xStp,yStp,doAllWaves)
            if (xStp == 999) then 
                isRoot(2) = .false.
            else
                isRoot(2) = .true.
            endif
            uStp = xStp/yStp
            k2 = h * dpda(aStp, pStp, m0, uStp) 

            aStp = a1+0.375*h
            pStp = p1+0.09375*k1+0.28125*k2
            kStp = sqrt(pStp**2.+m0**2.) - m0
            call Resonance(wave,wModel,m0,kStp,aStp,astar,xStp,yStp,doAllWaves)
            if (xStp == 999) then 
                isRoot(3) = .false.
            else
                isRoot(3) = .true.
            endif
            uStp = xStp/yStp
            k3 = h * dpda(aStp, pStp, m0, uStp) 

            aStp = a1+12*h/13
            pStp = p1 + (1932*k1-7200*k2+7296*k3)/2197
            kStp = sqrt(pStp**2.+m0**2.) - m0
            call Resonance(wave,wModel,m0,kStp,aStp,astar,xStp,yStp,doAllWaves)
            if (xStp == 999) then 
                isRoot(4) = .false.
            else
                isRoot(4) = .true.
            endif
            uStp = xStp/yStp
            k4 = h * dpda(aStp, pStp, m0, uStp)

            aStp = a1+h
            pStp = p1 + (8341*k1-32832*k2+29440*k3-845*k4)/4104
            kStp = sqrt(pStp**2.+m0**2.) - m0
            call Resonance(wave,wModel,m0,kStp,aStp,astar,xStp,yStp,doAllWaves)
            if (xStp == 999) then 
                isRoot(5) = .false.
            else
                isRoot(5) = .true.
            endif
            uStp = xStp/yStp 
            k5 = h * dpda(aStp, pStp, m0, uStp) 

            aStp = a1+0.5*h
            pStp = p1 + (-6080*k1+41040*k2-28352*k3-9295*k4-5643*k5)/20520
            kStp = sqrt(pStp**2.+m0**2.) - m0
            call Resonance(wave,wModel,m0,kStp,aStp,astar,xStp,yStp,doAllWaves)
            if (xStp == 999) then 
                isRoot(6) = .false.
            else
                isRoot(6) = .true.
            endif
            uStp = xStp/yStp
            k6 = h * dpda(aStp, pStp, m0, uStp) 

            ! Update next value of p 
            y1 = p1 + (3246625*k1+15397888*k3+15027480*k4-5610168*k5)/28050840
            p1 = p1 + (33440*k1+146432*k3+142805*k4-50787*k5+10260*k6)/282150

            if(.not.all(isRoot)) exit !if encounter no resonant roots break out & reduce da

            ! Update next value of a1 
            a1 = a1 + h 
            dah = dah + h

            ! Determine next step size
            if (p1 == y1) then 
                s = 1.0 
            else
                s = (tol/(2*abs(p1-y1)))**0.25
            end if

            h = s*h
        enddo

    end subroutine rkf45

end module wpicalc
