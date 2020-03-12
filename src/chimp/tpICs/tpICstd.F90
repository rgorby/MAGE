!Standard initialization routine for test particles
!Must specify initParticles / addIncoming
module usertpic
    use chmpdefs
    use tptypes
    use ebtypes
    use ebinterp
    use tputils
    use streamline
    use xml_input
    use math

    implicit none

    !Module variables to handle solar wind
    !zWind = value to bounce TPs back towards equator
    logical :: doWind = .false.
    real(rp) :: zWind = 0.0_rp

    !Values to handle outflow
    logical :: doOutflow = .false.

    contains

    !Standard TP init routine
    !Create equatorial particles w/ random L,phi,K,alpha
    subroutine initParticles(Model,ebState,tpState,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(inout)   :: ebState
        type(tpState_T), intent(inout)   :: tpState
        type(XML_Input_T), intent(inout) :: inpXML

        character(len=strLen) :: sStr
        integer :: Np,n,id0
        real(rp) :: rmin,rmax,kmin,kmax,pmin,pmax,amin,amax
        logical :: doLogK,doLogR
        real(rp), dimension(:), allocatable :: radius,energy,alpha,phi,psi,hgt,T0ps

        write(*,*) 'Initializing test particles ...'

        !Get species data
        call inpXML%Set_Val(sStr,"tps/species","H")
        call getSpecies(sStr,Model%m0,Model%q0)

        !Get number
        call inpXML%Set_Val(Np,"tps/Np",100)

        !Create basic setup
        tpState%Np = Np
        allocate(tpState%TPs(Np))

        tpState%TPs(:)%isIn = .false.
        tpState%TPs(:)%isInit = .false.
        
        if (Model%imeth == IFO .or. Model%imeth == IDYN) then
            !Start dynamic particles as full orbit
            tpState%TPs(:)%isGC = .false.
        else
            tpState%TPs(:)%isGC = .true.
        endif

        !Get starting ID
        id0 = 1+Np*(Model%Nblk-1)

    !Get phase space bounds
        !Radius (L0), Energy (keV), phi/alpha/psi (degrees)
        call getSample(Model,inpXML,"radius",Np,radius,5.0_rp, 25.0_rp)
        call getSample(Model,inpXML,"energy",Np,energy,1.0_rp,100.0_rp)
        call getSample(Model,inpXML,"phi"   ,Np,phi   ,0.0_rp,360.0_rp)
        call getSample(Model,inpXML,"alpha" ,Np,alpha ,0.0_rp,180.0_rp)
        call getSample(Model,inpXML,"psi"   ,Np,psi   ,0.0_rp,360.0_rp)
        call getSample(Model,inpXML,"height",Np,hgt   ,0.0_rp,  0.0_rp)
        call inpXML%Set_Val(Model%doEBInit,"energy/doEBInit",.false.)

        !Do TP birthdays
        if (Model%doStream) then
            !Generate random start times (specify defaults in unscaled time)
            call getSample(Model,inpXML,"stream",Np,T0ps,Model%T0/inTScl,Model%tFin/inTScl)
            tpState%TPs(:)%T0p = inTScl*T0ps
        else
            tpState%TPs(:)%T0p = Model%T0 !Already scaled
        endif

        !Handle solar wind TP stuff
        call inpXML%Set_Val(doWind,'height/doWind',.false.)
        if (doWind) then
            call inpXML%Set_Val(zWind,'height/zWind',zWind)
        endif


        !Convert to code units where necessary
        phi = phi/rad2deg
        psi = psi/rad2deg
        alpha = alpha/rad2deg

        !Handle outflow case
        !NOTE: For outflow, height is treated as latitude (in deg)
        call inpXML%Set_Val(doOutflow,'height/doOutflow',.false.)
        if (doOutflow) then
            hgt = hgt/rad2deg !Treating height as latitude
        endif

        do n=1,Np
            !Create particles and release if not streaming
            tpState%TPs(n)%id = id0-1+n
            call createParticle(Model,tpState%TPs(n),radius(n),phi(n),hgt(n),energy(n),alpha(n),psi(n))
            if (.not. Model%doStream) then
                !Always release
                call releaseParticle(Model,ebState,tpState%TPs(n),Model%T0)
            endif
        enddo
        if (Model%doStream) then
            call addIncoming(Model,ebState,tpState)
        endif

        !Count currently in
        tpState%NpT = count(tpState%TPs(:)%isIn)

        if (doOutflow .and. doWind) then
            write(*,*) "You can't have both wind and outflow, WTF?"
            stop
        endif
    end subroutine initParticles

    !Loop over particles, find uninitialized particles past their release date and release them
    !Uses releaseParticle defined in usertpic
    subroutine addIncoming(Model,ebState,tpState)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(tpState_T), intent(inout)   :: tpState

        integer :: n
        logical :: toDo

        !$OMP PARALLEL DO &
        !$OMP default(shared) &
        !$OMP private(n,toDo) &
        !$OMP schedule(guided)
        do n=1,tpState%Np
            toDo = (.not. tpState%TPs(n)%isInit) .and. (Model%t >= tpState%TPs(n)%T0p)
            if (toDo) then
                !Release particle
                call releaseParticle(Model,ebState,tpState%TPs(n),Model%t)
            endif
        enddo

    end subroutine addIncoming

    !Distinguish particle creation/particle release

    !Set time-independent particle quantities (others may have to wait until particle release)
    !At end of this, remaining variables to be set are
    !FO: PXFO,PYFO,PZFO
    !GC: P11GC,MUGC
    !ALL: ijk0, timestep
    subroutine createParticle(Model,prt,rad,phi,hgt,K0,alpha,psi)
        type(chmpModel_T), intent(in) :: Model
        type(prt_T), intent(inout) :: prt
        real(rp), intent(in) :: rad,phi,hgt,K0,alpha,psi

        real(rp) :: gamma

        !Set spatial position
        prt%Q(XPOS:ZPOS) = [rad*cos(phi), rad*sin(phi), hgt]
        !Note always setting initial EQ values regardless of initial z
        prt%Qeq(EQX:EQY) = prt%Q(XPOS:YPOS)

        if (doOutflow) then
            !In this case hgt is latitude (in radians at this point)
            prt%Q(XPOS:ZPOS) = [rad*cos(hgt)*cos(phi),rad*cos(hgt)*sin(phi),rad*sin(hgt)]
            prt%Qeq(EQX:EQY) = prt%Q(XPOS:YPOS)
        endif

        !Set particle energy
        !If GC particle set gamma
        !if FO store energy in EQKEV and use when particle is released
        gamma = (K0*1.0e-3)/(Model%m0*mec2) + 1.0
        if (prt%isGC) then
            prt%Q(GAMGC) = gamma
        endif
        !For now set equatorial energies (both) to K0, will reset when particle is released
        !Can't set final until know ExB
        prt%Qeq(EQKEV) = K0
        prt%Qeq(EQKEB) = K0 

        if (.not. prt%isGC) then
            prt%Q(PSITP) = psi
        endif

        !Set particle pitch angle (not a dynamic quantity)
        prt%alpha = alpha
        prt%Qeq(EQALP) = alpha

        !Initialize various quantities
        prt%ijk0 = 0
        prt%Ngc = 0
        prt%Nfo = 0
        prt%ddt = 0.0
    end subroutine createParticle

    !Release particle (finalize quantities) at a given time t
    subroutine releaseParticle(Model,ebState,prt,t)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(prt_T), intent(inout) :: prt
        real(rp), intent(in) :: t

        integer :: pSgn,OCb
        integer, dimension(NDIM) :: ijk
        logical :: inDom
        
        real(rp) :: alpha,p11,Mu,psi,ddt
        real(rp) :: gamma,ebGam
        real(rp) :: MagB,pMag
        real(rp), dimension(NDIM) :: req,r,p,pxy,E,B,vExB,xhat,yhat,bhat

        !Project to the real equator if necessary
        if (Model%doEQProj) then
            call getMagEQ(Model,ebState,prt%Q(XPOS:ZPOS),t,req,MagB,OCb)
            !If this is a closed line move TP to real equator
            if (OCb == 2) prt%Q(XPOS:ZPOS) = req
        endif
        r = prt%Q(XPOS:ZPOS)

        !Find particle
        call locate(r,ijk,Model,ebState%ebGr,inDom)
        if (.not. inDom) return
        prt%ijk0 = ijk

        !Calculate local field quantities and magnetic coordinate system
        call ebFields(r,t,Model,ebState,E,B,ijkO=ijk,vExB=vExB)
        call MagTriad(r,B,xhat,yhat,bhat)
        MagB = max(norm2(B),TINY)

        if (prt%isGC) then
            gamma = prt%Q(GAMGC)
        else
            gamma = (prt%Qeq(EQKEV)*1.0e-3)/(Model%m0*mec2) + 1.0
        endif

        !Get gamma in ExB frame
        ebGam = max(1.0,gamma - 0.5*dot_product(vExB,vExB))
        
        !Get pitch angle and set direction for p11
        alpha = prt%alpha
        if (alpha <= PI/2) then
            pSgn = +1
        else
            pSgn = -1
        endif

        if (prt%isGC) then
            !Need p11 & Mu
            pMag = Model%m0*sqrt(ebGam**2.0 - 1.0)
            p11 = pSgn*pMag*sqrt( 1 - sin(alpha)**2.0 )

            prt%Q(P11GC) = p11
            Mu = ( (Model%m0**2)*(ebGam**2.0 - 1.0) - p11*p11) / (2*Model%m0*MagB)
            prt%Q(MUGC) = Mu
            prt%Qeq(EQKEB) = (Model%m0*mec2*1.0e+3)*(ebGam-1.0)

            if (Model%do2D) then
                write(*,*) '2D Integration not implemented for GC yet'
                stop
            endif
            if (doWind) then
                write(*,*) 'Wind ICs not implemented for GC yet'
                stop
            endif
        else
            psi = prt%Q(PSITP)
            !Use either total energy = K, or in ExB = K
            pMag = Model%m0*sqrt(gamma**2.0 - 1.0)
            !Create parallel & gyro momentum components
            p11 = pSgn*pMag*sqrt( 1 - sin(alpha)**2.0 )
            pxy = pMag*sin(alpha)*( cos(psi)*xhat + sin(psi)*yhat )
            
            if (Model%doEBInit) then
                p = Model%m0*vExB
            else
                p = 0.0
            endif

            !Calculate total p (using already set frame momentum)
            p = p + p11*bhat + pxy
            if (doWind) then
                !Check if TP is beyond zWind and point it to the equator if so
                if ( (abs(r(ZPOS))>=zWind) .and. (r(ZPOS)*p(ZPOS)>=0) ) then
                    !Flip z momentum
                    p(ZPOS) = -p(ZPOS)
                endif
            endif !doWind

            prt%Q(PXFO:PZFO) = p

            !Force 2D if required
            if (Model%do2D) then
                call p2xy(prt%Q(PXFO:PZFO))
            endif
            
            !Reset certain quantities in case energy has changed (ie from going into ExB frame)
            gamma = pv2Gam(p,Model%m0) 
            prt%Qeq(EQKEV) = prt2kev(Model,prt)
            prt%Qeq(EQKEB) = p2kev(Model,p-Model%m0*vExB)
        endif

        !Calculate timestep using FO-style estimate w/ lazy over-estimate of momentum
        pMag = Model%m0*sqrt(gamma**2.0-1.0)
        p = pMag*(xhat+yhat+bhat)
        ddt = dtFO(Model,p,E,B)

        prt%ddt = min(ddt,Model%dt)

        !Particle is ready, let it go
        !Reset release time (may be off by dt)
        prt%T0p = t
        prt%Qeq(EQTIME) = prt%T0p
        prt%isIn = .true.
        prt%isInit = .true.

    end subroutine releaseParticle

end module usertpic