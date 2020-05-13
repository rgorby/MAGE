module psdutils
    use chmpdefs
    use ebtypes
    use streamline
    use psdtypes
    use chmpunits
    use ebinterp

    implicit none

    !Default density/temperatures
    real(rp) :: rho0,kT0
    real(rp) :: vc0=vc_cgs !Use CGS speed of light
    real(rp) :: MHDGamma = 5.0/3.0
    
    contains

    !Returns multi-index of TP state q(r,p,k,a) in psGr
    function psLoc(q,psGr) result(idx)
        real(rp), intent(in) :: q(NVARPS)
        type(PSEq_T), intent(in) :: psGr

        integer :: idx(NVARPS)
        logical :: inR,inP,inK,inA,inPS
        real(rp) :: r,p,K,a

        !Pull values and set signs for angles
        r = q(PSRAD  )
        p = q(PSPHI  )
        if (p<0) p=p+2*PI
        K = q(PSKINE )
        !Handle pitch angle, wrap around if necessary
        a = q(PSALPHA)
        if (a<0) a=a+PI
        !Wrap around pitch angle if necessary
        if (a>=psGr%dimBds(PSALPHA,2)) then
            !Alpha is larger than max pitch angle, wrap it back around
            !Ie, A=[0,90] instead of [0,180]
            a = PI-a
        endif

        !Check bounds
        inR = ( r >= psGr%dimBds(PSRAD   ,1) ) .and. ( r <= psGr%dimBds(PSRAD   ,2) )
        inP = ( p >= psGr%dimBds(PSPHI   ,1) ) .and. ( p <= psGr%dimBds(PSPHI   ,2) )
        inA = ( a >= psGr%dimBds(PSALPHA ,1) ) .and. ( a <= psGr%dimBds(PSALPHA ,2) )
        inK = ( K >= psGr%dimBds(PSKINE  ,1) ) .and. ( K <= psGr%dimBds(PSKINE  ,2) )
        inPS = inR .and. inP .and. inK .and. inA
        if (.not. inPS) then
            idx = 0
            return
        endif

        !Now we know we're in this PS grid
        idx(PSRAD  ) = maxloc(psGr%rI,dim=1,mask=r>=psGr%rI)
        idx(PSPHI  ) = maxloc(psGr%pI,dim=1,mask=p>=psGr%pI)
        idx(PSKINE ) = maxloc(psGr%kI,dim=1,mask=K>=psGr%kI)
        idx(PSALPHA) = maxloc(psGr%aI,dim=1,mask=a>=psGr%aI)
        idx(PSALPHA) = max(1,idx(PSALPHA))
        idx(PSALPHA) = min(psGr%Na,idx(PSALPHA))

    end function psLoc

    function psMagP(Model,psGr,idx) result(pMag)
        type(chmpModel_T), intent(in) :: Model
        type(PSEq_T), intent(in) :: psGr
        integer, intent(in) :: idx(NVARPS)
        real(rp) :: pMag
        real(rp) :: e0,W
        !P = sqrt(2mK  + (K/c)**2.0), momentum magnitude
        !  = sqrt(2K*mc^2 + K*K)/c

        e0 = Model%m0*mec2*1.0e+3 !Particle rest mass [keV]
        W = psGr%kC(idx(PSKINE))
        pMag = sqrt(2*e0*W + W*W)/vc0

    end function psMagP

    !Calculate metric term of a specific variable at PS cell idx
    function dGamma(Model,psGr,idx,dGVar) result(dl)
        type(chmpModel_T), intent(in) :: Model
        type(PSEq_T), intent(in) :: psGr
        integer, intent(in) :: dGVar,idx(NVARPS)

        real(rp) :: dl,W,e0,pMag
        real(rp) :: xScl,kScl,aScl
        integer :: ir,ip,ik,ia

        xScl = L0

        ir = idx(PSRAD)
        ip = idx(PSPHI)
        ik = idx(PSKINE)
        ia = idx(PSALPHA)

        W = psGr%kC(ik)
        e0 = Model%m0*mec2*1.0e+3 !Particle rest mass [keV]
        pMag = psMagP(Model,psGr,idx)
        
        select case(dGVar)
        !--------
        case(PSRAD)
            !dr
            dl = xScl*psGr%rD(ir)
        !--------
        case(PSPHI)
            !r*dphi
            dl = xScl*psGr%rC(ir)*psGr%pD(ip)
        !--------
        case(PSALPHA)
            !pMag*dA
            !Add scaling factor to account for a=[0,90] case
            aScl = PI/(psGr%dimBds(PSALPHA,2)-psGr%dimBds(PSALPHA,1))
            dl = aScl*pMag*psGr%aD(ia)
        !--------
        case(PSKINE)
            kScl = (W+e0)/(vc0*vc0*pMag)
            dl = kScl*psGr%kD(ik)
        !--------
        case(PSPSI)
            !Assuming full 2pi contribution
            !pMag*sin(alpha)*dPsi
            dl = pMag*sin(psGr%aC(ia))*2*PI
        !--------
        case default
            dl = 0.0

        end select
    end function dGamma

    !Update phase space (recalculate volume elements)
    subroutine updatePS(Model,psGr,ebState,t)
        type(chmpModel_T), intent(in) :: Model
        type(PSEq_T), intent(inout) :: psGr
        type(ebState_T), intent(in)  :: ebState
        real(rp), intent(in) :: t

        integer :: i,j,k,n,idx(NVARPS)
        real(rp) :: dGx,dGp,dk,da,ds,Vr
        real(rp) :: R,kT,x0(NDIM),Qmhd(NVARMHD)

        !Calculate flux tube volumes at each polar point
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(R,x0,Qmhd,kT,Vr)
        do i=1,psGr%Nr
            do j=1,psGr%Np
                !Get seed point, clean old stream, get new
                R = psGr%rC(i)
                x0(XDIR) = R*cos(psGr%pC(j))
                x0(YDIR) = R*sin(psGr%pC(j))
                x0(ZDIR) = 0.0

                call cleanStream(psGr%bLns(i,j))
                !Trace field line
                call genStream(Model,ebState,x0,t,psGr%bLns(i,j))

                !Turn field line into flux tube volume afa alpha
                call TubedV(Model,psGr,psGr%bLns(i,j),i,j)

                if (Model%doMHD) then
                    !Do flux-tube volume averaged quantities if possible
                    if (psGr%isClosed(i,j)) then
                        call dVFlow(Model,psGr,psGr%bLns(i,j),i,j,Qmhd,kT,Vr)
                    else
                        !Do MHD interpolation at equatorial cell center
                        call eqFlow(Model,ebState,x0,t,Qmhd,kT,Vr)
                    endif
                else
                    !Use default values
                    call eqFlow(Model,ebState,x0,t,Qmhd,kT,Vr)
                endif

                !----
                !Save "equatorial" values

                psGr%Qrp (i,j,:) = Qmhd
                psGr%kTeq(i,j)   = kT
                psGr%Vreq(i,j)   = Vr

            enddo
        enddo

        psGr%time = t
        !Now calculate full PS volume element
        !dGx ~ L0^3
        !dGp ~ keV^3/c^3

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,n,idx,dGx,dGp,dk,da,ds)
        do n=1,psGr%Na
            do k=1,psGr%Nk
                do j=1,psGr%Np
                    do i=1,psGr%Nr
                        idx = [i,j,k,n]
                        !Get flux tube volume [L0^3 -> cm3]
                        dGx = (L0**3.0)*psGr%dVb(i,j,n)

                        !Get volume element contributions from K,alpha,psi
                        dk = dGamma(Model,psGr,idx,PSKINE)
                        da = dGamma(Model,psGr,idx,PSALPHA)
                        ds = dGamma(Model,psGr,idx,PSPSI)
                        dGp = dk*da*ds

                        psGr%dG(i,j,k,n) = dGx*dGp

                    enddo
                enddo
            enddo
        enddo

    end subroutine updatePS

    !Get density/temperature at a cell (can also scale by flow speed)
    subroutine eqFlow(Model,ebState,x0,t,Qmhd,kT,Vr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)  :: ebState
        real(rp), intent(in) :: x0(NDIM), t
        real(rp), intent(out) :: Vr,kT,Qmhd(NVARMHD)

        real(rp) :: Cs,P,rho
        real(rp) :: x,y,Vx,Vy,r

        x = x0(XDIR)
        y = x0(YDIR)
        r = sqrt(x**2.0+y**2.0)

        if (Model%doMHD) then
            Qmhd = mhdInterp(x0,t,Model,ebState)
            rho = Qmhd(DEN)
            P = Qmhd(PRESSURE)

            !Qmhd units are [#/cc] and nPa
            ! (P/rho)*(1.0e-9)*(1.0e-2)**3 ,kT in J
            ! = (P/rho)*1.0e-15 [J]
            ! = (P/rho)(1.0e+16)*(1.0e-15)/1.6 [keV]
            kT = (P/rho)*(10.0/1.6)
            Vx = Qmhd(VELX)
            Vy = Qmhd(VELY)
        else
            rho = rho0
            Qmhd(DEN) = rho0
            Qmhd(VELX:VELZ) = 0.0
            Qmhd(PRESSURE) = 0.0
            kT  = kT0
            Vx = 0.0
            Vy = 0.0
        endif
        Vr = (x*Vx+y*Vy)/r
        
    end subroutine eqFlow

    subroutine dVFlow(Model,psGr,bTrc,i,j,Qmhd,kT,Vr)
        type(chmpModel_T), intent(in) :: Model
        type(PSEq_T), intent(in) :: psGr
        type(fLine_T), intent(in)  :: bTrc
        integer, intent(in) :: i,j
        real(rp), intent(out) :: Vr,kT,Qmhd(NVARMHD)

        real(rp) :: R,dR,dp,dA
        integer :: Nc,n,k
        real(rp), dimension(:), allocatable :: bAvg,dl,eD,eV,eP,dV
        real(rp), dimension(NDIM) :: Vxyz
        real(rp) :: VrP,VrM,dV0
        Qmhd = 0.0
        kT = 0.0
        Vr = 0.0

        !Get local geometry
        R = psGr%rC(i)
        dR = psGr%rI(i+1)-psGr%rI(i)
        dp = psGr%pI(j+1)-psGr%pI(j)
        dA = R*dR*dp

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)

        !Recenter to edges
        Nc = Nm+Np+1-1 !Centers
        allocate(bAvg(Nc)) !Edge field strength
        allocate(dl  (Nc)) !Edge length
        allocate(eD  (Nc)) !Edge-centered density
        allocate(eV  (Nc)) !Edge-centered Vr
        allocate(eP  (Nc)) !Edge-centered pressure
        allocate(dV  (Nc)) !Volume element
        n = 1
        do k=-Nm,Np-1
            dl(n)   = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            bAvg(n) = 0.5*(bTrc%lnVars(0)%V(k+1) + bTrc%lnVars(0)%V(k))
            eD(n)   = 0.5*(bTrc%lnVars(DEN)%V(k+1) + bTrc%lnVars(DEN)%V(k))
            eP(n)   = 0.5*(bTrc%lnVars(PRESSURE)%V(k+1) + bTrc%lnVars(PRESSURE)%V(k))
            
            Vxyz = [bTrc%lnVars(VELX)%V(k+1),bTrc%lnVars(VELY)%V(k+1),bTrc%lnVars(VELZ)%V(k+1)]
            VrP = dot_product(bTrc%xyz(k+1,:),Vxyz)/norm2(bTrc%xyz(k+1,:))
            Vxyz = [bTrc%lnVars(VELX)%V(k  ),bTrc%lnVars(VELY)%V(k  ),bTrc%lnVars(VELZ)%V(k  )]
            VrM = dot_product(bTrc%xyz(k  ,:),Vxyz)/norm2(bTrc%xyz(k  ,:))

            !VrP = dot_product(bTrc%xyz(k+1,:),bTrc%lnVars(VELX:VELZ)%V(k+1))/norm2(bTrc%xyz(k+1,:))
            !VrM = dot_product(bTrc%xyz(k  ,:),bTrc%lnVars(VELX:VELZ)%V(k  ))/norm2(bTrc%xyz(k  ,:))
            eV(n)   = 0.5*(VrP+VrM)
            n = n+1
        enddo

        dV = dA*minval(bAvg)*dl/bAvg
        dV0 = sum(dV)
        Qmhd(DEN)      = sum(eD*dV)/dV0
        Qmhd(PRESSURE) = sum(eP*dV)/dV0
        Qmhd(VELX) = bTrc%lnVars(VELX)%V(0)
        Qmhd(VELY) = bTrc%lnVars(VELY)%V(0)
        Qmhd(VELZ) = bTrc%lnVars(VELZ)%V(0)

        Vr = sum(eV*dV)/dV0

        !Qmhd units are [#/cc] and nPa
        ! (P/rho)*(1.0e-9)*(1.0e-2)**3 ,kT in J
        ! = (P/rho)*1.0e-15 [J]
        ! = (P/rho)(1.0e+16)*(1.0e-15)/1.6 [keV]
        kT = (Qmhd(PRESSURE)/Qmhd(DEN))*(10.0/1.6)

        end associate
    end subroutine dVFlow

    !Calculate flux tube @ i,j
    !Flux tube volume is in units of L0^3 (see chmpunits.F90)
    subroutine TubedV(Model,psGr,bTrc,i,j)
        type(chmpModel_T), intent(in) :: Model
        type(PSEq_T), intent(inout) :: psGr
        type(fLine_T), intent(in)  :: bTrc
        integer, intent(in) :: i,j

        real(rp) :: R,dr,dp,dA,Beq,bMin,aeq
        real(rp), dimension(:), allocatable :: dl,bAvg,bScl,aMir
        integer :: i0,Nc,n,k
        real(rp) :: rP,rM

        !Get local geometry
        R = psGr%rC(i)
        dR = psGr%rI(i+1)-psGr%rI(i)
        dp = psGr%pI(j+1)-psGr%pI(j)
        dA = R*dR*dp

        !Get field line quantities
        Beq = bTrc%lnVars(0)%V0

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
        !Find if this field line is closed
        rP = norm2(bTrc%xyz( Np,:))
        rM = norm2(bTrc%xyz(-Nm,:))

        if ( (rP <= psGr%rCrit) .and. (rM <= psGr%rCrit) ) then
            psGr%isClosed(i,j) = .true.
        else
            psGr%isClosed(i,j) = .false.
            !Don't bother with field line calculations
            psGr%dVb(i,j,:) = 0.0
            return
        endif

        !Recenter from nodes to centers
        !Track dl,bscl,aeq (which equatorial pitch angle mirrors here)
        !Working in degrees for mirror point calculation
        Nc = Nm+Np+1-1 !Centers
        allocate(dl  (Nc))
        allocate(bAvg(Nc))
        allocate(bScl(Nc))
        allocate(aMir(Nc))
        !Start with recentering
        n = 1
        do k=-Nm,Np-1
            dl(n)   = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            bAvg(n) = 0.5*(bTrc%lnVars(0)%V(k+1) + bTrc%lnVars(0)%V(k))
            n = n+1
        enddo

        !Calculate bmin from recentered values
        bMin = minval(bAvg)
        bScl = bMin/bAvg
        aMir = rad2deg*asin(sqrt(bScl))

        end associate

        !Now loop through alpha and calculate accessible flux tube volume
        !Assuming symmetric wrt alpha

        !Loop over cell center pitch angles
        do k=1,psGr%Na
            aeq = rad2deg*psGr%aC(k)
            !Map to 0,90
            if (aeq>=90) aeq = 180-aeq
            !Sum accessible flux tube volume
            if ( count(mask=aMir>=aeq)>0 ) then
                psGr%dVb(i,j,k) = dA*sum(bScl*dl,mask=aMir>=aeq)
            else
                psGr%dVb(i,j,k) = 0.0
            endif
        enddo

    end subroutine TubedV

    !Calculate diagonal pressure tensor components
    !Uses psGr (r,p) and psPop fPSD
    !Pxyz = |Nr,Np,NDIM|
    !Note, this calculation is explicitly non-relativistic

    subroutine CalcMoms(Model,psGr,psPop,Nk,Vxyz,Pxyz)
        type(chmpModel_T), intent(in) :: Model
        type(PSEq_T), intent(in) :: psGr
        type(psdPop_T), intent(in) :: psPop
        real(rp), intent(out) :: Nk(psGr%Nr,psGr%Np)
        real(rp), intent(out) :: Vxyz(psGr%Nr,psGr%Np,NDIM)
        real(rp), intent(out) :: Pxyz(psGr%Nr,psGr%Np,NDIM)

        integer :: Ns,ir,ip,ia,ik,is
        real(rp) :: PScl,e0,aScl
        !Inside loop variables
        real(rp) :: fC,dk,da,ds,K,psi,sA,cA,sS,cS,dG
        real(rp) :: Vb,vMag,Vx,Vy,Vz
        real(rp), dimension(NDIM) :: dvdv

        !Scaling factor for pressure
        !eb units (kev/cm3) -> nPa
        PScl = 1.60218e-1

        e0 = Model%m0*mec2*1.0e+3 !Particle rest mass [keV]
        !Add scaling factor to account for a=[0,90] case
        aScl = PI/(psGr%dimBds(PSALPHA,2)-psGr%dimBds(PSALPHA,1))

        !Lazily adding PSI binning here to deal with EB direction
        Ns = 16

        !Scaling factor for pressure
        !eb units (kev/cm3) -> nPa
        PScl = 1.60218e-1

        e0 = Model%m0*mec2*1.0e+3 !Particle rest mass [keV]
        !Add scaling factor to account for a=[0,90] case
        aScl = PI/(psGr%dimBds(PSALPHA,2)-psGr%dimBds(PSALPHA,1))

        !Lazily adding PSI binning here to deal with EB direction
        Ns = 16

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(is,ik,ia,ir,ip) &
        !$OMP private(dvdv,Vb,vMag,Vx,Vy,Vz) &
        !$OMP private(fC,dk,da,ds,K,psi,sA,cA,sS,cS,dG)
        do ip=1,psGr%Np
            do ir=1,psGr%Nr
                Nk  (ir,ip  ) = 0.0
                Vxyz(ir,ip,:) = 0.0
                Pxyz(ir,ip,:) = 0.0

                !Vb is MHD flow (assumed in x dir)
                Vb = 0 !Already subtracted by particle pusher
                !Now loop over K,alpha & psi* to integrate moment
                do ia=1,psGr%Na
                    do ik=1,psGr%Nk
                        do is=1,Ns

                            !Get various pieces
                            fC = psPop%fPSD(ir,ip,ik,ia) !Units of (keV*s)^-3
                            dk = psGr%kD(ik)
                            da = aScl*psGr%aD(ia)
                            ds = 2*PI/Ns
                            

                            K = psGr%kC(ik)
                            psi = 0.0 + (ds/2) + (is-1)*ds !Lazy psi bin centers

                            sA = sin(psGr%aC(ia))
                            cA = cos(psGr%aC(ia))
                            sS = sin(psi)
                            cS = cos(psi)

                            !Get velocities
                            !Vx,y,z = integration v
                            vMag = sqrt(2*K/e0)
                            Vx = vMag*sA*cS
                            Vy = vMag*sA*sS
                            Vz = vMag*cA
                            dvdv(XDIR) = (Vx-Vb)**2.0
                            dvdv(YDIR) = (Vy)**2.0
                            dvdv(ZDIR) = (Vz)**2.0
                            dG = (e0*sqrt(2*e0*K)*sA*dk*da*ds)/(vc_cgs**3.0)
                            
                            !Density
                            Nk(ir,ip) = Nk(ir,ip) + fC*dG
                            !Pressure
                            Pxyz(ir,ip,:) = Pxyz(ir,ip,:) + e0*dvdv(:)*fC*dG
                            !Velocity
                            Vxyz(ir,ip,:) = Vxyz(ir,ip,:) + 0.0

                        enddo
                    enddo
                enddo

                !End accumulation into pressure
                !At this point, P is in units of keV/cm3
                !Want nPa (Pa = J/m3)
                Pxyz(ir,ip,:) = PScl*Pxyz(ir,ip,:)

            enddo
        enddo

    end subroutine CalcMoms

end module psdutils
