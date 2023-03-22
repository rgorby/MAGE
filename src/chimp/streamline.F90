module streamline
    use chmpdefs
    use ebtypes
    use math
    use ebinterp
    use planethelper
    use streamutils
    use imaghelper

    implicit none

    real(rp), parameter, private :: ShueScl = 1.25 !Safety factor for Shue MP
    real(rp), parameter, private :: rShue   = 6.0  !Radius to start checking Shue
    integer , parameter, private :: NpChk   = 10   !Cadence for Shue checking

    contains

    !doNHO = T, assume doing RCM coupling
    subroutine genStream(Model,ebState,x0,t,fL,MaxStepsO,doShueO,doNHO)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(fLine_T), intent(inout) :: fL
        integer , intent(in), optional :: MaxStepsO
        logical , intent(in), optional :: doShueO,doNHO

        integer :: N1,N2,i,Np,Nm,n,MaxN
        real(rp) :: dx(NDIM)
        real(rp) :: Xn(0:MaxFL,NDIM,2),Vn(0:MaxFL,0:NumVFL,2)
        integer :: ijkn(0:MaxFL,NDIM,2)
        logical :: inDom,doShue,doNH

        !Start by emptying line
        call cleanStream(fL)
        fL%x0 = x0
        inDom = inDomain(x0,Model,ebState%ebGr)
        
        if (.not. inDom) then
            !Bad tube, seed point isn't in domain
            fL%isGood = .false.
            allocate(fL%xyz(0:0,NDIM))
            fL%xyz(0,:) = x0
            return
        endif
        
        !Get traced line
        fl%lnVars(0)%idStr = "B"

        if (Model%doMHD) then
            !Create holders for all MHD variables
            fl%lnVars(DEN)%idStr = "D"
            fl%lnVars(VELX)%idStr = "Vx"
            fl%lnVars(VELY)%idStr = "Vy"
            fl%lnVars(VELZ)%idStr = "Vz"
            fl%lnVars(PRESSURE)%idStr = "P"
        endif

        if (present(MaxStepsO)) then
            MaxN = MaxStepsO
        else
            MaxN = MaxFL
        endif
        if ((present(doNHO))) then
            doNH = doNHO
        else
            doNH = .false.
        endif

        if (present(doShueO)) then
            doShue = doShueO
        else
            doShue = .false.
        endif

        call genTrace(Model,ebState,x0,t,Xn(:,:,1),ijkn(:,:,1),Vn(:,:,1),N1,-1,MaxN,doShue)
        if (doNH) then
            N2 = 0
        else
            call genTrace(Model,ebState,x0,t,Xn(:,:,2),ijkn(:,:,2),Vn(:,:,2),N2,+1,MaxN,doShue)
        endif

        !Create field line
        fL%isGood = .true.
        fL%Nm = N1
        fL%Np = N2

        !Do allocations/set seed point value
        allocate(fL%xyz(-N1:N2,NDIM))
        allocate(fL%ijk(-N1:N2,NDIM))
        fL%xyz(0,:) = x0
        fL%ijk(0,:) = ijkn(0,1,1)
        do n=0,NumVFL
            allocate(fL%lnVars(n)%V(-N1:N2))
            !fL%lnVars(n)%V(0) = Vn(n,1,1)
            fL%lnVars(n)%V(0) = Vn(0,n,1)
            fL%lnVars(n)%V0   = Vn(0,n,1)
        enddo
        
        !Now load field line, positive dir then negative
        do i=1,N1 !Negative
            fL%xyz(-i,:)       = Xn(i,:,1)
            fL%ijk(-i,:)     = ijkn(i,:,1)
            !fL%lnVars(:)%V(-i) = Vn(i,:,1)
            do n=0,NumVFL
                fL%lnVars(n)%V(-i) = Vn(i,n,1)
            enddo
        enddo
        
        do i=1,N2 !Positive
            fL%xyz(i,:)       = Xn(i,:,2)
            fL%ijk(i,:)     = ijkn(i,:,2)
            !fL%lnVars(:)%V(i) = Vn(i,:,2)
            do n=0,NumVFL
                fL%lnVars(n)%V(i) = Vn(i,n,2)
            enddo

        enddo
        
    end subroutine genStream

    !!Gathers field line topology information
    subroutine SliceFL(Model,ebState,x0,t,ebTrc)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(ebTrc_T), intent(inout) :: ebTrc

        type(fLine_T) :: bTrc
        !Initialize the values
        ebTrc%OCb = 0.0
        ebTrc%dvB = 0.0
        ebTrc%bD  = 0.0
        ebTrc%bP  = 0.0
        ebTrc%bS  = 0.0
        ebTrc%bMin = 0.0

        ebTrc%MagEQ(:) = 0.0
        ebTrc%xEPm (:) = 0.0
        ebTrc%xEPp (:) = 0.0

        if (.not. inDomain(x0,Model,ebState%ebGr)) return
        !Trace field line
        call genStream(Model,ebState,x0,t,bTrc)

        !Get diagnostics
        ebTrc%OCb = 1.0*FLTop(Model,ebState%ebGr,bTrc)
        if (ebTrc%OCb > 0) then
            !Get flux-tube integrals
            call FLThermo(Model,ebState%ebGr,bTrc,ebTrc%bD,ebTrc%bP,ebTrc%dvB)
            ebTrc%bS   = FLEntropy(Model,ebState%ebGr,bTrc)

            !Get magnetic equator info
            call FLEq(Model,bTrc,ebTrc%MagEQ,ebTrc%bMin)

            !Get endpoints info
            associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
            ebTrc%xEPm = bTrc%xyz(-Nm,:)
            ebTrc%xEPp = bTrc%xyz(+Np,:)
            end associate

        endif

        !write(*,*) 'FL size = ', bTrc%Nm+bTrc%Np+1
    end subroutine SliceFL

!---------------------------------
!Field line diagnostics
    !Flux tube volume
    function FLVol(Model,ebGr,bTrc) result(dvB)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: bTrc
        real(rp) :: dvB

        integer :: OCb,k
        real(rp) :: Phi0,dl,bAvg,dA

        Phi0 = 1.0 !Fixed flux, could be B0*A0
        dvB = 0.0

        OCb = FLTop(Model,ebGr,bTrc)
        if (OCb<2) return
        
        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
        do k=-Nm,Np-1
            dl = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            bAvg = 0.5*(bTrc%lnVars(0)%V(k+1) + bTrc%lnVars(0)%V(k))
            dA = Phi0/bAvg
            dvB = dvB + dA*dl
        enddo
        end associate
    end function FLVol

    !Calculate arc length of field line
    function FLArc(Model,ebGr,bTrc) result(L)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: bTrc
        real(rp) :: L
        real(rp) :: dL
        integer :: k

        L = 0.0
        if (.not. bTrc%isGood) return

        do k=-bTrc%Nm,bTrc%Np-1
            dL = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            L = L + dL
        enddo

    end function FLArc

    !Calculate Alfven crossing time on line
    function FLAlfvenX(Model,ebGr,bTrc) result(dtX)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: bTrc
        real(rp) :: dtX

        integer :: k
        real(rp) :: dL,eD,bMag,Va

        dtX = 0.0
        if (.not. bTrc%isGood) return
        do k=-bTrc%Nm,bTrc%Np-1
            dL = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            dL = dL*L0*1.0e-5 !Corner units to km
            !Get egde-centered quantities
            eD = 0.5*(bTrc%lnVars(DEN)%V(k+1) + bTrc%lnVars(DEN)%V(k))
            bMag = 0.5*(bTrc%lnVars(0)%V(k+1) + bTrc%lnVars(0)%V(k))
            !Convert B to nT, eD in #/cc
            bMag = oBScl*bMag
            if (eD > TINY) then
                Va = 22.0*bMag/sqrt(eD) !Alfven speed in km/s, NRL formulary
                dtX = dtX + dL/Va
            endif
        enddo

    end function FLAlfvenX

    !Calculate Alfven+Sound crossing time on line
    function FLFastX(Model,ebGr,bTrc) result(dtX)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: bTrc
        real(rp) :: dtX

        integer :: k
        real(rp) :: dL,eD,eP,TiEV,bMag,Va,Cs

        dtX = 0.0
        if (.not. bTrc%isGood) return
        do k=-bTrc%Nm,bTrc%Np-1
            dL = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            dL = dL*L0*1.0e-5 !Corner units to km
            !Get egde-centered quantities
            eD   = 0.5*(bTrc%lnVars(DEN)     %V(k+1) + bTrc%lnVars(DEN)     %V(k))
            eP   = 0.5*(bTrc%lnVars(PRESSURE)%V(k+1) + bTrc%lnVars(PRESSURE)%V(k))
            bMag = 0.5*(bTrc%lnVars(0)       %V(k+1) + bTrc%lnVars(0)       %V(k))

            !Convert B to nT, eD in #/cc, eP in nPa
            bMag = oBScl*bMag
            if (eD > TINY) then
                Va = 22.0*bMag/sqrt(eD) !Alfven speed in km/s, NRL formulary
                !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
                TiEV = (1.0e+3)*DP2kT(eD,eP) !Temp in eV
                Cs = 9.79*sqrt((5.0/3)*TiEV)
                dtX = dtX + dL/(Va+Cs)
            endif
        enddo

    end function FLFastX

    !Averaged density/pressure
    subroutine FLThermo(Model,ebGr,bTrc,bD,bP,dvB,bBetaO)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: bTrc
        real(rp), intent(out) :: bD,bP,dvB
        real(rp), intent(out), optional :: bBetaO

        integer :: k
        real(rp) :: bMag,dl,eP,eD,ePb !Edge-centered values
        real(rp) :: bPb,bBeta
        
        if (.not. bTrc%isGood) return

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
        !Zero out accumulators
        bD = 0.0
        bP = 0.0
        dvB = 0.0
        bPb = 0.0
        bBeta = 0.0

        !Loop over edges
        do k=-Nm,Np-1
            !Get edge-centered quantities
            dl = norm2(bTrc%xyz(k+1,:) - bTrc%xyz(k,:)) !Edge length
            bMag = 0.5*(bTrc%lnVars(0)%V(k+1) + bTrc%lnVars(0)%V(k))

            eD = 0.5*(bTrc%lnVars(DEN)%V(k+1) + bTrc%lnVars(DEN)%V(k))
            eP = 0.5*(bTrc%lnVars(PRESSURE)%V(k+1) + bTrc%lnVars(PRESSURE)%V(k))

            !Get edge mag pressure, bmag=>nT(oBScl)=>T
            !(NRL Plasma formulary):
            !3.98x10^6 * (B/B0)^2 = Pb [dynes/cm2] = 0.1 Pa, x10^8 0.1 Pa => nPa
            !ePb = (1.0e+8)*(3.98*1.0e+6)*(bMag*oBScl*1.0e-9)**2.0
            ePb = 1.0e+14*(bMag*oBScl*1.0e-9/0.501)**2.0 !Edge mag pressure in nPa

            !Now accumulate into flux-tube integrals
            dvB = dvB +     dl/bMag
            bD  = bD  +  eD*dl/bMag
            bP  = bP  +  eP*dl/bMag
            bPb = bPb + ePb*dl/bMag
            bBeta = bBeta + (eP/ePb)*dl/bMag
        enddo

        !Now turn flux-tube integrals of quantities into flux-tube averages
        bD  = bD/dvB
        bP  = bP/dvB
        bPb = bPb/dvB

        !bBeta = bP/bPb
        bBeta = bBeta/dvB


        if (present(bBetaO)) then
            bBetaO = bBeta
        endif

        end associate
    end subroutine FLThermo
    
    !Flux tube entropy
    function FLEntropy(Model,ebGr,bTrc,GamO) result(S)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: bTrc
        real(rp), optional :: GamO

        real(rp) :: S

        integer :: OCb,k
        real(rp) :: Phi0,dl,bAvg,pAvg,dV,Gam
        if (present(GamO)) then
            Gam = GamO
        else
            Gam = 5.0/3.0
        endif

        Phi0 = 1.0
        S = 0.0

        OCb = FLTop(Model,ebGr,bTrc)
        if ( OCb<2 .or. (.not. Model%doMHD) ) return

        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
        do k=-Nm,Np-1
            dl = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            bAvg = 0.5*(bTrc%lnVars(0)%V(k+1) + bTrc%lnVars(0)%V(k))
            dV = dl*Phi0/bAvg
            pAvg = 0.5*(bTrc%lnVars(PRESSURE)%V(k+1) + bTrc%lnVars(PRESSURE)%V(k))
            if (pAvg > TINY) then
                S = S + dV*(pAvg**(1.0/Gam))
            endif
        enddo

        end associate

    end function FLEntropy

    function FLTop(Model,ebGr,bTrc) result(OCb)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: bTrc
        integer :: OCb

        logical :: isCP,isCM,isFin,isStart
        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)

        !Test topology
        !OCB =  0 (solar wind), 1 (half-closed), 2 (both ends closed)
        !OCB = -1 (plasmoid/timeout)

        if (.not. bTrc%isGood) then
            OCb = -1
            return
        endif

        isCP  = isClosed(bTrc%xyz(+Np,:),Model)
        isCM  = isClosed(bTrc%xyz(-Nm,:),Model)

        
        isFin = (Np<MaxFL-1) .and. (Nm<MaxFL-1) !Check if finished
        isStart = (Np>0) .and. (Nm>0) !Check if both sides went somewhere

        OCb = 0

        if ( isFin ) then
        !Both sides ended normally
            if (isCP .or. isCM) then
                !At least one side is closed
                if (isCP .and. isCM) then
                    OCb = 2
                else
                    OCb = 1
                endif
            else
                !Neither side closed and both sides ended normally
                !IMF
                OCb = 0
            endif
        else
        !At least one side ended badly
            OCb = -1 !Just set timeout flag
            !if ( (Np<MaxFL-1) ) OCb = OCb-1
            !if ( (Nm<MaxFL-1) ) OCb = OCb-1
        endif
        end associate
    end function FLTop

    !Get conjugate point from field line, assuming looking for southern hemisphere
    subroutine FLConj(Model,ebGr,bTrc,xyzC)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in) :: ebGr
        type(fLine_T), intent(in) :: bTrc
        real(rp), intent(out) :: xyzC(NDIM)

        real(rp), dimension(NDIM) :: xP,xM

        if (.not. bTrc%isGood) then
            xyzC = 0.0
            return
        endif
        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)
            
        !Get endpoints of field line
        xP = bTrc%xyz(+Np,:)
        xM = bTrc%xyz(-Nm,:)

        if (xP(ZDIR)<0.0) then
            xyzC = xP
        else if (xM(ZDIR)<0.0) then
            xyzC = xM
        else
            xyzC = 0.0
        endif

        end associate

    end subroutine FLConj

    !Get minimum field strength and location
    subroutine FLEq(Model,bTrc,xeq,Beq)
        type(chmpModel_T), intent(in) :: Model
        type(fLine_T), intent(in) :: bTrc
        real(rp), intent(out) :: xeq(NDIM),Beq

        integer :: i0,iMin,OCb
        associate(Np=>bTrc%Np,Nm=>bTrc%Nm)

        if (.not. bTrc%isGood) then
            !Bad field line
            xeq = 0.0
            Beq = 0.0
            return
        endif

        !Find minimum field and where it occurs
        Beq = minval(bTrc%lnVars(0)%V) !Min field strength
        i0  = minloc(bTrc%lnVars(0)%V,dim=1) !Note this is between 1:N
        !Need to offset i0, iMin = i0-Nm-1
        iMin = i0-Nm-1
        xeq = bTrc%xyz(iMin,:)

        end associate
    end subroutine FLEq

    !Get curvature radius at equator and ExB velocity [km/s]
    subroutine FLCurvRadius(Model,ebGr,ebState,bTrc,rCurv,vEB)
        type(chmpModel_T), intent(in)  :: Model
        type(ebGrid_T)   , intent(in)  :: ebGr
        type(ebState_T)  , intent(in)  :: ebState
        type(fLine_T)    , intent(in)  :: bTrc
        real(rp)         , intent(out) :: rCurv
        real(rp)         , intent(out) :: vEB

        type(gcFields_T) :: gcFields
        real(rp), dimension(NDIM,NDIM) :: Jacbhat
        real(rp), dimension(NDIM) :: xeq,E,B,vExB
        real(rp) :: t,Beq

        rCurv = 0.0
        if (.not. bTrc%isGood) return

        !Get equator
        call FLEq(Model,bTrc,xeq,Beq)

        !Now get field information there
        t = ebState%eb1%time
        call ebFields(xeq,t,Model,ebState,E,B,vExB=vExB,gcFields=gcFields)

        rCurv = getRCurv(B,gcFields%JacB)
        vEB = norm2(vExB)*oVScl
    end subroutine FLCurvRadius
!---------------------------------
!Projection routines
    !Project to SM EQ (Z=0)
    subroutine getEquatorProjection(Model,ebState,x0,t,xe)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(out) :: xe(NDIM) ! end point
        logical :: failEquator

        if (.not. inDomain(x0,Model,ebState%ebGr)) then
           xe = HUGE
           return
        endif
        
        if (x0(ZDIR)>=TINY) then
           ! assume first that northern hemisphere is always traced in -B direction
           call project(Model,ebState,x0,t,xe,-1,.true.,failEquator)
        else if (x0(ZDIR)<=-TINY) then 
           call project(Model,ebState,x0,t,xe,+1,.true.,failEquator)
        else
           ! if we're for some reason at equator don't do anything.
           xe = x0
           return
        endif

        ! at this point, we have either gotten to equator
        ! or to the domain boundary
        ! because we trace along -B for z>0 and along B for z<0
        ! can get to weird places on non-closed field lines or 
        ! even on closed if strongly tilted
        ! trap for those points here 
        if (failEquator) xe = HUGE

      end subroutine getEquatorProjection

    !Project to magnetic equator (min along field line)
    subroutine getMagEQ(Model,ebState,x0,t,xeq,Beq,tOpt)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(out) :: xeq(NDIM),Beq
        integer, intent(out), optional :: tOpt

        type(fLine_T) :: bTrc
        integer :: OCb
        
        !Initialize output
        xeq = 0.0
        Beq = 0.0
        if (present(tOpt)) tOpt = -1

        if (.not. inDomain(x0,Model,ebState%ebGr)) return

        !Trace field line
        call genStream(Model,ebState,x0,t,bTrc)

        !Get diagnostics
        call FLEq(Model,bTrc,xeq,Beq)
        OCb = FLTop(Model,ebState%ebGr,bTrc)

        if (present(tOpt)) then
            tOpt = OCb
        endif

    end subroutine getMagEQ

!---------------------------------
    !Project XYZ to lat-lon on ionosphere
    !TODO: Better merge this with voltron RCM code, too redundant now

    subroutine Map2NH(ebModel,ebState,xyz,t,x1,x2)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebState_T)  , intent(in) :: ebState
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(in) :: t
        real(rp), intent(out) :: x1,x2

        real(rp), dimension(NDIM) :: xE,xIon,xyz0
        real(rp) :: dX,rC,startEps,rEps
        logical :: isGood

        startEps = 0.05
        rEps = 0.125

        x1 = 0.0
        x2 = 0.0

        ! trap for when we're within epsilon of the inner boundary
        ! (really, it's probably only the first shell of nodes at R=Rinner_boundary that doesn't trace correctly)
        if ( (norm2(xyz)-rClosed)/rClosed < startEps ) then
           ! dipole-shift to startEps
           xyz0 = DipoleShift(xyz,norm2(xyz)+startEps)
        else
           xyz0 = xyz
        end if

        !Use one-sided projection routine from chimp
        !Trace along field line (i.e. to northern hemisphere)
        call project(ebModel,ebState,xyz0,t,xE,+1,toEquator=.false.)

        dX = norm2(xyz0-xE)
        rC = rClosed*(1.+rEps)
        isGood = (dX>TINY) .and. (norm2(xE) <= rC) .and. (xE(ZDIR) > 0)

        if (isGood) then
            !Get invariant lat/lon
            x1 = InvLatitude(xE)
            x2 = katan2(xE(YDIR),xE(XDIR))
            if (x2 < 0) x2 = x2 + 2*PI
        else
            !Set 0/0 for projection failure
            x1 = 0.0
            x2 = 0.0
        endif

    end subroutine Map2NH

!---------------------------------
!Tracing routines
    
    !Calculate one-sided trace (in sgn direction)
    subroutine genTrace(Model,ebState,x0,t,xyzn,ijkn,vM,Np,sgn,MaxStepsO,doShueO)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in) :: x0(NDIM),t
        real(rp), intent(inout) :: xyzn(0:MaxFL,NDIM),vM(0:MaxFL,0:NumVFL)
        integer , intent(inout) :: ijkn(0:MaxFL,NDIM)
        integer , intent(out) :: Np
        integer , intent(in) :: sgn
        integer , intent(in), optional :: MaxStepsO
        logical , intent(in), optional :: doShueO

        integer :: MaxSteps
        type(GridPoint_T) :: gPt
        real(rp) :: h
        real(rp), dimension(NDIM) :: B,oB,dx
        real(rp), dimension(NVARMHD) :: Q
        logical :: inDom,doShue

        if (present(MaxStepsO)) then
            MaxSteps = MaxStepsO
        else
            MaxSteps = MaxFL
        endif
        if (present(doShueO)) then
            doShue = doShueO
        else
            doShue = .false.
        endif

    !Initialize
        gPt%xyz = x0
        gPt%t   = t
        Np = 0
        inDom = inDomain(gPt%xyz,Model,ebState%ebGr)

    !Start main loop
        do while (inDom .and. Np <= MaxSteps)
            if (Np == 0) then
                !First time don't have guess
                call locate(gPt%xyz,gPt%ijkG,Model,ebState%ebGr,inDom)
                gPt%dl = getDiag(ebState%ebGr,gPt%ijkG)
                !Pick first h
                h = sgn*Model%epsds*gPt%dl
                !Get first B field
                B = FastMag(gPt%xyz,gPt%t,Model,ebState,inDom,gPt%ijkG)
            else
                !Locate w/ guess (last known cell)
                !call locate(gPt%xyz,gPt%ijkG,Model,ebState%ebGr,inDom,gPt%ijkG)
                !NOTE: Don't need to re-locate b/c last ijkG is correct from FSAL
                gPt%dl = getDiag(ebState%ebGr,gPt%ijkG)
            endif
            
        !Get values at this point on the tube
            !Already have B field from FSAL (first same as last) of previous iteration
            !B = fldInterp(gPt%xyz,gPt%t,Model,ebState,BFLD,inDom,gPt%ijkG)
            xyzn(Np,:) = gPt%xyz
            ijkn(Np,:) = gPt%ijkG
            vM  (Np,0) = norm2(B)
            if (Model%doMHD) then
                Q = mhdInterp(gPt%xyz,gPt%t,Model,ebState,gPt%ijkG)
                vM(Np,1:NVARMHD) = Q
            endif

        !Now do step
            !Call streamstep providing first B and getting last B
            call StreamStep(gPt,Model,ebState,Model%epsds,h,dx,B,oB)
            B = oB !Setup B field for next iteration
            gPt%xyz = gPt%xyz + dx
            !Verify sign of h to be sure
            h = sign(h,sgn*1.0_rp)

        !Check if we're done
            inDom = inDomain(gPt%xyz,Model,ebState%ebGr) .and. (norm2(dx)>TINY)
            if ( doShue .and. (modulo(Np,NpChk) == 0) .and. (norm2(gPt%xyz)>=rShue) ) then
                inDom = inDom .and. inShueMP_SM(gPt%xyz,ShueScl)
            endif

            if (inDom) Np = Np + 1
        enddo

        if (MaxFL == Np) then
            !$OMP CRITICAL
            write(*,*) ANSIRED
            write(*,*) "<WARNING! genTrace hit max tube size!>"
            write(*,*) "Seed: ", x0
            write(*,*) "End : ", xyzn(Np,:)
            write(*,'(a)',advance="no") ANSIRESET, ''
            !$OMP END CRITICAL
        endif

    end subroutine genTrace

    !Slimmed down projection to northern hemisphere for MAGE
    !RinO is optional cut-off inner radius when in northern hemisphere
    !epsO is optional epsilon (otherwise use Model default)
    subroutine mageproject(Model,ebState,x0,t,xyz,Np,isG,epsO,MaxStepsO)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in) :: x0(NDIM),t
        real(rp), intent(out) :: xyz(NDIM)
        integer, intent(out) :: Np
        logical, intent(out) :: isG
        real(rp), intent(in), optional :: epsO
        integer , intent(in), optional :: MaxStepsO

        type(GridPoint_T) :: gPt
        integer :: sgn,MaxSteps
        real(rp) :: eps,h
        real(rp), dimension(NDIM) :: dx,B,oB
        logical :: inDom,isSC,isDone

        if (present(epsO)) then
            eps = epsO
        else
            eps = Model%epsds
        endif

        if (present(MaxStepsO)) then
            MaxSteps = MaxStepsO
        else
            MaxSteps = MaxFL
        endif

        sgn = +1 !Step towards NH
        isG = .false.
    !Initialize
        gPt%xyz = x0
        gPt%t   = t
        Np = 0
        call CheckDone(Model,ebState,gPt%xyz,inDom,isSC,isDone,Np)

        do while ( (.not. isDone) .and. (Np <= MaxSteps) )
        !Locate and get fields
            !Get location in ijk using old ijk as guess if possible
            if (Np == 0) then
                !First step, no guess yet
                call locate(gPt%xyz,gPt%ijkG,Model,ebState%ebGr,inDom)
                gPt%dl = getDiag(ebState%ebGr,gPt%ijkG)
                !Pick first h
                h = sgn*eps*gPt%dl
                !Get first B field
                B = FastMag(gPt%xyz,gPt%t,Model,ebState,inDom,gPt%ijkG)
            else
                !Otherwise don't need to locate b/c last ijkG is correct
                gPt%dl = getDiag(ebState%ebGr,gPt%ijkG)                
            endif

        !Now do step
            call StreamStep(gPt,Model,ebState,eps,h,dx,B,oB)
            B = oB !Setup B field for next iteration
            gPt%xyz = gPt%xyz + dx
            !Verify sign of h to be sure
            h = sign(h,sgn*1.0_rp)

        !Check if we're done
            if (norm2(dx)<TINY) then
                !Multistep integrator had one out value
                isDone = .true.
                inDom  = .false.
                isSC   = .false.
            else
                call CheckDone(Model,ebState,gPt%xyz,inDom,isSC,isDone,Np)
            endif
            Np = Np + 1

        enddo !Main stepping loop

        xyz = gPt%xyz
    !Finished loop somehow, decide if it was good
        
        if (isSC .or. (Np>=MaxSteps)) then
            !Got short circuited or timed out, not good
            isG = .false.
            return
        endif

        if (inDom) then
            !This shouldn't happen
            !$OMP CRITICAL
            write(*,*) 'How did you get here? Bad thing in mageproject'
            !$OMP END CRITICAL
            stop
        else
            !Left domain, but not sure if left from inner boundary
            if (isClosed(xyz,Model)) then
                isG = .true.
            else
                isG = .false.
            endif
            return
        endif

        contains
             
            subroutine CheckDone(Model,ebState,xyz,inDom,isSC,isDone,Np)
                type(chmpModel_T), intent(in) :: Model
                type(ebState_T), intent(in)   :: ebState
                real(rp), intent(inout) :: xyz(NDIM)
                logical, intent(out) :: inDom,isSC,isDone
                integer, intent(in) :: Np
                logical :: inMP

                inDom = inDomain(xyz,Model,ebState%ebGr)
                if ( (modulo(Np,NpChk) == 0) .and. (norm2(gPt%xyz)>=rShue) ) then
                    inMP  = inShueMP_SM(xyz,ShueScl)
                else
                    inMP = .true.
                endif
                
                if (inDom) then
                    !In domain, check if in Shue MP w/ safety factor
                    if (inMP) then
                        isSC   = .false.
                        isDone = .false.
                    else
                        !Not in shue mp
                        isSC = .true.
                        isDone = .true.
                    endif
                else
                    !Straight up out of domain
                    isSC = .false.
                    isDone = .true.
                endif

            end subroutine CheckDone

    end subroutine mageproject

!====      
    !Calculate one-sided projection (in sgn direction)
    subroutine project(Model,ebState,x0,t,xn,sgn,toEquator,failEquator)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in) :: x0(NDIM),t
        integer :: Np
        integer, intent(in) :: sgn
        real(rp), intent(inout) :: xn(NDIM)
        logical, optional, intent(in) :: toEquator
        logical, optional, intent(out) :: failEquator

        logical :: inDom
        real(rp), dimension(NDIM) :: B,E,dx
        real(rp), dimension(NDIM) :: Jb,Jb2,Jb3,F1,F2,F3,F4
        real(rp), dimension(NDIM,NDIM) :: JacB
        real(rp) :: ds,dsmag,dl,MagB,MagJb,dzSgn
        real(rp), dimension(NVARMHD) :: Q
        integer, dimension(NDIM) :: ijk,ijkG
        type(gcFields_T) :: gcF

        !Prime pump
        call locate(x0,ijk,Model,ebState%ebGr,inDom)
        B = fldInterp(x0,t,Model,ebState,BFLD,inDom,ijk)
        
        !Prep for iteration
        Np = 0
        Xn = x0
        dl = getDiag(ebState%ebGr,ijk)
        ds = sgn*Model%epsds*dl

        ijkG = ijk

        if (present(toEquator)) then
           if ((toEquator).and.(.not.(present(failEquator)))) then
              write(*,*) 'Project operator called incorrectly.'
              stop
           endif
        end if

        if (present(failEquator)) failEquator = .true.

        !write(*,*) 'sgn/ds/X0 = ', sgn,ds,x0
        do while (inDom .and. Np <= MaxFL)
        !Locate and get fields
            !Get location in ijk using old ijk as guess
            call locate(Xn,ijk,Model,ebState%ebGr,inDom,ijkG)
            call ebFields(Xn,t,Model,ebState,E,B,ijk,gcFields=gcF)
            MagB = norm2(B)

            ! get the jacobian and new ds
            JacB = gcF%JacB
            MagJb = norm2(JacB)
            dl = getDiag(ebState%ebGr,ijk)

            if (MagJb <= TINY) then
                !Field is constant-ish, use local grid size
                dsmag = dl
            else
                dsmag = MagB/MagJb
            endif
            ds = sgn*Model%epsds*min(dl,dsmag)
        !Update position
            !Convert ds to streamline units
            ds = ds/max(MagB,TINY)
            
            !Get powers of jacobian
            Jb  = matmul(JacB,B  )
            Jb2 = matmul(JacB,Jb )
            Jb3 = matmul(JacB,Jb2)

            !Calculate steps
            F1 = ds*B
            F2 = F1 + (ds*ds/2)*Jb
            F3 = F2 + (ds*ds*ds/4)*Jb2
            F4 = ds*B + ds*ds*Jb + (ds**3.0/2.0)*Jb2 + (ds**4.0/4.0)*Jb3

            !Advance
            dx = (F1+2*F2+2*F3+F4)/6.0
            Xn = Xn + dx
        !Prep for next step, test exit criteria
            inDom = inDomain(xn,Model,ebState%ebGr)
            if (inDom) then
                Np = Np+1

                if (toEquator) then
                    !Check for x'ing

                   ! (trap for a rare weird case)
                   ! if we happened to be just at the equator before making the step, don't do anything
                   ! this can happen if we start tracing exactly from z=0
                   ! the calling getEquatorProjection function should capture this case
                   ! but still add this trap here, just in case
                    dzSgn = Xn(ZDIR)*( Xn(ZDIR)-dx(ZDIR) )
                    if ((dzSgn < 0).or.(Xn(ZDIR)-dx(ZDIR).eq.0.)) then
                        ! interpolate exactly to equator
                        Xn = Xn-dx*abs(Xn(ZDIR))/abs(dx(ZDIR))

                        failEquator = .false.
                        return
                    endif !dzSgn
                 end if
            endif !inDom
 
        enddo
        !write(*,*) 'Found sign/points/distance = ', sgn,Np,norm2(x0-Xn)

    end subroutine project
    

end module streamline
