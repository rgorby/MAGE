module streamline
    use chmpdefs
    use ebtypes
    use math
    use ebinterp

    implicit none

    contains

    subroutine genStream(Model,ebState,x0,t,fL)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(fLine_T), intent(inout) :: fL

        integer :: N1,N2,i,Np,Nm,n
        real(rp) :: dx(NDIM)
        real(rp) :: Xn(0:MaxFL,NDIM,2),Vn(0:MaxFL,0:NumVFL,2)
        integer :: ijkn(0:MaxFL,NDIM,2)
        logical :: inDom

        !Start by emptying line
        call cleanStream(fL)
        inDom = inDomain(x0,Model,ebState%ebGr)

        if (.not. inDom) then
            return
        endif
        
        
        !Get traced line
        fL%x0 = x0
        fl%lnVars(0)%idStr = "B"

        if (Model%doMHD) then
            !Create holders for all MHD variables
            fl%lnVars(DEN)%idStr = "D"
            fl%lnVars(VELX)%idStr = "Vx"
            fl%lnVars(VELY)%idStr = "Vy"
            fl%lnVars(VELZ)%idStr = "Vz"
            fl%lnVars(PRESSURE)%idStr = "P"
        endif

        call genTrace(Model,ebState,x0,t,Xn(:,:,1),ijkn(:,:,1),Vn(:,:,1),N1,-1)
        call genTrace(Model,ebState,x0,t,Xn(:,:,2),ijkn(:,:,2),Vn(:,:,2),N2,+1)
        
        !Create field line
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
        do k=-bTrc%Nm,bTrc%Np-1
            dL = norm2(bTrc%xyz(k+1,:)-bTrc%xyz(k,:))
            L = L + dL
        enddo

    end function FLArc

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

        ! !$OMP CRITICAL
        ! write(*,*) '---'
        ! write(*,*) 'dvB = ', dvB
        ! write(*,*) 'bP/bPb = ', bP,bPb

        ! write(*,*) 'Beta (avg,int) = ', bP/bPb,bBeta/dvB
        ! write(*,*) '---'
        ! !$OMP END CRITICAL

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
            S = S + dV*(pAvg**(1.0/Gam))
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

        !Find minimum field and where it occurs
        Beq = minval(bTrc%lnVars(0)%V) !Min field strength
        i0  = minloc(bTrc%lnVars(0)%V,dim=1) !Note this is between 1:N
        !Need to offset i0, iMin = i0-Nm-1
        iMin = i0-Nm-1
        xeq = bTrc%xyz(iMin,:)

        end associate
    end subroutine FLEq
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
!Tracing routines

    !Calculate one-sided trace (in sgn direction)
    subroutine genTrace(Model,ebState,x0,t,xyzn,ijkn,vM,Np,sgn,toEquatorO)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in) :: x0(NDIM),t
        real(rp), intent(inout) :: xyzn(0:MaxFL,NDIM),vM(0:MaxFL,0:NumVFL)
        integer, intent(inout) :: ijkn(0:MaxFL,NDIM)
        integer, intent(out) :: Np
        integer, intent(in) :: sgn
        logical, optional, intent(in) :: toEquatorO

        logical :: inDom,toEquator
        real(rp), dimension(NDIM) :: Xn,B,E,dx
        real(rp), dimension(NDIM) :: Jb,Jb2,Jb3,F1,F2,F3,F4
        real(rp), dimension(NDIM,NDIM) :: JacB
        real(rp) :: ds,dl,MagJb,dsmag
        real(rp), dimension(NVARMHD) :: Q
        integer, dimension(NDIM) :: ijk,ijkG
        type(gcFields_T) :: gcF

        if (present(toEquatorO)) then
            toEquator = toEquator
        else
            toEquator = .false.
        endif
        
        !Prime pump
        call locate(x0,ijk,Model,ebState%ebGr,inDom)
        B = fldInterp(x0,t,Model,ebState,BFLD,inDom,ijk)
        
        !Save initial values
        xyzn(0,:) = x0
        vM (0,0)  = norm2(B)
        ijkn(0,:) = ijk

        if (Model%doMHD) then
            Q = mhdInterp(x0,t,Model,ebState,ijk)
            vM(0,1:NVARMHD) = Q
        endif

        !Prep for iteration
        Np = 0
        Xn = x0
        dl = getDiag(ebState%ebGr,ijk)
        !Note: ds gets multipled by mag(B)
        ds = sgn*min( Model%epsds*dl/norm2(B), dl )
        
        ijkG = ijk

        !write(*,*) 'sgn/ds/X0 = ', sgn,ds,x0
        do while (inDom .and. Np <= MaxFL)

        !Locate and get fields
            !xyz(Np,:) = Xn
            !vM(Np,1) = norm2(B)
            !ijkG = ijk

            !Get location in ijk using old ijk as guess
            call locate(xN,ijk,Model,ebState%ebGr,inDom,ijkG)
            call ebFields(Xn,t,Model,ebState,E,B,ijk,gcFields=gcF)
        !Update position
            !Get powers of jacobian
            JacB = gcF%JacB
            Jb = matmul(JacB,B)
            Jb2 = matmul(JacB,Jb)
            Jb3 = matmul(JacB,Jb2)

            !Calculate steps
            F1 = ds*B
            F2 = F1 + (ds*ds/2)*Jb
            F3 = F2 + (ds*ds*ds/4)*Jb2
            F4 = ds*B + ds*ds*Jb + (ds**3.0/2.0)*Jb2 + (ds**4.0/4.0)*Jb3

            !Advance
            dx = (F1+2*F2+2*F3+F4)/6.0
            Xn = Xn + dx
            
            !Get field/inDom at new position
            ijkG = ijk !Use better guess
            call locate(xN,ijk,Model,ebState%ebGr,inDom,ijkG)
            ijkG = ijk !Update guess

            B = fldInterp(Xn,t,Model,ebState,BFLD,inDom,ijkG)
            
            if (inDom) then
                Np = Np+1

                xyzn(Np,:) = Xn
                ijkn(Np,:) = ijk
                vM(Np,0) = norm2(B)
                if (Model%doMHD) then
                    Q = mhdInterp(Xn,t,Model,ebState,ijk)
                    vM(Np,1:NVARMHD) = Q
                endif

                ! if tracing to equator and we just crossed it: quit
                if (toEquator) then
                   if ( Xn(ZDIR)*(Xn(ZDIR)-dx(ZDIR) ) < 0. ) return
                endif

                !Get new ds
                MagJb = sqrt(sum(JacB**2.0))
                if (MagJb <= TINY) then
                    !Field is constant-ish, use local grid size
                    dl = getDiag(ebState%ebGr,ijk)
                    dsmag = Model%epsds*dl/norm2(B)
                else
                    dsmag = Model%epsds/MagJb
                endif
                ds = sgn*min(dl,dsmag)
            endif
        enddo

    end subroutine genTrace

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
        real(rp) :: ds,dsmag,dl,MagJb,dzSgn
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
        ds = sgn*min( Model%epsds*dl/norm2(B), dl )
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

            ! get the jacobian
            JacB = gcF%JacB

            !Get new ds
            MagJb = sqrt(sum(JacB**2.0))
            if (MagJb <= TINY) then
                !Field is constant-ish, use local grid size
                dl = getDiag(ebState%ebGr,ijk)
                dsmag = Model%epsds*dl/norm2(B)
            else
                dsmag = Model%epsds/MagJb
            endif
            ds = sgn*min(dl,dsmag)     
        !Update position
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
    
    function getDiag(ebGr,ijk) result (dl)
        type(ebGrid_T), intent(in)   :: ebGr
        integer, intent(in) :: ijk(NDIM)
        real(rp) :: dl
        integer :: i,j,k
        i = ijk(IDIR) ; j = ijk(JDIR) ; k = ijk(KDIR)

        dl = norm2( ebGr%xyz(i+1,j+1,k+1,:)-ebGr%xyz(i,j,k,:) )

    end function getDiag

    subroutine cleanStream(fL)
        type(fLine_T), intent(inout) :: fL

        integer :: i
        if (allocated(fL%xyz)) deallocate(fL%xyz)
        if (allocated(fL%ijk)) deallocate(fL%ijk)
        do i=0,NumVFL
            if (allocated(fL%lnVars(i)%V)) deallocate(fL%lnVars(i)%V)
        enddo
        fL%x0 = 0.0
        fL%Nm = 0
        fL%Np = 0
    end subroutine cleanStream
end module streamline
