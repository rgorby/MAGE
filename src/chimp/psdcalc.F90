!Main computational routines for PSD calculations
module psdcalc
    use chmpdefs
    use ebtypes
    use psdtypes
    use psdutils
    use pdfuns

    implicit none

    logical :: doInitWgt = .true.

    abstract interface
        subroutine PSDShape_T(x0,xI,wX)
            import :: rp
            real(rp), intent(in)  :: x0
            real(rp), intent(in)  :: xI(-1:+2) !x0 is between xI(0) and xI(1)
            real(rp), intent(out) :: wX(-1:+1)
        end subroutine PSDShape_T
    end interface
    procedure(PSDShape_T), pointer, private :: ShapeWeight => ShapeWeight_CIC

    contains

    !Calculate weights for unweighted particles
    subroutine CalcWeights(Model,psGr,ebSt,psPop)
        type(chmpModel_T), intent(in)    :: Model
        type(psdPop_T)   , intent(inout) :: psPop
        type(PSEq_T)     , intent(in)    :: psGr
        type(ebState_T)  , intent(in)    :: ebSt

        integer, dimension(:,:,:,:), allocatable :: nPS
        integer :: NumW, n,ir,ip,ik,ia, NumOut
        integer, dimension(NVARPS) :: idx
        real(rp) :: n0,kT0,fC,kev,alpha,iScl,R,phi

        !Find how many particles need weighting
        NumW = count( (.not. psPop%isWgt) .and. (psPop%isIn) )
        if (NumW == 0) return

        !If not streaming particles, only weight once
        if ( (.not. Model%doStream) .and. (.not. doInitWgt) ) return

        !Otherwise, do some weighting
        if (doInitWgt) then
            write(*,*) '--- 1st Weighting ---'
            doInitWgt = .false.
        else
            write(*,*) '--- Re-Weighting  ---'
        endif
        
        write(*,*) '    Total TPs      = ', psPop%NumTP
        write(*,*) '    Active TPs     = ', count(psPop%isIn)
        write(*,*) '    Unweighted TPs = ', NumW, ' (before)'


        associate(Nr=>psGr%Nr,Np=>psGr%Np,Nk=>psGr%Nk,Na=>psGr%Na)
        allocate(nPS(Nr,Np,Nk,Na))
        nPS = 0
        NumOut = 0

    !Start by doing PS count
        !write(*,*) 'Doing PS count ...'
        !Localize all particles to PS grid and count per bin
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(idx,ir,ip,ik,ia)
        do n=1,psPop%NumTP
            if (.not. psPop%isIn(n)) cycle !Skip bad particles
            !Otherwise, find where you are in psGr
            idx = psLoc(psPop%TPs(n,:),psGr)
            !Test for good value
            if (dot_product(idx,idx) <= TINY) then
                !$OMP ATOMIC
                NumOut = NumOut + 1
                cycle
            endif
            ir = idx(PSRAD)
            ip = idx(PSPHI)
            ik = idx(PSKINE)
            ia = idx(PSALPHA)
            !Add to counter
            !$OMP ATOMIC
            nPS(ir,ip,ik,ia) = nPS(ir,ip,ik,ia) + 1
        enddo

    !Calc weights for new particles
        !Loop over all particles, for unweighted and in particles
        ! w ~ f*dG/ntp
        !write(*,*) 'Doing weighting loop ...'
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,idx,ir,ip,ik,ia)  &
        !$OMP private(n0,kT0,kev,alpha,fC,iScl,R,phi)
        do n=1,psPop%NumTP
            if ( (.not. psPop%isIn(n)) .or. (psPop%isWgt(n)) ) then
                !Don't bother if particle isn't in or has already been weighted
                cycle
            endif
            !Check if particle is in space
            idx = psLoc(psPop%TPs(n,:),psGr)
            if (dot_product(idx,idx) <= TINY) then
                !If not in PS grid, set weight to 0 and don't look back
                psPop%isWgt(n) = .true.
                psPop%wgt(n) = 0.0
                cycle
            endif

            !If we're still here, then let's get to work
            ir = idx(PSRAD)
            ip = idx(PSPHI)
            ik = idx(PSKINE)
            ia = idx(PSALPHA)

            !Get TP location
            R     = psGr%rC(ir)
            phi   = psGr%pC(ip)
            kev   = psGr%kC(ik)
            alpha = psGr%aC(ia)

            !Get flow density/temperature
            n0  = psGr%Qrp (ir,ip,DEN)
            kT0 = psGr%kTeq(ir,ip)*psPop%kTScl

            !Calculate PSD IC
            fC = fPSD0(Model,R,phi,kev,alpha,n0,kT0)
            
            !Rescale PSD IC if doing injection over time
            if (Model%doStream) then
                !Injection scaling, (tau*u)/dWr
                !Negative from turning Vr -> V-earthward
                iScl = -psPop%dTau*psGr%Vreq(ir,ip)/psPop%dShell
                iScl = max(iScl,0.0)
                fC = iScl*fC
            endif
            !Set weights
            psPop%wgt(n) = fC*psGr%dG(ir,ip,ik,ia)/nPS(ir,ip,ik,ia)
            psPop%isWgt(n) = .true.
            if (psPop%wgt(n) /= psPop%wgt(n)) then
                !$OMP CRITICAL
                write(*,*) '--------'
                write(*,*) 'n = ', n
                write(*,*) 'I = ',ir,ip,ik,ia
                write(*,*) 'fC = ', fC
                write(*,*) 'n0 = ',n0
                write(*,*) 'kT0 /kev = ', kT0,kev
                write(*,*) 'alpha = ', alpha
                write(*,*) 'dG = ', psGr%dG(ir,ip,ik,ia)
                write(*,*) 'nPS = ', nPS(ir,ip,ik,ia)
                write(*,*) '--------'
                !$OMP END CRITICAL
            endif
        enddo
        NumW = count( (.not. psPop%isWgt) .and. (psPop%isIn) )
        write(*,*) '    Unweighted TPs = ', NumW, ' (after)'
        write(*,*) '    Out of Dom TPs = ', NumOut
        

    !Finished
        end associate
        write(*,*) 'Total weight   = ', sum(psPop%wgt(:))
        write(*,*) 'Min/Max weight = ', minval(psPop%wgt),maxval(psPop%wgt)
    end subroutine CalcWeights

    !Given weights for all TPs, calculate PSD on given PS
    subroutine CalcPSD(Model,psGr,ebSt,psPop)
        type(chmpModel_T), intent(in)    :: Model
        type(psdPop_T)   , intent(inout) :: psPop
        type(PSEq_T)     , intent(in)    :: psGr
        type(ebState_T)  , intent(in)    :: ebSt

        integer :: n,ir,ip,ik,ia
        integer, dimension(NVARPS) :: idx

        real(rp) :: w,wTot,dG

        associate(Nr=>psGr%Nr,Np=>psGr%Np,Nk=>psGr%Nk,Na=>psGr%Na)

        !Just start fresh each time
        !TODO: Avoid reallocating if already allocated and correctly shaped
        if (allocated(psPop%fPSD)) deallocate(psPop%fPSD)
        if (allocated(psPop%nPSD)) deallocate(psPop%nPSD)
        allocate(psPop%fPSD(Nr,Np,Nk,Na))
        allocate(psPop%nPSD(Nr,Np,Nk,Na))
        psPop%fPSD = 0.0
        psPop%nPSD = 0


        !Start by looping over all particles and depositing weights
        !TODO: For now just using delta function, add shape-based weight deposition here
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(idx,ir,ip,ik,ia,w)  &
        !$OMP schedule(dynamic)
        do n=1,psPop%NumTP
            if (.not. psPop%isIn(n) .or. .not. psPop%isWgt(n) ) cycle !Skip bad particles
            !Otherwise, find where you are in psGr
            idx = psLoc(psPop%TPs(n,:),psGr)
            !Test for good value
            if (dot_product(idx,idx) <= TINY) then
                cycle
            endif
            ir = idx(PSRAD)
            ip = idx(PSPHI)
            ik = idx(PSKINE)
            ia = idx(PSALPHA)

            !Add weight and ps counter
            w = psPop%wgt(n)
            
            call depositWeight(Model,psGr,psPop,psPop%TPs(n,:),idx,w)
            
        enddo

        !Now, divide by dG to create PSD
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(wTot,dG,ia,ik,ip,ir)
        do ia=1,Na
            do ik=1,Nk
                do ip=1,Np
                    do ir=1,Nr
                        if (psGr%isClosed(ir,ip)) then
                            wTot = psPop%fPSD(ir,ip,ik,ia)
                            dG   = psGr%dG   (ir,ip,ik,ia)
                        else
                            wTot = 0.0
                            dG = 1.0
                        endif

                        if (dG>dGTINY) then
                            psPop%fPSD(ir,ip,ik,ia) = wTot/dG
                        else
                            psPop%fPSD(ir,ip,ik,ia) = 0.0
                        endif

                    enddo
                enddo
            enddo
        enddo


        end associate

    end subroutine CalcPSD

    subroutine depositWeight(Model,psGr,psPop,Q,idx,wgt)
        type(chmpModel_T), intent(in)    :: Model
        type(PSEq_T)     , intent(in)    :: psGr
        type(psdPop_T)   , intent(inout) :: psPop
        
        integer , intent(in) :: idx(NVARPS)
        real(rp), intent(in) :: wgt, Q(NVARPS)

        integer :: ir,ip,ik,ia
        real(rp), dimension(-1:1,-1:1,-1:1,-1:1) :: wX
        real(rp), dimension(-1:1) :: wR,wP,wK,wA

        real(rp) :: r0,p0,k0,a0,wIJK
        integer :: dr,dp,dk,da
        integer :: irp,ipp,ikp,iap
        logical :: inR,inP,inK,inA,inGrid

        ir = idx(PSRAD)
        ip = idx(PSPHI)
        ik = idx(PSKINE)
        ia = idx(PSALPHA)

        r0 = Q(PSRAD)
        p0 = Q(PSPHI)
        if (p0<0) p0=p0+2*PI

        k0 = Q(PSKINE)
        a0 = Q(PSALPHA)
        
        
        if (.not. psGr%doShape) then
            !$OMP ATOMIC
            psPop%fPSD(ir,ip,ik,ia) = psPop%fPSD(ir,ip,ik,ia) + wgt
            !$OMP ATOMIC
            psPop%nPSD(ir,ip,ik,ia) = psPop%nPSD(ir,ip,ik,ia) + 1

            return
        endif

        !If we're still here, do shaping
        wX = 0.0

        wR = 0.0
        wK = 0.0
        wA = 0.0
        wP = 0.0

        !Radius
        if ( (ir == 1) .or. (ir == psGr%Nr) ) then
            wR(0) = 1.0
        else
            !This is in middle
            call ShapeWeight(r0,psGr%rI(ir-1:ir+2),wR)
        endif

        !Energy
        if ( (ik == 1) .or. (ik == psGr%Nk) ) then
            wK(0) = 1.0
        else
            !This is in middle
            call ShapeWeight(k0,psGr%kI(ik-1:ik+2),wK)
        endif

        !Alpha
        if ( (ia == 1) .or. (ia == psGr%Na) ) then
            wA(0) = 1.0
        else
            !This is in middle
            call ShapeWeight(a0,psGr%aI(ia-1:ia+2),wA)
        endif

        !phi
        if ( (ip == 1) .or. (ip == psGr%Np) ) then
            wP(0) = 1.0
        else
            !This is in middle
            call ShapeWeight(p0,psGr%pI(ip-1:ip+2),wP)
        endif

        do dr=-1,1
            do dp=-1,1
                do dk=-1,1
                    do da=-1,1
                        irp = ir+dr
                        ipp = ip+dp
                        ikp = ik+dk
                        iap = ia+da

                        wIJK = wR(dr)*wP(dp)*wK(dk)*wA(da)
                        wX(dr,dp,dk,da) = wIJK

                        !Check if this is a good cell
                        inR = (irp >= 1) .and. (irp <= psGr%Nr)
                        inP = (ipp >= 1) .and. (ipp <= psGr%Np)
                        inK = (ikp >= 1) .and. (ikp <= psGr%Nk)
                        inA = (iap >= 1) .and. (iap <= psGr%Na)

                        inGrid = inR .and. inP .and. inK .and. inA
                        if (inGrid) then
                            !$OMP ATOMIC
                            psPop%fPSD(irp,ipp,ikp,iap) = psPop%fPSD(irp,ipp,ikp,iap) + wgt*wIJK

                        endif !Update

                    enddo
                enddo
            enddo
        enddo

        !$OMP ATOMIC
        psPop%nPSD(ir,ip,ik,ia) = psPop%nPSD(ir,ip,ik,ia) + 1

        !wIJK = sum(wX)-1.0
        !write(*,*) 'wIJK = ', wIJK

    end subroutine depositWeight

!--- Shape functions
!x0 is particle position, xI are the interfaces of the 3 closest cells (-1->0,0->1,1->2)
    subroutine ShapeWeight_CIC(x0,xI,wX)
        real(rp), intent(in)  :: x0
        real(rp), intent(in)  :: xI(-1:+2)        
        real(rp), intent(out) :: wX(-1:+1)

        real(rp) :: xCC,xM,xP,dX

        !Initialize
        wX = 0.0

        xM  = 0.5*(xI(-1)+xI( 0))
        xCC = 0.5*(xI( 0)+xI(+1))
        xP  = 0.5*(xI(+1)+xI(+2))

        if (x0>=xCC) then
            dX = xP-xCC
            wX( 0) = (xP-x0 )/dX
            wX(+1) = (x0-xCC)/dX
        else
            dX = xCC-xM
            wX(-1) = (xCC-x0)/dX
            wX( 0) = (x0 -xM)/dX
        endif
    end subroutine ShapeWeight_CIC

    !TODO: Fix bug in TSC shaping
    ! subroutine ShapeWeight_TSC(x0,xI,wX)
    !     real(rp), intent(in)  :: x0
    !     real(rp), intent(in)  :: xI(-1:+2)        
    !     real(rp), intent(out) :: wX(-1:+1)

    !     real(rp) :: dX,xPI,xMI,mArg,pArg
    !     !Initialize
    !     wX = 0.0

    !     dX = xI(+1)-xI( 0) !Width of center cell

    !     xPI = xI(+1)-x0
    !     xMI = x0-xI(0)

    !     mArg = 1.0-xMI
    !     pArg = 1.0-xPI

    !     wX(-1) = 0.5*mArg*mArg
    !     wX(+1) = 0.5*pArg*pArg
    !     wX( 0) = 1.0-wX(-1)-wX(+1)


    ! end subroutine ShapeWeight_TSC

end module psdcalc