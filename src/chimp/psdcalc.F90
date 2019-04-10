!Main computational routines for PSD calculations
module psdcalc
    use chmpdefs
    use ebtypes
    use psdtypes
    use psdutils
    use pdfuns

    implicit none

    logical :: doInitWgt = .true.

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
        real(rp) :: n0,kT0,fC,kev,alpha,iScl

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
        !$OMP private(n0,kT0,kev,alpha,fC,iScl)
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

            !Get flow density/temperature
            n0  = psGr%Qrp (ir,ip,DEN)
            kT0 = psGr%kTeq(ir,ip)*psPop%kTScl
            kev = psGr%kC(ik)
            alpha = psGr%aC(ia)

            !Calculate PSD IC
            fC = fPSD0(Model,n0,kT0,kev,alpha)
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
        write(*,*) 'Total weight = ', sum(psPop%wgt(:))
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
            
            !$OMP ATOMIC
            psPop%fPSD(ir,ip,ik,ia) = psPop%fPSD(ir,ip,ik,ia) + w
            !$OMP ATOMIC
            psPop%nPSD(ir,ip,ik,ia) = psPop%nPSD(ir,ip,ik,ia) + 1
            
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
        !TODO: Can check for under-sampled cells here and smooth based on neighbors


        end associate

    end subroutine CalcPSD


end module psdcalc