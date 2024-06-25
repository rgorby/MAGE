!Handle ring averaging about pole
!Notation
!nSi/nSe: Slices.  Cyl (Z), LFM (R)
!nR: Rings. Cyl (R), LFM (J)

module ringav
    use gamtypes
    use gamutils
    use fields
    use math
    use ringutils
    use ringrecon
    use multifluid
    
    implicit none

    !Enumerators for Fourier reduction coefficients
    enum, bind(C)
        enumerator :: FTCOS=1, FTSIN
    endenum
    
    contains

    !Takes state variable and performs ring average on both gas/mag variables
    subroutine RingAverage(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), target, intent(inout) :: State
        real(rp), dimension(:,:,:,:), allocatable, save :: SmoothE !Smoothing electric field
        
        !DIR$ attributes align : ALIGN :: SmoothE

    !----
    !Initialization
        !Start by making sure we should be here
        if (.not. (Model%Ring%doS .or. Model%Ring%doE)) then
            return
        endif
        call Tic("Ring-Init")
        !Initialize array for smoothing field
        call InitRAVec(Model,Gr,SmoothE)

        call Toc("Ring-Init")
    !----
    !Magnetic fields

        call Tic("RA-MAG")

        call Tic("SmoothE")
        !Start by getting smoothing field
        if (Model%Ring%doS) call CalcSmoothE(Model,Gr,State,SmoothE,SPOLE)
        if (Model%Ring%doE) call CalcSmoothE(Model,Gr,State,SmoothE,EPOLE)
        call Toc("SmoothE")

        call Tic("ApplyE")
        call PoleE(Model,Gr,SmoothE) !Singularity fix smoothing field
        call E2Flux(Model,Gr,State%magFlux,SmoothE) !Apply smoothing field
        
        !Clean field loops
        if (doCleanLoop) call CleanLoops(Model,Gr,State)
        call bFlux2Fld(Model,Gr,State%magFlux,State%Bxyz) !Recalculate cell-centered fields
        
        call Toc("ApplyE")

        call Toc("RA-MAG")

    !----
    !Hydro variables
        call Tic("RA-GAS")
        if (Model%Ring%doS) call SmoothGas(Model,Gr,State,SPOLE)
        if (Model%Ring%doE) call SmoothGas(Model,Gr,State,EPOLE)
        call Toc("RA-GAS")

    end subroutine RingAverage

    !Apply ring average to hydro variables
    subroutine SmoothGas(Model,Gr,State,XPOLE)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Gr
        type(State_T), intent(inout) :: State
        integer, intent(in) :: XPOLE

        integer :: nR,nS,nv,n
        integer :: NumCh, dC, lRs,lRe
        integer :: s,s0,sE
        integer :: jshift,jPhase(Model%Ring%NumR)
        real(rp), dimension(0:NFTMAX,FTCOS:FTSIN,NVAR) :: cFTxyz !Coefficients for Fourier reduction
        real(rp), dimension(Model%Ring%Np,NVAR) :: rW,raW
        logical , dimension(Model%Ring%Np)      :: gW

        !DIR$ attributes align : ALIGN :: rW,raW,gW

        associate(Np=>Model%Ring%Np)
        jPhase = 0
        if (doShift) then
            !Get random shift to apply to rings before chunking
            do n=1,Model%Ring%NumR
                jPhase(n) = nint( genRand(0.0_rp,1.0_rp*Np) )
            enddo
        endif

        s0 = BLK
        if (Model%doMultiF) then
            s0 = 1
            sE = Model%nSpc
        else
            sE = BLK
        endif

        !Loop over all species
        do s=s0,sE
            !Loop over slices, rings/slice
            !Do OMP on slices
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(nR,nv,n,lRs,lRe,NumCh,dC) &
            !$OMP private(rW,raW,gW,jshift,cFTxyz)
            do nS=Model%Ring%nSi,Model%Ring%nSe
                do nR=1,Model%Ring%NumR
                    NumCh = Model%Ring%Nch(nR) !Number of total chunks for this ring
                    dC = Np/NumCh !# of cells per chunk for this ring

                !Pull ring hydro vars and convert to ring vars
                    call PullRingCC(Model,Gr,State%Gas(:,:,:,1:NVAR,s),rW,nR,nS,NVAR,XPOLE)

                    !Shift if necessary
                    if (doShift) then
                        !Do positive shifts on pulled rings
                        !Different phase on each radius
                        jshift = jPhase(nR)
                        rW = cshift(rW,+jshift,dim=1)
                    endif

                    !Convert to ring variables
                    call Gas2Ring(Model,rW)

                    !Decide on good cells
                    if (Model%doMultiF) then
                        gW = (rW(:,DEN) > Spcs(s)%dVac)
                        if (.not. any(gW)) cycle
                    else
                        gW = .true.
                    endif

                !Create chunk averaged array
                    !Clean ring (Vxyz) if desired before chunk-averaging
                    if (doVClean) then
                        cFTxyz = 0.0
                        do nv=MOMX,MOMZ
                            call CleanRingWgt(Model,rW(:,nv),cFTxyz(:,:,nv),rW(:,DEN),gW,NFTVEL)
                        enddo
                    endif
                
                    !Loop over chunks, pull average over each chunk
                    do nv=1,NVAR
                        do n=1,NumCh
                            !Index bounds for this chunk in 1:Np ring
                            lRs = 1 + (n-1)*dC
                            lRe = lRs + dC - 1
                            raW(lRs:lRe,nv) = sum(rW(lRs:lRe,nv))/dC
                        enddo !Chunk loop
                    enddo !Var loop

                !Now reconstruct ring/s
                    !Reconstruct mass first w/ standard reconstruction
                    call ReconstructRing(Model,raW(:,DEN),NumCh,gW)

                    !Loop over remaining variables and either do pure recon or mass-recon
                    do nv=2,NVAR
                        if (Model%Ring%doMassRA) then
                           call WgtReconstructRing(Model,raW(:,nv),raW(:,DEN),NumCh,gW)
                        else
                           call ReconstructRing(Model,raW(:,nv),NumCh,gW)
                        endif  
                    enddo

                    !Dirty ring (Vxyz) if you cleaned before converting back
                    if (doVClean) then
                        do nv=MOMX,MOMZ
                            call DirtyRingWgt(Model,raW(:,nv),cFTxyz(:,:,nv),raW(:,DEN),gW,NFTVEL)
                        enddo
                    endif

                    !Convert back to gas variables
                    call Ring2Gas(Model,raW)

                    !Unshift if necessary
                    if (doShift) then
                        !Unshift before pushing back into state
                        !Doing negative of original shift
                        !Working on raW not rW
                        jshift = jPhase(nR)
                        raW = cshift(raW,-jshift,dim=1)
                    endif

                    !Push ring back
                    call PushRingCC(Model,Gr,State%Gas(:,:,:,1:NVAR,s),raW,nR,nS,NVAR,XPOLE)

                enddo !Ring loop
            enddo !Slice loop

        enddo !Species loop

        if (Model%doMultiF) then
            !Doing ALL cells, not just ring-avg ones (meh)
            call State2Bulk(Model,Gr,State)
        endif

        end associate
    end subroutine SmoothGas
    
    subroutine CalcSmoothE(Model,Gr,State,SmoothE,XPOLE)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp), intent(inout) :: SmoothE(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM)
        integer, intent(in) :: XPOLE

        integer :: nR,nS,d,n,lRs,lRe,m
        integer :: fD,eD,NumCh,dC

        integer, dimension(2) :: fDIRS,eDIRS !Face/E-field directions
        real(rp), dimension(2) :: eScls
        real(rp) :: eScl, cumFlx

        real(rp), dimension(Model%Ring%Np,NDIM) :: bFlx,dEr
        real(rp), dimension(Model%Ring%Np) :: rFlx,raFlx,dFlx
        real(rp), dimension(0:NFTMAX,FTCOS:FTSIN) :: cFT !Coefficients for Fourier reduction

        !DIR$ attributes align : ALIGN :: bFlx,dEr,rFlx,raFlx,dFlx

    !---
    !Start by getting directions
        select case (Model%Ring%GridID)
        !Set face flux->E field correlations
        case ("cyl")
            !Cylindrical
            fDIRS = [IDIR,KDIR]
            eDIRS = [KDIR,IDIR]
            eScls = [1.0,-1.0]
        case ("lfm")
            !LFM
            fDIRS = [JDIR,IDIR]
            eDIRS = [IDIR,JDIR]
            eScls = [1.0,-1.0]
        end select

    !---
    !Loop over slices,rings/slice
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(nR,nS,d,n,lRs,lRe,m,fD,eD,NumCh,dC) &
        !$OMP private(bFlx,cFT,dEr,eScl,cumFlx,rFlx,raFlx,dFlx)
        do nS=Model%Ring%nSi,Model%Ring%nSe+1
            do nR=1,Model%Ring%NumR
                NumCh = Model%Ring%Nch(nR) !Number of total chunks for this ring
                dC = Model%Ring%Np/NumCh !# of cells per chunk for this ring

                !Pull ring/s about pole, all flux components
                call PullRingI(Model,Gr,State%magFlux,bFlx,nR,nS,XPOLE)
            !---
            !Calculate components of electric field (2 directions)
                dEr = 0.0
                do d=1,2
                    fD   = fDIRS(d)
                    eD   = eDIRS(d)
                    eScl = eScls(d)

                    !Get this flux direction
                    rFlx = bFlx(:,fD)

                    !Clean ring and then chunk-average
                    call CleanRing(Model,rFlx,cFT,NFTMAG)

                    do n=1,NumCh
                        !Find indices for this chunk in the global ring
                        lRs = 1 + (n-1)*dC
                        lRe = lRs + dC - 1

                        raFlx(lRs:lRe) = sum(rFlx(lRs:lRe))/dC

                    enddo !Chunk loop

                    !Do reconstruction on filtered/chunk-averaged fluxes
                    call ReconstructRing(Model,raFlx,NumCh)

                    !Add back in Fourier modes
                    call DirtyRing(Model,raFlx,cFT,NFTMAG)

                    !Calculate changes in flux
                    dFlx = raFlx - bFlx(:,fD)

                !Accumulate into ring E field
                    !Loop over chunks and within chunks
                    do n=1,NumCh
                        !Find indices for this chunk in the global ring
                        lRs = 1 + (n-1)*dC
                        lRe = lRs + dC - 1

                        cumFlx = -dFlx(lRs)
                        do m=2,dC
                            dEr(lRs+m-1,eD) = eScl*cumFlx
                            cumFlx = cumFlx - dFlx(lRs+m-1)
                        enddo !Intra-chunk loop
                    enddo !Chunk loop

                enddo !E-field directions
            !---
            !Ring electric field complete, store into smoothing E field
                select case (Model%Ring%GridID)
                case ("cyl") !Cylindrical
                    !Ek starts @ is+1, rest at is
                    !Only one side of pole
                    if (XPOLE == SPOLE) then
                        SmoothE(Gr%is+nR  ,Gr%js:Gr%je,nS,KDIR) = dEr(:,KDIR)
                        SmoothE(Gr%is+nR-1,Gr%js:Gr%je,nS,IDIR) = dEr(:,IDIR)
                        SmoothE(Gr%is+nR-1,Gr%js:Gr%je,nS,JDIR) = dEr(:,JDIR)
                    endif
                case ("lfm")
                    !LFM-style singularity

                    if (XPOLE == SPOLE) then !+X pole
                        !Ei starts at js+1, rest at js (K-field is 0)
                        SmoothE(nS,Gr%js+nR  ,Gr%ks:Gr%ke,IDIR) = dEr(:,IDIR)
                        SmoothE(nS,Gr%js+nR-1,Gr%ks:gr%ke,JDIR) = dEr(:,JDIR)
                        SmoothE(nS,Gr%js+nR-1,Gr%ks:gr%ke,KDIR) = dEr(:,KDIR)

                        !Replicate for ke+1
                        SmoothE(nS,Gr%js+nR  ,Gr%ke+1,IDIR) = SmoothE(nS,Gr%js+nR  ,Gr%ks,IDIR)
                        SmoothE(nS,Gr%js+nR-1,Gr%ke+1,JDIR) = SmoothE(nS,Gr%js+nR-1,Gr%ks,JDIR)
                    endif

                    if (XPOLE == EPOLE) then !-X pole
                        SmoothE(nS,Gr%je-nR+1,Gr%ks:Gr%ke,IDIR:KDIR) = dEr(:,IDIR:KDIR)
                        !Replicate for ke+1
                        SmoothE(nS,Gr%je-nR+1,Gr%ke+1,IDIR) = SmoothE(nS,Gr%je-nR+1,Gr%ks,IDIR)
                        SmoothE(nS,Gr%je-nR+1,Gr%ke+1,JDIR) = SmoothE(nS,Gr%je-nR+1,Gr%ks,JDIR)
                    endif

                end select

            enddo !Rings, nR
        enddo !Slices, nS

    end subroutine CalcSmoothE

    !Clean loops around axis
    subroutine CleanLoops(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T) , intent(in) :: Gr
        type(State_T), intent(inout) :: State

        integer :: nS,nR,fD
        real(rp) :: tScl,avgFlx

        associate(bFlux=>State%magFlux)
        
        !Loop over slices
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(nS,nR,fD,tScl,avgFlx)
        do nS=Model%Ring%nSi,Model%Ring%nSe+1
            if (doFastLoop) then
                tScl = 1.0
            else
                tScl = 1.0/Model%Ring%Np
            endif
            
            select case (Model%Ring%GridID)
            case ("cyl")
                !Cylindrical
                fD = JDIR
                nR = 1
                bFlux(Gr%is,Gr%js:Gr%je+1,nS,fD) = bFlux(Gr%is,Gr%js:Gr%je+1,nS,fD) - sum(bFlux(Gr%is,Gr%js:Gr%je,nS,fD))/Model%Ring%Np
            case ("lfm")
                fD = KDIR
                !Diffuse over NumRings timesteps
                if (Model%Ring%doS) then
                    nR = Gr%js
                    avgFlx = sum(bFlux(nS,nR,Gr%ks:Gr%ke,fD))/Model%Ring%Np
                    bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) = bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) - tScl*avgFlx
                endif
                if (Model%Ring%doE) then
                    nR = Gr%je
                    avgFlx = sum(bFlux(nS,nR,Gr%ks:Gr%ke,fD))/Model%Ring%Np
                    bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) = bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) - tScl*avgFlx
                endif
            end select
        enddo

        end associate

    end subroutine CleanLoops

    !Set electric field values at pole prior to B-field update
    subroutine PoleE(Model,Gr,E)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), dimension(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM), intent(inout) :: E
        integer :: nS,nR
        
        select case (Model%Ring%GridID)
        case ("cyl")
            !Cylindrical
            !$OMP PARALLEL DO default(shared)
            do nS=Model%Ring%nSi,Model%Ring%nSe+1
                E(Gr%is,:,nS,JDIR) = 0.0 !Pole direction
                E(Gr%is,:,nS,KDIR) = sum(E(Gr%is,Gr%js:Gr%je,nS,KDIR))/Model%Ring%Np
                
            enddo
        case ("lfm")
            !Do LFM +/- pole
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(nS,nR)
            do nS=Model%Ring%nSi,Model%Ring%nSe+1
                !S pole
                if (Model%Ring%doS) then
                    E(nS,Gr%js,:,KDIR) = 0.0
                    E(nS,Gr%js,:,IDIR) = sum(E(nS,Gr%js,Gr%ks:Gr%ke,IDIR))/Model%Ring%Np

                endif
                if (Model%Ring%doE) then                
                    !E pole
                    E(nS,Gr%je+1,:,KDIR) = 0.0
                    E(nS,Gr%je+1,:,IDIR) = sum(E(nS,Gr%je+1,Gr%ks:Gr%ke,IDIR))/Model%Ring%Np
                endif

            enddo
            call FixEFieldLFM(Model,Gr,E)
        end select
        
    end subroutine PoleE

    !Given variaable defined over ring, remove 0-NFT modes and save coefficients
    subroutine CleanRing(Model,Q,cFT,numFT)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: Q(Model%Ring%Np)
        real(rp), intent(out) :: cFT(0:NFTMAX,FTCOS:FTSIN)
        integer, intent(in) :: numFT

        integer :: n,m
        real(rp) :: aScl,dp,phi

        !DIR$ ASSUME_ALIGNED Q: ALIGN

        cFT = 0.0
        aScl = 2.0/Model%Ring%Np

        !Calculate coefficients
        dp = 2.0*pi/Model%Ring%Np
        do n=1,Model%Ring%Np
            phi = dp*n-0.5*dp
            do m=0,numFT
                cFT(m,FTCOS) = cFT(m,FTCOS) + Q(n)*cos(1.0*m*phi)
                cFT(m,FTSIN) = cFT(m,FTSIN) + Q(n)*sin(1.0*m*phi)
            enddo
        enddo

        !Scale
        cFT(0,:)     = 0.5*aScl*cFT(0,:)
        if (numFT>0) then
            cFT(1:numFT,:) = 1.0*aScl*cFT(1:numFT,:)
        endif

        !Subtract modes
        do n=1,Model%Ring%Np
            phi = dp*n-0.5*dp
            do m=0,numFT
                Q(n) = Q(n) - cFT(m,FTCOS)*cos(1.0*m*phi) - cFT(m,FTSIN)*sin(1.0*m*phi)
            enddo !Loop over modes
        enddo !Loop over cells

    end subroutine CleanRing

    !Given variaable defined over ring, return 0-NFT modes
    subroutine DirtyRing(Model,Q,cFT,numFT)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: Q(Model%Ring%Np)
        real(rp), intent(in) :: cFT(0:NFTMAX,FTCOS:FTSIN)
        integer, intent(in) :: numFT

        integer :: n,m
        real(rp) :: dp,phi

        !DIR$ ASSUME_ALIGNED Q: ALIGN

        dp = 2.0*pi/Model%Ring%Np
        do n=1,Model%Ring%Np
            phi = dp*n-0.5*dp
            do m=0,numFT
                Q(n) = Q(n) + cFT(m,FTCOS)*cos(1.0*m*phi) + cFT(m,FTSIN)*sin(1.0*m*phi)
            enddo !Loop over modes
        enddo !Loop over cells

    end subroutine DirtyRing

!Clean and dirty w/ scaling
    !ie for momentum, want to remove modes from velocity but not density
    subroutine CleanRingWgt(Model,Q,cFT,w,isG,numFT)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: Q(Model%Ring%Np)
        real(rp), intent(out)   :: cFT(0:NFTMAX,FTCOS:FTSIN)
        real(rp), intent(in)    :: w(Model%Ring%Np)
        logical , intent(in)    :: isG(Model%Ring%Np)
        integer, intent(in) :: numFT
        integer :: n
        !DIR$ ASSUME_ALIGNED Q: ALIGN
        !DIR$ ASSUME_ALIGNED w: ALIGN

        cFT = 0.0
        if (.not. all(isG)) return !Don't clean if some cells in this ring are vacuum
        !If still here convert to Q/w and clean that
        do n=1,Model%Ring%Np
            Q(n) = Q(n)/w(n)
        enddo

        call CleanRing(Model,Q,cFT,numFT) !Done in place

        !Now go back to Q*w
        do n=1,Model%Ring%Np
            Q(n) = Q(n)*w(n)
        enddo

    end subroutine CleanRingWgt

    subroutine DirtyRingWgt(Model,Q,cFT,w,isG,numFT)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: Q(Model%Ring%Np)
        real(rp), intent(in)    :: cFT(0:NFTMAX,FTCOS:FTSIN)
        real(rp), intent(in)    :: w(Model%Ring%Np)
        logical , intent(in)    :: isG(Model%Ring%Np)
        integer, intent(in) :: numFT
        integer :: n
        !DIR$ ASSUME_ALIGNED Q: ALIGN
        !DIR$ ASSUME_ALIGNED w: ALIGN

        if (.not. all(isG)) return !Don't clean if some cells in this ring are vacuum
        !If still here convert to Q/w and dirty that
        do n=1,Model%Ring%Np
            Q(n) = Q(n)/w(n)
        enddo
        call DirtyRing(Model,Q,cFT,numFT)

        !Now go back to Q*w
        do n=1,Model%Ring%Np
            Q(n) = Q(n)*w(n)
        enddo

    end subroutine DirtyRingWgt

!Initialization stuff
    subroutine InitRAVec(Model,Gr,A)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: A

        logical :: doInit
        integer, dimension(4) :: flxDims
        integer :: dI

        if (.not. allocated(A)) then
            doInit = .true.
        else
            flxDims = [Gr%Ni,Gr%Nj,Gr%Nk,NDIM]
            dI = sum(abs(flxDims-shape(A)))
            if (dI>0) then
                doInit = .true.
            else
                doInit = .false.
            endif
        endif

        if (doInit) then
            if(allocated(A)) deallocate(A)
            allocate(A(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM))
            A = 0.0
        endif

    end subroutine InitRAVec

end module ringav
