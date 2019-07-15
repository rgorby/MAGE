!Handle ring averaging about pole
!Notation
!nSi/nSe: Slices.  Cyl (Z), LFM (R)
!nR: Rings. Cyl (R), LFM (J)

module ringav
    use types
    use gamutils
    use fields
    use math
    use ringutils
    use multifluid
    
    implicit none

    integer, parameter :: Ng = 2 !Number of ghost zones for ring reconstruction
    integer, parameter :: NFT = 2 !Number of Fourier modes (beyond 0th) to remove from signed quantities
    
    logical :: initRingMod = .true. !Do we need to initialize ring-av workspaces
    logical, parameter :: doShift = .true. !Whether to add random circular shift to ring chunking
    logical, parameter :: doRingH = .true. !Whether to thermalize kinetic energy destroyed by ring-avg

    !dE = Ring-avg E field, (Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM)
    real(rp), dimension(:,:,:,:), allocatable, private :: dE !Ring-avg E field
    logical :: doCleanLoop = .true. !Whether to remove magnetic field loops

    !RingLR_T
    !Generic reconstruction routine for ring-avg
    !5-point stencil -> L/R
    abstract interface
        subroutine RingLR_T(fm2,fm1,f,fp1,fp2,vm,vp)
            Import :: rp
            real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
            real(rp), intent(out) :: vm,vp
        end subroutine RingLR_T
    end interface

    !Set choice of ring reconstruction here
    procedure(RingLR_T), pointer :: RingLR

    !Enumerators for Fourier reduction coefficients
    enum, bind(C)
        enumerator :: FTCOS=1, FTSIN
    endenum

    contains    

    !Takes State variable, averages about pole
    subroutine RingAvgHydro(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), target, intent(inout) :: State

        integer :: nR,nS,nv,n
        integer :: NumCh, dC, lRs,lRe
        integer :: s,s0,sE
        integer :: jshift,jPhase(Model%Ring%NumR)

        !Holds data for ring (2 for sph/LFM case)
        real(rp), dimension(Np,NVAR) :: rWm,rWp,raWm,raWp
        !Holds data for "good"
        logical, dimension(Np) :: gWm,gWp

        !Holds initial kinetic energy that can be thermalized
        real(rp), dimension(Np) :: mKin0,pKin0,mKin1,pKin1
        real(rp) :: dKm,dKp
        real(rp), dimension(:,:,:,:), contiguous, pointer :: W
        !Holds mode information
        real(rp), dimension(0:NFT,FTCOS:FTSIN,NDIM) :: cFTm,cFTp !Coefficients for Fourier reduction (M-XYZ)

        !DIR$ attributes align : ALIGN :: rWm,rWp,raWm,raWp
        !DIR$ attributes align : ALIGN :: gWm,gWp
        !DIR$ attributes align : ALIGN :: mKin0,pKin0,mKin1,pKin1

        jPhase = 0
        if (doShift) then
            !Get random shift to apply to rings before chunking
            do n=1,Model%Ring%NumR
                jPhase(n) = nint( genRand(0.0_rp,1.0_rp*Np) )
            enddo
        endif

        s0 = BLK !0
        if (Model%doMultiF) then
            s0 = 1
            sE = Model%nSpc
        else
            sE = BLK
        endif
        do s=s0,sE
            W => State%Gas(:,:,:,:,s)

            !Loop over slices, rings/slice
            !Only do OMP on slices
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(nR,nv,n,lRs,lRe,NumCh,dC) &
            !$OMP private(rWm,rWp,raWm,raWp,gWp,gWm,jshift) &
            !$OMP private(mKin0,pKin0,mKin1,pKin1,dKm,dKp) &
            !$OMP private(cFTm,cFTp)
            do nS=Model%Ring%nSi,Model%Ring%nSe
                do nR=1,Model%Ring%NumR
                    NumCh = Model%Ring%Nch(nR) !Number of total chunks for this ring
                    dC = Np/NumCh !# of cells per chunk for this ring

                    !Pull ring/s about pole
                    call PullRing(Model,Grid,W,rWm,rWp,nR,nS,NVAR)

                    if (doShift) then
                        !Do positive shifts on pulled rings
                        !Different phase on each radius
                        jshift = jPhase(nR)
                        rWm = cshift(rWm,+jshift,dim=1)
                        rWp = cshift(rWp,+jshift,dim=1)
                    endif

                    !Convert to ring av variables
                    call Con2Ring(Model,rWm)
                    call Con2Ring(Model,rWp)

                    !Get initial kinetic energy over ring
                    if (doRingH) then
                        do n=1,Np
                            mKin0(n) = 0.5*( rWm(n,MOMX)**2.0 + rWm(n,MOMY)**2.0 + rWm(n,MOMZ)**2.0 )/rWm(n,DEN)
                            pKin0(n) = 0.5*( rWp(n,MOMX)**2.0 + rWp(n,MOMY)**2.0 + rWp(n,MOMZ)**2.0 )/rWp(n,DEN)
                        enddo
                    endif

                    !Clean (subtract modes) momentum (xyz) over ring
                    do nv=1,NDIM
                        n = MOMX + nv - 1 !Index in nvar, MOMX:MOMZ
                        call CleanRing(rWm(:,n),cFTm(:,:,nv))
                        call CleanRing(rWp(:,n),cFTp(:,:,nv))
                    enddo

                    !Create chunk averaged array
                    !Loop over chunks, pull average over each chunk
                    do nv=1,NVAR
                        do n=1,NumCh
                            !Index bounds for this chunk in 1:Np ring
                            lRs = 1 + (n-1)*dC
                            lRe = lRs + dC - 1
                            raWp(lRs:lRe,nv) = sum(rWp(lRs:lRe,nv))/dC
                            raWm(lRs:lRe,nv) = sum(rWm(lRs:lRe,nv))/dC
                        enddo !Chunk loop
                    enddo !Var loop

                    !Now reconstruct ring/s
                    if (Model%doMultiF) then
                        gWp = (rWp(:,DEN) > Spcs(s)%dFloor)
                        gWm = (rWm(:,DEN) > Spcs(s)%dFloor)
                    else
                        gWp = .true.
                        gWm = .true.
                    endif

                    !NOTE: Reconstructing both, but not necessary in cyl
                    do nv=1,NVAR
                        call ReconstructRing(Model,raWp(:,nv),NumCh,gWp)
                        call ReconstructRing(Model,raWm(:,nv),NumCh,gWm)
                    enddo

                    !Reconstruction done, now dirty (put back modes) into RECONSTRUCTED ring
                    do nv=1,NDIM
                        n = MOMX + nv - 1 !Index in nvar, MOMX:MOMZ
                        call DirtyRing(raWm(:,n),cFTm(:,:,nv))
                        call DirtyRing(raWp(:,n),cFTp(:,:,nv))
                    enddo

                    if (doRingH) then
                        !Get post-RA kinetic energy over ring
                        do n=1,Np
                            mKin1(n) = 0.5*( raWm(n,MOMX)**2.0 + raWm(n,MOMY)**2.0 + raWm(n,MOMZ)**2.0 )/raWm(n,DEN)
                            pKin1(n) = 0.5*( raWp(n,MOMX)**2.0 + raWp(n,MOMY)**2.0 + raWp(n,MOMZ)**2.0 )/raWp(n,DEN)
                        enddo

                        !Now loop over chunks and thermalize lost kinetic energy
                        do n=1,NumCh
                            !Index bounds for this chunk in 1:Np ring
                            lRs = 1 + (n-1)*dC
                            lRe = lRs + dC - 1
                            dKm = sum(mKin0(lRs:lRe)) - sum(mKin1(lRs:lRe))
                            dKp = sum(pKin0(lRs:lRe)) - sum(pKin1(lRs:lRe))
                            if (dKm>0) then
                                !Uniformly spread lost energy over chunk
                                raWm(lRs:lRe,ENERGY) = raWm(lRs:lRe,ENERGY) + dKm/dC
                            endif
                            if (dKp>0) then
                                !Uniformly spread lost energy over chunk
                                raWp(lRs:lRe,ENERGY) = raWp(lRs:lRe,ENERGY) + dKp/dC
                            endif                                
                        enddo !Chunk loop
                    endif !Ring thermalization

                    !Convert ring vars back to con
                    call Ring2Con(Model,raWm)
                    call Ring2Con(Model,raWp)

                    if (doShift) then
                        !Unshift before pushing back into state
                        !Doing negative of original shift
                        !Working on raWx not rWx
                        jshift = jPhase(nR)
                        raWm = cshift(raWm,-jshift,dim=1)
                        raWp = cshift(raWp,-jshift,dim=1)
                    endif
                    
                    !Push ring/s back
                    call PushRing(Model,Grid,W,raWm,raWp,nR,nS,NVAR)
                enddo !Ring/Slice loop
            enddo !Slice loop
        enddo !Species loop

        if (Model%doMultiF) then
            !Doing ALL cells, not just ring-avg ones (meh)
            call State2Bulk(Model,Grid,State)
        endif

    end subroutine RingAvgHydro

    !Set electric field values at pole prior to B-field update
    subroutine PoleE(Model,Gr,E)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), dimension(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM), intent(inout) :: E
        integer :: nS

        select case (Model%Ring%GridID)
        case ("cyl")
            !Cylindrical
            !$OMP PARALLEL DO default(shared)
            do nS=Model%Ring%nSi,Model%Ring%nSe+1
                E(Gr%is,:,nS,JDIR) = 0.0 !Pole direction
                E(Gr%is,:,nS,KDIR) = sum(E(Gr%is,Gr%js:Gr%je,nS,KDIR))/Np
                
            enddo
        case ("lfm")
            !Do LFM +/- pole
            !$OMP PARALLEL DO default(shared)            
            do nS=Model%Ring%nSi,Model%Ring%nSe+1
                !S pole
                if (Model%Ring%doS) then
                    E(nS,Gr%js,:,KDIR) = 0.0
                    E(nS,Gr%js,:,IDIR) = sum(E(nS,Gr%js,Gr%ks:Gr%ke,IDIR))/Np
                endif
                if (Model%Ring%doE) then                
                    !E pole
                    E(nS,Gr%je+1,:,KDIR) = 0.0
                    E(nS,Gr%je+1,:,IDIR) = sum(E(nS,Gr%je+1,Gr%ks:Gr%ke,IDIR))/Np
                endif                
            enddo
        end select
        
    end subroutine PoleE

    !Takes magFlux variable, averages about pole
    subroutine RingAvgMag(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), target, intent(inout) :: State

        integer :: nR,nS,d,n,lRs,lRe,m
        integer :: fD,eD,NumCh,dC

        integer, dimension(2) :: fDIRS,eDIRS !Face/E-field directions
        real(rp), dimension(2) :: eScls

        real(rp) :: eScl, cumFlx,tScl
        real(rp), dimension(:,:,:,:), contiguous, pointer :: bFlux

        !Holds data for ring (2 for sph/LFM case)
        real(rp), dimension(Np,NDIM) :: mFlx,pFlx, dErm,dErp
        real(rp), dimension(Np) :: rFlxM,rFlxP,raFlxM,raFlxP,dFlxM,dFlxP
        real(rp), dimension(0:NFT,FTCOS:FTSIN) :: cFTm,cFTp !Coefficients for Fourier reduction

        !DIR$ attributes align : ALIGN :: mFlx,pFlx, dErm,dErp
        !DIR$ attributes align : ALIGN :: rFlxM,rFlxP,raFlxM,raFlxP,dFlxM,dFlxP

        !Do initialization if needed
        if (initRingMod) then
            allocate(dE(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM))
            dE = 0.0
            initRingMod = .false.
        endif
        bFlux => State%magFlux
        
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

        call Tic("Smoothing-Field")
        !Loop over slices, rings/slice
        !Only do OMP on slices
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(nR,nS,d,n,lRs,lRe,m,fD,eD,NumCh,dC) &
        !$OMP private(mFlx,pFlx, cFTm,cFTp, dErm,dErp,eScl, cumFlx) &
        !$OMP private(rFlxM,rFlxP,raFlxM,raFlxP,dFlxM,dFlxP)
        do nS=Model%Ring%nSi,Model%Ring%nSe+1
            do nR=1,Model%Ring%NumR
                NumCh = Model%Ring%Nch(nR) !Number of total chunks for this ring
                dC = Np/NumCh !# of cells per chunk for this ring

                !Pull ring/s about pole, all flux components

                call PullRing(Model,Gr,bFlux,mFlx,pFlx,nR,nS,NDIM,.true.)                

                !Loop over dimensions (2 E field directions)
                dErm = 0.0
                dErp = 0.0
                do d=1,2
                    fD   = fDIRS(d)
                    eD   = eDIRS(d)
                    eScl = eScls(d)

                    !Pull this flux direction
                    rFlxM = mFlx(:,fD)
                    rFlxP = pFlx(:,fD)

                    
                    !Clean ring and then chunk-average
                    call CleanRing(rFlxM,cFTm)
                    call CleanRing(rFlxP,cFTp)

                    do n=1,NumCh
                        !Find indices for this chunk in the global ring
                        lRs = 1 + (n-1)*dC
                        lRe = lRs + dC - 1

                        raFlxM(lRs:lRe) = sum(rFlxM(lRs:lRe))/dC
                        raFlxP(lRs:lRe) = sum(rFlxP(lRs:lRe))/dC
                    enddo !Chunk

                    !Do reconstruction on filtered/chunk-averaged fluxes
                    call ReconstructRing(Model,raFlxM,NumCh)
                    call ReconstructRing(Model,raFlxP,NumCh)

                    !Add back in lowest modes
                    call DirtyRing(raFlxM,cFTm)
                    call DirtyRing(raFlxP,cFTp)

                    !Calculate changes in flux
                    dFlxM = raFlxM - mFlx(:,fD)
                    dFlxP = raFlxP - pFlx(:,fD)

                    !Accumulate into ring E field
                    !Loop over chunks and within chunks
                    do n=1,NumCh
                        !Find indices for this chunk in the global ring
                        lRs = 1 + (n-1)*dC
                        lRe = lRs + dC - 1

                        !P side
                        cumFlx = -dFlxP(lRs)
                        do m=2,dC
                            dErp(lRs+m-1,eD) = eScl*cumFlx
                            cumFlx = cumFlx - dFlxP(lRs+m-1)
                        enddo !Intra-chunk loop

                        !M side
                        cumFlx = -dFlxM(lRs)
                        do m=2,dC
                            dErm(lRs+m-1,eD) = eScl*cumFlx
                            cumFlx = cumFlx - dFlxM(lRs+m-1)
                        enddo !Intra-chunk loop
                        !write(*,*) 'cumFlx = ', cumFlx
                    enddo !Chunk loop

                enddo !Directions

                !Ring electric field complete, push into overall E field
                select case (Model%Ring%GridID)
                case ("cyl")
                    !Cylindrical
                    !Ek starts at is+1, rest at is
                    dE(Gr%is+nR  ,Gr%js:Gr%je,nS,KDIR) = dErm(:,KDIR)
                    dE(Gr%is+nR-1,Gr%js:Gr%je,nS,IDIR) = dErm(:,IDIR)
                    dE(Gr%is+nR-1,Gr%js:Gr%je,nS,JDIR) = dErm(:,JDIR)
                    
                case ("lfm")
                    !LFM
                    !+X pole
                    if (Model%Ring%doS) then
                        !Ei starts at js+1, rest at js (K-field is 0)
                        dE(nS,Gr%js+nR  ,Gr%ks:Gr%ke,IDIR) = dErp(:,IDIR)
                        dE(nS,Gr%js+nR-1,Gr%ks:Gr%ke,JDIR:KDIR) = dErp(:,JDIR:KDIR)
                    endif
                    if (Model%Ring%doE) then
                        !-X pole
                        dE(nS,Gr%je-nR+1,Gr%ks:Gr%ke,IDIR:KDIR) = dErm(:,IDIR:KDIR)
                    endif
                end select

            enddo !Rings
        enddo !Slices, done dE calculation
        call Toc("Smoothing-Field")

        call Tic("Ring-CT")
        call PoleE(Model,Gr,dE)

        !Use dE to do a CT update (fake dt=1)
        call E2Flux(Model,Gr,State%magFlux,dE)
        call Toc("Ring-CT")

        !Finally, correct for lingering field loops
        call Tic("Loop-Cleaning")
        if ( doCleanLoop .and. (Model%Ring%NumR > 0) ) then
            !Loop over slices
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(nR,nS,fD,NumCh,tScl) 
            do nS=Model%Ring%nSi,Model%Ring%nSe+1
                nR = 1 !Only need first ring
                select case (Model%Ring%GridID)
                case ("cyl")
                    !Cylindrical
                    fD = JDIR
                    bFlux(Gr%is,Gr%js:Gr%je+1,nS,fD) = bFlux(Gr%is,Gr%js:Gr%je+1,nS,fD) - sum(bFlux(Gr%is,Gr%js:Gr%je,nS,fD))/Np
                case ("lfm")
                    !LFM
                    fD = KDIR
                    !Diffuse net flux over ring with timescale NumCh timesteps
                    !Instantaneous diffusion of ring-flux is tScl=1.0
                    NumCh = Model%Ring%Nch(nR) !Number of total chunks for this ring
                    !tScl = 1.0/NumCh
                    tScl = 1.0
                    if (Model%Ring%doS) then
                        nR = Gr%js
                        bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) = bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) - tScl*sum(bFlux(nS,nR,Gr%ks:Gr%ke,fD))/Np             
                    endif
                    if (Model%Ring%doE) then
                        nR = Gr%je
                        bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) = bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) - tScl*sum(bFlux(nS,nR,Gr%ks:Gr%ke,fD))/Np
                    endif

                end select
            enddo !Loop over slices
            
        endif
        call Toc("Loop-Cleaning")
        
    end subroutine RingAvgMag


    !Reconstruct chunked ring values
    !m indices are over ring
    !n indices are over chunks
    !isGO (optional), which cells are good
    subroutine ReconstructRing(Model,rW,Nc,isGO)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: rW(Np)
        integer, intent(in) :: Nc !Number of chunks in ring
        logical, intent(in), optional :: isGO(Np)

        !Hold expanded ring (ie ghosts)
        real(rp) :: chW(1-Ng:Nc+Ng)
        real(rp) :: fm2,fm1,f,fp1,fp2,fL,fR
        real(rp) :: a,b,c,fI
        integer :: dJ,m,n,nv, mS,mE,ngood
        logical :: isG(Np)

        !DIR$ ASSUME_ALIGNED rW: ALIGN

        if (present(isGO)) then
            isG = isGO
        else
            isG = .true.
        endif

        dJ = Np/Nc

        !Fill inner region and add ghosts
        chW(1:Nc) = rW(1:Np:dJ)
        do n=1,Ng
            chW(1-n)  = chW(Nc-n+1)
            chW(Nc+n) = chW(n)
        enddo

        !Loop over each chunk, create interpolant for chunk
        do n=1,Nc
            !Indices of this chunk within 1,Np array
            mS = 1 + (n-1)*dJ
            mE = mS + dJ - 1
            if (Model%doMultiF) then
                ngood = count( isG(mS:mE) ) !Number of good elements in chunk
                if (ngood < dJ) then
                    cycle !Skip this chunk 
                endif
            endif

            !Grab stencil for interval LR's
            fm2 = chW(n-2)
            fm1 = chW(n-1)
            f   = chW(n  )
            fp1 = chW(n+1)
            fp2 = chW(n+2)

            call RingLR(fm2,fm1,f,fp1,fp2,fL,fR)
            
            !Calculate coefficients for parabolic interpolant
            !TODO: Check this against equation 2 of Bin's ring avg paper
            a = 3*(fL + fR - 2*f)
            b = 2*(3*f - fR - 2*fL)
            c = fL

            !Reconstruct within this chunk
            do m=1,dJ
                fI = (a/3.0)*(3*m*m-3*m+1)/(dJ*dJ) + 0.5*b*(2*m-1)/dJ + c
                rW(mS+m-1) = fI
            enddo
        enddo !Chunks

    end subroutine ReconstructRing

    !Given variaable defined over ring, remove 0-NFT modes and save coefficients
    subroutine CleanRing(Q,cFT)
        real(rp), intent(inout) :: Q(Np)
        real(rp), intent(out) :: cFT(0:NFT,FTCOS:FTSIN)

        integer :: n,m
        real(rp) :: aScl,dp,phi

        !DIR$ ASSUME_ALIGNED Q: ALIGN

        cFT = 0.0
        aScl = 2.0/Np

        !Calculate coefficients
        dp = 2.0*pi/Np
        do n=1,Np
            phi = dp*n-0.5*dp
            do m=0,NFT
                cFT(m,FTCOS) = cFT(m,FTCOS) + Q(n)*cos(1.0*m*phi)
                cFT(m,FTSIN) = cFT(m,FTSIN) + Q(n)*sin(1.0*m*phi)
            enddo
        enddo

        !Scale
        cFT(0,:)     = 0.5*aScl*cFT(0,:)
        cFT(1:NFT,:) = 1.0*aScl*cFT(1:NFT,:)

        !Subtract modes
        do n=1,Np
            phi = dp*n-0.5*dp
            do m=0,NFT
                Q(n) = Q(n) - cFT(m,FTCOS)*cos(1.0*m*phi) - cFT(m,FTSIN)*sin(1.0*m*phi)
            enddo !Loop over modes
        enddo !Loop over cells

    end subroutine CleanRing

    !Given variaable defined over ring, return 0-NFT modes
    subroutine DirtyRing(Q,cFT)
        real(rp), intent(inout) :: Q(Np)
        real(rp), intent(in) :: cFT(0:NFT,FTCOS:FTSIN)

        integer :: n,m
        real(rp) :: dp,phi

        !DIR$ ASSUME_ALIGNED Q: ALIGN

        dp = 2.0*pi/Np
        do n=1,Np
            phi = dp*n-0.5*dp
            do m=0,NFT
                Q(n) = Q(n) + cFT(m,FTCOS)*cos(1.0*m*phi) + cFT(m,FTSIN)*sin(1.0*m*phi)
            enddo !Loop over modes
        enddo !Loop over cells

    end subroutine DirtyRing

    !Lazy routine to make things equivalent to PCM
    !Both L & R are f
    subroutine PCM_IntLR(fm2,fm1,f,fp1,fp2,vm,vp)
        real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
        real(rp), intent(out) :: vm,vp

        vm = f
        vp = f

    end subroutine PCM_IntLR

    !Lazy routine to make things equivalent to PLM
    !Need that vp+vm=2*f to zero out quadratic term

    subroutine PLM_IntLR(fm2,fm1,f,fp1,fp2,vm,vp)
        real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
        real(rp), intent(out) :: vm,vp

        real(rp) :: dwL,dwR,dwC,dwM
        real(rp) :: min1,min2

        dwL = f - fm1
        dwR = fp1 - f
        dwC = ( fp1 - fm1 )/2.0

        min1 = min(2*abs(dwL),2*abs(dwR))
        min2 = min(min1,abs(dwC))
        !SIGN(A,B) returns the value of A with the sign of B
        dwM = sign(min2,dwC)

        vp = f + 0.5*dwM
        vm = f - 0.5*dwM

    end subroutine PLM_IntLR

    subroutine WENO_IntLR(fm2,fm1,f,fp1,fp2,vm,vp)
        real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
        real(rp), intent(out) :: vm,vp

        vm = Weno(fp2,fp1,f,fm1,fm2)
        vp = Weno(fm2,fm1,f,fp1,fp2)

    end subroutine WENO_IntLR

    !Calculates PPM *interval* L/R
    subroutine PPM_IntLR(fm2,fm1,f,fp1,fp2,vm,vp)
        real(rp), intent(in) :: fm2,fm1,f,fp1,fp2
        real(rp), intent(out) :: vm,vp

        real(rp) :: dvm1,dvm2,dvp1,dvp2
        real(rp) :: SM,Sm1,Sp1,Sp2
        real(rp) :: dvc,am,ap

        dvm2 = fm1 - fm2
        dvm1 = f   - fm1
        dvp1 = fp1 - f
        dvp2 = fp2 - fp1
 
        dvc = 0.5*(dvm2 + dvm1)
        SM  = 2.0*minmod(dvm2, dvm1)
        Sm1 = minmod(dvc, SM)
 
        dvc = 0.5*(dvm1 + dvp1)
        SM  = 2.0*minmod(dvm1, dvp1)
        Sp1 = minmod(dvc, SM)
 
        dvc = 0.5*(dvp2 + dvp1)
        SM  = 2.0*minmod(dvp2, dvp1)
        Sp2 = minmod(dvc, SM)
 
        vp = 0.5*(f + fp1) - (Sp2 - Sp1)/6.0
        vm = 0.5*(f + fm1) - (Sp1 - Sm1)/6.0
 
        ap = vp - f
        am = vm - f
 
        if (ap*am >= 0.0) then
            ap = 0.0
            am = 0.0
        else 
            if (abs(ap) >= 2.0*abs(am)) then
                ap = -2.0*am
            endif

            if (abs(am) >= 2.0*abs(ap)) then
                am = -2.0*ap
            endif
        endif
 
        vp = f + ap
        vm = f + am

    end subroutine PPM_IntLR

    !WENO-5 reconstruction
    function Weno(hm2,hm1,h,hp1,hp2) result(hI)
        real(rp), intent(in) :: hm2,hm1,h,hp1,hp2
        real(rp) :: hI

        real(rp) :: h0,h1,h2, is0,is1,is2
        real(rp) :: a0,a1,a2, w0,w1,w2

        h0 =  (1.0/3)*hm2 - (7.0/6)*hm1 + (11.0/6)*h
        h1 = -(1.0/6)*hm1 + (5.0/6)*h   +  (1.0/3)*hp1
        h2 =  (1.0/3)*h   + (5.0/6)*hp1 -  (1.0/6)*hp2

        is0 = (13.0/12)*(hm2 - 2.0*hm1 +   h)**2.0 + (1.0/4)*(hm2 - 4.0*hm1 + 3.0*h)**2.0
        is1 = (13.0/12)*(hm1 - 2.0*h   + hp1)**2.0 + (1.0/4)*(hm1 - hp1)**2.0
        is2 = (13.0/12)*(h   - 2.0*hp1 + hp2)**2.0 + (1.0/4)*(3.0*h - 4.0*hp1 + hp2)**2.0

        a0 = (1.0/10)/(is0+TINY)**2.0
        a1 = (6.0/10)/(is1+TINY)**2.0
        a2 = (3.0/10)/(is2+TINY)**2.0

        w0 = a0/(a0+a1+a2)
        w1 = a1/(a0+a1+a2)
        w2 = a2/(a0+a1+a2)
 
        hI = w0*h0 + w1*h1 + w2*h2
    end function Weno
    
end module ringav
