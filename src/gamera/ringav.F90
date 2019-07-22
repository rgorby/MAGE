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
    use ringrecon
    use multifluid
    
    implicit none

    logical, parameter, private :: doCleanLoop = .true. !Whether to remove magnetic field loops
    !Information for Fourier reductions
    integer, parameter, private :: NFT = 1 !Number of Fourier modes (beyond 0th) to remove from signed quantities
    
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

        real(rp), dimension(:,:,:,:), allocatable :: SmoothE !Smoothing electric field
        
        !DIR$ attributes align : ALIGN :: SmoothE
    !----
    !Initialization
        !Start by making sure we should be here
        if (.not. (Model%Ring%doS .or. Model%Ring%doE)) then
            return
        endif
        call Tic("Ring-Init")
        !Initialize array for smoothing field
        allocate(SmoothE(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM))
        SmoothE = 0.0
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
        call bFlux2Fld(Model,Gr,State%magFlux,State%Bxyz) !Recalculate cell-centered fields
        !Clean field loops
        if (doCleanLoop) call CleanLoops(Model,Gr,State)
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
        
        real(rp), dimension(Np,NVAR) :: rW,raW
        real(rp), dimension(Np,NDIM) :: rdB,rB0,rB
        real(rp), dimension(0:NFT,FTCOS:FTSIN,NDIM) :: cFT !Holds Fourier coefficients
        logical, dimension(Np) :: gW

        !DIR$ attributes align : ALIGN :: rW,raW,rdB,rB0,rB,gW

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
            !$OMP private(rW,raW,rdB,rB0,rB,gW,cFT)
            do nS=Model%Ring%nSi,Model%Ring%nSe
                do nR=1,Model%Ring%NumR
                    NumCh = Model%Ring%Nch(nR) !Number of total chunks for this ring
                    dC = Np/NumCh !# of cells per chunk for this ring

                    !Pull ring hydro vars
                    call PullRingCC(Model,Gr,State%Gas(:,:,:,1:NVAR,s),rW,nR,nS,NVAR,XPOLE)

                    !Pull ring cell-centered fields
                    if (Model%doMHD) then
                        call PullRingCC(Model,Gr,State%Bxyz(:,:,:,1:NDIM),rdB,nR,nS,NDIM,XPOLE)
                        if (Model%doBackground) then
                            call PullRingCC(Model,Gr,Gr%B0(:,:,:,1:NDIM),rB0,nR,nS,NDIM,XPOLE)
                        else
                            rB0 = 0.0
                        endif
                    else
                        !Just zero out fields
                        rdB = 0.0
                        rB0 = 0.0
                    endif
                    rB = rdB+rB0

                    !Convert to ring variables
                    call Gas2Ring(Model,rW,rB)

                    !Clean (subtract modes) momentum (xyz) over ring
                    do nv=1,NDIM
                        n = MOMX+nv-1 !Index in nvar, MOMX:MOMZ
                        call CleanRing(rW(:,n),cFT(:,:,nv))
                    enddo !nv loop

                !Create chunk averaged array
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
                    !Decide on good cells
                    if (Model%doMultiF) then
                        gW = (rW(:,DEN) > Spcs(s)%dVac)
                    else
                        gW = .true.
                    endif

                    do nv=1,NVAR
                        call ReconstructRing(Model,raW(:,nv),NumCh,gW)
                    enddo !nv loop

                !Dirty rings and put back
                    do nv=1,NDIM
                        n = MOMX+nv-1
                        call DirtyRing(raW(:,n),cFT(:,:,nv))
                    enddo

                    !Convert back to gas variables
                    call Ring2Gas(Model,raW,rB)

                    !Push ring back
                    call PushRingCC(Model,Gr,State%Gas(:,:,:,1:NVAR,s),raW,nR,nS,NVAR,XPOLE)
                enddo !Ring loop
            enddo !Slice loop

        enddo !Species loop

        if (Model%doMultiF) then
            !Doing ALL cells, not just ring-avg ones (meh)
            call State2Bulk(Model,Gr,State)
        endif

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

        real(rp), dimension(Np,NDIM) :: bFlx,dEr
        real(rp), dimension(Np) :: rFlx,raFlx,dFlx
        real(rp), dimension(0:NFT,FTCOS:FTSIN) :: cFT !Coefficients for Fourier reduction

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
                dC = Np/NumCh !# of cells per chunk for this ring

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
                    call CleanRing(rFlx,cFT)

                    do n=1,NumCh
                        !Find indices for this chunk in the global ring
                        lRs = 1 + (n-1)*dC
                        lRe = lRs + dC - 1

                        raFlx(lRs:lRe) = sum(rFlx(lRs:lRe))/dC

                    enddo !Chunk loop

                    !Do reconstruction on filtered/chunk-averaged fluxes
                    call ReconstructRing(Model,raFlx,NumCh)

                    !Add back in Fourier modes
                    call DirtyRing(raFlx,cFT)

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
                        SmoothE(nS,Gr%js+nR  ,Gr%ks:Gr%ke,IDIR) = dEr(:,IDIR)
                        SmoothE(nS,Gr%js+nR-1,Gr%ks:gr%ke,JDIR) = dEr(:,JDIR)
                        SmoothE(nS,Gr%js+nR-1,Gr%ks:gr%ke,KDIR) = dEr(:,KDIR)
                    endif

                    if (XPOLE == EPOLE) then !-X pole
                        SmoothE(nS,Gr%je-nR+1,Gr%ks:Gr%ke,IDIR:KDIR) = dEr(:,IDIR:KDIR)
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
        real(rp) :: tScl

        tScl = 1.0
        associate(bFlux=>State%magFlux)
        
        !Loop over slices
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(nS,nR,fD)
        do nS=Model%Ring%nSi,Model%Ring%nSe+1
            select case (Model%Ring%GridID)
            case ("cyl")
                !Cylindrical
                fD = JDIR
                nR = 1
                bFlux(Gr%is,Gr%js:Gr%je+1,nS,fD) = bFlux(Gr%is,Gr%js:Gr%je+1,nS,fD) - sum(bFlux(Gr%is,Gr%js:Gr%je,nS,fD))/Np
            case ("lfm")
                fD = KDIR
                if (Model%Ring%doS) then
                    nR = Gr%js
                    bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) = bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) - tScl*sum(bFlux(nS,nR,Gr%ks:Gr%ke,fD))/Np
                endif
                if (Model%Ring%doE) then
                    nR = Gr%je
                    bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) = bFlux(nS,nR,Gr%ks:Gr%ke+1,fD) - tScl*sum(bFlux(nS,nR,Gr%ks:Gr%ke,fD))/Np
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
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(nS)
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
    
end module ringav
