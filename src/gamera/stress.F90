module stress
    use gamtypes
    use clocks
    use gamutils
    use gridutils
    use recon
    use ringav
    use ringutils
    use multifluid
    
    implicit none

    !Enumerators for L/R sides
    enum, bind(C)
        enumerator :: LEFT=1, RIGHT
    endenum
    
    !Signs for left/right going fluxes @ interfaces
    integer, parameter, dimension(2), private :: SgnLR=[-1,1]

    logical, parameter, private :: doRingFlux = .true.
    logical, parameter, private :: doNuke = .true. !Do nuclear option
    logical, parameter, private :: doHogs11 = .true. !Do // magnetic hogs diffusion

    !cLim: Vile magic number, when to apply nuclear option (v>cLim*Ca)
    !LFM uses 1.5
    real(rp), parameter, private :: cLim = 1.5

    contains

!Calculates plasma deltas (flux updates) based on Reynolds/Maxwell stresses
    !Note: Allocating MHD arrays for hydro, but unused
    !Overview:
    !Perform initialization
    !Loop over flux directions, dimensions and call Flux routine on i-blocks for vectorization
    !Loop-ordering is direction dep. for cache optimization
        !Calculate flux at "vecLen" interfaces at a time, (iBlk:iBlk+vecLen,j,k)
        !Vary loop order to increase cache locality
        !i-fluxes: (k,j) & (iBlk,i)
        !j-fluxes: (k,iBlk) & (j,i)
        !k-fluxes: (j,iBlk) & (k,i)
        !Inner loop is always i=1,vecLen/iMax for vectorization
    subroutine CalcStress(Model,Gr,State,gFlx,mFlx,dGasH,dGasM)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp), intent(out) :: gFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(out) :: mFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM,1:NDIM)
        real(rp), intent(out) :: dGasH(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,BLK:Model%nSpc)
        real(rp), optional, intent(out) :: dGasM(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM)

        integer :: i,iB,j,k,d,nFlx
        integer :: ks,ke !For 2.5D handling
        logical :: doMaxwell
        real(rp) :: dV

        !DIR$ ASSUME_ALIGNED dGasH: ALIGN
        !DIR$ ASSUME_ALIGNED dGasM: ALIGN
        !DIR$ ASSUME_ALIGNED gFlx: ALIGN
        !DIR$ ASSUME_ALIGNED mFlx: ALIGN

    !Do init stuff
    !-------------
        !Set # of flux directions to calculate
        nFlx = NDIM
        ks = Gr%ks
        ke = Gr%ke
        if (Model%do25D) then
            nFlx = 2
            !Only do 1 k slice if doing 2.5D
            ks = Gr%ks
            ke = Gr%ks
        endif

        !Set whether to do magnetic stresses
        doMaxwell = .false.
        if ( present(dGasM) .and. Model%doMHD  ) doMaxwell = .true.

    !Main parallel region
        call Tic("Fluxes")
    !Open one big parallel region, attach to individual loops
        !$OMP PARALLEL default(shared) private(i,iB,j,k,d,dV)
        
        !Flux loop, calculate brick of fluxes for all species
        !---------------------------
        do d=1,nFlx
            !Structure loops based on direction
            select case(d)
            case(IDIR)
                !Fluxes in i direction
                !---------------------
                !$OMP DO collapse(2)
                do k=ks,ke
                    do j=Gr%js,Gr%je
                        do iB=Gr%is,Gr%ie+1,vecLen
                            !Do block of i-interface fluxes
                            call Fluxes(Model,Gr,State,iB,j,k,d,gFlx,mFlx)
                        enddo
                    enddo
                enddo !K loop
                !$OMP ENDDO NOWAIT
            case(JDIR)
                !Fluxes in j direction
                !---------------------
                !$OMP DO collapse(2)
                do k=ks,ke
                    do iB=Gr%is,Gr%ie,vecLen
                        do j=Gr%js,Gr%je+1
                            !Do block of j-interface fluxes
                            call Fluxes(Model,Gr,State,iB,j,k,d,gFlx,mFlx)
                        enddo
                    enddo
                enddo !K loop
                !$OMP ENDDO NOWAIT
            case(KDIR)
                !Fluxes in k direction
                !---------------------
                !$OMP DO collapse(2)
                do j=Gr%js,Gr%je
                    do iB=Gr%is,Gr%ie,vecLen
                        do k=ks,ke+1
                            !Do block of k-interface fluxes
                            call Fluxes(Model,Gr,State,iB,j,k,d,gFlx,mFlx)
                        enddo
                    enddo
                enddo
                !$OMP ENDDO NOWAIT
            end select !Flux direction

        enddo !Flux directions

        !Close big parallel region
        !$OMP END PARALLEL

        !Do various hacks to fluxes before conversion to deltas
        !Fix fluxes on ring if necessary
        if (doRingFlux) call RingFlux(Model,Gr,gFlx,mFlx)

        if (associated(Model%HackFlux)) then
            call Tic("HackFlux")
            call Model%HackFlux(Model,Gr,gFlx,mFlx)
            call Toc("HackFlux")
        endif

        call Toc("Fluxes")

    !Turn fluxes into deltas
    !---------------------------
        call Tic("Flux2Deltas")
                
        !$OMP PARALLEL DO default(shared) collapse (2) &
        !$OMP private(i,j,k,dV)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie
                    dV = Gr%volume(i,j,k)
                    !Do all species here
                    !Fluxes have already been scaled by face areas!
                    dGasH(i,j,k,1:NVAR,:) = (  gFlx(i,j,k,1:NVAR,IDIR,:) - gFlx(i+1,j,k,1:NVAR,IDIR,:) &
                                             + gFlx(i,j,k,1:NVAR,JDIR,:) - gFlx(i,j+1,k,1:NVAR,JDIR,:) &
                                             + gFlx(i,j,k,1:NVAR,KDIR,:) - gFlx(i,j,k+1,1:NVAR,KDIR,:) )/dV

                    if (doMaxwell) then
                        dGasM(i,j,k,XDIR:ZDIR) = (  mFlx(i,j,k,XDIR:ZDIR,IDIR) - mFlx(i+1,j,k,XDIR:ZDIR,IDIR) &
                                                  + mFlx(i,j,k,XDIR:ZDIR,JDIR) - mFlx(i,j+1,k,XDIR:ZDIR,JDIR) &
                                                  + mFlx(i,j,k,XDIR:ZDIR,KDIR) - mFlx(i,j,k+1,XDIR:ZDIR,KDIR) )/dV
                        
                        if (Model%doBackground) then
                            !Add background field force terms
                            dGasM(i,j,k,XDIR:ZDIR) = dGasM(i,j,k,XDIR:ZDIR) + Gr%dpB0(i,j,k,XDIR:ZDIR)
                        endif
                    endif !doMax
                enddo ! i loop
            enddo !J loop
        enddo !K loop
        
        call Toc("Flux2Deltas")

    end subroutine CalcStress


    !Calculate vecLen block of fluxes (iB:iB+vecLen,j,k) in direction d
    !dGasM passed but unused in hydro case
    subroutine Fluxes(Model,Gr,State,iB,j,k,dN,gFlx,mFlx)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        integer, intent(in) :: iB,j,k,dN
        real(rp), intent(inout) :: gFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(inout) :: mFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM,1:NDIM)
        

        !Block variables, local to thread
        !Reconstruction stencils
        real(rp) :: VolB(vecLen,recLen)
        real(rp) :: ConB(vecLen,recLen,NVAR)
        real(rp) :: MagB(vecLen,recLen,NDIM)
        real(rp) :: TfB(vecLen,NDIM*NDIM) !Face transform
        real(rp) :: faB(vecLen) !Face area in this direction

        !Reconstructed values, L/Rs
        real(rp) :: PrimLRB(vecLen,NVAR,2), MagLRB(vecLen,NDIM,2)
        !Misc.
        real(rp) :: gFlxB(vecLen,NVAR), mFlxB(vecLen,NDIM) !Gas/Mag-kinetic fluxes
        real(rp) :: VnB(vecLen,2) !L/R normal velocities @ interfaces
        real(rp) :: dW(vecLen,NVAR) !R-L values for HOGs flux
        real(rp) :: VaD(vecLen) !Diffusive Alfven speed @ interface for hogs
        real(rp) :: B0(vecLen,NDIM), B0n(vecLen), Bn(vecLen)
        real(rp) :: bbD(vecLen), dVel(vecLen,NDIM)
        logical :: doFlx(vecLen,2) !Whether L/R wrt interface are good

        !Scalar holders
        real(rp) :: vL,vR,Vn
        real(rp) :: RhoML,RhoMR
        real(rp), dimension(NDIM) :: deeP,bAvg,bhat

        !Indices
        integer :: ie,iMax,isB,ieB,iG !Vector direction indices
        integer :: i,nv,q,s
        logical :: isBulk,doI

        !DIR$ attributes align : ALIGN :: VolB,ConB,MagB,TfB,faB
        !DIR$ attributes align : ALIGN :: PrimLRB,MagLRB,gFlxB,mFlxB
        !DIR$ attributes align : ALIGN :: VnB,dW,VaD,B0,B0n,Bn
        !DIR$ attributes align : ALIGN :: bbD,dVel
        !DIR$ ASSUME_ALIGNED gFlx: ALIGN
        !DIR$ ASSUME_ALIGNED mFlx: ALIGN

    !Initialize bounds/zero out accumulators
    !---------------------------
        ie = Gr%ie + ijkD(dN,IDIR) !Last active interface, ie if j/k or ie+1 if i
        !Check for max if not full vecLen of data
        iMax = min(vecLen,ie-iB+1)

        !Set beginning/end i's of this block
        isB = iB
        ieB = iB+iMax-1

        !Zero out accumulators
        gFlxB = 0.0
        mFlxB = 0.0
        dW = 0.0
        VaD = 0.0
        bbD = 0.0
        dVel = 0.0
        VolB = 0.0
        MagB = 0.0
        ConB = 0.0
        TfB = 0.0
        B0 = 0.0
        Bn = 0.0
        B0n = 0.0


        doFlx(:,:) = .true.
    !Load non-species data into local work arrays
    !---------------------------
        !Get geometric information, Need recLen/2 radius about each i,j,k in brickette
        call LoadBlock(Model,Gr,VolB,Gr%volume,iB,j,k,iMax,dN)

        !Get relevant face transforms/area (no stencil needed)
        TfB(1:iMax,:) = Gr%Tf  (isB:ieB,j,k,:,dN)
        faB(1:iMax  ) = Gr%face(isB:ieB,j,k,  dN)

        !If MHD pull field info and do L/R's and interface calculations
        if (Model%doMHD) then
            do nv=1,NDIM
                call LoadBlock(Model,Gr,MagB(:,:,nv),State%Bxyz(:,:,:,nv),iB,j,k,iMax,dN)
            enddo
            call BlockLRs(VolB,MagB,MagLRB(:,:,LEFT),MagLRB(:,:,RIGHT),NDIM)

            Bn(1:iMax) = State%magFlux(isB:ieB,j,k,dN)/faB(1:iMax)
            if (Model%doBackground) then
                B0(1:iMax,:) = Gr%fcB0(isB:ieB,j,k,:,dN)
                B0n(1:iMax) = Gr%bFlux0(isB:ieB,j,k,dN)/faB(1:iMax)
            endif

        endif

    !Loop over all species (including bulk) and calculate gas flux.  Calculate mag flux when doing bulk species
    !---------------------------
        do s=BLK,Model%nSpc
            !Reset accumulators
            gFlxB = 0.0
            dW = 0.0
            dVel = 0.0

            doFlx(:,:) = .true.

            if (s == BLK) then
                isBulk = .true.
            else
                isBulk = .false.
            endif
        !Load conserved quantities for this species and get LR's
        !---------------------------
            do nv=1,NVAR
                call LoadBlock(Model,Gr,ConB(:,:,nv),State%Gas(:,:,:,nv,s),iB,j,k,iMax,dN)
            enddo
            call BlockStateLRs(Model,VolB,ConB,PrimLRB(:,:,LEFT),PrimLRB(:,:,RIGHT))
            if (.not. isBulk) then
                !For non-bulk species, test both sides of interface
                doFlx(:,:) = ( PrimLRB(:,DEN,:) >= Spcs(s)%dVac )
            endif

            !TODO: Maybe move this check up after first loadblock?
            if (.not. isBulk) then
                !For non-bulk species, test both sides of interface
                doFlx(:,:) = ( PrimLRB(:,DEN,:) >= Spcs(s)%dVac )

                !Bail out of this species if no good interfaces
                if (.not. any(doFlx)) cycle
            endif

        !Reynolds stress fluxes (all species)
        !---------------------------
            do q=1,2 !L/R directions
                call GasKinFlux(Model,q,PrimLRB(:,:,q),TfB,gFlxB,dW,dVel,VnB(:,q),doFlx(:,q))
            enddo

            if (isBulk .and. Model%doMHD) then
            !Maxwell stress fluxes (only bulk)
            !---------------------------
                !Don't worry about doFlx because bulk is always good
                do q=1,2 !L/R directions
                    call MagKinFlux(Model,q,PrimLRB(:,:,q),MagLRB(:,:,q),Bn,VnB(:,q), &
                                   & TfB,mFlxB,B0,B0n,bbD,VaD)
                enddo
            endif !Maxwell fluxes

        !HOGS diffusive terms (all species), must be after calculation of diffusive alfven speeds
        !---------------------------
            if (Model%doHogs .and. Model%doMHD) then
                !Hydro hogs
                do nv=1,NVAR
                    do i=1,iMax
                        if ( doFlx(i,LEFT) .or. doFlx(i,RIGHT) ) then
                            gFlxB(i,nv) = gFlxB(i,nv) - Model%cHogH*VaD(i)*dW(i,nv)
                        endif
                    enddo
                enddo

                !Boris hogs (only do for bulk species)
                if (Model%doBoris .and. isBulk) then
                    do i=1,iMax

                        deeP = bbD(i)*dVel(i,XDIR:ZDIR)/(Model%Ca**2.0)

                        if (.not. doHogs11) then
                        !Limit to perp component
                            !Get average XYZ field @ interface
                            bAvg = B0(i,XDIR:ZDIR) + 0.5*( MagLRB(i,XDIR:ZDIR,LEFT) + MagLRB(i,XDIR:ZDIR,RIGHT) )
                            bhat = normVec(bAvg)
                            !Pull out parallel component
                            deeP = Vec2Perp(deeP,bhat)
                        endif

                        !Now apply mag hogs
                        mFlxB(i,XDIR:ZDIR) = mFlxB(i,XDIR:ZDIR) - Model%cHogM*VaD(i)*deeP

                        ! !Calculate del(rho_m v), the magnetic momentum
                        ! !Two choices here: 
                        ! !1) Delta mag momentum or Magmass x delta-vee
                        ! !2) Total direction or perp to field

                        ! RhoML = norm2( B0(i,XDIR:ZDIR) + MagLRB(i,XDIR:ZDIR,LEFT ) )**2.0/(Model%Ca**2.0)
                        ! RhoMR = norm2( B0(i,XDIR:ZDIR) + MagLRB(i,XDIR:ZDIR,RIGHT) )**2.0/(Model%Ca**2.0)
                        ! deeP = RhoMR*PrimLRB(i,VELX:VELZ,RIGHT) - RhoML*PrimLRB(i,VELX:VELZ,LEFT)



                        ! !Now apply mag hogs
                        ! mFlxB(i,XDIR:ZDIR) = mFlxB(i,XDIR:ZDIR) - Model%cHogM*VaD(i)*deeP

                    enddo
                endif !Boris HOGS

                !Do super diffusion if interface speed is too much faster than "light" (all species)
                if (doNuke .and. Model%doBoris) then
                    !Loop over interfaces, check if average interface velocity is > Ca
                    do i=1,iMax
                        !Vn = sqrt( (vL^2 + vR^2 )/2 )
                        vL = norm2(PrimLRB(i,VELX:VELZ,LEFT))
                        vR = norm2(PrimLRB(i,VELX:VELZ,RIGHT))
                        Vn = sqrt(vL**2.0 + vR**2.0)/sqrt(2.0)
                        doI = (Vn >= cLim*Model%Ca) .and. ( doFlx(i,LEFT) .and. doFlx(i,RIGHT) )

                        if (doI) then
                            !Go nuts and add a shit-ton of diffusion
                            !Use diffusive speed of Ca, and Delta of cell values i,i-1 (instead of L/R's)
                            gFlxB(i,DEN   ) = gFlxB(i,DEN   ) - Model%Ca*( ConB(i,Nr2+1,DEN   ) - ConB(i,Nr2,DEN   ) )
                            gFlxB(i,MOMX  ) = gFlxB(i,MOMX  ) - Model%Ca*( ConB(i,Nr2+1,MOMX  ) - ConB(i,Nr2,MOMX  ) )
                            gFlxB(i,MOMY  ) = gFlxB(i,MOMY  ) - Model%Ca*( ConB(i,Nr2+1,MOMY  ) - ConB(i,Nr2,MOMY  ) )
                            gFlxB(i,MOMZ  ) = gFlxB(i,MOMZ  ) - Model%Ca*( ConB(i,Nr2+1,MOMZ  ) - ConB(i,Nr2,MOMZ  ) )
                            gFlxB(i,ENERGY) = gFlxB(i,ENERGY) - Model%Ca*( ConB(i,Nr2+1,ENERGY) - ConB(i,Nr2,ENERGY) )

                        endif !Going nuts
                    enddo
                endif !Nuclear option

            endif !HOGS

        !Store fluxes into main arrays, multiply by face area *HERE*
        !---------------------------
            do nv=1,NVAR
                do i=1,iMax
                    iG = isB+i-1 !Global index
                    gFlx(iG,j,k,nv,dN,s) = faB(i)*gFlxB(i,nv)
                enddo
            enddo

            if (isBulk .and. Model%doMHD) then
                do nv=1,NDIM
                    do i=1,iMax
                        iG = isB+i-1 !Global index
                        mFlx(iG,j,k,nv,dN) = faB(i)*mFlxB(i,nv)
                    enddo
                enddo
            endif !Mag fluxes

        !Done with this species, move on
        enddo !Main species loop

    end subroutine Fluxes

    !Calculate 1-sided Reynolds flux
    subroutine GasKinFlux(Model,q,PrimLRB,TfB,gFlxB,dW,dVel,VnB,doFlx)
        type(Model_T), intent(in) :: Model
        integer, intent(in) :: q
        real(rp), intent(in) :: PrimLRB(vecLen,NVAR), TfB(vecLen,NDIM*NDIM)
        real(rp), intent(inout) :: gFlxB(vecLen,NVAR), dW(vecLen,NVAR), dVel(vecLen,NDIM)
        real(rp), intent(out) :: VnB(vecLen)
        logical, intent(in) :: doFlx(vecLen)

        integer :: i
        real(rp) :: D,P,Vx,Vy,Vz,E,lambda
        real(rp) :: Vn,Vt1,Vt2,Vn0,Vn1,fMn,fMt1,fMt2,fMx,fMy,fMz

        !DIR$ ASSUME_ALIGNED PrimLRB: ALIGN
        !DIR$ ASSUME_ALIGNED TfB: ALIGN
        !DIR$ ASSUME_ALIGNED gFlxB: ALIGN
        !DIR$ ASSUME_ALIGNED dW: ALIGN
        !DIR$ ASSUME_ALIGNED dVel: ALIGN
        !DIR$ ASSUME_ALIGNED VnB: ALIGN

        !Bail out if none of these cells have "real" fluid in this species
        if (.not. any(doFlx)) return

        do i=1,vecLen
            if ( doFlx(i) ) then
                !Get primitive values, calculate lambda
                D = PrimLRB(i,DEN)
                P = PrimLRB(i,PRESSURE)
                Vx = PrimLRB(i,VELX)
                Vy = PrimLRB(i,VELY)
                Vz = PrimLRB(i,VELZ)
                E = 0.5*D*(Vx**2.0 + Vy**2.0 + Vz**2.0) + P/(Model%gamma-1)
                lambda = D/(2*P)

                !Rotate into face normal, save normal
                Vn  =  TfB(i,NORMX)*Vx + &
                &      TfB(i,NORMY)*Vy + &
                &      TfB(i,NORMZ)*Vz
                Vt1 =  TfB(i,TAN1X)*Vx + &
                &      TfB(i,TAN1Y)*Vy + &
                &      TfB(i,TAN1Z)*Vz
                Vt2 =  TfB(i,TAN2X)*Vx + &
                &      TfB(i,TAN2Y)*Vy + &
                &      TfB(i,TAN2Z)*Vz

                VnB(i) = Vn

                !Calculate first/second moments
                Vn0 = 0.5*erfc(SgnLR(q)*sqrt(lambda)*Vn)
                Vn1 = Vn*Vn0 - 0.5*SgnLR(q)*exp(-lambda*Vn*Vn)/sqrt(PI*lambda)

                !Calculate fluxes in face-normal system
                gFlxB(i,DEN)    = gFlxB(i,DEN) + D*Vn1
                gFlxB(i,ENERGY) = gFlxB(i,ENERGY) + (E+0.5*P)*Vn1 + 0.5*P*Vn0*Vn

                fMn  = D*Vn *Vn1 + P*Vn0
                fMt1 = D*Vt1*Vn1
                fMt2 = D*Vt2*Vn1

                !Rotate back into x,y,z system
                fMx = TfB(i,NORMX)*fMn  + &
                    & TfB(i,TAN1X)*fMt1 + &
                    & TfB(i,TAN2X)*fMt2
                fMy = TfB(i,NORMY)*fMn  + &
                    & TfB(i,TAN1Y)*fMt1 + &
                    & TfB(i,TAN2Y)*fMt2
                fMz = TfB(i,NORMZ)*fMn  + &
                    & TfB(i,TAN1Z)*fMt1 + &
                    & TfB(i,TAN2Z)*fMt2

                !Accumulate for L/R-moving fluxes
                gFlxB(i,MOMX) = gFlxB(i,MOMX) + fMx
                gFlxB(i,MOMY) = gFlxB(i,MOMY) + fMy
                gFlxB(i,MOMZ) = gFlxB(i,MOMZ) + fMz

                !Accumulate dW (R-L) for hogs
                dW(i,DEN)    = dW(i,DEN)    + SgnLR(q)*D
                dW(i,MOMX)   = dW(i,MOMX)   + SgnLR(q)*D*Vx
                dW(i,MOMY)   = dW(i,MOMY)   + SgnLR(q)*D*Vy
                dW(i,MOMZ)   = dW(i,MOMZ)   + SgnLR(q)*D*Vz
                dW(i,ENERGY) = dW(i,ENERGY) + SgnLR(q)*E

                !Accumulate dVel for mag hogs
                dVel(i,XDIR) = dVel(i,XDIR) + SgnLR(q)*Vx
                dVel(i,YDIR) = dVel(i,YDIR) + SgnLR(q)*Vy
                dVel(i,ZDIR) = dVel(i,ZDIR) + SgnLR(q)*Vz
            endif !doFlx(i)
        enddo !i loop
        
    end subroutine GasKinFlux

    subroutine MagKinFlux(Model,q,PrimLRB,MagLRB,Bn,VnB,TfB,mFlxB,B0,B0n,bbD,VaD)
        type(Model_T), intent(in) :: Model
        integer, intent(in) :: q
        real(rp), intent(in) :: PrimLRB(vecLen,NVAR), MagLRB(vecLen,NDIM), TfB(vecLen,NDIM*NDIM)
        real(rp), dimension(vecLen), intent(in) :: Bn,B0n,VnB
        real(rp), intent(in) :: B0(vecLen,NDIM)
        real(rp), intent(inout) :: mFlxB(vecLen,NDIM), bbD(vecLen), VaD(vecLen)

        integer :: i
        real(rp) :: D,P,Bx,By,Bz,dPb,Va2,Vac2,Nx,Ny,Nz
        real(rp) :: B0x,B0y,B0z,BdB0,Pb0
        real(rp) :: lambda, Vn0

        !DIR$ ASSUME_ALIGNED PrimLRB: ALIGN
        !DIR$ ASSUME_ALIGNED MagLRB: ALIGN
        !DIR$ ASSUME_ALIGNED Bn: ALIGN
        !DIR$ ASSUME_ALIGNED VnB: ALIGN
        !DIR$ ASSUME_ALIGNED TfB: ALIGN
        !DIR$ ASSUME_ALIGNED mFlxB: ALIGN
        !DIR$ ASSUME_ALIGNED B0: ALIGN
        !DIR$ ASSUME_ALIGNED B0n: ALIGN
        !DIR$ ASSUME_ALIGNED bbD: ALIGN
        !DIR$ ASSUME_ALIGNED VaD: ALIGN
        
        do i=1,vecLen
            !Calculate lambda
            D = PrimLRB(i,DEN)
            P = PrimLRB(i,PRESSURE)
            Bx = MagLRB(i,XDIR)
            By = MagLRB(i,YDIR)
            Bz = MagLRB(i,ZDIR)

            B0x = B0(i,XDIR)
            B0y = B0(i,YDIR)
            B0z = B0(i,ZDIR)

            dPb = 0.5*(Bx**2.0 + By**2.0 + Bz**2.0) !Pressure in residual field
            if (Model%doBackground) then
                Va2 = ( (Bx+B0x)**2.0 + (By+B0y)**2.0 + (Bz+B0z)**2.0 )/D
                bbD(i) = bbD(i) + 0.5*( (Bx+B0x)**2.0 + (By+B0y)**2.0 + (Bz+B0z)**2.0)
            else
                Va2 = 2*dPb/D
                bbD(i) = bbD(i) + 0.5*(Bx**2.0 + By**2.0 + Bz**2.0)
            endif
            
            if (Model%doBoris) then
                Vac2 = Va2/(1 + (Va2/Model%Ca**2.0) ) !Boris-corrected Va^2
            else
                Vac2 = Va2
            endif

            !Calculate 1st moment
            !lambda = D/2*P_tot, Va2 = B.B/D = 2*Pb/D
            lambda = 1/(2*P/D + Vac2)
            Vn0 = 0.5*erfc(SgnLR(q)*sqrt(lambda)*VnB(i)) !Expensive

            !Accumulate Va for diffusive Alfven speed for hogs
            VaD(i) = VaD(i) + 0.5*sqrt(Vac2)

            !Now calculate fluxes
            Nx = TfB(i,NORMX)
            Ny = TfB(i,NORMY)
            Nz = TfB(i,NORMZ)

            !Accumulate for L/R-moving
            mFlxB(i,XDIR) = mFlxB(i,XDIR) + Vn0*( dPb*Nx - Bx*Bn(i) )
            mFlxB(i,YDIR) = mFlxB(i,YDIR) + Vn0*( dPb*Ny - By*Bn(i) )
            mFlxB(i,ZDIR) = mFlxB(i,ZDIR) + Vn0*( dPb*Nz - Bz*Bn(i) )

            if (Model%doBackground) then
                !Start by adding cross-terms
                BdB0 = B0x*Bx + B0y*By + B0z*Bz
                mFlxB(i,XDIR) = mFlxB(i,XDIR) + Vn0*( -B0x*Bn(i) - Bx*B0n(i) + BdB0*Nx) !+ 0.5*(- B0x*B0n(i) + Pb0(i)*Nx)
                mFlxB(i,YDIR) = mFlxB(i,YDIR) + Vn0*( -B0y*Bn(i) - By*B0n(i) + BdB0*Ny) !+ 0.5*(- B0y*B0n(i) + Pb0(i)*Ny)
                mFlxB(i,ZDIR) = mFlxB(i,ZDIR) + Vn0*( -B0z*Bn(i) - Bz*B0n(i) + BdB0*Nz) !+ 0.5*(- B0z*B0n(i) + Pb0(i)*Nz)
            endif

        enddo
    end subroutine MagKinFlux

end module stress
