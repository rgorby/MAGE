module mhdgroup
    use gamtypes
    use clocks
    use gamutils
    use gridutils
    use stress
    use fields
    use ringav
    use multifluid

    implicit none

    contains

    subroutine AdvanceMHD(Model,Grid,State,oState,Solver,dt)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State,oState
        type(Solver_T), intent(inout) :: Solver
        real(rp), intent(in) :: dt

        integer :: n

        !Use predictor to create half state
        call Tic("Predictor")
        call Predictor(Model,Grid,oState,State,Solver%StateHf,0.5*dt)
        call Toc("Predictor")

        !Get electric field for MHD case
        if (Model%doMHD) then
            !NOTE, CalcElecField and CalcStress are independent of each other
            call Tic("E-Field")
            call CalcElecField(Model,Grid,Solver%StateHf,Solver%Vf,State%Efld)
            !Call user hack function if defined
            if (associated(Model%HackE)) then
                call Tic("HackE")
                call Model%HackE(Model,Grid,State)
                call Toc("HackE")
            endif
            call Toc("E-Field")
            
        endif
        
        !Get plasma stresses
        !Send dGasM in hydro/MHD, but unallocated w/ hydro
        !Calculate stresses for all species
        call Tic("Stress")
        call CalcStress(Model,Grid,Solver%StateHF,Solver%gFlx,Solver%mFlx,Solver%dGasH,Solver%dGasM)
        call Toc("Stress")

        !Finalize, apply stresses and save State->oState for next predictor step
        call Tic("Update")
        call applyUpdate(Model,Grid,State,oState,dt,Solver%dGasH,Solver%dGasM,State%Efld)
        call Toc("Update")

        !Apply gravity if necessary (before ring avg)
        if (Model%doGrav) then
            call Tic("Gravity")
            call applyGrav(Model,Grid,State,oState,dt)
            call Toc("Gravity")
        endif

        !Ring average if necessary
        !Ring average configuration and this grid contains part of singularity
        call Tic("Ring-Average")
        if ( Model%doRing .and. (Model%Ring%doS .or. Model%Ring%doE) ) then
            call RingAverage(Model,Grid,State)
        endif
        call Toc("Ring-Average")

    end subroutine AdvanceMHD
    
    !Updates plasma state using plasma deltas and electric field
    subroutine applyUpdate(Model,Grid,State,oState,dt,dGasH,dGasM,E)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State
        type(State_T), intent(inout) :: oState
        real(rp), intent(in) :: dt
        real(rp), intent(in) :: dGasH(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NVAR,BLK:Model%nSpc)
        real(rp), optional, intent(in) :: dGasM(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), optional, intent(inout) :: E(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM)

        integer :: i,j,k,nv,s
        real(rp), dimension(NVAR,0:Model%nSpc) :: U,oU,dPg
        real(rp), dimension(NDIM) :: B,oB,dPm
        logical :: BorisUpdate
        real(rp) :: Ca
        !DIR$ ASSUME_ALIGNED dGasH: ALIGN
        !DIR$ ASSUME_ALIGNED dGasM: ALIGN
        !DIR$ ASSUME_ALIGNED E: ALIGN
    !--------------------
    !Copy current->old
        call Tic("Copy2Old")
        !Start by saving State->oState for predictor on next step
        !Only need magFlux/State

        !NOTE, doing this with threads to avoid poor scaling
        !$OMP PARALLEL default(shared)

        !$OMP DO collapse(2)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg
                    oState%Gas(i,j,k,:,:) = State%Gas(i,j,k,:,:)
                    oState%Bxyz(i,j,k,:)  = State%Bxyz(i,j,k,:)
                enddo
            enddo
        enddo
        !$OMP END DO NOWAIT !Rely on barrier at end of parallel region

        
        !$OMP DO collapse(2)
        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                do i=Grid%isg,Grid%ieg+1
                    oState%magFlux(i,j,k,:) = State%magFlux(i,j,k,:)
                enddo
            enddo
        enddo
        !$OMP END PARALLEL

        oState%time = State%time
        State%time = State%time + dt

        call Toc("Copy2Old")

    !--------------------
    !Apply hydro stresses
        !State%Gas = State%Gas + dt*dGasH

        call Tic("Reynolds")
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(U,dPg)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=Grid%is,Grid%ie
                    U   = State%Gas(i,j,k,:,:)
                    dPg = dGasH    (i,j,k,:,:)
                    call CellReynolds(Model,U,dPg,Model%dt)
                    State%Gas(i,j,k,:,:) = U
                enddo
            enddo
        enddo
        call Toc("Reynolds")

        !If we're only doing hydro then we're done here
        if (.not. Model%doMHD) then
            return
        endif

    !--------------------
    !Update magnetic fields
        !Only here if we're doing MHD

        call Tic("CT-Update")
        !Now apply E-field update to calculate updated MAG FLUXES
        !Average pole field if doing ring average
        if ( Model%doRing .and. (Model%Ring%doS .or. Model%Ring%doE) ) then
            call PoleE(Model,Grid,E)
        endif

        call E2Flux(Model,Grid,State%magFlux,E,dt)
        !Calculate Bxyz's
        call bFlux2Fld(Model,Grid,State%magFlux,State%Bxyz)
        call Toc("CT-Update")
        
    !--------------------
    !Apply Maxwell/Boris stresses
        !If multifluid then always do Boris (use HUGE if no value specified)
        if (Model%doMultiF) then
            BorisUpdate = .true.
            if (Model%doBoris) then
                Ca = Model%Ca
            else
                Ca = HUGE
            endif
        else
            !Not multifluid
            if (Model%doBoris) then
                BorisUpdate = .true.
                Ca = Model%Ca
            else
                BorisUpdate = .false.
                Ca = 0.0
            endif
        endif

        call Tic("Maxwell")
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(U,oU,B,oB,dPg,dPm) 
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=Grid%is,Grid%ie
                    U  =  State%Gas(i,j,k,:,:)
                    oU = oState%Gas(i,j,k,:,:)
                    B  =  State%Bxyz(i,j,k,:)
                    oB = oState%Bxyz(i,j,k,:)
                    if (Model%doBackground) then
                        B  =  B + Grid%B0(i,j,k,:)
                        oB = oB + Grid%B0(i,j,k,:) 
                    endif
                    dPg = dGasH(i,j,k,:,:)
                    dPm = dGasM(i,j,k,:)
                    if (BorisUpdate) then
                        call CellBoris(Model,U,oU,B,oB,dPg,dPm,Ca,dt)
                    else
                        call CellMaxwell(Model,U,dPm,dt)
                    endif
                    State%Gas(i,j,k,:,:) = U

                enddo
            enddo
        enddo

        call Toc("Maxwell")

    end subroutine applyUpdate

    !Apply reynolds update to cell over all species
    !In multifluid case do accumulation of species->bulk
    subroutine CellReynolds(Model,U,dGasH,dt)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U(NVAR,BLK:Model%nSpc)
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(inout) :: dGasH
        real(rp), intent(in) :: dt

        integer :: s,n
        logical, dimension(Model%nSpc) :: isGood

        if (Model%doMultiF) then
            isGood = ( U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
            !Update individual species then accumulate
            do s=1,Model%nSpc
                if (isGood(s)) then
                    U(:,s) = U(:,s) + dt*dGasH(:,s)
                else
                    !Always apply update (even to cells that aren't above floor)
                    !Otherwise they'll never get above floor
                    !TODO: Test only updating "bad" fluids if delta-rho > 0, i.e. don't take from vacuum
                    U(DEN,s) = U(DEN,s) + dt*dGasH(DEN,s)
                endif
            enddo
            call MultiF2Bulk(Model,U)
            !Reset bulk delta to sum of component species
            do n=1,NVAR
                dGasH(n,BLK) = sum(dGasH(n,1:Model%nSpc),mask=isGood)
            enddo
        else
            !Just do bulk flow and get outta here
            U(:,BLK) = U(:,BLK) + dt*dGasH(:,BLK)
        endif

    end subroutine CellReynolds

    !Apply Maxwell update to cell (only used in single species case)
    subroutine CellMaxwell(Model,U,dGasM,dt)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U(NVAR,BLK:Model%nSpc)
        real(rp), dimension(NDIM), intent(in) :: dGasM
        real(rp), intent(in) :: dt

        real(rp) :: D,Mx,My,Mz,TotE,KinE,P

        D = max( U(DEN,BLK),dFloor )
        Mx = U(MOMX,BLK)
        My = U(MOMY,BLK)
        Mz = U(MOMZ,BLK)
        TotE = U(ENERGY,BLK)

        KinE = 0.5*(Mx**2.0 + My**2.0 + Mz**2.0)/D
        P = max((Model%gamma-1)*( TotE - KinE ),pFloor)

        !Perform momentum update
        Mx = Mx + dt*dGasM(XDIR)
        My = My + dt*dGasM(YDIR)
        Mz = Mz + dt*dGasM(ZDIR)
        KinE = 0.5*(Mx**2.0 + My**2.0 + Mz**2.0)/D
    
        !Update variables
        U(DEN   ,BLK) = D !Incorporate floor
        U(MOMX  ,BLK) = Mx
        U(MOMY  ,BLK) = My
        U(MOMZ  ,BLK) = Mz
        U(ENERGY,BLK) = KinE + P/(Model%gamma-1)

    end subroutine CellMaxwell


    !Apply Boris update to cell over all species
    !Assuming already have updated hydro w/ dGasH
    !U = hydro-updated state, do update in place
    !Assuming B2/B0 are *TOTAL* fields (ie includes background)
    subroutine CellBoris(Model,U2,U0,B2,B0,dGasH,dGasM,Ca,dt)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U2(NVAR,BLK:Model%nSpc)
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(in) :: U0,dGasH
        real(rp), dimension(NDIM), intent(in) :: B2, B0, dGasM
        real(rp), intent(in) :: Ca,dt

        integer :: s
        real(rp) :: D0 ,D2 ,P2 ,alpha ,pCon ,dRho
        real(rp) :: D0s,D2s,P2s,alphaS,pConS,dRhoS,Rs,Rsbar
        real(rp), dimension(NDIM) :: B1,bhat1,bhat2,dPMag
        real(rp), dimension(NDIM) :: V0,V2,V0s,Vtmp,Vperp
        real(rp), dimension(NDIM) :: dPGas,dPGasS,dMom
        real(rp), dimension(NVAR) :: W2

    !Start w/ bulk update even in multifluid case
        !Get D, pressure from hydro update
        call CellC2P(Model,U2(:,BLK),W2)
        D2 = max(W2(DEN),dFloor)
        P2 = max(W2(PRESSURE),pFloor)

        !Get old D/V
        D0 = U0(DEN      ,BLK)
        V0 = U0(MOMX:MOMZ,BLK)/max(D0,dFloor)

        !B0 = n, B1 = n+1/2, B2 = B
        B1 = 0.5*(B2+B0) !Half-step field
        !Get terms for updating
        alpha = dot_product(B1,B1)/(D2*Ca*Ca) !Ratio of magnetic/plasma mass
        dRho  = dt*dGasH(DEN      ,BLK)
        dPGas = dt*dGasH(MOMX:MOMZ,BLK)
        dPMag = dt*dGasM(XDIR:ZDIR)

        pCon = 1.0/(1.0+alpha) !Dg/(Db+Dg), ratio of plasma mass and total (w/ magnetic) mass
        bhat1 = normVec(B1)
        bhat2 = normVec(B2) !Final field direction

        if (alpha > TINY) then
            !Note, only using dPMag in perp 
            dMom = Vec2Para(dPGas,bhat1) + pCon*Vec2Perp(dPGas+dPMag,bhat1) &
                 + alpha*pCon*dRho*Vec2Perp(V0,bhat1) 
        else 
            !Null field region, just do Maxwell update
            dMom = dPGas + dPMag
        endif !alpha

        !Finish bulk fluid update using dMom, reset primitives then convert to conserved
        V2 = (D0*V0 + dMom)/D2
        W2(DEN      ) = D2
        W2(VELX:VELZ) = V2
        W2(PRESSURE ) = P2
        call CellP2C(Model,W2,U2(:,BLK))

        !Now handle multifluid update
        if (Model%doMultiF) then
            !Do species independent things
            Vperp = Vec2Perp(V2,bhat2)

            !Loop over species
            do s=1,Model%nSpc
                !Get D, pressure from hydro update
                call CellC2P(Model,U2(:,s),W2)
                D2s = max(W2(DEN),dFloor)
                P2s = max(W2(PRESSURE),pFloor)
                
                !Get old D/V
                D0s = U0(DEN      ,s)
                V0s = U0(MOMX:MOMZ,s)/max(D0s,dFloor)

                if (D2s >= Spcs(s)%dVac) then
                    alphaS = alpha*0.5*(D2s+D0s)/D2s
                    pConS = 1.0/(1.0+alphaS)
                    Rs = D2s/D2
                    Rsbar = 0.5*( D2s/D2 + D0s/D0 )
                    dRhoS = D2s-D0s
                    dPGasS = dt*dGasH(MOMX:MOMZ,s)

                    !Note, dPMag depends on pCon (not pConS)
                    dMom = pCon*Rs*dPMag                          &
                         + alphaS*pConS*dRhoS*Vec2Perp(V0s,bhat1) &
                         + pConS*Rsbar*Vec2Perp(dPGas,bhat1)      &
                         + Vec2Para(dPGasS,bhat1)
                         
                    Vtmp = (D0s*V0s + dMom)/D2s
                else 
                    !Not a good species, just use bulk flow
                    Vtmp = V2
                endif
                W2(DEN      ) = D2s
                W2(VELX:VELZ) = Vperp + Vec2Para(Vtmp,bhat2)
                W2(PRESSURE ) = P2s
                call CellP2C(Model,W2,U2(:,s))
            enddo

            !Done calculating species conserved quantities, recalculate bulk
            call MultiF2Bulk(Model,U2)

        endif !Multifluid

    end subroutine CellBoris

    !Predictor on single cell, [oU,U] -> pU
    subroutine CellPredictor(Model,ht,oU,U,pU)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: ht
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(in)  :: oU,U
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(out) :: pU

        real(rp), dimension(NVAR) :: oW,W,pW
        integer :: s,s0,sE
        logical, dimension(0:Model%nSpc) :: isGood,isGood1,isGood2

        isGood  = .true.
        isGood1 = .false.
        isGood2 = .false.

        !Figure out which fluids/times are good
        if (Model%doMultiF) then
            isGood1(1:Model%nSpc) = ( oU(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
            isGood2(1:Model%nSpc) = (  U(DEN,1:Model%nSpc) >= Spcs(:)%dVac )
            isGood (1:Model%nSpc) = isGood1(1:Model%nSpc) .and. isGood2(1:Model%nSpc)
            s0 = 1
            sE = Model%nSpc            
        else
            s0 = BLK
            sE = BLK
        endif

        do s=s0,sE
            if ( isGood(s) ) then
                !Both states are good, do standard predictor
                !Convert both states to primitive
                call CellC2P(Model,oU(:,s),oW)
                call CellC2P(Model, U(:,s), W)
                !Do predictor on primitives
                pW = W + ht*(W - oW)

                !Apply global (not species) clamps
                pW(DEN) = max(pW(DEN),dFloor)
                pW(PRESSURE) = max(pW(PRESSURE),pFloor)

                !Return to conserved for final predicted state
                call CellP2C(Model,pW,pU(:,s))
            !Otherwise at least 1 state is bad, just use most recent
            else
                pU(:,s) = U(:,s)
            endif
            
        enddo

        if (Model%doMultiF) then
            !Accumulate to bulk
            call MultiF2Bulk(Model,pU)
        endif

    end subroutine CellPredictor

    subroutine Predictor(Model,Grid,oState,State,pState,pdt)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(in) :: State, oState
        type(State_T), intent(inout) :: pState
        real(rp), intent(in) :: pdt

        integer :: i,j,k
        real(rp) :: odt,ht
        logical :: isCC

        !Do timing stuff
        odt = State%time-oState%time !Sep. between states
        pState%time = State%time + pdt
        ht = pdt/odt

        !One big // region
        !$OMP PARALLEL default(shared) private(i,j,k,isCC)

        !Loop over grid and perform predictor on each cell of fluid/s
        !XYZ fields and interface fluxes
        !$OMP DO collapse(2)
        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                do i=Grid%isg,Grid%ieg+1
                    !Check if loop iteration is in interior
                    isCC = (k<=Grid%keg) .and. (j<=Grid%jeg) .and. (i<=Grid%ieg)
                    if (isCC) then
                        call CellPredictor(Model,ht,oState%Gas(i,j,k,:,:),State%Gas(i,j,k,:,:),pState%Gas(i,j,k,:,:))
                    endif
                    if (Model%doMHD) then
                        !Do interface fluxes
                        pState%magFlux(i,j,k,:) = State%magFlux(i,j,k,:) + (pdt/odt)*(State%magFlux(i,j,k,:) - oState%magFlux(i,j,k,:))
                    endif !MHD
                enddo !I loop
            enddo
        enddo !K loop (implicit barrier)

        !Now do flux->field where necessary
        if (Model%doMHD) then
            !Loop through grid and replace Bxyz w/ flux-> field 
            !$OMP DO collapse(2)
            do k=Grid%ksg,Grid%keg
                do j=Grid%jsg,Grid%jeg
                    do i=Grid%isg,Grid%ieg                        
                        pState%Bxyz(i,j,k,:) = CellBxyz(Model,Grid,pState%magFlux,i,j,k)
                    enddo
                enddo
            enddo
        endif

        !Close big parallel region
        !$OMP END PARALLEL

        if (Model%doRing) call RingPredictorFix(Model,Grid,pState)
        
        if (associated(Model%HackPredictor)) then
            call Model%HackPredictor(Model,Grid,pState)
        endif

    end subroutine Predictor

    !Updates plasma state using static gravitational potential
    !Use time-averaged density from n,n+1 and update momentum
    subroutine applyGrav(Model,Gr,State,oState,dt)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State
        type(State_T), intent(inout) :: oState
        real(rp), intent(in) :: dt

        integer :: i,j,k,s,s0,sE
        
        real(rp) :: Dhf,Dblk
        real(rp), dimension(NDIM) :: g,B,bMxyz,cMxyz
        real(rp), dimension(NDIM,NDIM) :: Lam,Laminv
        real(rp), dimension(NVAR) :: pW,pCon
        real(rp), dimension(0:Model%nSpc) :: RhoMin


        RhoMin(BLK) = 0.0
        if (Model%doMultiF) then
            !Don't do bulk
            s0 = 1
            sE = Model%nSpc
            RhoMin(1:Model%nSpc) = Spcs(:)%dVac
        else
            !Only do bulk
            s0 = BLK
            sE = BLK
        endif

        !Add grav forces
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(s,i,j,k,Dhf,Dblk,bMxyz,cMxyz) &
        !$OMP private(g,B,Lam,Laminv,pW,pCon)
        do s=s0,sE
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%is,Gr%ie
                        if ( State%Gas(i,j,k,DEN,s) >= RhoMin(s) ) then
                        !Get necessary state info
                            pCon = State%Gas(i,j,k,1:NVAR,s) !Conserved
                            call CellC2P(Model,pCon,pW) !Primitive
                            Dblk = State%Gas(i,j,k,DEN,BLK) !Bulk density

                            !Get average density
                            Dhf = 0.5*( State%Gas(i,j,k,DEN,s) + oState%Gas(i,j,k,DEN,s) )
                        !Map classical momentum to boris momentum
                            !NOTE: transform is I if Boris isn't on
                            cMxyz = pCon(MOMX:MOMZ) !Classical momentum
                            B = State%Bxyz(i,j,k,:)
                            call Mom2Rel(Model,Dblk,B,Lam)
                            bMxyz = matmul(Lam,cMxyz)
                        !Apply gravitational force to boris momentum
                            g = Gr%gxyz(i,j,k,XDIR:ZDIR)
                            bMxyz = bMxyz + Model%dt*Dhf*g
                        !Map updated Boris momentum back to classical
                            call Rel2Mom(Model,Dblk,B,Laminv)
                            cMxyz = matmul(Laminv,bMxyz) !Classical momentum
                        !Set new primitive variables, convert and store
                            !Only V has changed
                            pW(VELX:VELZ) = cMxyz/max(pW(DEN),dFloor)
                            call CellP2C(Model,pW,pCon)
                            State%Gas(i,j,k,:,s) = pCon
                        endif
                    enddo !i
                enddo !j
            enddo !k
        enddo !s

        if (Model%doMultiF) then
            call State2Bulk(Model,Gr,State)
        endif

    end subroutine applyGrav

end module mhdgroup
