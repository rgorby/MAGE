module mhdgroup
    use types
    use clocks
    use gamutils
    use gridutils
    use stress
    use fields
    use ringav
    use multifluid

    implicit none

    !TODO: Look into removing half-state holder
    type(State_T) :: StateHf
    logical :: doInit = .true.
    real(rp), dimension(:,:,:,:,:), allocatable :: dGasH
    real(rp), dimension(:,:,:,:), allocatable :: dGasM

    contains

    subroutine AdvanceMHD(Model,Grid,State,oState,dt)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State,oState
        real(rp), intent(in) :: dt

        integer :: n

        !Initialize integrator if necessary
        if (doInit) then
            call allocState(Model,Grid,StateHf)
            call allocGridHydroSpc(Model,Grid,dGasH)  
            if (Model%doMHD) call allocGridVec(Model,Grid,dGasM)       
            doInit = .false.
        endif

        !Use predictor to create half state
        call Tic("Predictor")
        call Predictor(Model,Grid,oState,State,StateHf,0.5*dt)
        call Toc("Predictor")

        !Get electric field for MHD case
        if (Model%doMHD) then
            !NOTE, CalcElecField and CalcStress are independent of each other
            call Tic("E-Field")
            call CalcElecField(Model,Grid,StateHf,State%Efld)
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
        call CalcStress(Model,Grid,StateHF,dGasH,dGasM)
        call Toc("Stress")

        !Finalize, apply stresses and save State->oState for next predictor step
        call Tic("Update")
        call applyUpdate(Model,Grid,State,oState,dt,dGasH,dGasM,State%Efld)
        call Toc("Update")

        !Apply gravity if necessary (before ring avg)
        if (Model%doGrav) then
            call Tic("Gravity")
            call applyGrav(Model,Grid,State,oState,dt)
            call Toc("Gravity")
        endif

        !Ring average if necessary
        !Ring average configuration and this grid contains part of singularity
        if ( Model%doRing .and. (Model%Ring%doS .or. Model%Ring%doE) ) then
            call Tic("Ring-Average")
            call Tic("RA-Hydro")
            call RingAvgHydro(Model,Grid,State)
            call Toc("RA-Hydro")
            if (Model%doMHD) then
                call Tic("RA-Mag")
                call RingAvgMag(Model,Grid,State)
                call Toc("RA-Mag")
            endif
            call Toc("Ring-Average")
        endif

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

        integer :: s
        if (Model%doMultiF) then

            !Update individual species then accumulate
            do s=1,Model%nSpc
                !Always apply update (even to cells that aren't above floor)
                !Otherwise they'll never get above floor
                U(:,s) = U(:,s) + dt*dGasH(:,s)
            enddo
            call MultiF2Bulk(Model,U)
            !Reset bulk delta to sum of component species
            dGasH(:,BLK) = sum(dGasH(:,1:Model%nSpc),dim=2)
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
    !Assuming B/oB are *TOTAL* fields (ie includes background)
    subroutine CellBoris(Model,U,oU,B,oB,dGasH,dGasM,Ca,dt)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: U(NVAR,BLK:Model%nSpc)
        real(rp), dimension(NVAR,BLK:Model%nSpc), intent(in) :: oU,dGasH
        real(rp), dimension(NDIM), intent(in) :: B, oB, dGasM
        real(rp), intent(in) :: Ca,dt

        real(rp) :: D,Mx,My,Mz,TotE,KinE,P,dRho
        real(rp) :: alfnRatio,perpRatio,dvAlf,BdotV
        real(rp) :: oD,Dblk,oDblk,spcR,dRhoS,alfStar,perpCons
        real(rp), dimension(NDIM) :: B1,dPg,dPm,V0,Vtmp,Pnew,dPs
        real(rp), dimension(NDIM) :: Vblk,Vpara,Vperp
        real(rp), dimension(NDIM) :: bhat1,bhat2,Fperp

        integer :: s

        !Start with bulk update even in multifluid case
        Mx = U(MOMX,BLK)
        My = U(MOMY,BLK)
        Mz = U(MOMZ,BLK)
        D  = max(U(DEN,BLK),dFloor)
        TotE  = U(ENERGY,BLK)
        KinE = 0.5*(Mx**2.0 + My**2.0 + Mz**2.0)/D
        P = max((Model%gamma-1)*( TotE - KinE ),pFloor)

        !B0 = n, B1 = n+1/2, B2 = B
        B1 = 0.5*(B+oB) !Half-step field
        !Get ratios for Boris update
        alfnRatio = dot_product(B1,B1)/(D*Ca*Ca)
        perpRatio = 1.0/(1+alfnRatio)

        !Get deltas
        dRho = dt*dGasH(DEN      ,BLK)
        dPg  = dt*dGasH(MOMX:MOMZ,BLK)
        dPm  = dt*dGasM(XDIR:ZDIR)

        if (alfnRatio>TINY) then
            !Get old velocity
            V0 = oU(MOMX:MOMZ,BLK)/max(oU(DEN,BLK),dFloor)
            dvAlf = alfnRatio*dRho
            Vtmp  = alfnRatio*(dPg-dRho*V0)

            !Calculate B dot V
            BdotV = dot_product(B1,Vtmp)/dot_product(B1,B1)

            !Perform FINAL momentum update, updating from state n->n+1
            Pnew = oU(MOMX:MOMZ,BLK) + perpRatio*( dPg + dPm + dvAlf*V0 + BdotV*B1 )
        else
            !Null field region, just do Maxwell update
            Pnew = U(MOMX:MOMZ,BLK) + dPm
        endif

        !Apply bulk update
        KinE = 0.5*dot_product(Pnew,Pnew)/D
        U(DEN,      BLK) = D
        U(MOMX:MOMZ,BLK) = Pnew
        U(ENERGY   ,BLK) = KinE + P/(Model%gamma-1)

        !Now handle multifluid case
        if (Model%doMultiF) then
            !Do species independent things
            !Bulk flow
            Dblk  = max( U(DEN,BLK),dFloor)
            oDblk = max(oU(DEN,BLK),dFloor)
            Vblk = U(MOMX:MOMZ,BLK)/Dblk
            !Magnetic field directions
            bhat2 = normVec(B)
            bhat1 = normVec(B1)

            Vperp = Vblk - dot_product(Vblk,bhat2)*bhat2
            Fperp = dPm  - dot_product(dPm ,bhat1)*bhat1

            !Loop over species
            !Need species dependent V-parallel
            do s=1,Model%nSpc
                
                !Start by getting species state
                D = max(U(DEN,s),dFloor)
                Pnew = U(MOMX:MOMZ,s)
                TotE = U(ENERGY,s)
                KinE = 0.5*dot_product(Pnew,Pnew)/D
                P = max((Model%gamma-1)*( TotE - KinE ),pFloor)

                if (D >= Spcs(s)%dFloor) then
                    !This is a good species, calculate Vpara
                    oD = max(oU(DEN,s),dFloor)
                    spcR = 0.5*( oD/oDblk + D/Dblk )
                    dRhoS = D - oD
                    alfStar = alfnRatio*( 1 - 0.5*dRhoS/D )
                    perpCons = 1.0/(1+alfStar)

                    !Rescale dRhoS to species
                    dRhoS = dRhoS*alfStar
                    dPs = (1+alfStar)*dt*dGasH(MOMX:MOMZ,s) - dRhoS*oU(MOMX:MOMZ,s)/oD
                    Pnew = oU(MOMX:MOMZ,s) + perpCons*dRhoS*oU(MOMX:MOMZ,s)/oD &
                       & + perpCons*spcR*Fperp                   &
                       & + perpCons*dot_product(dPs,bhat1)*bhat1
                    Vtmp = (Pnew/D) + perpRatio*dPm/Dblk
                else
                    !This isn't a good species, use bulk
                    Vtmp = Vblk
                endif
                Vpara = dot_product(Vtmp,bhat2)*bhat2
                !Finish species update
                U(DEN      ,s) = D
                U(MOMX:MOMZ,s) = D*(Vpara + Vperp)
                KinE = 0.5*D*dot_product(Vpara + Vperp,Vpara + Vperp)
                U(ENERGY   ,s) = KinE + P/(Model%gamma-1)
            enddo

            !Done calculating species conserved quantities, recalculate bulk
            call MultiF2Bulk(Model,U)

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
            isGood1(1:Model%nSpc) = ( oU(DEN,1:Model%nSpc) >= Spcs(:)%dFloor )
            isGood2(1:Model%nSpc) = (  U(DEN,1:Model%nSpc) >= Spcs(:)%dFloor )
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
        type(Grid_T), intent(in) :: Grid
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

        !Loop over grid and perform predictor on each cell of fluid/s
        !XYZ fields and interface fluxes
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,isCC)
        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                do i=Grid%isg,Grid%ieg+1
                    !Check if loop iteration is in interior
                    isCC = (k<=Grid%keg) .and. (j<=Grid%jeg) .and. (i<=Grid%ieg)
                    if (isCC) then
                        call CellPredictor(Model,ht,oState%Gas(i,j,k,:,:),State%Gas(i,j,k,:,:),pState%Gas(i,j,k,:,:))
                        !Do XYZ fields
                        if (Model%doMHD) then
                            pState%Bxyz(i,j,k,:) = State%Bxyz(i,j,k,:) + (pdt/odt)*(State%Bxyz(i,j,k,:) - oState%Bxyz(i,j,k,:))
                        endif !MHD
                    endif
                    if (Model%doMHD) then
                        !Do interface fluxes
                        pState%magFlux(i,j,k,:) = State%magFlux(i,j,k,:) + (pdt/odt)*(State%magFlux(i,j,k,:) - oState%magFlux(i,j,k,:))
                    endif !MHD
                enddo !I loop
            enddo
        enddo !K loop

        !Now do flux->field where necessary
        if (Model%doMHD) then
            !Replace Bxyz w/ flux->field calculations in xxMG region
            !$OMP PARALLEL DO default(shared) collapse(2)
            do k=Grid%ksMG,Grid%keMG
                do j=Grid%jsMG,Grid%jeMG
                    do i=Grid%isMG,Grid%ieMG
                        pState%Bxyz(i,j,k,:) = CellBxyz(Model,Grid,pState%magFlux,i,j,k)
                    enddo
                enddo
            enddo
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
        real(rp) :: Dhf,IntE
        real(rp), dimension(NDIM) :: p,g
        real(rp), dimension(NVAR) :: pW,pCon
        real(rp), dimension(0:Model%nSpc) :: RhoMin

        RhoMin(BLK) = 0.0
        if (Model%doMultiF) then
            !Don't do bulk
            s0 = 1
            sE = Model%nSpc
            RhoMin(1:Model%nSpc) = Spcs(:)%dFloor
        else
            !Only do bulk
            s0 = BLK
            sE = BLK
        endif

        !Add grav forces
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(s,i,j,k,Dhf,IntE,p,g,pW,pCon)
        do s=s0,sE
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%is,Gr%ie
                        if ( State%Gas(i,j,k,DEN,s) >= RhoMin(s) ) then
                            !Save thermal energy of state
                            pCon = State%Gas(i,j,k,1:NVAR,s)
                            call CellC2P(Model,pCon,pW)
                            IntE = pW(PRESSURE)/(Model%gamma-1)

                            !Get average density
                            Dhf = 0.5*( State%Gas(i,j,k,DEN,s) + oState%Gas(i,j,k,DEN,s) )

                            !Update momentum
                            g = Gr%gxyz(i,j,k,XDIR:ZDIR)
                            p = pCon(MOMX:MOMZ) + Model%dt*Dhf*g

                            !Reset conserved state
                            pCon(MOMX:MOMZ) = p
                            pCon(ENERGY) = IntE + 0.5*dot_product(p,p)/pCon(DEN)
                            State%Gas(i,j,k,:,s) = pCon
                        endif
                    enddo
                enddo
            enddo
        enddo
        
        if (Model%doMultiF) then
            call State2Bulk(Model,Gr,State)
        endif

    end subroutine applyGrav

end module mhdgroup
