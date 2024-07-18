module raijuPreAdvancer
    use clocks
    use planethelper
    use kai2geo

    use raijudefs
    use raijutypes
    use raijuetautils
    use raijuBCs
    use raijulosses
    use raijuRecon

    use shellGrid
    use shellInterp

    implicit none

    contains


!------
! Main high-level functions
!------

    subroutine raijuPreAdvance(Model, Grid, State, isFirstCplO)
        !! Takes a state and calculates what is needed in order to advance
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        logical, optional, intent(in) :: isFirstCplO

        logical :: isFirstCpl
        integer :: k

        if (present(isFirstCplO)) then
            isFirstCpl = isFirstCplO
        else
            isFirstCpl = .false.
        endif

        ! Clear things that will be accumulated over the advance
        State%precipType_ele = 0.0
        State%precipNFlux = 0.0
        State%precipEFlux = 0.0

        ! Moments to etas, initial active shell calculation
        call Tic("BCs")
        call applyRaijuBCs(Model, Grid, State, doWholeDomainO=isFirstCpl) ! If fullEtaMap=True, mom2eta map is applied to the whole domain
        call Toc("BCs")

        ! Handle edge cases that may effect the validity of information carried over from last coupling period
        ! TODO: do this in predictor function
        call prepEtaLast(Grid%shGrid, State, isFirstCpl)

        ! Calc cell velocities
        call Tic("Calc face velocities")
        call calcPotGrads(Model, Grid, State)
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call calcVelocity(Model, Grid, State, k, State%iVel(:,:,k,:))  ! Get velocity at cell interfaces
            ! Calc sub-time step. Each channel will do this on its own, but this way we can output the step sizes everyone is using
            State%dtk(k) = activeDt(Model, Grid, State, k)
        enddo
        call Toc("Calc face velocities")

        call Tic("Calc cell-centered velocities")
        call calcPotGrads_cc(Model, Grid, State)
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call calcVelocityCC_gg(Model, Grid, State, k, State%cVel(:,:,k,:))  ! Get velocity at cell interfaces
            ! Calc sub-time step. Each channel will do this on its own, but this way we can output the step sizes everyone is using
            !State%dtk(k) = activeDt(Model, Grid, State, k)
            call reconVelocityLRs(Model, Grid, State, k, State%iVelL(:,:,k,:), State%iVelR(:,:,k,:))
        enddo
        call Toc("Calc cell-centered velocities")

        !if (Model%doDebugOutput) then
        !    call Tic("Calc CC velocities")
        !    call velFace2CC(Model, Grid, State)
        !    call Toc("Calc CC velocities")
        !endif

        ! Loss rate calc depends on up-to-date densities, so we should run EvalMoments first
        call Tic("Moments Eval PreAdvance")
        call EvalMoments(Grid, State)
        call Toc("Moments Eval PreAdvance")

        call Tic("Calc loss rates")
        if (Model%doLosses) then
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(k)
            do k=1,Grid%Nk
                call calcChannelLossRates(Model, Grid, State, k)
            enddo
        endif
        call Toc("Calc loss rates")

    end subroutine raijuPreAdvance


!------
! Transitions between coupling chunks
!------

    subroutine prepEtaLast(sh, State, isFirstCpl)
        ! To be called before main advance loop
        ! Edits the initial eta_last state to account for edge cases
        type(ShellGrid_T), intent(in) :: sh
        type(raijuState_T), intent(inout) :: State
        logical, intent(in) :: isFirstCpl

        integer :: i,j

        ! If its the first coupling instance, we have no good previous information
        ! Set previous cpl info to current state
        if (isFirstCpl) then
            State%active_last = State%active
            State%eta_last    = State%eta
            return
        endif

        ! For a given cell, if it is currently active, but wasn't active last step, 
        !   we need to set its eta_last to something usable in the first eta_half calculation
        ! We do this in the buffer region as well, because this was freshly populated with MHD moments and any existing eta_last is stale
        ! After the first dt, eta_last will be set with valid information and no further adjustments are necessay
        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg
                if ( (State%active(i,j) .eq. RAIJUACTIVE .and. State%active_last(i,j) .eq. RAIJUINACTIVE) &
                .or. (State%active(i,j) .eq. RAIJUBUFFER) ) then
                    State%eta_last(i,j,:) = State%eta(i,j,:)
                endif
            enddo
        enddo

        ! NOTE: Any case not addressed here should mean that eta_last and active_last are valid in those cases

    end subroutine prepEtaLast


!------
! Cell Potentials and their gradients
!------
    !! Potential variables come in pre-allocated, since they may be subsets of larger arrays
    subroutine potExB(sh, State, pExB)
        ! Trivial, but putting it here in case we need extra options for it later
        type(ShellGrid_T), intent(in) :: sh
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(inout) :: pExB
        
        pExB = State%espot * 1.e3  ! [kV -> V]

    end subroutine potExB

    subroutine potCorot(planet, sh, pCorot, doGeoCorotO)
        ! Simple corotation potential [V]
        type(planet_T), intent(in) :: planet
        type(ShellGrid_T), intent(in) :: sh
        real(rp), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(inout) :: pCorot
        logical, intent(in), optional :: doGeoCorotO
        
        integer :: i,j

        if (present(doGeoCorotO)) then
            if (doGeoCorotO) then
                pCorot = 0.0
                do j=sh%jsg,sh%jeg+1
                    do i=sh%isg,sh%ieg+1
                        call geocorotation_from_SMTP(sh%th(i), sh%ph(j), pCorot(i,j))
                    enddo
                enddo
                ! geopack corotation returns in kV, we need V
                pCorot = pCorot * 1.e3  ! [kV -> V]
                return
            endif
        endif

        ! IfdoGeoCorotO was not provided, or it was false, we default to corotation on aligned dipole and rotational axis
        do j=sh%jsg,sh%jeg+1
            pCorot(:,j) = -planet%psiCorot*(planet%rp_m/planet%ri_m)*sin(sh%th)**2 * 1.e3  ! [kV -> V]
        enddo

    end subroutine potCorot


    subroutine potGC(shGrid, alamc, bVol, pGC)
        ! Simple gradient curvature potential [V]
        type(ShellGrid_T), intent(in) :: shGrid
        real(rp), intent(in) :: alamc
        real(rp), dimension(shGrid%isg:shGrid%ieg+1,&
                            shGrid%jsg:shGrid%jeg+1), intent(in) :: bVol

        real(rp), dimension(shGrid%isg:shGrid%ieg+1,&
                            shGrid%jsg:shGrid%jeg+1), intent(inout) :: pGC

        pGC = alamc*bVol**(-2./3.)

    end subroutine potGC


    subroutine calcEffectivePotential(Model, Grid, State, pEff)
        !! Note: This is not used to calculate velocities in the solver itself
        !! Calculates effective potential [V] for all lambda channels as the sum of pExB, pCorot, and pGC
        type(raijuModel_T) , intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1,&
                            Grid%Nk), intent(inout) :: pEff
        
        integer :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: pExB, pCorot, pGC
        
        pEff = 0.0

        ! Grab 2D potentials
        call potExB(Grid%shGrid, State, pExB)
        call potCorot(Model%planet, Grid%shGrid, pCorot, Model%doGeoCorot)
        
        ! Build 3D effective potential
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(k, pGC)
        do k=1,Grid%Nk
            call potGC(Grid%shGrid, Grid%alamc(k), State%bVol, pGC)
            pEff(:,:,k) = pExB + pCorot + pGC
        enddo

    end subroutine calcEffectivePotential


    subroutine calcPotGrads(Model, Grid, State)
        !! Calcs gradient of ionospheric potential, corotation potential, and (flux tube volume ^ -2/3)
        !! Only needs to be called when any of the above quantities are updated (e.g. beginning of every coupling timestep)
        !! Note: FTV is not a potential, alamc*bVol**(-2./3.) is
        !! Units (gradE and gradCorot): V / m
        !! Units (gradVM): Vol^(-2/3) / m
        type(raijuModel_T), intent(in)    :: Model
        type(raijuGrid_T ), intent(in)    :: Grid
        type(raijuState_T), intent(inout) :: State


        integer :: i,j
        ! Cell corners
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: isGCorner
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: pExB, pCorot
        ! Cell faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: gradVM
        
        ! Determine which cell corners we consider good enough to calculate gradients from
        where (State%topo .eq. RAIJUCLOSED)
            isGCorner = .true.
        elsewhere
            isGCorner = .false.
        end where

        State%gradVM = 0.0
        
        ! Ionospheric and corotation potentials are just simple derivatives
        call potExB(Grid%shGrid, State, pExB)  ! [V]
        call potCorot(Model%planet, Grid%shGrid, pCorot, Model%doGeoCorot)  ! [V]
        call calcGradIJ(Model%planet%rp_m, Grid, isGCorner, pExB  , State%gradPotE    )  ! [V/m]
        call calcGradIJ(Model%planet%rp_m, Grid, isGCorner, pCorot, State%gradPotCorot)  ! [V/m]

        ! GC drifts depend on grad(lambda * V^(-2/3))
        ! lambda is constant, so just need grad(V^(-2/3) )
        ! grad(V^(-2/3)) = -2/3*V^(-5/3) * grad(V)
        !call calcGradFTV(Model%planet%rp_m, Model%planet%ri_m, Model%planet%magMoment, Grid, isGCorner, State%bvol, gradVM)
        !State%gradVM(:,:,RAI_TH) = (-2./3.) * State%bvol**(-5./3.) * gradVM(:,:,RAI_TH)  ! [Vol^(-2/3)/m]
        !State%gradVM(:,:,RAI_PH) = (-2./3.) * State%bvol**(-5./3.) * gradVM(:,:,RAI_PH)  ! [Vol^(-2/3)/m]
        call calcGradFTV(Model%planet%rp_m, Model%planet%ri_m, Model%planet%magMoment, Grid, isGCorner, State%bvol, gradVM)

        do i=Grid%shGrid%isg,Grid%shGrid%ieg+1
            do j=Grid%shGrid%jsg,Grid%shGrid%jeg+1
                if (State%bvol(i,j) > 0.0) then
                    State%gradVM(i,j,RAI_TH) = (-2./3.) * State%bvol(i,j)**(-5./3.) * gradVM(i,j,RAI_TH)  ! [Vol^(-2/3)/m]
                    State%gradVM(i,j,RAI_PH) = (-2./3.) * State%bvol(i,j)**(-5./3.) * gradVM(i,j,RAI_PH)  ! [Vol^(-2/3)/m]
                endif
            enddo
        enddo

    end subroutine calcPotGrads


    subroutine calcPotGrads_cc(Model, Grid, State)
        !! Calculate cell-centered potential gradients using Green-Gauss method
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State


        integer :: i,j
        real(rp) :: bvol_cc
        ! Cell corners
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: isGCorner
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: pExB, pCorot
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradVM

        State%cVel = 0.0

        where (State%topo .eq. RAIJUCLOSED)
            isGCorner = .true.
        elsewhere
            isGCorner = .false.
        end where

        associate(sh=>Grid%shGrid)
        ! Gauss-Green calculation of cell-averaged phi gradients
        ! NOTE: These are not across cell faces where gradPot(th) is phi gradient stored on theta edge.
        !       grad(phi)_th is the gradient in the theta direction
        call potExB(Grid%shGrid, State, pExB)  ! [V]
        call potCorot(Model%planet, Grid%shGrid, pCorot, Model%doGeoCorot)  ! [V]
        call calcGradIJ_cc(Model%planet%rp_m, Grid, isGCorner, pExB  , State%gradPotE_cc    )  ! [V/m
        call calcGradIJ_cc(Model%planet%rp_m, Grid, isGCorner, pCorot, State%gradPotCorot_cc)  ! [V/m]

        ! GC drifts depend on grad(lambda * V^(-2/3))
        ! lambda is constant, so just need grad(V^(-2/3) )
        ! grad(V^(-2/3)) = -2/3*V^(-5/3) * grad(V)
        call calcGradFTV_cc(Model%planet%rp_m, Model%planet%ri_m, Model%planet%magMoment, Grid, isGCorner, State%bvol, gradVM)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg    
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (all(isGCorner(i:i+1,j:j+1))) then
                    State%gradVM_cc(i,j,RAI_TH) = (-2./3.) * State%bvol_cc(i,j)**(-5./3.) * gradVM(i,j,RAI_TH)  ! [Vol^(-2/3)/m]
                    State%gradVM_cc(i,j,RAI_PH) = (-2./3.) * State%bvol_cc(i,j)**(-5./3.) * gradVM(i,j,RAI_PH)  ! [Vol^(-2/3)/m]
                endif
            enddo
        enddo
        end associate
        
    end subroutine calcPotGrads_cc


    subroutine calcGradIJ(Rp_m, Grid, isG, Q, gradQ)
        !! Calc gradient in spherical coordinates of corner variable Q across cell edges/faces
        real(rp), intent(in) :: Rp_m
            !! Planet radius in m
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: isG
            !! Mask for corners that are safe to use in calculating the gradient across the attached faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: Q
            !! Variable we are taking the gradient of across faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: gradQ
            !! [units(Q)/m] we return
        
        integer :: i,j

        gradQ = 0.0

        associate(sh=>Grid%shGrid)

        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg+1
                if ( isG(i,j) .and. isG(i,j+1) ) then
                    gradQ(i,j,RAI_TH) = (Q(i,j+1) - Q(i,j)) / ( Grid%lenFace(i,j,RAI_TH) * Rp_m )
                endif
            enddo
        enddo

        do j=sh%jsg,sh%jeg+1
            do i=sh%isg,sh%ieg
                if ( isG(i,j) .and. isG(i+1, j) ) then
                    gradQ(i,j,RAI_PH) = (Q(i+1,j) - Q(i,j)) / ( Grid%lenFace(i,j,RAI_PH) * Rp_m )
                endif
            enddo
        enddo

        end associate

    end subroutine calcGradIJ


    subroutine calcGradIJ_cc(Rp_m, Grid, isG, Q, gradQ)
        !! Uses Green-Gauss theorem to get cell-averaged gradient, using corner-located values
        real(rp), intent(in) :: Rp_m
            !! Planet radius in m
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: isG
            !! Mask for corners that are safe to use in calculating the gradient across the attached faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: Q
            !! Variable we are taking the gradient of across faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradQ

        integer :: i,j
        real(rp) :: qLow, qHigh 
        
        gradQ = 0.0

        associate(sh=>Grid%shGrid)

        do j=sh%jsg,sh%jeg
            do i=sh%isg,sh%ieg
                if (all(isG(i:i+1, j:j+1))) then
                    ! Theta-direction gradient
                    !qLow  = Grid%lenFace(i  ,j,RAI_TH)/2.0 * (Q(i  ,j+1) + Q(i  ,j))
                    !qHigh = Grid%lenFace(i+1,j,RAI_TH)/2.0 * (Q(i+1,j+1) + Q(i+1,j))
                    !gradQ(i,j,RAI_TH) = (qHigh - qLow) / (Grid%areaCC(i,j) * Rp_m)  ! [Q/m]
                    qLow  = 0.5 * (Q(i  ,j+1) + Q(i  ,j))
                    qHigh = 0.5 * (Q(i+1,j+1) + Q(i+1,j))
                    gradQ(i,j,RAI_TH) = (qHigh - qLow) / (Grid%lenFace(i,j,RAI_PH) * Rp_m)  ! [Q/m]

                    ! Phi-direction gradient
                    !qLow  = Grid%lenFace(i,j  ,RAI_PH)/2.0 * (Q(i+1,j  ) + Q(i, j  ))
                    !qHigh = Grid%lenFace(i,j+1,RAI_PH)/2.0 * (Q(i+1,j+1) + Q(i, j+1))
                    !gradQ(i,j,RAI_PH) = (qhigh - qLow) / (Grid%areaCC(i,j) * Rp_m)
                    qLow  = 0.5 * (Q(i+1,j  ) + Q(i, j  ))
                    qHigh = 0.5 * (Q(i+1,j+1) + Q(i, j+1))
                    gradQ(i,j,RAI_PH) = (qhigh - qLow) / (0.5*(Grid%lenFace(i,j,RAI_TH)+Grid%lenFace(i+1,j,RAI_TH)) * Rp_m)
                endif
            enddo
        enddo

        end associate

    end subroutine calcGradIJ_cc


    subroutine calcGradFTV(Rp_m, RIon_m, B0, Grid, isG, V, gradV)
        !! Calc gradient in spherical coordinates of corner variable V (flux tube volume) across cell edges/faces
        !! Returns gradient in [bVol / m]
        real(rp), intent(in) :: Rp_m
            !! Planet radius [m]]
        real(rp), intent(in) :: RIon_m
            !! Ionosphere radius [m]
        real(rp), intent(in) :: B0
            !! Planet's surface field strength [Gauss]
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: isG
            !! Mask for corners that are safe to use in calculating the gradient across the attached faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: V
            !! Flux tube volume (FTV)
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: gradV
            !! grad(FTV) we return [units(FTV)/m]

        integer :: i
        real(rp) :: dcl_dm
            !! d colat / d meters, used to convert dipole derivative w.r.t. colat to meters
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: V0, dV, dV0_dth
            !! V0 = dipole FTV, dV = V - V0

        associate(sh=>Grid%shGrid)
            
            ! Calculate dipole FTV
            do i=sh%isg, sh%ieg+1
                !V0(i,:) = DipFTV_colat(sh%th(i), B0)
                ! DipFTV_colat takes a colat at 1 Rp, make sure we use that value instead of colat in ionosphere (shGrid%th)
                V0(i,:) = DipFTV_colat(Grid%thRp(i), B0)
            enddo


            ! Break up into background + perturbation
            dV = 0.0
            where (isG)
                dV =  V - V0
            end where
            
            !where (dV < TINY)
            !    dV = 0.0
            !end where

            ! Take gradients of each
            ! Analytic gradient of dipole FTV
            dV0_dth = 0.0
            dcl_dm = 1.0/RIon_m
            do i=sh%isg, sh%ieg+1
                !! DerivDipFTV takes gradient w.r.t. theta
                !! We need to convert to be w.r.t. arc len in meters
                !dV0_dth(i,:) = DerivDipFTV(sh%th(i), B0) * dcl_dm
                dV0_dth(i,:) = DerivDipFTV(Grid%thRp(i), B0) * dcl_dm
            enddo

            ! Gradient of perturbation
            call calcGradIJ(Rp_m, Grid, isG, dV, gradV)

            ! Add on grad from dipole (Gradient of the dipole is in the theta direction, stored on the phi face)
            gradV(:,:,RAI_PH) = gradV(:,:,RAI_PH) + dV0_dth
            ! No need to do theta faces cause dipole doesn't have gradient in phi

        end associate

    end subroutine calcGradFTV


    subroutine calcGradFTV_cc(Rp_m, RIon_m, B0, Grid, isG, V, gradV, doSmoothO)
        !! Same as calcGradFTV, but calculating cell-averaged gradient instead
        real(rp), intent(in) :: Rp_m
            !! Planet radius [m]]
        real(rp), intent(in) :: RIon_m
            !! Ionosphere radius [m]
        real(rp), intent(in) :: B0
            !! Planet's surface field strength [Gauss]
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: isG
            !! Mask for corners that are safe to use in calculating the gradient across the attached faces
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: V
            !! Flux tube volume (FTV)
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradV
            !! grad(FTV) we return [units(FTV)/m]
        logical, optional, intent(in) :: doSmoothO

        integer :: i
        real(rp) :: dcl_dm
            !! d colat / d meters, used to convert dipole derivative w.r.t. colat to meters
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1) :: V0, dV
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: dV0_dth
            !! V0 = dipole FTV, dV = V - V0
        logical :: doSmooth

        if (present(doSmoothO)) then
            doSmooth = doSmoothO
        else
            doSmooth = .true.
        endif

        associate(sh=>Grid%shGrid)

        do i=sh%isg, sh%ieg+1
            ! DipFTV_colat takes a colat at 1 Rp, make sure we use that value instead of colat in ionosphere (shGrid%th)
            V0(i,:) = DipFTV_colat(Grid%thRp(i), B0)
        enddo

        dV = 0.0
        where(isG)
            dV = V - V0
        end where

        ! Analytic gradient of dipole FTV at cell centers
        dV0_dth = 0.0
        dcl_dm = 1.0/RIon_m
        do i=sh%isg,sh%ieg
            dV0_dth(i,:) = DerivDipFTV(Grid%thcRp(i), B0) * dcl_dm
        enddo

        ! Gradient of perturbation
        call calcGradIJ_cc(Rp_m, Grid, isG, dV, gradV)

        gradV(:,:,RAI_TH) = gradV(:,:,RAI_TH) + dV0_dth

        end associate

    end subroutine calcGradFTV_cc

!------
! Velocity calculations
!------

    subroutine calcVelocity(Model, Grid, State, k, Vtp)
        !! Uses gradPots stored in State to calculate interface velocity [m/s] for lambda channel k
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: Vtp
                            !! Our interface velocity that we return

        integer :: i,j
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: gradPot

        ! alamc [eV * (Rp/nT)^(2/3)]
        ! gradVM [(Rp/nT)^(-2/3) / m]
        ! alamc*gradVM [eV/m]
        ! v_driftGC = gradPot / (q * B)
        ! alamc*gradVM [eV/m] / q = [V/m]
        gradPot = State%gradPotE + State%gradPotCorot + Grid%alamc(k)*State%gradVM  ! [V/m]        

        ! Vel = ( Bvec x grad(pot) ) / B^2  = ( bhat x grad(pot) ) / B
        ! [gradPot] = [V/m] = [T*m/s]
        ! Note, 1/B term used is actually Br
        ! gradPot / B [T] = [m/s]

        ! NOTE: Physically, this is where we would take the cross product
        ! But: we have stored the gradient of the potential across Theta face, which is in the phi direction, in (i,j,RAI_TH)
        !  In other words, the cross product was already taken (except for the sign) when doing the gradients along faces
        !  So just need to include the sign here to complete the cross product
        Vtp(:,:,RAI_TH) =      gradPot(:,:,RAI_TH) / (Grid%BrFace(:,:,RAI_TH)*1.0e-9)  ! [m/s]
        Vtp(:,:,RAI_PH) = -1.0*gradPot(:,:,RAI_PH) / (Grid%BrFace(:,:,RAI_PH)*1.0e-9)  ! [m/s]

    end subroutine calcVelocity


    subroutine calcVelocityCC_gg(Model, Grid, State, k, Vtp)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: Vtp
        
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradPot

        gradPot = State%gradPotE_cc + State%gradPotCorot_cc + Grid%alamc(k)*State%gradVM_cc

        Vtp(:,:,RAI_TH) =      gradPot(:,:,RAI_PH) / (Grid%Brcc*1.0e-9)  ! [m/s]
        Vtp(:,:,RAI_PH) = -1.0*gradPot(:,:,RAI_TH) / (Grid%Brcc*1.0e-9)  ! [m/s]
    end subroutine calcVelocityCC_gg


    subroutine reconVelocityLRs(Model, Grid, State, k, iVelL, iVelR)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2), intent(inout) :: iVelL, iVelR

        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGCC
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg+1,2) :: tmpVelL, tmpVelR
        
        iVelL = 0.0
        iVelR = 0.0
        isGCC = .false.
        ! The good mask here indicates if a cell can be used in reconstruction, which for us includes ACTIVE and BUFFER zones
        where (State%active .ne. RAIJUINACTIVE)
            isGCC = .true.
        end where

        ! ReconFaces is gonna take one cc value and reconstruct at cell faces
        ! But phi velocity at theta edge is meaningless for us
        ! So we make two temporary arrays, and then just pack the meaningful components into iVelL and iVelR at the end
        call ReconFaces(Model, Grid, isGCC, State%cVel(:,:,k,RAI_TH), tmpVelL, tmpVelR)
        iVelL(:,:,RAI_TH) = tmpVelL(:,:,RAI_TH)
        iVelR(:,:,RAI_TH) = tmpVelR(:,:,RAI_TH)
        call ReconFaces(Model, Grid, isGCC, State%cVel(:,:,k,RAI_PH), tmpVelL, tmpVelR)
        iVelL(:,:,RAI_PH) = tmpVelL(:,:,RAI_PH)
        iVelR(:,:,RAI_PH) = tmpVelR(:,:,RAI_PH)


    end subroutine reconVelocityLRs

!------
! time handling
!------

    function activeDt(Model, Grid, State, k) result(dt)
        !! Calculates min dt needed to safele evolve active domain for given energy invariant channel k
        !! TODO: Consider dynamic CFL factor
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        
        integer :: i,j
        real(rp) :: dl, velij_th, velij_ph
        real(rp), dimension(Grid%shGrid%is:Grid%shGrid%ie+1,Grid%shGrid%js:Grid%shGrid%je+1) :: dtArr
        real(rp) :: dt

        associate (sh => Grid%shGrid)

        dtArr = HUGE

        do j=sh%js,sh%je+1
            do i=sh%is,sh%ie+1
                velij_th = TINY
                velij_ph = TINY
                ! NOTE: We are only checking faces bordering non-ghost cells because those are the only ones we use for evolution
                ! TODO: Strictly speaking, there are some faces that are included here that shouldn't be because they're never used to evolve anything
                !  so we should make the loop js:je, is:ie, and handle the last row and column afterwards

                ! We are responsible for face (i,j,TH) and (i,j,PH)
                ! Theta faces first

                ! This is a weird pattern. Basically, we can't cycle if the overall desired condition isn't met, because we are doing theta and phi directions in the same loop
                ! At the same time, we don't want to write one massive if condition with all options covered, or a bunch of nested if statements
                ! So instead, we break them up, such that if a single if condition is true it means we have failed the physical condition for us to calculate a valid timestep
                ! If all if conditions fail then the physical conditions are met and we can do our calculations in the else block
                if (Model%doActiveShell .and. ( .not. State%activeShells(i-1,k) .or. .not. State%activeShells(i,k)) ) then
                    ! In order for a theta-dir face to be usable, we need both sides to be active
                    continue
                else if ( State%active(i-1,j) .ne. RAIJUACTIVE .or. State%active(i,j) .ne. RAIJUACTIVE ) then
                    continue
                else
                    ! We are good lets calculate a dt for this face
                    velij_th = max(abs(State%iVel(i,j,k,RAI_TH)), TINY)  ! [m/s]
                    !dtArr(i,j,RAI_TH) = ( Grid%delTh(i) * Model%planet%ri_m ) / velij  ! [s]
                endif
                
                ! Phi faces
                if (Model%doActiveShell .and. .not. State%activeShells(i,k)) then
                    ! In order for a phi-dir face to be usable, we just need this i shell to be active
                    continue
                else if (State%active(i,j-1) .ne. RAIJUACTIVE  .or. State%active(i,j) .ne. RAIJUACTIVE ) then 
                    continue
                else
                    velij_ph = max(abs(State%iVel(i,j,k,RAI_PH)), TINY)
                    !dtArr(i,j,RAI_PH) = ( Grid%delPh(j) * Model%planet%ri_m ) / velij  ! [s]
                endif

                dl = min(Grid%delTh(i), Grid%delPh(j)*sin(sh%thc(i)))
                dtArr(i,j) = (dl*Model%planet%ri_m) / sqrt(velij_th**2 + velij_ph**2)

            enddo
        enddo

        dt = Model%CFL*minval(dtArr)

        end associate
            
    end function activeDt


!------
! Extras
!------

    subroutine velFace2CC(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: i,j,k

        associate(sh=>Grid%shGrid)

        !State%cVel(:,:,:,RAI_TH) = 0.5*(State%iVel(Grid%shGrid%isg+1:,:,:,RAI_TH) + State%iVel(:Grid%shGrid%ieg,:,:,RAI_TH))
        !State%cVel(:,:,:,RAI_PH) = 0.5*(State%iVel(:,Grid%shGrid%jsg+1:,:,RAI_PH) + State%iVel(:,:Grid%shGrid%jeg,:,RAI_PH))
        do i=sh%isg,sh%ieg
            do j=sh%jsg,sh%jeg
                do k=1,Grid%Nk
                    State%cVel(i,j,k,RAI_TH) = 0.5*(State%iVel(i+1,j  ,k,RAI_TH)+State%iVel(i,j,k,RAI_TH))
					State%cVel(i,j,k,RAI_PH) = 0.5*(State%iVel(i  ,j+1,k,RAI_PH)+State%iVel(i,j,k,RAI_PH))
                enddo
            enddo
        enddo
        end associate

    end subroutine velFace2CC

!------
! Graveyard
!------

!    subroutine calcGradIJ_cc(RIon, Grid, isG, Q, gradQ)
!        !! Calc gradint in spherical coordinates across entire grid, including ghosts
!        !! Up to someone else to overwrite ghosts
!        real(rp), intent(in) :: RIon
!            !! Ionosphere radius in m
!        type(raijuGrid_T), intent(in) :: Grid
!        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isG
!        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: Q
!        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradQ
!            !! I think [unitsQ/rad]
!
!        integer :: i,j
!        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg) :: sinTh
!        
!        gradQ = 0.0
!        
!        associate(sh=>Grid%shGrid)
!        ! Calc sin(theta) for everyone to use
!        sinTh = sin(sh%thc)
!
!        ! Leave out the end rows and columns
!        !$OMP PARALLEL DO default(shared) collapse(1) &
!        !$OMP schedule(dynamic) &
!        !$OMP private(i,j) &
!        !$OMP IF(.false.)
!        do j=sh%jsg+1,sh%jeg-1
!            do i=sh%isg+1,sh%ieg-1
!                if (.not. isG(i,j)) then
!                    cycle
!                end if
!
!                ! Theta dir
!                if(isG(i-1,j) .and. isG(i+1,j)) then
!                    gradQ(i,j,RAI_TH) = (Q(i+1,j) - Q(i-1,j))/RIon/(Grid%delTh(i)+Grid%delTh(i+1))
!                else if (isG(i-1,j)) then
!                    gradQ(i,j,RAI_TH) = (Q(i,j)-Q(i-1,j))/RIon/Grid%delTh(i)
!                else if (isG(i+1,j)) then
!                    gradQ(i,j,RAI_TH) = (Q(i+1,j)-Q(i,j))/RIon/Grid%delTh(i+1)
!                end if
!
!                ! Phi dir
!                if(isG(i,j-1) .and. isG(i,j+1)) then
!                    gradQ(i,j,RAI_PH) = (Q(i,j+1) - Q(i,j-1))/RIon/sinTh(i)/(Grid%delPh(j)+Grid%delPh(j+1))
!                else if (isG(i,j-1)) then
!                    gradQ(i,j,RAI_PH) = (Q(i,j)-Q(i,j-1))/RIon/sinTh(i)/Grid%delPh(j)
!                else if (isG(i,j-1)) then
!                    gradQ(i,j,RAI_PH) = (Q(i,j+1)-Q(i,j))/RIon/sinTh(i)/Grid%delPh(j+1)
!                end if
!            enddo
!        enddo
!
!        ! Edge rows and columns
!        ! Set to same gradient so that resulting velocity is the same
!        ! Doesn't really matter since they won't be included in reconstructions anyways
!        gradQ(sh%isg,:,:) = gradQ(sh%isg+1,:,:)
!        gradQ(sh%ieg,:,:) = gradQ(sh%ieg-1,:,:)
!        gradQ(:,sh%jsg,:) = gradQ(:,sh%jsg+1,:)
!        gradQ(:,sh%jeg,:) = gradQ(:,sh%jeg-1,:)
!
!        end associate
!    end subroutine calcGradIJ_cc


!    subroutine calcVelocityCC(Model, Grid, State, k, Vtp)
!        !! Uses gradPots stored in State to calculate cell-centered velocity [m/s] for lambda channel k
!        type(raijuModel_T), intent(in) :: Model
!        type(raijuGrid_T ), intent(in) :: Grid
!        type(raijuState_T), intent(in) :: State
!        integer, intent(in) :: k
!        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
!                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: Vtp
!
!        integer :: i,j
!        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg) :: cosdip
!        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradPot
!
!
!        
!        gradPot = State%gradPotE + State%gradPotCorot + Grid%alamc(k)*State%gradVM  ! [V/m]        
!
!        ! Vel = ( Bvec x grad(pot) ) / B^2  = ( bhat x grad(pot) ) / B
!        ! [gradPot] = [V/m] = [T*m/s]
!        ! Note, 1/B term includes dip angle
!        !associate (colat => Grid%shGrid%thc)
!        !do j=Grid%shGrid%jsg, Grid%shGrid%jeg
!        !    !!NOTE: Assuming dipole field dip angle
!        !    cosdip(:,j) = 2.0*cos(colat)/sqrt(1.0 + 3.0*cos(colat)**2.0)
!        !enddo
!        !end associate
!
!        Vtp(:,:,RAI_TH) =      gradPot(:,:,RAI_PH) / Grid%cosdip / (Grid%Bmag*1.0e-9)  ! [m/s]
!        Vtp(:,:,RAI_PH) = -1.0*gradPot(:,:,RAI_TH) / Grid%cosdip / (Grid%Bmag*1.0e-9)  ! [m/s]
!
!    end subroutine calcVelocityCC



end module raijuPreAdvancer