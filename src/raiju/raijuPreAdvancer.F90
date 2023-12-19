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
    use raijuAdvancer, only : activeDt

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
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG_vFaces
            !! We calculate which cells are good for the purpose of reconstructing velocity at all faces

        if (present(isFirstCplO)) then
            isFirstCpl = isFirstCplO
        else
            isFirstCpl = .false.
        endif

        ! Calc isG for velocity reconstruction
        where (State%active .ne. RAIJUINACTIVE)
            isG_vFaces = .true.
        elsewhere
            isG_vFaces = .false.
        end where

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
        call Tic("Calc cell-center velocities")
        call calcPotGrads(Model, Grid, State)
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call calcVelocityCC(Model, Grid, State, k, State%cVel(:,:,k,:))
            ! Calc sub-time step
            State%dtk(k) = activeDt(Model, Grid, State, k)
            ! We can also calculate velocity at faces here because it won't change during sub-stepping
            call ReconFaces(Grid, isG_vFaces, State%cVel(:,:,k,1), State%iVel(:,:,k,:), Qcc_phO=State%cVel(:,:,k,2))
        enddo
        call Toc("Calc cell-center velocities")

        ! Loss rate calc depends on up-to-date densities, so we should run EvalMoments first
        call Tic("Moments Eval PreAdvance")
        call EvalMoments(Grid, State)
        call Toc("Moments Eval PreAdvance")

        call Tic("Calc loss rates")
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call calcChannelLossRates(Model, Grid, State, k)
        enddo
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
    !! Potential variables come in pre-allocated, since they may be subests of larger arrays
    subroutine potExB(sh, State, pExB)
        ! Trivial, but putting it here in case we need extra options for it later
        type(ShellGrid_T), intent(in) :: sh
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(inout) :: pExB
        
        pExB = State%espot * 1.e3  ! [kV -> V]

    end subroutine potExB

    subroutine potCorot(planet, sh, pCorot, doGeoCorotO)
        ! Simple corotation potential [V]
        type(planet_T), intent(in) :: planet
        type(ShellGrid_T), intent(in) :: sh
        logical, intent(in), optional :: doGeoCorotO

        real(rp), dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(inout) :: pCorot
        integer :: i,j

        if (present(doGeoCorotO)) then
            if (doGeoCorotO) then
                pCorot = 0.0
                do j=sh%jsg,sh%jeg
                    do i=sh%isg,sh%ieg
                        call geocorotation_from_SMTP(sh%thc(i), sh%phc(j), pCorot(i,j))
                    enddo
                enddo
                ! geopack corotation returns in kV, we need V
                pCorot = pCorot * 1.e3  ! [kV -> V]
                return
            endif
        endif

        ! IfdoGeoCorotO was not provided, or it was false, we default to corotation on aligned dipole and rotational axis

        do j=sh%jsg,sh%jeg
            pCorot(:,j) = -planet%psiCorot*(planet%rp_m/planet%ri_m)*sin(sh%thc)**2 * 1.e3  ! [kV -> V]
        enddo

    end subroutine potCorot


    subroutine potGC(shGrid, alamc, bVol, pGC)
        ! Simple gradient curvature potential [V]
        type(ShellGrid_T), intent(in) :: shGrid
        real(rp), intent(in) :: alamc
        real(rp), dimension(shGrid%isg:shGrid%ieg,&
                            shGrid%jsg:shGrid%jeg), intent(in) :: bVol

        real(rp), dimension(shGrid%isg:shGrid%ieg,&
                            shGrid%jsg:shGrid%jeg), intent(inout) :: pGC

        pGC = alamc*bVol**(-2./3.)

    end subroutine potGC


    subroutine calcEffectivePotential(Model, Grid, State, pEff)
        !! Calculates effective potential [V] for all lambda channels as the sum of pExB, pCorot, and pGC
        !! Note: This is not used to calculate velocities
        type(raijuModel_T) , intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg,&
                            Grid%Nk), intent(inout) :: pEff
        
        integer :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: pExB, pCorot, pGC
        
        pEff = 0.0

        ! Grab 2D potentials
        call potExB(Grid%shGrid, State, pExB)
        call potCorot(Model%planet, Grid%shGrid, pCorot, Model%doGeoCorot)
        
        ! Build 3D effective potential
        !$OMP PARALLEL DO default(shared) collapse(1) &
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
        !! Units (gradVM): Vol / m / lambda
        type(raijuModel_T), intent(in)    :: Model
        type(raijuGrid_T ), intent(in)    :: Grid
        type(raijuState_T), intent(inout) :: State

        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: pExB, pCorot
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradVM
        ! What's good
        where (State%active .ne. RAIJUINACTIVE)
            isGood = .true.
        elsewhere
            isGood = .false.
        end where
        
        ! Ionospheric and corotation potentials are just simple derivatives
        call potExB(Grid%shGrid, State, pExB)  ! [V]
        call potCorot(Model%planet, Grid%shGrid, pCorot, Model%doGeoCorot)  ! [V]
        call calcGradIJ(Model%planet%ri_m, Grid, isGood, pExB  , State%gradPotE    )  ! [V/m]
        call calcGradIJ(Model%planet%ri_m, Grid, isGood, pCorot, State%gradPotCorot)  ! [V/m]

        ! GC drifts depend on grad(lambda * V^(-2/3))
        ! lambda is constant, so just need grad(V^(-2/3) )
        ! grad(V^(-2/3)) = -2/3*V^(-5/3) * grad(V)
        call calcGradFTV(Model%planet%ri_m, Model%planet%magMoment, Grid, isGood, State%bvol, gradVM)
        State%gradVM(:,:,RAI_TH) = (-2./3.) * State%bvol**(-5./3.) * gradVM(:,:,RAI_TH)  ! [Vol/m/lambda]
        State%gradVM(:,:,RAI_PH) = (-2./3.) * State%bvol**(-5./3.) * gradVM(:,:,RAI_PH)  ! [Vol/m/lambda]

    end subroutine calcPotGrads


    subroutine calcGradIJ(RIon, Grid, isG, Q, gradQ)
        !! Calc gradint in spherical coordinates across entire grid, including ghosts
        !! Up to someone else to overwrite ghosts
        real(rp), intent(in) :: RIon
            !! Ionosphere radius in m
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isG
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: Q
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradQ
            !! I think [unitsQ/rad]

        integer :: i,j
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg) :: sinTh
        
        gradQ = 0.0
        
        associate(sh=>Grid%shGrid)
        ! Calc sin(theta) for everyone to use
        sinTh = sin(sh%thc)

        ! Leave out the end rows and columns
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j) &
        !$OMP IF(.false.)
        do j=sh%jsg+1,sh%jeg-1
            do i=sh%isg+1,sh%ieg-1
                if (.not. isG(i,j)) then
                    cycle
                end if

                ! Theta dir
                if(isG(i-1,j) .and. isG(i+1,j)) then
                    gradQ(i,j,RAI_TH) = (Q(i+1,j) - Q(i-1,j))/RIon/(Grid%delTh(i)+Grid%delTh(i+1))
                else if (isG(i-1,j)) then
                    gradQ(i,j,RAI_TH) = (Q(i,j)-Q(i-1,j))/RIon/Grid%delTh(i)
                else if (isG(i+1,j)) then
                    gradQ(i,j,RAI_TH) = (Q(i+1,j)-Q(i,j))/RIon/Grid%delTh(i+1)
                end if

                ! Phi dir
                if(isG(i,j-1) .and. isG(i,j+1)) then
                    gradQ(i,j,RAI_PH) = (Q(i,j+1) - Q(i,j-1))/RIon/sinTh(i)/(Grid%delPh(j)+Grid%delPh(j+1))
                else if (isG(i,j-1)) then
                    gradQ(i,j,RAI_PH) = (Q(i,j)-Q(i,j-1))/RIon/sinTh(i)/Grid%delPh(j)
                else if (isG(i,j-1)) then
                    gradQ(i,j,RAI_PH) = (Q(i,j+1)-Q(i,j))/RIon/sinTh(i)/Grid%delPh(j+1)
                end if
            enddo
        enddo

        ! Edge rows and columns
        ! Set to same gradient so that resulting velocity is the same
        ! Doesn't really matter since they won't be included in reconstructions anyways
        gradQ(sh%isg,:,:) = gradQ(sh%isg+1,:,:)
        gradQ(sh%ieg,:,:) = gradQ(sh%ieg-1,:,:)
        gradQ(:,sh%jsg,:) = gradQ(:,sh%jsg+1,:)
        gradQ(:,sh%jeg,:) = gradQ(:,sh%jeg-1,:)

        end associate
    end subroutine calcGradIJ


    subroutine calcGradFTV(RIon, B0, Grid, isG, V, gradV)
        !! Calculates derivative of flux tube volume (bvol) in units [bVol / m]
        real(rp), intent(in) :: RIon
            !! Iono radius in meters
        real(rp), intent(in) :: B0
            !! Planet's surface field strength
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isG
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: V
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradV

        integer :: i
        real(rp) :: dcl_dm
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: V0, dV, dV0_dth
            !! V0 = dipole V, dV = V - V0

        associate(sh=>Grid%shGrid)
            
            ! Calculate dipole FTV
            do i=sh%isg, sh%ieg
                V0(i,:) = DipFTV_colat(sh%thc(i), B0)
            enddo

            ! Break up into background + perturbation
            dV = 0.0
            where (isG)
                dV =  V - V0
            end where

            ! Take gradients of each
            ! Analytic gradient of dipole FTV
            dV0_dth = 0.0
            do i=sh%isg+1, sh%ieg-1
                dcl_dm = 0.5*(sh%thc(i-1) - sh%thc(i+1))/RIon
                dV0_dth(i,:) = DerivDipFTV(sh%thc(i), B0) * dcl_dm
            enddo

            ! Gradient of perturbation
            call calcGradIJ(RIon, Grid, isG, dV, gradV)

            ! Add on grad from dipole
            gradV(:,:,RAI_TH) = gradV(:,:,RAI_TH) + dV0_dth
            ! No need to do phi direction cause dipole doesn't have gradient in phi

        end associate

    end subroutine calcGradFTV



!------
! Velocity calculations
!------
    subroutine calcVelocityCC(Model, Grid, State, k, Vtp)
        !! Uses gradPots stored in State to calculate cell-centered velocity [m/s] for lambda channel k
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: Vtp

        integer :: i,j
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg) :: cosdip
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradPot


        
        gradPot = State%gradPotE + State%gradPotCorot + Grid%alamc(k)*State%gradVM  ! [V/m]        

        ! Vel = ( Bvec x grad(pot) ) / B^2  = ( bhat x grad(pot) ) / B
        ! [gradPot] = [V/m] = [T*m/s]
        ! Note, 1/B term includes dip angle
        !associate (colat => Grid%shGrid%thc)
        !do j=Grid%shGrid%jsg, Grid%shGrid%jeg
        !    !!NOTE: Assuming dipole field dip angle
        !    cosdip(:,j) = 2.0*cos(colat)/sqrt(1.0 + 3.0*cos(colat)**2.0)
        !enddo
        !end associate

        Vtp(:,:,RAI_TH) =      gradPot(:,:,RAI_PH) / Grid%cosdip / (Grid%Bmag*1.0e-9)  ! [m/s]
        Vtp(:,:,RAI_PH) = -1.0*gradPot(:,:,RAI_TH) / Grid%cosdip / (Grid%Bmag*1.0e-9)  ! [m/s]

    end subroutine calcVelocityCC



end module raijuPreAdvancer