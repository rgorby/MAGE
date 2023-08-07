module raijuadvancer
    use clocks
    use planethelper
    use kai2geo

    use raijudefs
    use raijutypes
    use raijuetautils
    use raijuBCs
    use raijulosses

    implicit none

    contains

!------
! Main high-level functions
!------

    subroutine raijuPreAdvance(Model, Grid, State, fullEtaMapO)
        !! Takes a state and calculates what is needed in order to advance
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        logical, optional, intent(in) :: fullEtaMapO

        logical :: fullEtaMap
        integer :: k

        if (present(fullEtaMapO)) then
            fullEtaMap = fullEtaMapO
        else
            fullEtaMap = .false.
        endif

        ! Moments to etas
        call Tic("BCs")
        call applyRaijuBCs(Model, Grid, State, fullEtaMap) ! If fullEtaMap=True, mom2eta map is applied to the whole domain
        call Toc("BCs")

        ! Calc cell velocities
        call Tic("Calc cell-center velocities")
        call calcPotGrads(Model, Grid, State)
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call calcVelocityCC(Model, Grid, State, k, State%cVel(:,:,k,:))
            ! Calc sub-time step
            State%dtk(k) = activeDt(Grid%shGrid, Grid, State, k)
        enddo
        call Toc("Calc cell-center velocities")

    end subroutine raijuPreAdvance

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
        !! Units (gradVM): V / m / lambda
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
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
        call calcGradFTV(Model%planet%ri_m, Grid, isGood, State%bvol, gradVM)
        State%gradVM(:,:,1) = (-2./3.) * State%bvol**(-5./3.) * gradVM(:,:,1)  ! [V/m/lambda]
        State%gradVM(:,:,2) = (-2./3.) * State%bvol**(-5./3.) * gradVM(:,:,2)  ! [V/m/lambda]


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
        !$OMP PARALLEL DO default(shared) collapse(2) &
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
                    gradQ(i,j,1) = (Q(i+1,j) - Q(i-1,j))/RIon/(Grid%delTh(i)+Grid%delTh(i+1))
                else if (isG(i-1,j)) then
                    gradQ(i,j,1) = (Q(i,j)-Q(i-1,j))/RIon/Grid%delTh(i)
                else if (isG(i+1,j)) then
                    gradQ(i,j,1) = (Q(i+1,j)-Q(i,j))/RIon/Grid%delTh(i+1)
                end if

                ! Phi dir
                if(isG(i,j-1) .and. isG(i,j+1)) then
                    gradQ(i,j,2) = (Q(i,j+1) - Q(i,j-1))/RIon/sinTh(i)/(Grid%delPh(j)+Grid%delPh(j+1))
                else if (isG(i,j-1)) then
                    gradQ(i,j,2) = (Q(i,j)-Q(i,j-1))/RIon/sinTh(i)/Grid%delPh(j)
                else if (isG(i,j-1)) then
                    gradQ(i,j,2) = (Q(i,j+1)-Q(i,j))/RIon/sinTh(i)/Grid%delPh(j+1)
                end if
            enddo
        enddo

        ! Edge rows and columns
        ! Set to same gradient so that resulting velocity is the same
        gradQ(sh%isg,:,:) = gradQ(sh%isg+1,:,:)
        gradQ(sh%ieg,:,:) = gradQ(sh%ieg-1,:,:)
        gradQ(:,sh%jsg,:) = gradQ(:,sh%jsg+1,:)
        gradQ(:,sh%jeg,:) = gradQ(:,sh%jeg-1,:)

        end associate
    end subroutine calcGradIJ

    subroutine calcGradFTV(RIon, Grid, isG, V, gradV)
        real(rp), intent(in) :: RIon
            !! Iono radius in Rp
        type(raijuGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isG
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: V
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradV

        integer :: i
        real(rp) :: dcl_dth
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: V0, dV, dV0_dth
            !! V0 = dipole V, dV = V - V0

        associate(sh=>Grid%shGrid)
            
            ! Calculate dipole FTV
            do i=sh%isg, sh%ieg
                V0(i,:) = DipFTV_colat(sh%thc(i))
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
                dcl_dth = 0.5*(sh%thc(i-1) - sh%thc(i+1))
                dV0_dth(i,:) = DerivDipFTV(sh%thc(i)) * dcl_dth
            enddo

            ! Gradient of perturbation
            call calcGradIJ(RIon, Grid, isG, dV, gradV)

            ! Add on grad from dipole
            gradV(:,:,1) = gradV(:,:,1) + dV0_dth
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
        associate (colat => Grid%shGrid%thc)
        do j=Grid%shGrid%jsg, Grid%shGrid%jeg
            !!NOTE: Assuming dipole field dip angle
            cosdip(:,j) = 2.0*cos(colat)/sqrt(1.0 + 3.0*cos(colat)**2.0)
        enddo
        end associate

        Vtp(:,:,1) =      gradPot(:,:,2) / cosdip / (Grid%Bmag*1.0e-9)  ! [m/s]
        Vtp(:,:,2) = -1.0*gradPot(:,:,1) / cosdip / (Grid%Bmag*1.0e-9)  ! [m/s]

    end subroutine calcVelocityCC


    function activeDt(sh, Grid, State, k) result(dt)
        !! TODO: eventually this should use iVel instead
        !! TODO: Consider dynamic CFL factor
        type(ShellGrid_T), intent(in) :: sh
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        
        integer :: i,j
        real(rp), dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg,2) :: vel, dtArr
        logical, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg) :: asIJ, isGood
        real(rp) :: dt

        ! make activeShell 2D
        do i=sh%isg,sh%ieg
            asIJ(i,:) = State%activeShells(i,k)
        enddo
        
        !vMag = norm2(State%cVel(:,:,k,:), 3)
        ! Only consider points that are in activeShell and not RAIJUINACTIVE
        where (State%active .ne. RAIJUINACTIVE .and. asIJ)
            isGood = .true.
        elsewhere
            isGood = .false.
        end where

        vel(:,:,1) = merge(abs(State%cVel(:,:,k,1)), TINY, isGood)
        vel(:,:,2) = merge(abs(State%cVel(:,:,k,2)), TINY, isGood)

        ! Min dt in theta direction
        do j=sh%jsg,sh%jeg
            dtArr(:,j,1) = 0.5*(Grid%delTh(:sh%ieg)+Grid%delTh(sh%isg+1:)) / vel(:,j,1)
        enddo
        ! In Phi direction
        do i=sh%isg,sh%ieg
            dtArr(i,:,2) = 0.5*(Grid%delPh(:sh%jeg)+Grid%delPh(sh%jsg+1:)) / vel(i,:,2)
        enddo

        dt = minval(dtArr)

            
    end function activeDt


!------
! Advance entry point
!------

    subroutine AdvanceState(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: k

        ! Clear things that will be accumulated over the advance
        State%precipNFlux = 0.0
        State%precipEFlux = 0.0

        ! Make sure moments are up to date, some things (like coulomb losses) need them
        call EvalMoments(Grid, State)

        ! Send each channel off on their own
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call AdvanceLambda(Model, Grid, State, k)
        enddo
    end subroutine AdvanceState


    subroutine AdvanceLambda(Model, Grid, State, k)
        !! Advances a single lambda channel from State time t to t+dt
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k

        integer :: n, s , Nmax
            !! n = step counter
            !! s = species index
        real(rp) :: t, dt, tEnd
            !! t, dt = Current time and delta-t for this channel
            !! tEnd = Time we stop at

        
        ! Initial settings
        s = Grid%k2spc(k)
        t = State%t
        tEnd = State%t + State%dt

        !Nsteps = int(State%dt / State%dtk(k))+1
        !dt = State%dt / (1.0_rp*Nsteps)
        Nmax = 200
        associate(sh=>Grid%shGrid, spc=>Grid%spc(s))

            ! Here we go!
            n = 1  ! counter
            do while ( tEnd-t > TINY)
                !write(*,*)k,n,t,dt
                ! Calc new active shells
                
                !! Just mute for now, we're not actually using it
                !State%activeShells(:,k) = setLambdaActiveShells(sh, spc, State%bVol, &
                !                            State%eta(:,:,spc%kStart:spc%kEnd), k, worthyFracO = Model%worthyFrac)
                

                ! Calc next time step
                !! Also muting
                !dt = activeDt(sh, Grid, State, k)

                !! BAD: Boost dt to be around 200 iters max
                !dt = max(dt, (tEnd-t)/(Nmax-n))
                dt = (tEnd-t)/(Nmax-n)
                ! If needed, reduce dt to fit within remaining time
                if (t + dt > tEnd) then
                    dt = tEnd - t
                endif

                ! Advection
                
                ! Losses
                if (Model%doLosses) then
                    call calcStepLosses(Model, Grid, State, k, dt)
                endif

                t = t + dt
                n = n+1
            enddo

            if (Model%doLosses) then
                ! Divide precip fluxes by big dt to turn them into proper rates
                State%precipNFlux(:,:,k) = State%precipNFlux(:,:,k)/State%dt
                State%precipEFlux(:,:,k) = State%precipEFlux(:,:,k)/State%dt
            endif

        end associate

    end subroutine AdvanceLambda



end module raijuadvancer