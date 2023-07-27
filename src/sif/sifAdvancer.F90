module sifadvancer
    use clocks
    use planethelper

    use sifdefs
    use siftypes
    use sifetautils
    use siflosses

    implicit none

    contains


!------
! Active I Shell calculations
!------

    subroutine setActiveShellsByContribution(shGrid, spc, bVol, etas, activeShellsOut, nSpacesO, worthyFracO)
        !! For each lambda channel, calculate active shells based on how much they contribute to the total pressure and/or density
        type(ShellGrid_T), intent(in) :: shGrid
        type(SIFSpecies_T), intent(in) :: spc
        real(rp), dimension(shGrid%isg:shGrid%ieg,shGrid%jsg:shGrid%jeg), intent(in) :: bVol
        real(rp), dimension(shGrid%isg:shGrid%ieg,shGrid%jsg:shGrid%jeg,spc%kStart:spc%kEnd), intent(in) :: etas
        logical, dimension(shGrid%isg:shGrid%ieg, spc%kStart:spc%kEnd), intent(inout) :: activeShellsOut
        integer, optional, intent(in) :: nSpacesO
            !! Number of i spaces between last good value and active i for species
        real(rp), optional, intent(in) :: worthyFracO

        integer :: i,j,k, iL, iU
        integer :: nSpaces
        real(rp) :: worthyFrac
        real(rp) :: alamc, kDen, kPress
        real(rp), dimension(shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg) :: spcDen, spcPress
        logical, dimension(shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg) :: as2D

        if(present(nSpacesO)) then
            nSpaces = nSpacesO
        else
            nSpaces = nSpacesDef
        endif

        if(present(worthyFracO)) then
            worthyFrac = worthyFracO
        else
            worthyFrac = fracWorthyDef
        endif
        
        as2D = .false.
        activeShellsOut = .false.

        ! Start by getting total density and pressure for each i,j point
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j)
        do j=shGrid%jsg,shGrid%jeg
            do i=shGrid%isg,shGrid%ieg
                spcDen  (i,j) = SpcEta2Den  (spc, etas(i,j,:), bVol(i,j))
                spcPress(i,j) = spcEta2Press(spc, etas(i,j,:), bVol(i,j))
            enddo
        enddo

        ! Then calculate active shells for each k
        do k=spc%kStart,spc%kEnd

            ! Setup for this k
            as2D = .false.
            alamc = 0.5*abs(spc%alami(k) + spc%alami(k+1))

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,kDen,kPress)
            do j=shGrid%jsg,shGrid%jeg
                do i=shGrid%isg,shGrid%ieg

                    kDen   = etak2Den  (etas(i,j,k), bVol(i,j))
                    kPress = etak2Press(etas(i,j,k), alamc, bVol(i,j))
            
                    if (kDen/spcDen(i,j) > worthyFrac .or. kPress/spcPress(i,j) > worthyFrac) then
                        as2D(i,j) = .true.
                    endif
                enddo
            enddo

            ! Now collapse j
            !! Re-use as2D, first j element 
            do i=shGrid%isg,shGrid%ieg
                as2D(i,shGrid%jsg) = any(as2D(i,:))
            enddo

            ! Finally, assign to activeShells, including buffer spaces
            do i=shGrid%isg,shGrid%ieg
                iL = max(i-nSpaces, shGrid%isg)
                iU = min(i+nSpaces, shGrid%ieg)
                activeShellsOut(i,k) = any(as2D(iL:iU, shGrid%jsg))
            enddo

        enddo
    end subroutine setActiveShellsByContribution


!------
! Cell Potential and Velocity calculations
!------
    !! Potential variables come in pre-allocated, since they may be subests of larger arrays
    subroutine potExB(sh, State, pExB)
        ! Trivial, but putting it here in case we need extra options for it later
        type(ShellGrid_T), intent(in) :: sh
        type(sifState_T), intent(in) :: State
        real(rp), dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(inout) :: pExB
        
        pExB = State%espot * 1.e3  ! [kV -> V]

    end subroutine potExB

    subroutine potCorot(planet, sh, pCorot)
        ! Simple corotation potential [V]
        type(planet_T), intent(in) :: planet
        type(ShellGrid_T), intent(in) :: sh

        real(rp), dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(inout) :: pCorot
        integer :: j

        do j=sh%jsg,sh%jeg
            pCorot(:,j) = -planet%psiCorot*(planet%rp_m/planet%ri_m)*sin(sh%th)**2 * 1.e3  ! [kV -> V]
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


    subroutine calcEffectivePotential(planet, Grid, State)
        ! Calculates effective potential [V] for all lambda channels as the sum of pExB, pCorot, and pGC
        type(planet_T), intent(in) :: planet
        type(sifGrid_T) , intent(in) :: Grid
        type(sifState_T), intent(inout) :: State
        
        integer :: k
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg,&
                            Grid%Nk) :: pEff
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: pExB, pCorot, pGC
        
        State%pEff = 0.0

        ! Grab 2D potentials
        call potExB(Grid%shGrid, State, pExB)
        call potCorot(planet, Grid%shGrid, pCorot)
        
        ! Build 3D effective potential
        do k=1,Grid%Nk
            call potGC(Grid%shGrid, Grid%alamc(k), State%bVol, pGC)
            State%pEff(:,:,k) = pExB + pCorot + pGC
        enddo

    end subroutine calcEffectivePotential

    subroutine calcPotGrads(Model, Grid, State)
        !! Calcs gradient of ionospheric potential, corotation potential, and flux tube volume raised to -2/3
        !! Only needs to be called when any of the above quantities are updated (e.g. beginning of every coupling timestep)
        !! Note: FTV is not a potential yet, alamc*bVol**(-2./3.) is
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T ), intent(in) :: Grid
        type(sifState_T), intent(inout) :: State

        real(rp) :: RIon
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: pExB, pCorot
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradV
        ! What's good
        where (State%active .ne. SIFINACTIVE)
            isGood = .true.
        elsewhere
            isGood = .false.
        end where
        RIon = Model%planet%ri_m/Model%planet%Rp_m  ! Iono radius [Rp]

        
        ! Ionospheric and corotation potentials are just simple derivatives
        call potExB(Grid%shGrid, State, pExB)
        call potCorot(Model%planet, Grid%shGrid, pCorot)
        call calcGradIJ(RIon, Grid, isGood, pExB  , State%gradPotE    )
        call calcGradIJ(RIon, Grid, isGood, pCorot, State%gradPotCorot)

        ! GC drifts depend on grad(lambda * V^(-2/3))
        ! lambda is constant, so just need grad(V^(-2/3) )
        ! grad(V^(-2/3)) = -2/3*V^(-5/3) * grad(V)
        call calcGradFTV(RIon, Grid, isGood, State%bvol, gradV)
        State%gradVM(:,:,1) = (-2./3.) * State%bvol**(-5./3.) * gradV(:,:,1)
        State%gradVM(:,:,2) = (-2./3.) * State%bvol**(-5./3.) * gradV(:,:,2)


    end subroutine calcPotGrads



    function calcVelocityCC_DIP(Model, Grid, State, alamc) result (Vtp)
        !! Note: takes a 2D potential for an arbitrary lambda channel
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T ), intent(in) :: Grid
        type(sifState_T), intent(in) :: State
        real(rp), intent(in) :: alamc
            !! [ev*(Rp/nT)^(2/3)]

        integer :: i,j
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg) :: rmincc
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: Vtp

       

        ! But for now we will cheat and just start with typical drift velocity in a dipole
        Vtp = 0.0

        ! Drift velocity in eq plane
        associate(sh=>Grid%shGrid, active=>State%active)
            rmincc = norm2(State%xyzMincc, dim=3)  ! [Rp]

            ! No velocity in theta, only phi
            Vtp(:,:,2) = 0.016 * (State%bvol**(-2./3.)*alamc) * rmincc**2 / Model%planet%rp_m ! [m/s -> Rp/s]


            where (active .ne. SIFINACTIVE)
                isGood = .true.
            elsewhere
                isGood = .false.
            end where

            ! Map to ionosphere
            do i=sh%isg,sh%ieg
                do j=sh%js,sh%je
                    ! Loop over just active j cells, will wrap later
                    if (isGood(i,j) .eq. .false.) then
                        Vtp(i,j,2) = 0.0
                    else
                        ! Doesn't matter if adjacent points aren't good, still use their location for mapping
                        ! Not true for open adjacent neighbors though
                        Vtp(i,j,2) = Vtp(i,j,2) * (Grid%delPh(j)+Grid%delPh(j+1)) / abs(rmincc(i,j-1)-rmincc(i,j+1))  ! [rad/s]? Divide by 2 drops out of both
                    end if
                enddo
            enddo

        end associate

    end function calcVelocityCC_DIP

    function calcVelocityCC(Model, Grid, State, k) result (Vtp)
        !! Note: takes a 2D potential for an arbitrary lambda channel
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T ), intent(in) :: Grid
        type(sifState_T), intent(in) :: State
        integer, intent(in) :: k

        integer :: i,j
        real(rp) :: RIon
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg) :: pE, pExB, pCorot, pGC
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: gradPotE, gradPotGC
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg, 2) :: Vtp
         ! For cell centers, just set velocities based on cells' active mode
        ! Later when we calc interfaces, consider inactive2buffer, buffer2active, etc
        
        !! Note: Need to calc ghost cells in i direction
        !!  but not in j, because they will be overwritten with active cell info when wrapping

        ! Vel = ( Bvec x grad(pot) ) / B^2

        ! Get individual effective potentials

        ! Take grad of pExB and pCorot
        ! Take grad of pGC
        ! Calc total grad(pot) in ionosphere
        ! Calc velocity using B_iono
        Vtp = 0.0

        ! What's good
        where (State%active .ne. SIFINACTIVE)
            isGood = .true.
        elsewhere
            isGood = .false.
        end where
        RIon = Model%planet%ri_m/Model%planet%Rp_m  ! Iono radius [Rp]
        
        call potExB(Grid%shGrid, State, pExB)
        call potCorot(Model%planet, Grid%shGrid, pCorot)
        pE = pExB + pCorot ! Iono + corotation fields
        call potGC(Grid%shGrid, Grid%alamc(k), State%bVol, pGC) ! Gradient-curvature for this lambda
        

        ! Do E-field grad first
        call calcGradIJ(RIon, Grid, isGood, pE, gradPotE)
        !!FIXME: Should implement grad specifically for GC, skipping for now
        call calcGradIJ(RIon, Grid, isGood, pGC, gradPotGC)
        ! Just combine all into gradPotE
        gradPotE = gradPotE + gradPotGC

        ! [gradPot] = [V/rad] = [T*m^2/s/rad]
        ! Multiply by (rad/m)^2 to get [T*rad/s], so we can just divide by B[T] and we done
        gradPotE = gradPotE / Model%planet%ri_m**2  ! I think

        ! Calc velocity
        Vtp(:,:,1) =      gradPotE(:,:,2) / (Grid%Bmag * 1.0e-9)  ! [rad/s]
        Vtp(:,:,2) = -1.0*gradPotE(:,:,1) / (Grid%Bmag * 1.0e-9)  ! [rad/s]


    end function calcVelocityCC


    subroutine calcGradFTV(RIon, Grid, isG, V, gradV)
        real(rp), intent(in) :: RIon
            !! Iono radius in Rp
        type(sifGrid_T), intent(in) :: Grid
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


    subroutine calcGradIJ(RIon, Grid, isG, Q, gradQ)
        ! Calc gradint in spherical coordinates across entire grid, including ghosts
        ! Up to someone else to overwrite ghosts
        real(rp), intent(in) :: RIon
            !! Iono radius in Rp
        type(sifGrid_T), intent(in) :: Grid
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isG
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: Q
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg, 2), intent(inout) :: gradQ
            !! I think [unitsQ/rad]

        integer :: i,j
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg) :: sinTh
        
        gradQ = 0.0
        
        associate(sh=>Grid%shGrid)
        ! Cal sin(theta) for everyone to use
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
                    gradQ(i,j,1) = 0.5*(Q(i+1,j) - Q(i-1,j))/RIon/(Grid%delTh(i)+Grid%delTh(i+1))
                else if (isG(i-1,j)) then
                    gradQ(i,j,1) = (Q(i,j)-Q(i-1,j))/RIon/Grid%delTh(i)
                else if (isG(i+1,j)) then
                    gradQ(i,j,1) = (Q(i+1,j)-Q(i,j))/RIon/Grid%delTh(i+1)
                end if

                ! Phi dir
                if(isG(i,j-1) .and. isG(i,j+1)) then
                    gradQ(i,j,2) = 0.5*(Q(i,j+1) - Q(i,j-1))/RIon/sinTh(i)/(Grid%delPh(j)+Grid%delPh(j+1))
                else if (isG(i,j-1)) then
                    gradQ(i,j,2) = (Q(i,j)-Q(i,j-1))/RIon/sinTh(i)/Grid%delPh(j)
                else if (isG(i,j-1)) then
                    gradQ(i,j,2) = (Q(i,j+1)-Q(i,j))/RIon/sinTh(i)/Grid%delPh(j+1)
                end if
            enddo
        enddo

        ! Edge rows and columns
        !! TODO: Think more carefully about gradTh vs gradPhi here, and corner cells
        gradQ(sh%isg,:,:) = gradQ(sh%isg+1,:,:)
        gradQ(sh%ieg,:,:) = gradQ(sh%ieg-1,:,:)
        gradQ(:,sh%jsg,:) = gradQ(:,sh%jsg+1,:)
        gradQ(:,sh%jeg,:) = gradQ(:,sh%jeg-1,:)

        end associate
    end subroutine calcGradIJ


    function activeDt(sh, Grid, State, k) result(dt)
        !! TODO: eventually this should use iVel instead
        !! TODO: Consider dynamic CFL factor
        type(ShellGrid_T), intent(in) :: sh
        type(sifGrid_T), intent(in) :: Grid
        type(sifState_T), intent(in) :: State
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
        ! Only consider points that are in activeShell and not SIFINACTIVE
        where (State%active .ne. SIFINACTIVE .and. asIJ)
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
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T), intent(in) :: Grid
        type(sifState_T), intent(inout) :: State

        integer :: k

        ! Clear things that will be accumulated over the advance
        State%precipNFlux = 0.0
        State%precipEFlux = 0.0

        ! Make sure moments are up to date, some things (like coulomb losses) need them
        call EvalMoments(Grid, State)

        ! Send everyone off
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call AdvanceLambda(Model, Grid, State, k)
        enddo
    end subroutine AdvanceState


    subroutine AdvanceLambda(Model, Grid, State, k)
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T), intent(in) :: Grid
        type(sifState_T), intent(inout) :: State
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
                call calcStepLosses(Model, Grid, State, k, dt)

                t = t + dt
                n = n+1
            enddo

            ! Divide precip fluxes by big dt to turn them into proper rates
            State%precipNFlux(:,:,k) = State%precipNFlux(:,:,k)/State%dt
            State%precipEFlux(:,:,k) = State%precipEFlux(:,:,k)/State%dt

        end associate

    end subroutine AdvanceLambda


end module sifadvancer