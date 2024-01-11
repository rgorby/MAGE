module raijuAdvancer
    use clocks
    use planethelper
    use kai2geo

    use raijudefs
    use raijutypes
    use raijugrids
    use raijuetautils
    use raijuBCs
    use raijulosses
    use raijuRecon

    implicit none

    contains

    function activeDt(Model, Grid, State, k) result(dt)
        !! TODO: Consider dynamic CFL factor
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        
        integer :: i,j
        !real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1,2) :: vel, dtArr
        !logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: asIJ, isGood
        real(rp) :: velij
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: dtArr
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: isGood
        real(rp) :: dt

        associate (sh => Grid%shGrid)
        ! make activeShell 2D
        !do i=sh%isg,sh%ieg
        !    asIJ(i,:) = State%activeShells(i,k)
        !enddo
        !
        !! Only consider points that are in activeShell and not RAIJUINACTIVE
        !where (State%topo .ne. RAIJUCLOSED .and. asIJ)
        !    isGood = .true.
        !elsewhere
        !    isGood = .false.
        !end where

        dtArr = HUGE

        ! Determine which faces can be considered good
        isGood = .false.
        do j=sh%js,sh%je+1
            do i=sh%is,sh%ie+1
                ! NOTE: We are only checking faces bordering non-ghost cells because those are the only ones we use for evolution

                ! We are responsible for face (i,j,TH) and (i,j,PH)
                ! Theta faces first

                ! This is a weird pattern. Basically, we can't cycle if the overall desired condition isn't met, because we are doing theta and phi directions in the same loop
                ! At the time, we don't want to write one massive if condition with all options covered, or a bunch of nested if statements
                ! So instead, we break them up, such that if a single if condition is true it means we have failed the physical condition for us to calculate a valid timestep for
                ! If all if conditions fail then the physical conditions are met and we can do our calculations in the else block
                if (Model%doActiveShell .and. (State%activeShells(i-1,k) .eq. .false. .or. State%activeShells(i,k) .eq. .false.) ) then
                    continue
                else if ( State%active(i-1,j) .eq. .false. .or. State%active(i,j) .eq. .false.) then
                    continue
                else
                    ! We are good lets calculate a dt for this face
                    velij = max(abs(State%iVel(i,j,k,RAI_TH)), TINY)
                    dtArr(i,j,RAI_TH) = Grid%delTh(i) / ( velij / Model%planet%ri_m)  ! [s]
                endif
                
                ! Phi faces
                if (Model%doActiveShell .and. State%activeShells(i,k) .eq. .false. ) then
                    continue
                else if (State%active(i,j-1) .eq. .false. .or. State%active(i,j) .eq. .false. ) then 
                    continue
                else
                    velij = max(abs(State%iVel(i,j,k,RAI_PH)), TINY)
                    dtArr(i,j,RAI_PH) = Grid%delPh(j) / ( velij / Model%planet%ri_m)  ! [s]
                endif

            enddo
        enddo

        !vel(:,:,RAI_TH) = merge(abs(State%iVel(:,:,k,RAI_TH)), TINY, isGood(:,:,RAI_TH))
        !vel(:,:,RAI_PH) = merge(abs(State%iVel(:,:,k,RAI_PH)), TINY, isGood(:,:,RAI_PH))
!
        !! dt in theta direction
        !do j=sh%js,sh%je+1
        !    dtArr(:,j,RAI_TH) = Grid%delTh(:) / (vel(:,j,RAI_TH) / Model%planet%ri_m)
        !    !dtArr(:,j,RAI_TH) = 0.5*(Grid%delTh(:sh%ie)+Grid%delTh(sh%is+1:)) / (vel(:,j,RAI_TH) / Model%planet%ri_m)
        !    !dtArr(:,j,RAI_TH) = 0.5*(Grid%lenFace(:sh%ieg, j, RAI_TH)+Grid%lenFace(sh%isg+1:, j, RAI_TH)) / (vel(:,j,RAI_TH) / Model%planet%ri_m)
        !enddo
        !! In phi direction
        !do i=sh%is,sh%ie+1
        !    dtArr(:,j,RAI_PH) = Grid%delPh(:) / (vel(i,:,RAI_PH) / Model%planet%ri_m)
        !    !dtArr(i,:,RAI_PH) = 0.5*(Grid%delPh(:sh%je)+Grid%delPh(sh%js+1:)) / (vel(i,:,RAI_PH) / Model%planet%ri_m)
        !    !dtArr(i,:,RAI_PH) = 0.5*(Grid%lenFace(i, :sh%jeg, RAI_PH)+Grid%lenFace(i, sh%jsg+1:, RAI_PH)) / (vel(i,:,2) / Model%planet%ri_m)
        !enddo

        dt = Model%CFL*minval(dtArr)

        end associate

            
    end function activeDt


!------
! Advance entry point
!------

    subroutine AdvanceState(Model, Grid, State)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: k

        ! Send each channel off on their own
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(k)
        do k=1,Grid%Nk
            call AdvanceLambda(Model, Grid, State, k)
        enddo

        ! Wrap up, save new time and cpl step number
        ! Save current active domain to active_last
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


        associate(sh=>Grid%shGrid, spc=>Grid%spc(s))

            ! Here we go!
            n = 0  ! counter
            do while ( tEnd-t > TINY)
                !write(*,*)k,n,t,dt
                
                ! Calc new active shells
                !if (Model%doActiveShell) then
                !    State%activeShells(:,k) = setLambdaActiveShells(sh, spc, State%bVol, &
                !                                State%eta(:,:,spc%kStart:spc%kEnd), k, worthyFracO = Model%worthyFrac)
                !endif
                

                ! Calc next time step
                dt = activeDt(Model, Grid, State, k)

                ! If needed, reduce dt to fit within remaining time
                if (t + dt > tEnd) then
                    dt = max(tEnd - t, TINY)  ! Make sure we never go back in time and advance at least a little bit so we will eventually end
                endif

                ! Advection
                call stepLambda(Model, Grid, State, k, dt)
                
                ! Losses
                if (Model%doLosses) then
                    call applyLosses(Model, Grid, State, k, dt)
                endif

                t = t + dt
                n = n+1

                if (n/State%dt > Model%maxItersPerSec) then
                    write(*,*)"ERROR: Too many advance steps. Dying"
                    write(*,*)" k, dt, nSteps= ", k, dt, n
                    stop
                endif
            enddo

            if (Model%doLosses) then
                ! Divide precip fluxes by big dt to turn them into proper rates
                State%precipNFlux(:,:,k) = State%precipNFlux(:,:,k)/State%dt
                State%precipEFlux(:,:,k) = State%precipEFlux(:,:,k)/State%dt
            endif

            State%nStepk(k) = State%nStepk(k) + n

        end associate

    end subroutine AdvanceLambda


    subroutine stepLambda(Model, Grid, State, k, dt)
        !! Advances a single lambda channel a single time step
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt

        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1, &
                            Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: Qflux
            !! TODO: May not want to keep this at this scope

        integer :: i,j

        ! Calculate eta at half-step
        ! For now, just implement simple AB sub-stepping
        State%eta_half(:,:,k) = 1.5_rp*State%eta(:,:,k) - 0.5_rp*State%eta_last(:,:,k)

        ! Calculate eta face fluxes
        call calcFluxes(Model, Grid, State, k, State%eta_half(:,:,k), Qflux)  ! [eta * Rp^2/s]

        ! Save eta to eta_last
        State%eta_last(:,:,k) = State%eta(:,:,k)
        
        ! Calc new eta
        do j=Grid%shGrid%js,Grid%shGrid%je
            do i=Grid%shGrid%is,Grid%shGrid%ie
                State%eta(i,j,k) = State%eta(i,j,k) + dt/Grid%areaCC(i,j) &
                                                      * ( Qflux(i,j,RAI_TH) - Qflux(i+1,j,RAI_TH) + Qflux(i,j,RAI_PH) - Qflux(i,j+1,RAI_PH) )
            enddo
        enddo

        call wrapJcc(Grid%shGrid, State%eta(:,:,k))

    end subroutine stepLambda


end module raijuAdvancer
