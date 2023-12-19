module raijuAdvancer
    use clocks
    use planethelper
    use kai2geo

    use raijudefs
    use raijutypes
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
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg,2) :: vel, dtArr
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: asIJ, isGood
        real(rp) :: dt

        associate (sh => Grid%shGrid)
        ! make activeShell 2D
        do i=sh%isg,sh%ieg
            asIJ(i,:) = State%activeShells(i,k)
        enddo
        
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
            !dtArr(:,j,1) = 0.5*(Grid%delTh(:sh%ieg)+Grid%delTh(sh%isg+1:)) / (vel(:,j,1) / Model%planet%ri_m)
            dtArr(:,j,1) = 0.5*(Grid%lenFace(:sh%ieg, j, 1)+Grid%lenFace(sh%isg+1:, j, 1)) / (vel(:,j,1) / Model%planet%ri_m)
            !if (k == 20 .and. any(vel(:,j,1) > 400)) then
            !    write(*,*)"dt_theta------"
            !    write(*,*) 0.5*(Grid%lenFace(:sh%ieg, j, 1)+Grid%lenFace(sh%isg+1:, j, 1))
            !    write(*,*)"------"
            !endif
        enddo
        ! In Phi direction
        do i=sh%isg,sh%ieg
            !dtArr(i,:,2) = 0.5*(Grid%delPh(:sh%jeg)+Grid%delPh(sh%jsg+1:)) / (vel(i,:,2) / Model%planet%ri_m)
            dtArr(i,:,2) = 0.5*(Grid%lenFace(i, :sh%jeg, 2)+Grid%lenFace(i, sh%jsg+1:, 2)) / (vel(i,:,2) / Model%planet%ri_m)
            !if (k == 20 .and. any(vel(i,:,2) > 400)) then
            !    write(*,*)"dt_phi------"
            !    write(*,*) 0.5*(Grid%lenFace(i, :sh%jeg, 1)+Grid%lenFace(i, sh%jsg+1:, 1))
            !    write(*,*)"------"
            !endif
        enddo

        dt = Model%CFL*minval(dtArr)

        !if (k == 20) then
        !    write(*,*)"------"
        !    write(*,*) minval(Model%CFL*dtArr), maxval(Model%CFL*dtArr)
        !    write(*,*)"------"
        !endif

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
        !$OMP PARALLEL DO default(shared) collapse(1) &
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


        !Nsteps = int(State%dt / State%dtk(k))+1
        !dt = State%dt / (1.0_rp*Nsteps)
        associate(sh=>Grid%shGrid, spc=>Grid%spc(s))

            ! Here we go!
            n = 1  ! counter
            do while ( tEnd-t > TINY)
                !write(*,*)k,n,t,dt
                
                ! Calc new active shells
                if (Model%doActiveShell) then
                    State%activeShells(:,k) = setLambdaActiveShells(sh, spc, State%bVol, &
                                                State%eta(:,:,spc%kStart:spc%kEnd), k, worthyFracO = Model%worthyFrac)
                endif
                

                ! Calc next time step
                !! Also muting
                dt = activeDt(Model, Grid, State, k)

                !! BAD: Boost dt to be around 200 iters max
                !dt = (tEnd-t)/(Nmax-n)
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

    end subroutine stepLambda


end module raijuAdvancer