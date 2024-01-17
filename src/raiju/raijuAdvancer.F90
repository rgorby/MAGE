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
        !! Calculates min dt needed to safele evolve active domain for given energy invariant channel k
        !! TODO: Consider dynamic CFL factor
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        integer, intent(in) :: k
        
        integer :: i,j
        real(rp) :: velij
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: dtArr
        logical , dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 2) :: isGood
        real(rp) :: dt

        associate (sh => Grid%shGrid)

        dtArr = HUGE

        ! Determine which faces can be considered good
        isGood = .false.
        do j=sh%js,sh%je+1
            do i=sh%is,sh%ie+1
                ! NOTE: We are only checking faces bordering non-ghost cells because those are the only ones we use for evolution
                ! TODO: Strictly speaking, there are some faces that are included here that shouldn't be because they're never used to evolve anything
                !  so we should make the loop js:je, is:ie, and handle the last row and column afterwards

                ! We are responsible for face (i,j,TH) and (i,j,PH)
                ! Theta faces first

                ! This is a weird pattern. Basically, we can't cycle if the overall desired condition isn't met, because we are doing theta and phi directions in the same loop
                ! At the same time, we don't want to write one massive if condition with all options covered, or a bunch of nested if statements
                ! So instead, we break them up, such that if a single if condition is true it means we have failed the physical condition for us to calculate a valid timestep
                ! If all if conditions fail then the physical conditions are met and we can do our calculations in the else block
                if (Model%doActiveShell .and. (State%activeShells(i-1,k) .eq. .false. .or. State%activeShells(i,k) .eq. .false.) ) then
                    ! In order for a theta-dir face to be usable, we need both sides to be active
                    continue
                else if ( State%active(i-1,j) .eq. .false. .or. State%active(i,j) .eq. .false.) then
                    continue
                else
                    ! We are good lets calculate a dt for this face
                    velij = max(abs(State%iVel(i,j,k,RAI_TH)), TINY)  ! [m/s]
                    dtArr(i,j,RAI_TH) = ( Grid%delTh(i) * Model%planet%ri_m ) / velij  ! [s]
                endif
                
                ! Phi faces
                if (Model%doActiveShell .and. State%activeShells(i,k) .eq. .false. ) then
                    ! In order for a phi-dir face to be usable, we just need this i shell to be active
                    continue
                else if (State%active(i,j-1) .eq. .false. .or. State%active(i,j) .eq. .false. ) then 
                    continue
                else
                    velij = max(abs(State%iVel(i,j,k,RAI_PH)), TINY)
                    dtArr(i,j,RAI_PH) = ( Grid%delPh(j) * Model%planet%ri_m ) / velij  ! [s]
                endif

            enddo
        enddo

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
                ! This way, all channels will end exactly at tEnd
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
        ! Save eta to eta_last
        State%eta_last(:,:,k) = State%eta(:,:,k)

        ! Calculate eta face fluxes
        call calcFluxes(Model, Grid, State, k, State%eta_half(:,:,k), Qflux)  ! [eta * Rp/s]

        
        ! Calc new eta
        do j=Grid%shGrid%js,Grid%shGrid%je
            do i=Grid%shGrid%is,Grid%shGrid%ie
                State%eta(i,j,k) = State%eta(i,j,k) + dt/Grid%areaCC(i,j)/Grid%BrCC(i,j) &
                                                    * ( Qflux(i  ,j  ,RAI_TH)*Grid%lenFace(i  ,j  ,RAI_TH) &
                                                      - Qflux(i+1,j  ,RAI_TH)*Grid%lenFace(i+1,j  ,RAI_TH) &
                                                      + Qflux(i  ,j  ,RAI_PH)*Grid%lenFace(i  ,j  ,RAI_PH) &
                                                      - Qflux(i  ,j+1,RAI_PH)*Grid%lenFace(i  ,j+1,RAI_PH) )
            enddo
        enddo

        call wrapJcc(Grid%shGrid, State%eta(:,:,k))

    end subroutine stepLambda


end module raijuAdvancer
