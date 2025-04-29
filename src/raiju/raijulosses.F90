module raijulosses
    
    use kdefs
    use XML_Input

    use raijudefs
    use raijutypes
    use raijuspecieshelper, only : spcIdx

    use raijuLoss_CX
    use raijuLoss_CC
    use raijuLoss_SS

    implicit none

    contains


    subroutine initRaiLosses(Model, Grid, State, xmlInp)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        integer :: iLP
            !! iterator
        integer :: numLPs
            !! number of loss proccces we're gonna have
        integer :: initCX, initCC, initSS, initFLC
            !! Flag for if we need to allocate and init this LP

        State%lossRates = 0.0
        State%lossRatesPrecip = 0.0
        State%dEta_dt = 0.0
        State%CCHeatFlux = 0.0

        if (.not. Model%doLosses) then
            return
        endif

        initCX  = merge(1, 0, Model%doCX)
        initCC  = merge(1, 0, Model%doCC)
        initSS  = merge(1, 0, Model%doSS)
        initFLC = merge(1, 0, Model%doFLC)
        
        numLPs = initCX + initCC + initSS + initFLC
        
        allocate(State%lps(numLPs))

        do iLP=1,numLPs
            ! Determine which loss process is next in line for initting
            if (initCX==1) then
                if (allocated(State%lps(iLP)%p)) deallocate(State%lps(iLP)%p)
                allocate( raiLoss_CX_T :: State%lps(iLP)%p )
                initCX = 0
            elseif(initCC==1) then
                if (allocated(State%lps(iLP)%p)) deallocate(State%lps(iLP)%p)
                allocate( raiLoss_CC_T :: State%lps(iLP)%p )
                initCC = 0
                State%lp_cc_idx = iLP
            elseif(initSS==1) then
                if (allocated(State%lps(iLP)%p)) deallocate(State%lps(iLP)%p)
                allocate( raiLoss_SS_T :: State%lps(iLP)%p )
                initSS = 0
            endif

            ! Init newly-allocated LP
            call State%lps(iLP)%p%doInit(Model, Grid, xmlInp)
        enddo

    end subroutine initRaiLosses


    subroutine updateRaiLosses(Model, Grid, State)
        !! Update loss processes with step-specific information
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: nLP, iLP
        integer :: k
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood

        where (State%active .eq. RAIJUACTIVE)
            isGood = .true.
        elsewhere
            isGood = .false.
        end where

        if (allocated(State%lps)) then
            nLP = size(State%lps)
        else
            return
        endif

        do iLP=1,nLP
            call State%lps(iLP)%p%doUpdate(Model, Grid, State)
        enddo

        ! Prep for accumulation this coupling step
        State%dEta_dt = 0.0
        State%precipType_ele = 0.0
        !State%precipNFlux = 0.0
        !State%precipEFlux = 0.0
        State%CCHeatFlux = 0.0
        ! initialize all precip fluxes to zero and masks to false.
        do k=1,Grid%Nk
            State%precipNFlux(k)%data = 0.0
            State%precipEFlux(k)%data = 0.0
            State%precipNFlux(k)%mask = isGood
            State%precipEFlux(k)%mask = isGood
        enddo
    end subroutine updateRaiLosses


    subroutine calcChannelLossRates(Model, Grid, State, k)
        !! Calculate 2D loss rates for channel k
        !! Usually this will stay constant over a coupling period, so it can be called during pre-advance and not touched afterwards
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k

        integer :: i,j,l
        real(rp) :: rate, rateSS, ratePrecip
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG


        ! Electrons are special, handle them on their own
        !if (Grid%spc(Grid%k2spc(k))%spcType .eq. RAIJUELE) then
        !    !call calcElectronLossRate(Model, Grid, State, k)
        !    return
        !endif

        ! Otherwise, do default loss behavior
        State%lossRates(:,:,k)       = 0.0  ! 1/s, so 0 means we lose nothing no matter the dt
        State%lossRatesPrecip(:,:,k) = 0.0

        ! Mask inactive, go ahead and calc losses in buffer just in case anyone wants it
        where (State%active .ne. RAIJUINACTIVE)
            isG = .true.
        elsewhere
            isG = .false.
        end where

        associate(spc=>Grid%spc(Grid%k2spc(k)), lps=>State%lps)

        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (.not. isG(i,j)) then
                    cycle
                endif
                do l=1,size(lps)
                    if ( .not. lps(l)%p%isValidSpc(spc) ) then
                        cycle
                    endif

                    rate = 1.0_rp/lps(l)%p%calcTau(Model, Grid, State, i,j,k)
                    if (rate < 0.0) then
                        write(*,*)"ERROR in raijulosses: lossRate<0 for i,j,k=",i,j,k
                        stop
                    endif
                    if (rate .le. TINY) then
                        rate = 0.0
                    endif
                    
                    State%lossRates(i,j,k) = State%lossRates(i,j,k) + rate
                    if (lps(l)%p%isPrecip) then
                        State%lossRatesPrecip(i,j,k) = State%lossRatesPrecip(i,j,k) + rate
                    endif
                enddo
                
            enddo
        enddo

        end associate

    end subroutine calcChannelLossRates


!------
! High-level calc
!------

    subroutine calcElectronLossRate(Model, Grid, State, k)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k

    end subroutine calcElectronLossRate

!------
! Apply losses and calc useful info
!------

    subroutine applyLosses(Model, Grid, State, k, dt)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt
            !! Time delta [s]

        integer :: i,j
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood
        real(rp) :: deleta, eta0, pNFlux, tau

        where (State%active .eq. RAIJUACTIVE)
            isGood = .true.
        elsewhere
            isGood = .false.
        end where
        
        ! ! !$OMP PARALLEL DO default(shared) collapse(1) &
        ! ! !$OMP schedule(dynamic) &
        ! ! !$OMP private(j,i,eta0, deleta)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (.not. isGood(i,j)) then
                    cycle
                endif
                eta0 = State%eta(i,j,k)
                ! First update eta using total lossRates over dt
                deleta = eta0*(1.0-exp(-dt*State%lossRates(i,j,k)))
                State%eta(i,j,k) = max(0.0, eta0 - deleta)
                ! Note: Not dividing by dt here, we are accumulating over 
                ! coupling step and dividing by step dt at the end
                State%dEta_dt(i,j,k) = State%dEta_dt(i,j,k) + deleta


                ! Then calculate precipitation flux using lossRatesPrecip
                deleta = eta0*(1.0-exp(-dt*State%lossRatesPrecip(i,j,k)))
                pNFlux = deleta2NFlux(deleta, Model%planet%rp_m, Grid%Brcc(i,j), dt)
                ! Just accumulate total #/cm2 and erg/cm2, we divide by coupling dt at the end of advance
                !State%precipNFlux(i,j,k) = State%precipNFlux(i,j,k) + pNFlux*dt
                !State%precipEFlux(i,j,k) = State%precipEFlux(i,j,k) + nFlux2EFlux(pNFlux, Grid%alamc(k), State%bVol_cc(i,j))
                State%precipNFlux(k)%data(i,j) = State%precipNFlux(k)%data(i,j) + pNFlux*dt
                State%precipEFlux(k)%data(i,j) = State%precipEFlux(k)%data(i,j) + nFlux2EFlux(pNFlux, Grid%alamc(k), State%bVol_cc(i,j))

                ! Do special stuff for Coulomb collision effects
                if (Model%doCC .and. State%lps(State%lp_cc_idx)%p%isValidSpc(Grid%spc(Grid%k2spc(k)))) then
                    ! We can estimate heat transfer to plasmasphere electrons by energy lost from RC species to CC
                    ! So we can follow same prodecure as above, by using just CC tau and dividing later by coupling dt to get average heat flux
                    ! Treating this separately from precipication since its not actually precipitating ions
                    tau = max(TINY, State%lps(State%lp_cc_idx)%p%calcTau(Model, Grid, State, i,j,k))
                    deleta = eta0*(1.0 - exp(-dt/tau))
                    pNFlux = deleta2NFlux(deleta, Model%planet%rp_m, Grid%Brcc(i,j), dt)*dt
                    ! Just accumulate total #/cm2 and erg/cm2, we divide by coupling dt at the end of advance
                    !State%CCHeatFlux(i,j,k) = State%CCHeatFlux(i,j,k) + nFlux2EFlux(pNFlux, Grid%alamc(k), State%bvol_cc(i,j))
                    State%CCHeatFlux(i,j,k) = State%CCHeatFlux(i,j,k) + pNFlux*abs(Grid%alamc(k))*State%bVol_cc(i,j)**(-2./3.)  ! [eV/cm2]
                endif
            enddo
        enddo

        

    end subroutine applyLosses


    function deleta2NFlux(eta, Rp_m, Bmag, dt)
        !! Convert precipitating eta to precipitating number flux
        !! Breaking out here so unit conversion is clear
        real(rp), intent(in) :: eta, Rp_m, Bmag, dt
            !! [#/cc * Rp/T], [m], [nT], [s]
        real(rp) :: deleta2NFlux

        deleta2NFlux = (eta/sclEta) * (Rp_m*1.e2) * Bmag / dt
            !! ((#/cm^3 * Rp/T) * T/nT) * cm * nT / s  = #/cm^2/s
        
    end function deleta2NFlux


    function nFlux2EFlux(nFlux, alamc, bVol)
        !! Convert precipitating number flux to precipitating energy flux
        !! Breaking out here so unit conversion is clear
        real(rp), intent(in) :: nflux, alamc, bVol
            !! #/cm^2/s , eV*(Rp/nT)^(2/3) , Rp/nT

        real(rp) :: nFlux2EFlux

        nFlux2EFlux = nFlux*abs(alamc)*bVol**(-2./3.)*1.0e-3*kev2erg  ! [erg/cm^2/s]

    end function nFlux2EFlux




!    subroutine calcElectronLossRate(Model, Grid, State, k)
!        type(raijuModel_T), intent(inout) :: Model
!        type(raijuGrid_T), intent(in) :: Grid
!        type(raijuState_T), intent(inout) :: State
!        integer, intent(in) :: k
!        
!        integer :: i,j
!        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
!                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG
!        real(rp), dimension(2) :: lossRate2
!
!        State%lossRates(:,:,k) = 0.0
!        State%lossRatesPrecip(:,:,k) = 0.0
!
!        ! Calc regions where we actually need to evaluate
!        where (State%active .eq. RAIJUACTIVE)
!            isG = .true.
!        elsewhere
!            isG = .false.
!        end where
!
!        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
!            do i=Grid%shGrid%isg,Grid%shGrid%ieg
!                if (.not. isG(i,j)) then
!                    cycle
!                endif
!
!                ! Calc loss rate
!                lossRate2 = Model%eLossRateFn(Model, Grid, State, i, j, k)  ! [1/s]
!                    !! First is actual rate, second is the loss type
!                State%lossRates(i,j,k)      = lossRate2(1)
!                State%precipType_ele(i,j,k) = lossRate2(2)
!                ! All lost electrons assumed to precipitate
!                State%lossRatesPrecip(i,j,k) = State%lossRates(i,j,k)
!            enddo
!        enddo
!
!    end subroutine calcElectronLossRate
end module raijulosses