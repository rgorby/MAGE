module raijulosses
    
    use kdefs

    use raijudefs
    use raijutypes
    use raijuspecieshelper, only : spcIdx

    implicit none

    contains

!------
! Loss entry-point for arbitrary, individual lambda channel
!------

    subroutine calcStepLosses(Model, Grid, State, k, dt)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt
            !! Time delta [s]
        
        if (Grid%spc(Grid%k2spc(k))%spcType .eq. RAIJUHPLUS) then
            call protonLosses(Model, Grid, State, k, dt)
        endif

        if (Grid%spc(Grid%k2spc(k))%spcType .eq. RAIJUELE) then
            call electronLosses(Model, Grid, State, k, dt)
        endif

    end subroutine calcStepLosses


    subroutine protonLosses(Model, Grid, State, k, dt)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt
            !! Time delta [s]
        
        integer :: i,j, psphIdx
        real(rp) :: deleta, pNFlux, pEFlux, lossRate
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: rateSS, rateCC, rateCX, rateFLC
        

        associate(spc=>Grid%spc(Grid%k2spc(k)), Rp_m=>Model%planet%rp_m)
            ! We have nothing to do for plasmasphere
            if (spc%flav .eq. F_PSPH) then
                return
            end if

            ! Otherwise we have some protons to lose

            ! Calc regions where we actually need to evaluate
            where (State%active .eq. RAIJUACTIVE)
                isG = .true.
            elsewhere
                isG = .false.
            end where

            
            rateSS  = HUGE
            rateCC  = HUGE
            rateCX  = HUGE
            rateFLC = HUGE

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,deleta, pNFlux, pEFlux, lossRate, psphIdx) &
            !$OMP IF(.false.)
            do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                do i=Grid%shGrid%isg,Grid%shGrid%ieg
                    if (.not. isG(i,j)) then
                        cycle
                    endif

                    ! Precipitation first

                    ! Always do SS for baseline
                    rateSS(i,j) = IonSSRate(Rp_m, spc%amu, Grid%alamc(k), State%bVol(i,j), Grid%Bmag(i,j))

                    if (Model%doCC .and. Model%doPlasmasphere) then
                        ! Can only do coulomb collisions if we have a cold plasmasphere
                        ! If we do, use the last evaluated Density-plasmasphere
                        psphIdx = spcIdx(Grid, F_PSPH)  ! Add 1 cause we're grabbing from density, which has bulk as first element
                        if (psphIdx == -1) then
                            write(*,*) "ERROR: doPlasmasphere=t but can't find psph species"
                            write(*,*)"  goodbye"
                            stop
                        endif
                        rateCC(i,j) = CCRate(spc%spcType, Grid%alamc(k), &
                                        State%bVol(i,j)**(-2./3.), State%Den(i,j,1+psphIdx))
                    endif

                    ! Calc our precip loss rate
                    ! Nothing should be faster than strong scattering
                    lossRate = min(rateCC(i,j), rateSS(i,j))


                    deleta = State%eta(i,j,k)*(1.0-exp(-dt*lossRate))
                        !! Total eta lost over dt
                    State%eta(i,j,k) = max(0.0, State%eta(i,j,k) - deleta)

                    ! Assuming everything in deleta precipitates, calc precip fluxes
                    pNFlux = deleta2NFlux(deleta, Rp_m, Grid%Bmag(i,j), dt)
                    pEFlux = nFlux2EFlux(pNFlux, Grid%alamc(k), State%bVol(i,j))
                    State%precipNFlux(i,j,k) = State%precipNFlux(i,j,k) + pNFlux
                    State%precipEFlux(i,j,k) = State%precipEFlux(i,j,k) + pEFlux

                    ! TODO: CX loss rates
                enddo
            enddo

        end associate

    end subroutine protonLosses


    subroutine electronLosses(Model, Grid, State, k, dt)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt
            !! Time delta [s]
        
        integer :: i,j, psphIdx
        real(rp) :: deleta, pNFlux, pEFlux, lossRate
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG
        !real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
        !            Grid%shGrid%jsg:Grid%shGrid%jeg) :: rateLoss

        ! Calc regions where we actually need to evaluate
        where (State%active .eq. RAIJUACTIVE)
            isG = .true.
        elsewhere
            isG = .false.
        end where

        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (.not. isG(i,j)) then
                    cycle
                endif

                ! Calc loss rate
                lossRate = Model%calcELossRate_WM(Model, Grid, state, k)  ! [1/s]

                deleta = State%eta(i,j,k)*(1.0-exp(-dt*lossRate))
                    !! Total eta lost over dt
                State%eta(i,j,k) = max(0.0, State%eta(i,j,k) - deleta)

                ! Assuming everything in deleta precipitates, calc precip fluxes
                pNFlux = deleta2NFlux(deleta, Model%planet%rp_m, Grid%Bmag(i,j), dt)
                pEFlux = nFlux2EFlux(pNFlux, Grid%alamc(k), State%bVol(i,j))
                State%precipNFlux(i,j,k) = State%precipNFlux(i,j,k) + pNFlux
                State%precipEFlux(i,j,k) = State%precipEFlux(i,j,k) + pEFlux

            enddo
        enddo

    end subroutine electronLosses

!------
! Loss rates
!------

    function IonSSRate(Rp_m, amu, alam, bVol, Bfp) result(rateSS)
        !! Calculates strong scattering rate, according to Schulz 1998
        !! tau ~ [2*FTV*Bfp/(1-eta)](gamma*m0/p)
        !! FTV = flux tube volume, Bfp = B-field at foot point, eta - back-scattering rate
        !! Note: Assuming we don't have any relativistic protons, so implemented equation is:
        !! tau ~ [2*FTV*Bfp/(1-eta)]/V
        !! (This function should work for electrons too, if relativistic factor is included)
        real(rp), intent(in) :: Rp_m
            !! Planet radius [m]
        real(rp), intent(in) :: amu
            !! Ion mass in amu
        real(rp), intent(in) :: alam
            !! Lambda in eV*(Rp/nt)^(2/3)
        real(rp), intent(in) :: bVol
            !! FTV [(Rp/nT)]
        real(rp), intent(in) :: Bfp
            !! B-field strength at foot point / ionosphere [nT]

        real(rp) :: eta = 2./3.
            !! Back-scatter rate
        real(rp) :: vm, K, V
            !! vm = FTV^(-2/3) , K = KE [J] , V = velocity [Rp/s]
        real(rp) :: tauSS, rateSS
            !! Loss timescale [s], rate [1/s]

        
        vm = bVol**(-2./3.)  
        K = abs(alam)*vm*1e-3*kev2J  ! [J]
        V = sqrt(2*K/(amu*dalton)) / Rp_m  ! m/s -> Rp/s

        tauSS = 2.0*bVol*Bfp/(1.0-eta) / V  ! [s]
        if(tauSS > TINY) then
            rateSS = 1.0/tauSS
        else
            rateSS = 0.0
        endif

    end function IonSSRate


    ! Simple Coulomb collision losses, using fit to Ebihara+ 98 Fig #5
    function CCRate(spcType,alam,vm,Dpp) result(rateCC)
        integer, intent(in) :: spcType
        real(rp), intent(in) :: alam,vm,Dpp !Dpp is plasmasphere density in #/cc
        real(rp) :: rateCC

        real(rp), parameter :: a3 = -0.113288
        real(rp), parameter :: a2 = +0.659057
        real(rp), parameter :: a1 = +0.319542
        real(rp), parameter :: a0 = +2.16253        
        real(rp), parameter :: day2s = 24.0*60.0*60

        real(rp) :: K,x,y,nTau,Tau

        K = abs( alam*vm*(1.0e-3) )!Energy [keV]
        x = log10(K)
        
        rateCC = 0.0
        if (Dpp < TINY) return

        if (spcType == RAIJUHPLUS) then
            y = a3*(x**3.0) + a2*(x**2.0) + a1*x + a0
            nTau = 10.0**y !Normalized lifetime, days/cm3
            Tau = nTau*day2s/Dpp !Lifetime, [s]
            
            if (Tau>TINY) then
                rateCC = 1.0/Tau !1/s
            else
                rateCC = 0.0
            endif
        else
            rateCC = 0.0
            return
        endif

    end function CCRate



!------
! Conversions
!------

    function deleta2NFlux(eta, Rp_m, Bmag, dt)
        !! Convert precipitating eta to precipitating number flux
        !! Breaking out here so unit conversion is clear
        real(rp), intent(in) :: eta, Rp_m, Bmag, dt
            !! [#/cc * Rp/T], [m], [nT], [s]
        real(rp) :: deleta2NFlux

        deleta2NFlux = (eta/sclEta) * (Rp_m*1.e2) * Bmag / dt
            !! (#/cm^3 * Rp/T * T/nT) * cm * nT / s  = #/cm^2/s
        
    end function deleta2NFlux


    function nFlux2EFlux(nFlux, alamc, bVol)
        !! Convert precipitating number flux to precipitating energy flux
        !! Breaking out here so unit conversion is clear
        real(rp), intent(in) :: nflux, alamc, bVol
            !! #/cm^2/s , eV*(Rp/nT)^(2/3) , Rp/nT

        real(rp) :: nFlux2EFlux

        nFlux2EFlux = nFlux*abs(alamc)*bVol**(-2./3.)*1.0e-3*kev2erg  ! [erg/cm^2/s]

    end function nFlux2EFlux

end module raijulosses