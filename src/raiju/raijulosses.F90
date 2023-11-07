module raijulosses
    
    use kdefs

    use raijudefs
    use raijutypes
    use raijuspecieshelper, only : spcIdx

    implicit none

    contains

!------
! High-level calc
!------

    subroutine calcChannelLossRates(Model, Grid, State, k)
        !! Calculate 2D loss rates for channel k
        !! Usually this will stay constant over a coupling period, so it can be called during pre-advance and not touched afterwards
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        
        if (Grid%spc(Grid%k2spc(k))%spcType .eq. RAIJUHPLUS) then
            call calcProtonLossRate(Model, Grid, State, k)
        endif

        if (Grid%spc(Grid%k2spc(k))%spcType .eq. RAIJUELE) then
            call calcElectronLossRate(Model, Grid, State, k)
        endif

    end subroutine calcChannelLossRates


    subroutine calcProtonLossRate(Model, Grid, State, k)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        
        integer :: i,j, psphIdx
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: rateSS, rateCC, rateCX, rateFLC
        

        State%lossRates(:,:,k)       = 0.0  ! 1/s, so 0 means we lose nothing no matter the dt
        State%lossRatesPrecip(:,:,k) = 0.0

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

            psphIdx = spcIdx(Grid, F_PSPH)

            rateSS  = 0.0
            rateCC  = 0.0
            rateCX  = 0.0
            rateFLC = 0.0

            do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                do i=Grid%shGrid%isg,Grid%shGrid%ieg
                    if (.not. isG(i,j)) then
                        cycle
                    endif

                    ! Always do SS for baseline
                    rateSS(i,j) = StrongScatterRate(Rp_m, spc%amu, Grid%alamc(k), State%bVol(i,j), Grid%Bmag(i,j))

                    if (Model%doCC .and. Model%doPlasmasphere) then
                        ! Can only do coulomb collisions if we have a cold plasmasphere
                        ! If we do, use the last evaluated Density-plasmasphere
                        if (psphIdx == -1) then
                            write(*,*) "ERROR: doPlasmasphere=t but can't find psph species"
                            write(*,*)"  goodbye"
                            stop
                        endif
                        rateCC(i,j) = CCRate(spc%spcType, Grid%alamc(k), &
                                        State%bVol(i,j)**(-2./3.), State%Den(i,j,1+psphIdx)) ! Add 1 cause we're grabbing from density, which has bulk as first element
                    endif

                    ! TODO: CX loss rates

                    ! Calc our total loss rate and how much goes into precipitation
                    ! Nothing should be faster than strong scattering
                    State%lossRates      (i,j,k) = min(rateCC(i,j), rateSS(i,j))
                    State%lossRatesPrecip(i,j,k) = min(rateCC(i,j), rateSS(i,j))
                enddo
            enddo

        end associate

    end subroutine calcProtonLossRate


    subroutine calcElectronLossRate(Model, Grid, State, k)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        integer, intent(in) :: k
        
        integer :: i,j
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG
        real(rp), dimension(2) :: lossRate2

        State%lossRates(:,:,k) = 0.0
        State%lossRatesPrecip(:,:,k) = 0.0

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
                lossRate2 = Model%eLossRateFn(Model, Grid, State, i, j, k)  ! [1/s]
                    !! First is actual rate, second is the loss type
                State%lossRates(i,j,k)      = lossRate2(1)
                State%precipType_ele(i,j,k) = lossRate2(2)
                ! All lost electrons assumed to precipitate
                State%lossRatesPrecip(i,j,k) = State%lossRates(i,j,k)
            enddo
        enddo

    end subroutine calcElectronLossRate

!------
! Loss rate calculations
!------

    function StrongScatterRate(Rp_m, amu, alam, bVol, Bfp) result(rateSS)
        !! Calculates strong scattering rate, according to Schulz 1998
        !! tau ~ [2*FTV*Bfp/(1-eta)](gamma*m0/p)
        !! FTV = flux tube volume, Bfp = B-field at foot point, eta - back-scattering rate
        !! eta is backscatter rate at alitude h, here eta=2/3.
        !! gamma = m/m0 is relativisitc factor, p is particle momentum.
        !!       = mc2/m0c2 = (m0c2+K)/m0c2 = 1+K/mec2 ! mec2=0.511 is me*c^2 in MeV
        !! m = m0/sqrt(1-v^2/c^2)
        !! V = c*1/sqrt(1-1/gammar2)
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
        real(rp) :: gammar
            !! Relativistic gamma
        real(rp) :: tauSS, rateSS
            !! Loss timescale [s], rate [1/s]

        
        vm = bVol**(-2./3.)  
        K = abs(alam)*vm*1e-3  ! Energy [keV]
        gammar = 1.0+(K*1e-3)/mec2  ! gamma with 1 + MeV/MeV
        !V = sqrt(2*(K*keV2J)/(amu*dalton)) / Rp_m  ! m/s -> Rp/s
        !tauSS = 2.0*bVol*Bfp/(1.0-eta) / V  ! [s]
        V = (vc_cgs*1e-2)*sqrt(1.0-1.0/gammar**2)/Rp_m ! Rp/s
        tauSS = 2.0*bVol*Bfp/(1.0-eta) / V * gammar ! [s]
        
        if(tauSS > TINY) then
            rateSS = 1.0/tauSS
        else
            rateSS = 0.0
        endif

    end function StrongScatterRate


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
        real(rp) :: deleta, pNFlux, pEFlux
        
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(j,i,deleta)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                ! First update eta using total lossRates over dt
                deleta = State%eta(i,j,k)*(1.0-exp(-dt*State%lossRates(i,j,k)))
                State%eta(i,j,k) = max(0.0, State%eta(i,j,k) - deleta)

                ! Then calculate precipitation flux using lossRatesPrecip
                deleta = State%eta(i,j,k)*(1.0-exp(-dt*State%lossRatesPrecip(i,j,k)))
                State%precipNFlux(i,j,k) = State%precipNFlux(i,j,k) + deleta2NFlux(deleta, Model%planet%rp_m, Grid%Bmag(i,j), dt)
                State%precipEFlux(i,j,k) = State%precipEFlux(i,j,k) + nFlux2EFlux(pNFlux, Grid%alamc(k), State%bVol(i,j))
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