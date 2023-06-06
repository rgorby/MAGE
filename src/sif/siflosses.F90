module siflosses
    
    use kdefs

    use sifdefs
    use siftypes

    implicit none

    contains

!------
! Loss entry-point for arbitrary, individual lambda channel
!------

    subroutine calcStepLosses(Model, Grid, State, k, dt)
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T), intent(in) :: Grid
        type(sifState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt
            !! Time delta [s]
        
        if (Grid%spc(Grid%k2spc(k))%spcType .eq. SIFHPLUS) then
            call protonLosses(Model, Grid, State, k, dt)
        endif

        
    end subroutine calcStepLosses


    subroutine protonLosses(Model, Grid, State, k, dt)
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T), intent(in) :: Grid
        type(sifState_T), intent(inout) :: State
        integer, intent(in) :: k
        real(rp), intent(in) :: dt
            !! Time delta [s]
        
        integer :: i,j
        real(rp) :: deleta, pNFlux, pEFlux
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: isG
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg, &
                    Grid%shGrid%jsg:Grid%shGrid%jeg) :: tauSS, tauCC, tauCX, tauFLC
        

        associate(spc=>Grid%spc(Grid%k2spc(k)), Rp_m=>Model%planet%rp_m)
            ! We have nothing to do for plasmasphere
            if (spc%flav .eq. F_PSPH) then
                return
            end if

            ! Otherwise we have some protons to lose

            ! Calc regions where we actually need to evaluate
            where (State%active .eq. SIFACTIVE)
                isG = .true.
            elsewhere
                isG = .false.
            end where

            
            tauSS = HUGE
            tauCC = HUGE
            tauCX = HUGE
            tauFLC = HUGE

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,deleta, pNFlux, pEFlux)
            do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                do i=Grid%shGrid%isg,Grid%shGrid%ieg
                    if (.not. isG(i,j)) then
                        cycle
                    end if


                    if (Model%doSS) then
                        tauSS(i,j) = IonSSTau(Rp_m, spc%amu, Grid%alamc(k), State%bVol(i,j), Grid%Bmag(i,j))
                    end if

                    deleta = State%eta(i,j,k)*(1.0-exp(-dt/tauSS(i,j)))

                    ! Assuming everything in deleta precipitates, calc precip fluxes
                    pNFlux = deleta2NFlux(deleta, Rp_m, Grid%Bmag(i,j), dt)
                    pEFlux = nFlux2EFlux(pNFlux, Grid%alamc(k), State%bVol(i,j))
                    State%precipNFlux(i,j,k) = State%precipNFlux(i,j,k) + pNFlux
                    State%precipEFlux(i,j,k) = State%precipEFlux(i,j,k) + pEFlux

                    ! Finally, update eta
                    State%eta(i,j,k) = max(0.0, State%eta(i,j,k) - deleta)

                    if(i==60 .and. j==60) then
                        write(*,*) k,dt,tauSS(i,j),pNFlux,pEFlux
                    endif
                enddo
            enddo

        end associate

    end subroutine protonLosses


    function IonSSTau(Rp_m, amu, alam, bVol, Bfp) result(tauSS)
        !! Calculates strong scattering rate, according to Schulz 1998
        !! tau ~ [2*FTV*Bfp/(1-eta)](gamma*m/p)
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
        real(rp) :: tauSS
            !! Loss timescale [s]

        
        vm = bVol**(-2./3.)  
        K = abs(alam)*vm*1e-3*kev2J  ! [J]
        V = sqrt(2*K/(amu*dalton)) / Rp_m  ! m/s -> Rp/s

        tauSS = 2.0*bVol*Bfp/(1.0-eta) / V  ! [s]
    end function IonSSTau


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

end module siflosses