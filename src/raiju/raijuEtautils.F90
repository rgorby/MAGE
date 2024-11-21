module raijuetautils

    use planethelper
    use earthhelper

    use raijudefs
    use raijutypes
    use raijuspecieshelper

    implicit none

    real(rp), private, parameter :: kapDefault = 6.0

    !Quad prec. parameters for erf difference
    real(qp), parameter, private :: p  =  0.3275911  , &
                                    a1 =  0.254829592, &
                                    a2 = -0.284496736, &
                                    a3 =  1.421413741, &
                                    a4 = -1.453152027, &
                                    a5 =  1.061405429

    contains

!------
! High-level control
!------
    subroutine EvalMoments(Grid, State)
        !! Calculate D,P, and vAvg for all species across entire grid
        type(raijuGrid_T) , intent(in)    :: Grid
        type(raijuState_T), intent(inout) :: State

        integer :: i,j,s  ! i,j,species iterators

        do s=0,Grid%nSpc
            State%Den  (s)%data = 0.0
            State%Press(s)%data = 0.0
            State%vAvg (s)%data = 0.0
            State%Den  (s)%mask = .true.
            State%Press(s)%mask = .true.
            State%vAvg (s)%mask = .false.
        enddo

        associate (shG => Grid%shGrid, spc => Grid%spc)
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,s)
            do j=shG%jsg,shG%jeg
                do i=shG%isg,shG%ieg
                    do s=1,Grid%nSpc

                        ! isGood determination
                        if (State%active(i,j) .eq. RAIJUINACTIVE) then
                            State%Den  (s)%mask(i,j) = .false.
                            State%Press(s)%mask(i,j) = .false.
                            State%vAvg (s)%mask(i,j) = .false.
                            cycle
                        endif

                        State%Den(s)%data(i,j) = SpcEta2Den(spc(s), &  ! Species details
                            State%eta(i,j,spc(s)%kStart:spc(s)%kEnd), &  ! Etas for this species
                            State%bvol_cc(i,j)) &
                            * spc(s)%amu  ! [#/cc -> amu/cc]

                        State%Press(s)%data(i,j) = SpcEta2Press(spc(s), &  ! Species details
                            State%eta(i,j,spc(s)%kStart:spc(s)%kEnd), &  ! Etas for this species
                            State%bvol_cc(i,j))                        
                    enddo  ! s
                enddo  ! j
            enddo  ! i
            ! Then add each species moment to the bulk
            do s=1,Grid%nSpc
                ! Don't include electrons to total number density
                if(Grid%spc(s)%spcType .ne. RAIJUELE) then
                    State%Den(0)%data = State%Den(0)%data + State%Den(s)%data
                endif
                State%Press(0)%data = State%Press(0)%data + State%Press(s)%data
            enddo
        end associate

    end subroutine EvalMoments



!------
! Species-level control
!------
    function etak2Den(etak, bVol) result (Dk)
        !! Takes a single eta value and converts to density
        real(rp), intent(in) :: etak
        real(rp), intent(in) ::  bVol
            !! Flux tube volume [Rx/nT]

        real(rp) :: Dk
            !! Density [#/cc]

        Dk = (etak/sclEta)/bVol
    end function etak2Den


    function etak2Press(etak, alamc, bVol) result (Pk)
        !! Takes a single eta value and converts to pressure
        real(rp), intent(in) :: etak, alamc
        real(rp), intent(in) ::  bVol
            !! Flux tube volume [Rx/nT]

        real(rp) :: Pk
            !! Pressure [nPa]

        !! Note: 10^-9 from 1/sclEta cancels with 10^9 from Pa -> nPa
        Pk = 2./3.*etak*alamc*bVol**(-5./3.) * ev2J * 1.e6
    end function etak2Press


    function SpcEta2Den(spc, eta, bVol) result(D)
        !! Take a species' eta at a specific point, and sum moments to get its density and pressure
        type(raijuSpecies_T), intent(in) :: spc
            !! Species info
        real(rp), dimension(spc%kStart:spc%kEnd), intent(in) :: eta
            !! Etas we are summing
        real(rp), intent(in) ::  bVol
            !! Flux tube volume [Rx/nT]
        
        integer :: k
        real(rp) :: D
            !! Density [#/cc]

        D = 0.0

        if (bVol <= 0) return

        do k=spc%kStart,spc%kEnd
            D = D + etak2Den(eta(k), bVol)
        enddo

    end function SpcEta2Den


    function SpcEta2Press(spc, eta, bVol) result(P)
        !! Take a species' eta at a specific point, and sum moments to get its density and pressure
        type(raijuSpecies_T), intent(in) :: spc
            !! Species info
        real(rp), dimension(spc%kStart:spc%kEnd), intent(in) :: eta
            !! Etas we are summing
        real(rp), intent(in) ::  bVol
            !! Flux tube volume [Rx/nT]
        
        integer :: k
        real(rp) :: alamc
            !! Cell-center lambda value
        real(rp) :: P
            !! Pressure [nPa]

        P = 0.0

        if (bVol <= 0) return

        do k=spc%kStart,spc%kEnd
            alamc = 0.5*abs(spc%alami(k) + spc%alami(k+1))
            !P = P + eta(k)*alamc*vm**2.5 * ev2J * 1.e6
            !! Note: 10^-9 from 1/sclEta cancels with 10^9 from Pa -> nPa
            P = P + etak2Press(eta(k), alamc, bVol)
        enddo

    end function SpcEta2Press


    function spcEta2DPS(Model, Grid, State, spc, isGood) result(dpsdst)
        !! Calculate total DPS-Dst for given species within the defined isGood domain
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        type(raijuSpecies_T), intent(in) :: spc
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg, Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood
            !! Eval mask, true = point is included in calculation

        real(rp) :: dpsdst
        integer :: i,j,k
        real(rp) :: press, bVol, energyDen, energy

        dpsdst = 0.0

        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (.not. isGood(i,j)) cycle
                bVol = State%bvol_cc(i,j)
                press = SpcEta2Press(spc, State%eta(i,j,spc%kStart:spc%kEnd), bVol)  ! [nPa]
                energyDen = (press*1.0D-9) * (bVol*Model%planet%ri_m*1.0D9) * (Grid%Brcc(i,j)*1.0D-9)/kev2J  ! p[J/m^3] * bVol[m/T] * B[T]  = [J/m^2] * keV/J = [keV/m^2]
                energy = energyDen*(Grid%areaCC(i,j)*Model%planet%ri_m**2) !  [keV/m^2]* Re^2[m^2] = [keV]
                dpsdst = dpsdst - 4.2*(1.0D-30)*energy  ! [nT]
            enddo
        enddo
    end function spcEta2DPS


    subroutine DkT2SpcEta(Model, spc, eta, D, kT, vm, doAccumulateO, etaBelowO)
        !! Take a density and pressure, and map it to RAIJU eta channels for given species
        type(raijuModel_T), intent(in) :: Model
        type(raijuSpecies_T), intent(in) :: spc
            !! Species info
        real(rp), dimension(spc%kStart:spc%kEnd), intent(inout) :: eta
            !! len(spc%N) etas we need to populate
        real(rp), intent(in) :: D, kT, vm
            !! Density [#/cc], Energy [keV], bVol^-2/3 [(Rx/nT)^(-2/3)]
        real(rp), optional, intent(out) :: etaBelowO
            !! If provided, we will return the eta that's below the lowest channel of this species
        logical, optional, intent(in) :: doAccumulateO
            !! Whether or not to zero out etas before putting in new stuff

        
        integer :: k
        logical :: doAccumulate

        if (present(doAccumulateO)) then
            doAccumulate = doAccumulateO
        else
            doAccumulate = .false.
        endif

        if (.not. doAccumulate) then
            eta = 0.0
        endif

        ! Trap for zero density
        if (D < TINY) then
            return
        endif
        
        if (size(eta) .eq. 1) then
            ! Just throw all the density into the one channel
            eta(spc%kStart) = eta(spc%kStart) + D/(vm**1.5)*sclEta ! #/cc * Rx/nT * nT/T -> #/cc*Rx/T
            return
        endif

        ! Trap for zero energy. If doing zero-energy plasmasphere, would be handled above
        if (kT < TINY .and. abs(spc%alami(spc%kStart+1)) > TINY) then
                            ! If upper bound of lowest cell isn't above TINY, there's nowhere for zero-energy stuff to go
            return
        endif


        do k=spc%kStart,spc%kEnd
            eta(k) = eta(k) + Model%dp2etaMap(Model,D,kT,vm,abs(spc%alami(k)),abs(spc%alami(k+1)))
                !! Model%dp2etaMap may point to one of the functions below, like Maxwell2Eta
                !! or it may point to a user-defined mapping function
        enddo

        if (present(etaBelowO)) then
            etaBelowO = Model%dp2etaMap(Model,D,kT,vm,0.0_rp,abs(spc%alami(spc%kStart)))
        endif

        
    end subroutine DkT2SpcEta

    
    function Kappa2Eta(Model,D,kT,vm,amin,amax) result(etaK)
        !! Convert density and temperature to eta at specific lambda value
        !! Adapted from eqn 3.12 from 10.1007/s11214-013-9982-9 ?
        type(raijuModel_T), intent(in) :: Model
        real(rp), intent(in) :: D,kT,vm,amin,amax
            !! Density [#/cc], kT [keV], vm [(Rx/nT)^(-2/3)],
            !! min and max lambda vals [eV * (Rx/nT)^(2/3)]
            !! Optional Kappa value

        real(rp) :: kap,kap15,E_ev,E0_ev
        real(rp) :: A0,kapgam,kapbar,kArg,delscl
        real(rp) :: etaK
            !! Eta in units of [#/cc * Rx/T]
        etaK = 0.0

        kap = Model%kappa

        kap15 = kap-1.5
        kapgam = gamma(kap+1.0)/gamma(kap-0.5)
        kapbar = kap15 !Should be kap-3/2 or kap
            ! TODO: check 10.1002/2015JA021166

        E0_ev = kT*kap15/kap*1.e3  ! [eV]
        E_ev = abs(amin+amax)/2.0*vm  ! Lambda center energy [eV]

        kArg = 1.0 + (E_ev/E0_ev)/kapbar
        delscl = vm*abs(amax-amin)/E0_ev  ! Channel width / temp
            ! TODO: why?

        A0 = (2.0/sqrt(PI)) * D/(vm**1.5)*sclEta  ! #/cc * Rx/nT * (1/nT -> 1/T)
        etak = A0*kapgam/(kapbar**1.5) * sqrt(E_ev/E0_ev)*delscl*((kArg)**(-kap-1.0))

    end function Kappa2Eta


    function Maxwell2Eta(Model,D,kT,vm,amin,amax) result(etaK)
        !! This calculates Maxwellian DP2Eta mapping using experfdiff
        !! Original equation comes from B5 from Pembroke+ 2012
        !! But has been adapted by Kareem to calculate the erfdiff as the diff between two expansions (high vs low cell interface)
        !! This way, the 1 term cancels and you don't have a huge value minus another huge value
        !! From Kareem:
        !!   Difference of erf's using Abramowitz & Stegun, 7.1.26
        !!   erfdiff(qp,qm) = erf(qp)-erf(qm)
        !!   erf(x) ~ 1 - (a1.t + a2.t^2 + a3.t^3 + a4.t^5 + a5.t^5)*exp(-x^2) + eps(x)
        !!   t = 1/(1+px)
        !!   |eps(x)| <= 1.5e-7
        type(raijuModel_T), intent(in) :: Model
        real(rp), intent(in) :: D,kT,vm,amin,amax
            !! Density [#/cc], kT [keV], vm [(Rx/nT)^(-2/3)],
            !! min and max lambda vals [eV * (Rx/nT)^(2/3)]
        real(rp) :: A0
            !! Flux tube content
        real(qp) :: xp,xm,tp,tm,ep,em,erfdiff,expdiff,etaq
            !! (Quad precision) Variables for experfdiff calculation
        real(rp) :: etaK
            !! Eta in units of [#/cc * Rx/T]
        etaK = 0.0

        A0 = D/(vm**1.5)*sclEta ! #/cc * Rx/nT * nT/T -> #/cc*Rx/T
        xp = sqrt(abs(amax)*vm / (kT*1.e3))  ! [eV/eV]
        xm = sqrt(abs(amin)*vm / (kT*1.e3))

        tp = 1.0/(1.0+p*xp)
        tm = 1.0/(1.0+p*xm)
        ep = exp(-(xp**2.0))
        em = exp(-(xm**2.0))

        erfdiff = - (a1*tp + a2*(tp**2.0) + a3*(tp**3.0) + a4*(tp**4.0) + a5*(tp**5.0))*ep &
                + (a1*tm + a2*(tm**2.0) + a3*(tm**3.0) + a4*(tm**4.0) + a5*(tm**5.0))*em

        expdiff = 2.0*(xp*ep - xm*em)/sqrt(QPI) !Using quad prec PI
        etaq = A0*(erfdiff-expdiff)
        etaK = etaq !Cast back down to rp

    end function Maxwell2Eta


! General helpers
    function eta2intensity(spc, bVol, eta) result (intensity)
        type(raijuSpecies_T), intent(in) :: spc
        real(rp), intent(in) :: bVol
            !! bVol = flux tube volume [Rx/nT]
        real(rp), dimension(spc%kStart:spc%kEnd), intent(in) :: eta
            !! eta   = etas for a single species

        real(rp), dimension(:), allocatable :: intensity
        integer :: N
        real(rp), dimension(:), allocatable :: alamc, lamdiff   

        N = size(eta)
        allocate(alamc(N))
        allocate(lamdiff(N))
        allocate(intensity(N))

        alamc   = abs(spc%alami(spc%kStart+1:)+spc%alami(spc%kStart:spc%kEnd))/2.0
        lamdiff = abs(spc%alami(spc%kStart+1:)-spc%alami(spc%kStart:spc%kEnd))

        intensity = sclIntens/sqrt(spc%amu)  &
                  * sqrt(alamc)/lamdiff  &
                  * bVol**(-2./3.) * eta

    end function eta2intensity



!------
! Plasmasphere initialization
!------

    subroutine setRaijuInitPsphere(Model, Grid, State, Kp)
        type(raijuModel_T) , intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        real(rp) :: Kp

        integer :: psphIdx

        if (Model%doPlasmasphere .and. spcExists(Grid, F_PSPH)) then
            psphIdx = spcIdx(Grid, F_PSPH)
            State%eta     (:,:,Grid%spc(psphIdx)%kStart) = getInitPsphere(Grid, State, Kp)
            State%eta_last(:,:,Grid%spc(psphIdx)%kStart) = State%eta(:,:,Grid%spc(psphIdx)%kStart)
        endif
    
    end subroutine


    function getInitPsphere(Grid, State, Kp) result(etaPsph)
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        real(rp) :: Kp

        integer :: i,j
        real(rp) :: den, vm
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: etaPsph

        write(*,*) "RAIJU initializing plasmasphere with Kp =",Kp

        etaPsph = 0.0

        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,den,vm)
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                den = GallagherXY(State%xyzMincc(i,j,XDIR), State%xyzMincc(i,j,YDIR), Kp)  ! [#/cc]
                !vm = State%bvol_cc(i,j)**(-2./3.)
                etaPsph(i,j) = den*State%bvol_cc(i,j)*sclEta
            enddo
        enddo

    end function getInitPsphere

end module raijuetautils
