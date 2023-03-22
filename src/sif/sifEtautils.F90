module sifetautils

    use planethelper

    use sifdefs
    use siftypes
    use sifspecieshelper

    implicit none

    real(rp), private, parameter :: kapDefault = 6.0

    contains

! High-level control
    subroutine EvalMoments(Grid, State)
        !! Calculate D,P, and vAvg for all species across entire grid
        type(sifGrid_T) , intent(in)    :: Grid
        type(sifState_T), intent(inout) :: State

        integer :: i,j,s  ! i,j,species iterators

        State%Den = 0.0
        State%Press = 0.0
        State%vAvg = 0.0

        associate (shG => Grid%shGrid, spc => Grid%spc)

            do i=shG%is,shG%ie
                do j=shG%js,shG%je
                    do s=1,Grid%nSpc
                        ! TODO: handle isGood regions?

                        ! Calc moments for this species
                        State%Den(i,j,s+1) = SpcEta2Den(spc(s), &  ! Species details
                            State%eta(i,j,spc(s)%kStart:spc(s)%kEnd), &  ! Etas for this species
                            State%bvol(i,j)) &
                            * spc(s)%amu  ! [#/cc -> amu/cc]

                        State%Press(i,j,s+1) = SpcEta2Press(spc(s), &  ! Species details
                            State%eta(i,j,spc(s)%kStart:spc(s)%kEnd), &  ! Etas for this species
                            State%bvol(i,j))

                        ! TODO: vBlk
                        
                    enddo  ! s
                enddo  ! j
            enddo  ! i

            ! Then add each species moment to the bulk
            do s=2,Grid%nSpc+1
                State%Den  (:,:,1) = State%Den  (:,:,1) + State%Den  (:,:,s)
                State%Press(:,:,1) = State%Press(:,:,1) + State%Press(:,:,s)
                ! State%vAvg (:,:,1) = State%vAvg (:,:,1) + State%Den(:,:,s)*State%vAvg (:,:,s)
            enddo
            ! State%vAvg (:,:,1) = State%vAvg (:,:,1) / State%Den(:,:,1)
        end associate

    end subroutine EvalMoments




! Species-level control

    function SpcEta2Den(spc, eta, bVol) result(D)
        !! Take a species' eta at a specific point, and sum moments to get its density and pressure
        type(SIFSpecies_T), intent(in) :: spc
            !! Species info
        real(rp), dimension(:), intent(in) :: eta
            !! Etas we are summing
        real(rp), intent(in) ::  bVol
            !! Flux tube volume [Rx/nT]
        
        integer :: k
        real(rp) :: D
            !! Density [#/cc]

        D = 0.0

        if (bVol <= 0) return

        do k=1,spc%N
            ! Di = Di + density_factor*sclmass(ikflavc(k))*eta(k)*vm**1.5
            D = D + (eta(k)/sclEta)/bVol
        enddo

    end function SpcEta2Den


    function SpcEta2Press(spc, eta, bVol) result(P)
        !! Take a species' eta at a specific point, and sum moments to get its density and pressure
        type(SIFSpecies_T), intent(in) :: spc
            !! Species info
        real(rp), dimension(:), intent(in) :: eta
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

        do k=1,spc%N
            alamc = abs(spc%alami(k) + spc%alami(k+1))/2.0
            !P = P + eta(k)*alamc*vm**2.5 * ev2J * 1.e6
            !! Note: 10^-9 from 1/sclEta cancels with 10^9 from Pa -> nPa
            P = P + eta(k)*alamc*bVol**(-5./3.) * ev2J * 1.e6
        enddo

    end function SpcEta2Press


    subroutine DkT2SpcEta(spc, eta, D, kT, vm)
        !! Take a density and pressure, and map it to SIF eta channels for given flavor
        type(SIFSpecies_T), intent(in) :: spc
            !! Species info
        real(rp), dimension(:), intent(inout) :: eta
            !! Len(spc%N) etas we need to populate
        real(rp), intent(in) :: D, kT, vm
            !! Density [#/cc], Energy [keV], bVol^-2/3 [(Rx/nT)^(-2/3)]

        integer :: k
            !! loop variable

        ! TODO: Logical check to determine which distribution we should do.
        !       Hard coding as Kappa for now

        do k=1,spc%N
            eta(k) = Kappa2Eta(D,kT,vm, spc%alami(k),spc%alami(k+1))
        enddo

        
    end subroutine DkT2SpcEta

    function Kappa2Eta(D,kT,vm,amin,amax,kapO) result(etaK)
        !! Convert density and temperature to eta at specific lambda value
        !! Adapted from eqn 3.12 from 10.1007/s11214-013-9982-9 ?

        real(rp), intent(in) :: D,kT,vm,amin,amax
            !! Density [#/cc], kT [keV], vm [(Rx/nT)^(-2/3)],
            !! min and max lambda vals [eV * (Rx/nT)^(2/3)]
        real(rp), intent(in), optional :: kapO
            !! Optional Kappa value

        real(rp) :: kap,kap15,E_ev,E0_ev
        real(rp) :: A0,kapgam,kapbar,kArg,delscl, etaK

        etaK = 0.0

        if (present(kapO)) then
            kap = kapO
        else
            kap = kapDefault
        endif

        kap15 = kap-1.5
        kapgam = gamma(kap+1.0)/gamma(kap-0.5)
        kapbar = kap15 !Should be kap-3/2 or kap
            ! TODO: check 10.1002/2015JA021166

        E0_ev = kT*kap15/kap*1.e3  ! [eV]
        E_ev = abs(amin+amax)/2.0*vm  ! Lambda center energy [eV]

        kArg = 1.0 + (E_ev/E0_ev)/kapbar
        delscl = vm*abs(amax-amin)/E0_ev  ! Channel width / temp
            ! TODO: why?

        A0 = (2.0/sqrt(PI)) * D/(vm**1.5)*sclEta  ! #/cc * Rx/nT * 1/nT -> 1/T
        etak = A0*kapgam/(kapbar**1.5) * sqrt(E_ev/E0_ev)*delscl*((kArg)**(-kap-1.0))
        !write(*,*)"D=",D,"vm=",vm
        !write(*,*)"  A0=",A0,"delscl=",delscl
        !write(*,*)"  etaK=",etaK

    end function Kappa2Eta

end module sifetautils