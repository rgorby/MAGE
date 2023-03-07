module sifetautils

    use planethelper

    use siftypes
    use sifspecieshelper

    implicit none

    real(rp), private, parameter :: kapDefault = 6.0

    contains

    subroutine DkT2SpcEta(spc, eta, D, kT, vm)
        !! Take a density and pressure, and map it to SIF eta channels for given flavor
        type(SIFSpecies_T), intent(in) :: spc
            !! Species info
        real(rp), dimension(:), intent(inout) :: eta
            !! Len(spc%N) etas we need to populate
        real(rp), intent(in) :: D, kT, vm
            !! Density [#/cc], Energy [keV], bVol^-2/3 [(Rx/nT)^(-2/3)]


        integer :: idx
            !! Index of species with flav in Grid%spc
        integer :: k
            !! loop variable

        ! TODO: Logical check to determine which distribution we should do.
        !       Hard coding as Kappa for now

        write(*,*)"flav=",spc%flav,"N=",spc%N
        !do k=1,spc%N-1
        do k=1,1
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

        A0 = (2/sqrt(PI)) * D/(vm**1.5)*1e9  ! [#/cc * Rx/T]
        etak = A0*kapgam/(kapbar**1.5) * sqrt(E_ev/E0_ev)*delscl*((kArg)**(-kap-1.0))
        write(*,*)"D=",D,"vm=",vm
        write(*,*)"  A0=",A0,"delscl=",delscl
        write(*,*)"  etaK=",etaK

    end function Kappa2Eta

end module sifetautils