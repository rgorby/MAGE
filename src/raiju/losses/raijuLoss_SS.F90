module raijuLoss_SS
    !! Strong scattering for RAIJU electrons

    use kdefs, only : mec2, vc_cgs
    use raijudefs
    use raijutypes
    use raijuSpeciesHelper

    implicit none

    type, extends(baseRaijuLoss_T) :: raiLoss_SS_T
        logical :: reqsGood = .false.
        contains
        
        procedure :: doInit => SSLossInit
        !procedure :: doUpdate
        procedure :: isValidSpc => SSLossValidSpc
        procedure :: calcTau    => SSLossCalcTau

    end type raiLoss_SS_T

    contains

    subroutine SSLossInit(this, Model, Grid, xmlInp)
        class(raiLoss_SS_T), intent(inout) :: this
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(XML_Input_T) , intent(in) :: xmlInp

        ! Idk what planet you are, electrons gonna scatter
        this%reqsGood = .true.

    end subroutine SSLossInit

    function SSLossValidSpc(this, spc) result(isValid)
        class(raiLoss_SS_T), intent(in) :: this
        type(raijuSpecies_T), intent(in) :: spc
        logical :: isValid
        
        isValid = .false.
        if (.not. this%reqsGood) return

        if ( (spc%spcType .eq. RAIJUELE) ) then
            isValid = .true.
        endif
    end function SSLossValidSpc


    function SSLossCalcTau(this, Model, Grid, State, i, j, k) result(tau)
        !! Calculates strong scattering rate, according to Schulz 1998
        !! tau ~ [2*FTV*Bfp/(1-eta)](gamma*m0/p)
        !! FTV = flux tube volume, Bfp = B-field at foot point, eta - back-scattering rate
        !! eta is backscatter rate at alitude h, here eta=2/3.
        !! gamma = m/m0 is relativisitc factor, p is particle momentum.
        !!       = mc2/m0c2 = (m0c2+K)/m0c2 = 1+K/mec2 ! mec2=0.511 is me*c^2 in MeV
        !! m = m0/sqrt(1-v^2/c^2)
        !! V = c*1/sqrt(1-1/gammar2)
        class(raiLoss_SS_T), intent(in) :: this
        type(raijuModel_T ), intent(in) :: Model
        type(raijuGrid_T  ), intent(in) :: Grid
        type(raijuState_T ), intent(in) :: State
        integer, intent(in) :: i, j, k
        real(rp) :: tau

        real(rp) :: KE, gammar, V
        real(rp) :: eta_scatter = 2./3.

        tau = HUGE

        KE = abs(Grid%alamc(k))*State%bvol_cc(i,j)**(-2./3.) * 1.0D-3  ! Energy [keV]
        gammar = 1.0 + (KE*1.0D-3)/mec2  ! Gamma with 1 + MeV/MeV
        V = (vc_cgs*1e-2)*sqrt(1.0 - 1.0/gammar**2)/Model%planet%rp_m  ! [Rp/s]
        tau = 2.0*State%bvol_cc(i,j)*Grid%Bmag(i,j)/(1.0 - eta_scatter) / V*gammar  ! [Rp/nT * nT / (Rp/s) = s]

        !if (i==22 .and. j==180 ) then
        !    write(*,*)"ijk=",i,j,k
        !    write(*,*)" ",KE,gammar,V,tau
        !endif

    end function SSLossCalcTau


end module raijuLoss_SS