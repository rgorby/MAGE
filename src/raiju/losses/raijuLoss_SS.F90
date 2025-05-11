module raijuLoss_SS
    !! Strong scattering for RAIJU electrons

    use kdefs, only : mec2, vc_cgs
    use raijudefs
    use raijutypes
    use raijuSpeciesHelper
    use raijuEleLossHelper

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

        ! Idc what planet you are, electrons gonna scatter
        this%reqsGood = .true.
        ! And they gonna precipitate
        this%isPrecip = .true.

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

        tau = CalcTau_StrongScattering(Model, Grid, State, i, j, k)

    end function SSLossCalcTau
    

end module raijuLoss_SS