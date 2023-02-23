module sifetautils

    use siftypes

    implicit none

    contains

    subroutine DkT2SpcEta(Model, State, flav, D, kT)
        type(sifModel_T), intent(in) :: Model
            !! Used to get species info and Moments2Eta settings
        type(sifState_T), intent(inout) :: State
            !! Where we gotta put the etas
        integer, intent(in) :: flav
            !! Flavor species identifier
        real(rp), intent(in) :: D, kT
            !! Density [#/cc] and energy [keV]

        write(*,*)"TODO: DkT2SpcEta"
        
    end subroutine DkT2SpcEta

end module sifetautils