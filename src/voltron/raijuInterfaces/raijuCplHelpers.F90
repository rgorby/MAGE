module raijuCplHelpers

    use imagtubes

    use raijutypes
    use raijuSpeciesHelper

    implicit none

    contains

    subroutine defaultMHD2SpcMap(Model, Grid, State, ijTubes)
        !! Assumes:
        !!  MHD: single fluid
        !!  RAIJU: zero-energy psphere, hot electrons, hot protons
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        type(IMAGTube_T),  dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                                     Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: ijTubes

        integer :: i, j, k, sIdx
        real(rp) :: P, D
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood

        ! First clear out our previous moments input state
        State%Pavg = 0.0
        State%Davg = 0.0

        where (State%active .eq. RAIJUBUFFER .or. State%active .eq. RAIJUACTIVE)
            isGood = .true.
        elsewhere
            isGood = .false.
        end where

        associate(sh=>Grid%shGrid)

            do i=sh%isg,sh%ieg
                do j=sh%jsg,sh%jeg
                    if (isGood(i,j)) then
                        ! This means all 4 corners are good, can do cell centered stuff
                        P = toCenter2D(ijTubes(i:i+1,j:j+1)%Pave) * 1.0e+9  ! Pa -> nPa
                        D = toCenter2D(ijTubes(i:i+1,j:j+1)%Nave) * 1.0e-6  ! #/m^3 -> #/cc                      

                        ! First do hot protons
                        sIdx = spcIdx(Grid, F_HOTP)
                        State%Pavg(i,j,sIdx) = P / (1.0 + 1.0/Model%tiote)
                        State%Davg(i,j,sIdx) = D

                        ! Electrons
                        sIdx = spcIdx(Grid, F_HOTE)
                        State%Pavg(i,j,sIdx) = P / (1.0 + Model%tiote)
                        State%Davg(i,j,sIdx) = D

                        !! Note: Plasmasphere input density is zero because we're doing a single fluid
                    endif
                enddo
            enddo
            sIdx = spcIdx(Grid, F_HOTP)
            write(*,*)"Max ",sIdx," Davg_in=",maxval(State%Davg(:,:,sIdx))
            sIdx = spcIdx(Grid, F_HOTE)
            write(*,*)"Max ",sIdx," Davg_in=",maxval(State%Davg(:,:,sIdx))
        end associate

    end subroutine defaultMHD2SpcMap

end module raijuCplHelpers