module raijuCplHelpers

    ! Base
    use voltCplTypes

    ! Raiju
    use raijutypes
    use raijuDomain
    use raijuCplTypes
    use raijuSpeciesHelper

    implicit none

    contains

    subroutine imagTubes2RAIJU(Model, Grid, State, ijTubes, f_MHD2SpcMap)
        !! Map 2D array of IMAGTubes to RAIJU State
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        type(IMAGTube_T), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                                    Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: ijTubes
        procedure(raijuMHD2SpcMap_T), pointer, intent(in) :: f_MHD2SpcMap
        

        integer :: i,j


        associate(sh=>Grid%shGrid)

            ! Map ijTube's definition of topology to RAIJU's
            where (ijTubes%topo == 2)
                State%topo = RAIJUCLOSED
            elsewhere
                State%topo = RAIJUOPEN
            end where

            ! Now that topo is set, we can calculate active domain
            call setActiveDomain(Model, Grid, State)

            ! Assign corner quantities
            !$OMP PARALLEL DO default(shared) collapse(1) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%isg,sh%ieg+1
                do j=sh%jsg,sh%jeg+1
                    State%xyzMin(i,j,:)  = ijTubes(i,j)%X_bmin / Model%planet%rp_m  ! xyzMin in Rp
                    State%thcon(i,j)     = PI/2-ijTubes(i,j)%latc
                    State%phcon(i,j)     = ijTubes(i,j)%lonc
                    State%Bmin(i,j,ZDIR) = ijTubes(i,j)%bmin * 1.0e+9  ! Tesla -> nT
                    State%bvol(i,j)      = ijTubes(i,j)%Vol  * 1.0e-9  ! Rp/T -> Rp/nT
                enddo
            enddo

            ! Assign cell-centered quantities
            !$OMP PARALLEL DO default(shared) collapse(1) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%isg,sh%ieg
                do j=sh%jsg,sh%jeg
                    ! Note: we are doing this for all cells regardless of their goodness
                    State%xyzMincc(i,j,XDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,XDIR))  ! [Rp]
                    State%xyzMincc(i,j,YDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,YDIR))  ! [Rp]
                    State%xyzMincc(i,j,ZDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,ZDIR))  ! [Rp]
                enddo
            enddo

            ! Use provided definition to map moments
            call f_MHD2SpcMap(Model, Grid, State, ijTubes)

        end associate
        
    end subroutine imagTubes2RAIJU


    subroutine defaultMHD2SpcMap(Model, Grid, State, ijTubes)
        !! Assumes:
        !!  MHD: single fluid
        !!  RAIJU: zero-energy psphere, hot electrons, hot protons
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        type(IMAGTube_T),  dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                                     Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: ijTubes

        integer :: i, j, k, s, sIdx
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
                        do s=0,Model%nFluidIn
                            ! This means all 4 corners are good, can do cell centered stuff
                            P = 0.25*(ijTubes(i  ,j)%Pave(s) + ijTubes(i  ,j+1)%Pave(s) &
                                    + ijTubes(i+1,j)%Pave(s) + ijTubes(i+1,j+1)%Pave(s)) * 1.0D+9 ! [Pa -> nPa]
                            D = 0.25*(ijTubes(i  ,j)%Nave(s) + ijTubes(i  ,j+1)%Nave(s) &
                                    + ijTubes(i+1,j)%Nave(s) + ijTubes(i+1,j+1)%Nave(s)) * 1.0D-6 ! [#/m^3 --> #/cc]
                                                 

                            State%Pavg(i,j,s) = P
                            State%Davg(i,j,s) = D

                        enddo
                    endif
                enddo
            enddo
        end associate

    end subroutine defaultMHD2SpcMap

end module raijuCplHelpers