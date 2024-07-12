module raijuCplHelpers

    ! Base
    use voltCplTypes
    use planethelper

    ! Raiju
    use raijutypes
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
        

        integer :: i,j,s
        real(rp) :: P, D
        !real(rp) :: VaMKS, Tiev, csMKS


        associate(sh=>Grid%shGrid)

            ! Copy over all the tube info we want to have available to us

            ! Map ijTube's definition of topology to RAIJU's
            where (ijTubes%topo == 2)
                State%topo = RAIJUCLOSED
            elsewhere
                State%topo = RAIJUOPEN
            end where

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
                    State%vaFrac(i,j)    = ijTubes(i,j)%wIMAG
                enddo
            enddo

            State%bvol_cc = 0.0
            ! Assign cell-centered quantities
            !$OMP PARALLEL DO default(shared) collapse(1) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,s,P,D)
            do i=sh%isg,sh%ieg
                do j=sh%jsg,sh%jeg
                    ! Note: we are doing this for all cells regardless of their goodness
                    State%xyzMincc(i,j,XDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,XDIR))  ! [Rp]
                    State%xyzMincc(i,j,YDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,YDIR))  ! [Rp]
                    State%xyzMincc(i,j,ZDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,ZDIR))  ! [Rp]


                    ! Now we only calculate values for good cells
                    if (all(State%topo(i:i+1,j:j+1) .eq. RAIJUCLOSED)) then

                        do s=0,Model%nFluidIn
                            ! This means all 4 corners are good, can do cell centered stuff
                            P = 0.25*(ijTubes(i  ,j)%Pave(s) + ijTubes(i  ,j+1)%Pave(s) &
                                    + ijTubes(i+1,j)%Pave(s) + ijTubes(i+1,j+1)%Pave(s)) * 1.0D+9 ! [Pa -> nPa]
                            D = 0.25*(ijTubes(i  ,j)%Nave(s) + ijTubes(i  ,j+1)%Nave(s) &
                                    + ijTubes(i+1,j)%Nave(s) + ijTubes(i+1,j+1)%Nave(s)) * 1.0D-6 ! [#/m^3 --> #/cc]
                                                 

                            State%Pavg(i,j,s) = P
                            State%Davg(i,j,s) = D

                            State%bvol_cc(i,j) = toCenter2D(State%bvol(i:i+1,j:j+1))
                        enddo

                        ! Do our own "wIMAG" calculation here so we ensure we use the fluid we want to (Bulk)
                        ! Calculate the fraction of Alfven speed to total velocity
                        !VaMKS = flux tube arc length [km] / Alfven crossing time [s]
                        !VaMKS = (ijTube%Lb*planet%rp_m*1.0e-3)/ijTube%Tb 
                        !!CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
                        !TiEV = (1.0e+3)*DP2kT(State%Davg(i,j,0),State%Pavg(i,j,0)) !Temp in eV
                        !CsMKS = 9.79*sqrt((5.0/3)*TiEV)
                        !State%vaFrac(i,j) = VaMKS/( sqrt(VaMKS**2.0 + CsMKS**2.0) + VebMKS)

                        ! Never mind, we should let MHD decide since its more aware of what the right sound speed should be
                        
                    endif
                enddo
            enddo

            ! Use provided definition to map moments
            !call f_MHD2SpcMap(Model, Grid, State, ijTubes)

        end associate
        
    end subroutine imagTubes2RAIJU


    subroutine defaultMHD2SpcMap(Model, Grid, State, ijTubes)
        !! Sort of out of date
        !! Keeping it here in case we do eventually want to map some non-trivial way
        !! to our Pave, Dave. But I don't think we will
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        type(IMAGTube_T),  dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                                     Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: ijTubes

        integer :: i, j, k, s, sIdx
        real(rp) :: P, D

        ! First clear out our previous moments input state
        State%Pavg = 0.0
        State%Davg = 0.0

        associate(sh=>Grid%shGrid)
            do i=sh%isg,sh%ieg
                do j=sh%jsg,sh%jeg
                    if (all(State%topo(i:i+1,j:j+1) .eq. RAIJUCLOSED)) then

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