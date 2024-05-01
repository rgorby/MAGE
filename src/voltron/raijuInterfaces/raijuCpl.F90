module raijuCpl
    !! Contains functions used to pass data between voltron and RAIJU

    use math
    use imagtubes
    use planethelper
    use shellGrid

    use raijudefs
    use raijuTypes
    use raijuCplTypes
    use raijuSpeciesHelper
    use raijuDomain

    implicit none

    contains

    subroutine raijuCpl_init(vApp, raiApp, cplBase)
        type(voltApp_T), intent(in) :: vApp
        type(raijuApp_T), intent(in) :: raiApp
        type(raiju_cplBase_T), intent(inout) :: cplBase

        integer, dimension(4) :: chGhosts

        associate(fromV=>cplBase%fromV, toV=>cplBase%toV, &
            sh=>raiApp%Grid%shGrid)

        ! Init fromV first
            ! Allocations
            allocate(fromV%fLines (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))
            allocate(fromV%ijTubes(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))

            ! Shell Grid inits
            !chGhosts(NORTH) = sh%Ngn
            !chGhosts(SOUTH) = sh%Ngs
            !chGhosts(EAST) = sh%Nge
            !chGhosts(WEST) = sh%Ngw
            chGhosts(1) = sh%Ngn
            chGhosts(2) = sh%Ngs
            chGhosts(3) = sh%Nge
            chGhosts(4) = sh%Ngw
            call GenChildShellGrid(sh, fromV%shGr, "raijuCpl", nGhosts=chGhosts)
            call initShellVar(fromV%shGr, SHCORNER, fromV%pot)
            !allocate(fromV%pot    (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))

            ! Initial values
            fromV%tLastUpdate = -1.0*HUGE
            fromV%pot%data = 0.0
            fromV%pot%mask = .true.

        ! Init toV next
            ! Allocations
            !Initial values
            toV%tLastUpdate = -1.0*HUGE

        end associate

        ! If using user IC, let user determine coupling
        !  (this assumes icStr was already set by raijuInitState)
        if (trim(raiApp%Model%icStr) .eq. "USER") then
            write(*,*) "Can't do user initCpl yet"
            stop
        else
            ! Point to default coupling functions
            cplBase%convertToRAIJU => raijuCpl_Volt2RAIJU
            cplBase%convertToVoltron => raijuCpl_RAIJU2Volt
            cplBase%fromV%mhd2spcMap => defaultMHD2SpcMap
        endif

    end subroutine raijuCpl_init

    
    subroutine raijuCpl_Volt2RAIJU(cplBase, vApp, raiApp)
        !! Take info from cplBase%fromV and put it into RAIJU state
        class(raiju_cplBase_T), intent(in) :: cplBase
        type(voltApp_T)  , intent(in   ) :: vApp
        type(raijuApp_T ), intent(inout) :: raiApp

        ! Populate raiju state with coupling info
        ! Tubes
        call imagTubes2RAIJU(raiApp%Model, raiApp%Grid, raiApp%State, &
                cplBase%fromV%ijTubes, &
                cplBase%fromV%mhd2spcMap)
        ! Potential
        raiApp%State%espot(:,:) = cplBase%fromV%pot%data(:,:)

    end subroutine raijuCpl_Volt2RAIJU


    subroutine raijuCpl_RAIJU2Volt(cplBase, vApp, raiApp)
        class(raiju_cplBase_T), intent(inout) :: cplBase
        type(voltApp_T), intent(inout) :: vApp
        type(raijuApp_T) , intent(in   ) :: raiApp
    end subroutine raijuCpl_RAIJU2Volt

!------
! Helpers and defaults
!------
    
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


end module raijuCpl