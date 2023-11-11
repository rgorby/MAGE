module raijuCpl
    !! Contains functions used to pass data between voltron and RAIJU

    use math
    use imagtubes
    use planethelper

    use raijudefs
    use raijuTypes
    use raijuCplTypes
    use raijuSpeciesHelper
    implicit none

    contains

    subroutine raijuCpl_init(vApp, sApp, cplBase)
        type(voltApp_T), intent(in) :: vApp
        type(raijuApp_T), intent(in) :: sApp
        type(raiju_cplBase_T), intent(inout) :: cplBase

        associate(fromV=>cplBase%fromV, toV=>cplBase%toV, &
            sh=>sApp%Grid%shGrid)

        ! Init fromV first
            ! Allocations
            allocate(fromV%fLines (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))
            allocate(fromV%ijTubes(sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))
            allocate(fromV%pot    (sh%isg:sh%ieg  , sh%jsg:sh%jeg  ))

            ! Initial values
            fromV%tLastUpdate = -1.0*HUGE
            fromV%pot = 0.0

        ! Init toV next
            ! Allocations
            !Initial values
            toV%tLastUpdate = -1.0*HUGE

        end associate

        ! If using user IC, let user determine coupling
        !  (this assumes icStr was already set by raijuInitState)
        if (trim(sApp%Model%icStr) .eq. "USER") then
            !call RAIJUinitCplUserIC(sApp%Model, sApp%Grid, sApp%State, cplBase)
            !call userInitCplFunc(vApp, sApp, cplBase)
            write(*,*) "Can't do user initCpl yet"
            stop
        else
            ! Point to default coupling functions
            cplBase%convertToRAIJU => raijuCpl_Volt2RAIJU
            cplBase%convertToVoltron => raijuCpl_RAIJU2Volt
            cplBase%fromV%mhd2spcMap => defaultMHD2SpcMap
        endif

    end subroutine raijuCpl_init

    
    subroutine raijuCpl_Volt2RAIJU(cplBase, vApp, sApp)
        !! Take info from cplBase%fromV and put it into RAIJU state
        class(raiju_cplBase_T), intent(in) :: cplBase
        type(voltApp_T), intent(in   ) :: vApp
        type(raijuApp_T ), intent(inout) :: sApp
        
        
        ! Start with calculating the active domain so that ingestion can make decisions based on it
        call setActiveDomain(sApp%Grid%shGrid, sApp%Grid%nB, cplBase%fromV%ijTubes, sApp%State)

        ! Populate raiju state with coupling info
        ! Tubes
        call imagTubes2RAIJU(sApp%Model, sApp%Grid, sApp%State, &
                cplBase%fromV%ijTubes, &
                cplBase%fromV%mhd2spcMap)
        ! Potential
        sApp%State%espot(:,:) = cplBase%fromV%pot(:,:)

    end subroutine raijuCpl_Volt2RAIJU


    subroutine raijuCpl_RAIJU2Volt(cplBase, vApp, sApp)
        type(raiju_cplBase_T), intent(inout) :: cplBase
        type(voltApp_T), intent(inout) :: vApp
        type(raijuApp_T) , intent(in   ) :: sApp
    end subroutine raijuCpl_RAIJU2Volt

!------
! Helpers and defaults
!------

    subroutine setActiveDomain(sh, nB, ijTubes, State)
        type(ShellGrid_T), intent(in) :: sh
        integer, intent(in) :: nB
            !! Number of cells between open boundary and active domain
        type(IMAGTube_T), dimension(sh%isg:sh%ieg+1,sh%jsg:sh%jeg+1), intent(in) :: ijTubes
        type(raijuState_T), intent(inout) :: State

        integer :: i,j
        logical, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg) :: closedCC


        ! We make sure the whole thing is initialized later, but just in case
        State%active = RAIJUINACTIVE

        ! Mark any cell center with an open corner as open
        closedCC = .true.
        do i=sh%isg, sh%ieg
            do j=sh%jsg, sh%jeg
                !!! NOTE: At this point, we are still using ijTube's definition of topo,
                !!!  where topo = 0 (solar wind), 1 (half-closed), 2 (both ends closed)
                if (any(ijTubes(i:i+1, j:j+1)%topo < 2)) then
                    closedCC(i,j) = .false.
                endif
            enddo
        enddo

        State%OCBDist = CalcOCBDist(sh, closedCC, nB)

        ! Set zones
        where (State%OCBDist .eq. 0)
            State%active = RAIJUINACTIVE
        else where (State%OCBDist .le. nB)
            State%active = RAIJUBUFFER
        elsewhere
            State%active = RAIJUACTIVE
        end where


        contains

        function CalcOCBDist(sh, closedCC, nBnd) result(ocbDist)
            type(ShellGrid_T), intent(in) :: sh
                !! RAIJU shell grid
            logical, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(in) :: closedCC
                !! Whether cell centers are closed or open
            integer, intent(in) :: nBnd
                !! Number of desired layers between open/closed boundary and active domain

            integer :: iLayer, i, j, iL, iU, jL, jU
            integer, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg) :: ocbDist

            where (closedCC)
                ocbDist = nBnd + 1
            elsewhere
                ocbDist = 0
            end where


            ! Grow out from open/closed boundary, set proper distance for closed points
            do iLayer=1,nBnd
                !! Note: We can safely parallelize within eah iLayer loop if needed. But this isn't a bottleneck right now
                do i=sh%isg, sh%ieg
                    do j=sh%jsg, sh%jeg
                        iL = max(i-1, sh%isg)
                        iU = min(i+1, sh%ieg)
                        jL = max(j-1, sh%jsg)
                        jU = min(j+1, sh%jeg)

                        if (closedCC(i,j) .eq. .false.) then
                            cycle
                        else if ( (ocbDist(i,j) .eq. nBnd+1) .and. any(ocbDist(iL:iU,jL:jU) .eq. iLayer-1) ) then
                            !! If current point's distance hasn't been decided and its bordering cell with iL-1, this point is distance iL from ocb
                            ocbDist(i,j) = iLayer
                        endif
                    enddo
                enddo
            enddo
            
        end function CalcOCBDist

    end subroutine setActiveDomain

    subroutine imagTubes2RAIJU(Model, Grid, State, ijTubes, f_MHD2SpcMap)
        !! Map 2D array of IMAGTubes to RAIJU State
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        type(IMAGTube_T), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                                    Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: ijTubes
        procedure(raijuMHD2SpcMap_T), pointer, intent(in) :: f_MHD2SpcMap
        

        integer :: i,j
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                            Grid%shGrid%jsg:Grid%shGrid%jeg) :: isGood

        associate(sh=>Grid%shGrid)
            !! TODO: Actual checks and cleaning of data
            !!       e.g. choosing what happens when not all 4 corners are good
            !! TODO: Handle multispecies

            ! Corner quantities
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%isg,sh%ieg+1
                do j=sh%jsg,sh%jeg+1
                    State%xyzMin(i,j,:)   = ijTubes(i,j)%X_bmin / Model%planet%rp_m  ! xyzMin in Rp
                    State%topo(i,j)       = ijTubes(i,j)%topo
                    State%thcon(i,j) = PI/2-ijTubes(i,j)%latc
                    State%phcon(i,j)      = ijTubes(i,j)%lonc
                enddo
            enddo

            ! Map ijTube's definition of topology to RAIJU's
            where (ijTubes%topo == 2)
                State%topo = RAIJUCLOSED
            elsewhere
                State%topo = RAIJUOPEN
            end where


            ! Make sure we can safely cell-average (all 4 corners are closed field lines)
            where (State%active .eq. RAIJUBUFFER .or. State%active .eq. RAIJUACTIVE)
                isGood = .true.
            elsewhere
                isGood = .false.
            end where

            ! Cell-centered quantities
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%isg,sh%ieg
                do j=sh%jsg,sh%jeg
                    if (isGood(i,j)) then
                        State%Bmin(i,j,ZDIR) = toCenter2D(ijTubes(i:i+1,j:j+1)%bmin) * 1.0e+9  ! Tesla -> nT
                        State%bvol(i,j)      = toCenter2D(ijTubes(i:i+1,j:j+1)%Vol) * 1.0e-9  ! Rp/T -> Rp/nT
                    endif

                    ! Always do xyzmin
                    State%xyzMincc(i,j,XDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,XDIR))
                    State%xyzMincc(i,j,YDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,YDIR))
                    State%xyzMincc(i,j,ZDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,ZDIR))
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