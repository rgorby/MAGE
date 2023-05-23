module sifCpl
    !! Contains functions used to pass data between voltron and SIF

    use math
    use imagtubes

    use sifdefs
    use sifTypes
    use sifCplTypes

    implicit none

    contains

    subroutine sifCpl_init(vApp, sApp, cplBase)
        type(voltApp_T), intent(in) :: vApp
        type(sifApp_T), intent(in) :: sApp
        type(sif_cplBase_T), intent(inout) :: cplBase

        associate(fromV=>cplBase%fromV, toV=>cplBase%toV, &
            sh=>sApp%Grid%shGrid)

        ! Init fromV first
            ! Allocations
            allocate(fromV%fLines (sh%is:sh%ie+1, sh%js:sh%je+1))
            allocate(fromV%ijTubes(sh%is:sh%ie+1, sh%js:sh%je+1))
            allocate(fromV%pot    (sh%Nt, sh%Np))

            ! Initial values
            fromV%tLastUpdate = -1.0*HUGE
            fromV%pot = 0.0

        ! Init toV next
            ! Allocations
            !Initial values
            toV%tLastUpdate = -1.0*HUGE

        end associate

        ! Point to default couplers

        cplBase%convertToSIF => sifCpl_Volt2SIF
        cplBase%convertToVoltron => sifCpl_SIF2Volt

    end subroutine sifCpl_init

    
    subroutine sifCpl_Volt2SIF(cplBase, vApp, sApp)
        class(sif_cplBase_T), intent(inout) :: cplBase
        type(voltApp_T), intent(in   ) :: vApp
        type(sifApp_T ), intent(inout) :: sApp
        
        integer :: i,j
        real(rp) :: seedR
        type(fLine_T) :: fLine
        ! Take info from cplBase%fromV and put it into SIF state

        ! Populate sif state with coupling info
        ! Tubes
        call imagTubes2SIF(cplBase%fromV%ijTubes, sApp%Model, sApp%Grid, sApp%State)
        ! Potential
        sApp%State%espot = cplBase%fromV%pot

        ! Determine active domain
        call setActiveDomain(sApp%Grid%shGrid, sApp%Grid%nB, sApp%State)


    end subroutine sifCpl_Volt2SIF

    subroutine sifCpl_SIF2Volt(cplBase, vApp, sApp)
        type(sif_cplBase_T), intent(inout) :: cplBase
        type(voltApp_T), intent(inout) :: vApp
        type(sifApp_T) , intent(in   ) :: sApp
    end subroutine sifCpl_SIF2Volt

!------
! Helpers
!------

    subroutine imagTubes2SIF(ijTubes, Model, Grid, State)
        !! Map 2D array of IMAGTubes to SIF State
        type(IMAGTube_T), dimension(:,:), intent(in) :: ijTubes
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T ), intent(in) :: Grid
        type(sifState_T), intent(inout) :: State

        integer :: i,j

        associate(sh=>Grid%shGrid)
            !! TODO: Actual checks and cleaning of data
            !!       e.g. choosing what happens when not all 4 corners are good
            !! TODO: Handle multispecies

            ! Corner quantities
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%is,sh%ie+1
                do j=sh%js,sh%je+1
                    State%xyzMin(i,j,:)  = ijTubes(i,j)%X_bmin / Model%planet%rp_m  ! xyzMin in Rp
                    State%topo(i,j)      = ijTubes(i,j)%topo
                enddo
            enddo
            ! Cell-centered quantities
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%is,sh%ie
                do j=sh%js,sh%je
                    State%Pavg(i,j,1)    = toCenter2D(ijTubes(i:i+1,j:j+1)%Pave) * 1.0e+9  ! Pa -> nPa
                    State%Davg(i,j,1)    = toCenter2D(ijTubes(i:i+1,j:j+1)%Nave) * 1.0e-6  ! #/m^3 -> #/cc
                    State%Bmin(i,j,ZDIR) = toCenter2D(ijTubes(i:i+1,j:j+1)%bmin) * 1.0e+9  ! Tesla -> nT
                    State%thc(i,j)  = toCenter2D(PI/2-ijTubes(i:i+1,j:j+1)%latc)
                    State%phc(i,j)       = toCenter2D(ijTubes(i:i+1,j:j+1)%lonc)
                    State%bvol(i,j)      = toCenter2D(ijTubes(i:i+1,j:j+1)%Vol) * 1.0e-9  ! Rp/T -> Rp/nT
                enddo
            enddo
        end associate

        ! TODO: Save calculating active domain for last, once we have all other info
        !State%active(i,j)   = ??

        
    end subroutine imagTubes2SIF


    subroutine setActiveDomain(shGrid, nB, State)
        type(ShellGrid_T), intent(in) :: shGrid
        integer, intent(in) :: nB
            !! Number of cells between open boundary and active domain
        type(sifState_T), intent(inout) :: State

        integer :: i,j

        logical, dimension(shGrid%is:shGrid%ie, &
                            shGrid%js:shGrid%je) :: closedCC


        State%active = SIFINACTIVE

        ! Mark any cell center with an open corner as open
        closedCC = .true.
        do i=shGrid%is, shGrid%ie
            do j=shGrid%js, shGrid%je
                if (any(State%topo(i:i+1, j:j+1) < 2)) then
                    closedCC(i,j) = .false.
                endif
            enddo
        enddo

        State%OCBDist = CalcOCBDist(shGrid, closedCC, nB)
        write(*,*) State%OCBDist
        ! Set zones
        where (State%OCBDist .eq. 0)
            State%active = SIFINACTIVE
        else where (State%OCBDist .le. nB)
            State%active = SIFBUFFER
        elsewhere
            State%active = SIFACTIVE
        end where


        contains

        function CalcOCBDist(shGrid, closedCC, nBnd) result(ocbDist)
            type(ShellGrid_T), intent(in) :: shGrid
                !! SIF shell grid
            logical, dimension(shGrid%is:shGrid%ie, &
                    shGrid%js:shGrid%je), intent(in) :: closedCC
                !! Whether cell centers are closed or open
            integer, intent(in) :: nBnd
                !! Number of desired layers between open/closed boundary and active domain

            integer :: iL, i, j
            integer, dimension(shGrid%is:shGrid%ie, &
                    shGrid%js:shGrid%je) :: ocbDist

            where (closedCC)
                ocbDist = nBnd + 1
            elsewhere
                ocbDist = 0
            end where


            ! Grow out from open/closed boundary, set proper distance for closed points
            do iL=1,nBnd
                do i=shGrid%is, shGrid%ie
                    do j=shGrid%js, shGrid%je
                        if (closedCC(i,j) .eq. .false.) then
                            cycle
                        else if ( (ocbDist(i,j) .eq. nBnd+1) .and. any(ocbDist(i-1:i+1,j-1:j+1) .eq. iL-1) ) then
                            !! If current point's distance hasn't been decided and its bordering cell with iL-1, this point is distance iL from ocb
                            ocbDist(i,j) = iL
                        endif
                    enddo
                enddo
            enddo
            
        end function CalcOCBDist

    end subroutine setActiveDomain

end module sifCpl