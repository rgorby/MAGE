module raijuDomain

    use raijudefs
    use raijuTypes

    implicit none

    contains


    subroutine setActiveDomain(Model, Grid, State)
        !! Sets State%active array
        !! Determines Inactive -> buffer -> active cells
        type(raijuModel_T), intent(in   ) :: Model
        type(raijuGrid_T ), intent(in   ) :: Grid
        type(raijuState_T), intent(inout) :: State

        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: isInactive
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: isBuffer


        call getInactiveCells(Model, Grid, State, isInactive)
        call getBufferCells  (Model, Grid, State, isInactive, isBuffer, State%OCBDist)

        ! Set zones
        where (isInactive)
            State%active = RAIJUINACTIVE
        else where (isBuffer)
            State%active = RAIJUBUFFER
        elsewhere
            State%active = RAIJUACTIVE
        end where

    end subroutine setActiveDomain


    subroutine getInactiveCells(Model, Grid, State, isInactive)
        !! Applies series of criteria to determine which cells should be marked as inactive
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: isInactive

        integer :: i,j
        real(rp), dimension(2,2) :: bminSquare

        isInactive = .false.
        do j=Grid%shGrid%jsg, Grid%shGrid%jeg
            !do i=Grid%shGrid%ieg, Grid%shGrid%isg,-1
            do i=Grid%shGrid%isg, Grid%shGrid%ieg
                !! NOTE: idir in descending order
                !! This means we are going from inner R boundary to outer R

                ! Any cell with an open corner is bad
                if (any(State%topo(i:i+1,j:j+1) .eq. RAIJUOPEN)) then
                    isInactive(i,j) = .true.
                endif

                if (any(State%vaFrac(i:i+1,j:j+1) .le. Model%vaFracThresh)) then
                    isInactive(i,j) = .true.
                endif

                bminSquare(1,1) = norm2(State%Bmin(i  ,j  ,:))
                bminSquare(2,1) = norm2(State%Bmin(i+1,j  ,:))
                bminSquare(1,2) = norm2(State%Bmin(i  ,j+1,:))
                bminSquare(2,2) = norm2(State%Bmin(i+1,j+1,:))
                if (any( bminSquare .le. Model%bminThresh) ) then
                !if (any( norm2(State%Bmin(i:i+1,j:j+1,:),dim=3) .le. Model%bminThresh) ) then
                    isInactive(i,j) = .true.
                endif

            enddo
        enddo
    end subroutine getInactiveCells


    subroutine getBufferCells(Model, Grid, State, isInactive, isBuffer, ocbDistO)
        !! Determines which cells should be marked as inactive
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(in) :: State
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(in) :: isInactive
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: isBuffer
        integer, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), optional, intent(inout) :: ocbDistO
            !! Optional array to get ocbDist out of this function if desired

        integer, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg) :: ocbDist
            !! Distance each cell is to open-closed boundary, up to nB + 1
        
        ! As a first pass, we will blanket enforce that there must always be nB buffer between any inactive and active cells
        ! We can do more complex stuff later

        call CalcOCBDist(Grid%shGrid, isInactive, Grid%nB, ocbDist)

        where (ocbDist > 0 .and. ocbDist .le. Grid%nB)
            isBuffer = .true.
        elsewhere
            isBuffer = .false.
        end where

        if (present(ocbDistO)) then
            ocbDistO = ocbDist
        endif

    end subroutine getBufferCells


    subroutine CalcOCBDist(sh, isInactive, nBnd, ocbDist)
        !! Calculates distance each cell is from open/closed boundary, up to nBnd cells away
        type(ShellGrid_T), intent(in) :: sh
            !! RAIJU shell grid
        logical, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(in) :: isInactive
            !! Whether cell centers are inactive
        integer, intent(in) :: nBnd
            !! Number of desired layers between open/closed boundary and active domain
        integer, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(out) :: ocbDist


        integer :: iLayer, i, j, iL, iU, jL, jU


        where ( isInactive )
            ocbDist = 0
        elsewhere
            ocbDist = nBnd + 1
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

                    if (isInactive(i,j)) then
                        cycle
                    !else if ( (ocbDist(i,j) .eq. nBnd+1) .and. any(ocbDist(iL:iU,jL:jU) .eq. iLayer-1) ) then  ! This version includes the diagonal cells, which are technically not adjacent to our cell
                    else if ( (ocbDist(i,j) .eq. nBnd+1) .and. ( any(ocbDist(iL:iU,j) .eq. iLayer-1) .or. any(ocbDist(i,jL:jU) .eq. iLayer-1) ) ) then  ! Doesn't include diagonals
                        !! If current point's distance hasn't been decided and its bordering cell with iL-1, this point is distance iL from ocb
                        ocbDist(i,j) = iLayer
                    endif
                enddo
            enddo
        enddo
        
    end subroutine CalcOCBDist



end module raijuDomain