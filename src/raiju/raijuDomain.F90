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
        real(rp) :: xyMin
        real(rp), dimension(2,2) :: bminSquare

        isInactive = .false.
        do j=Grid%shGrid%jsg, Grid%shGrid%jeg
            !do i=Grid%shGrid%ieg, Grid%shGrid%isg,-1
                !! NOTE: idir in descending order
                !! This means we are going from inner R boundary to outer R
            do i=Grid%shGrid%isg, Grid%shGrid%ieg
                

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

                xyMin = sqrt(State%xyzMin(i,j,XDIR)**2 + State%xyzMin(i,j,YDIR)**2)
                ! Simple circle limit
                if ( (xyMin > Model%maxTail_buffer) .or. xyMin > Model%maxSun_buffer) then
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

        integer :: i,j
        real(rp) :: xyMin
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

        ! Also limit to max allowed sun/tail extent
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg
                if (ocbDist(i,j) > 0) then  !! OCBdist will be zero in inactive zone, we want to avoid editing those cells
                    xyMin = sqrt(State%xyzMin(i,j,XDIR)**2 + State%xyzMin(i,j,YDIR)**2)
                    ! Simple circle limit
                    if ( (xyMin > Model%maxTail_active) .or. xyMin > Model%maxSun_active) then
                        isBuffer(i,j) = .true.
                    endif
                endif
            enddo
        enddo

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


    subroutine calcMapJacNorm(Grid, bMin, jacNorm)
        type(raijuGrid_T), intent(in) :: Grid
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,Grid%shGrid%jsg:Grid%shGrid%jeg+1, 3), intent(in) :: bMin
            !! xyz coordinates of bMin surface
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg,Grid%shGrid%jsg:Grid%shGrid%jeg), intent(inout) :: jacNorm
            !! Cell-centered jacobian norm we calculate

        integer :: i,j, dim
        real(rp) :: ddTheta, ddPhi, aij
        jacNorm = 0.0

        ! Iterate through cell centers
        do j=Grid%shGrid%jsg,Grid%shGrid%jeg
            do i=Grid%shGrid%isg,Grid%shGrid%ieg

                ! Doing L 2,1 norm (https://en.wikipedia.org/wiki/Matrix_norm)
                ! Derivative by mapping cell corner values to edges, gradient across cell
                ! Deriv w.r.t. theta first:
                ddTheta = 0.0
                do dim = 1,3
                    aij = 0.5*(bMin(i+1,j+1,dim) + bMin(i+1,j,dim) \
                            - (bMin(i,j+1,dim) + bMin(i,j,dim))) /  Grid%lenFace(i,j,RAI_PH)
                    ddTheta = ddtheta + aij**2
                enddo

                ddPhi = 0.0
                do dim=1,3
                    aij = 0.5*(bMin(i+1,j+1,dim) + bMin(i,j+1,dim) \
                            - (bMin(i+1,j,dim) - bMin(i,j,dim))) / (0.5*(Grid%lenFace(i+1,j,RAI_TH) + Grid%lenFace(i,j,RAI_TH)))
                    ddPhi = ddPhi + aij**2
                enddo
                
                jacNorm(i,j) = sqrt(ddTheta) + sqrt(ddPhi)
            enddo
        enddo
    end subroutine calcMapJacNorm

end module raijuDomain