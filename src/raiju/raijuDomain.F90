module raijuDomain

    use raijudefs
    use raijuTypes

    implicit none

    contains


    subroutine setActiveDomain(sh, nB, State)
        type(ShellGrid_T), intent(in) :: sh
        integer, intent(in) :: nB
            !! Number of cells between open boundary and active domain
        type(raijuState_T), intent(inout) :: State

        integer :: i,j
        logical, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg) :: closedCC


        ! We make sure the whole thing is initialized later, but just in case
        State%active = RAIJUINACTIVE

        ! Mark any cell center with an open corner as open
        closedCC = .true.
        do i=sh%isg, sh%ieg
            do j=sh%jsg, sh%jeg
                if (any(State%topo .eq. RAIJUOPEN)) then
                    closedCC(i,j) = .false.
                endif
            enddo
        enddo

        call CalcOCBDist(sh, closedCC, nB, State%OCBDist)

        ! Set zones
        where (State%OCBDist .eq. 0)
            State%active = RAIJUINACTIVE
        else where (State%OCBDist .le. nB)
            State%active = RAIJUBUFFER
        elsewhere
            State%active = RAIJUACTIVE
        end where

    end subroutine setActiveDomain


    subroutine CalcOCBDist(sh, closedCC, nBnd, ocbDist)
        !! Calculates distance each cell is from open/closed boundary, up to nBnd cells away
        type(ShellGrid_T), intent(in) :: sh
            !! RAIJU shell grid
        logical, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(in) :: closedCC
            !! Whether cell centers are closed or open
        integer, intent(in) :: nBnd
            !! Number of desired layers between open/closed boundary and active domain
        integer, dimension(sh%isg:sh%ieg,sh%jsg:sh%jeg), intent(out) :: ocbDist


        integer :: iLayer, i, j, iL, iU, jL, jU


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
        
    end subroutine CalcOCBDist



end module raijuDomain