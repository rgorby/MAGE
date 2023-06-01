module sifBCs

    use planethelper

    use sifdefs
    use siftypes
    use sifetautils
    use sifadvancer

    implicit none

    contains

    subroutine applySifBCs(Model, Grid, State, doWholeDomainO)
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T) , intent(in) :: Grid
        type(sifState_T), intent(inout) :: State
        logical, optional, intent(in) :: doWholeDomainO

        integer :: i,j,s,k
        real(rp) :: kT,vm
        logical :: doWholeDomain
            !! Whether we should apply certain BCs (moments2eta mapping) to entire active domain or not
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: doBC

        if (present(doWholeDomainO)) then
            doWholeDomain = doWholeDomainO
        else
            doWholeDomain = .false.
        endif

        ! Determine where to do BCs
        if(doWholeDomain) then
            where (State%active .eq. SIFBUFFER .or. State%active .eq. SIFACTIVE)
                doBC = .true.
            elsewhere
                doBC = .false.
            end where
        else
            where (State%active .eq. SIFBUFFER)
                doBC = .true.
            elsewhere
                doBC = .false.
            end where
        endif

        
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,vm,kT)
        do i=Grid%shGrid%isg,Grid%shGrid%ieg
            do j=Grid%shGrid%jsg,Grid%shGrid%jeg
                ! Skip if we should leave point alone
                if(.not. doBC(i,j)) then
                    cycle
                endif

                vm = State%bvol(i,j)**(-2./3.)
                do s=1,Grid%nSpc
                    kT = DP2kT(State%Davg(i,j,s), State%Pavg(i,j,s))  ! [keV]
                    call DkT2SpcEta(Model,Grid%spc(s), &
                        State%eta(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd),&
                        State%Davg(i,j,s), kT, vm)
                    !write(*,*)s,i,j,kT,maxval(State%eta(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd))
                enddo  ! s
            enddo  ! j
        enddo  ! i

        ! Do first round of determining active shells for each k
        do s=1,Grid%nSpc
            do k=Grid%spc(s)%kStart,Grid%spc(s)%kEnd
                State%activeShells(:,k) = setLambdaActiveShells(Grid%shGrid, Grid%spc(s), State%bVol, &
                        State%eta(:,:,Grid%spc(s)%kStart:Grid%spc(s)%kEnd), &
                        k, &
                        worthyFracO=Model%worthyFrac)
            enddo  ! k
        enddo  ! s

    end subroutine applySifBCs


end module sifBCs