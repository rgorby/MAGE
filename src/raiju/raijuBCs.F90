module raijuBCs

    use planethelper
    use raijudefs
    use raijutypes
    use raijuetautils

    implicit none

    contains

    subroutine applyRaijuBCs(Model, Grid, State, doWholeDomainO)
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        logical, optional, intent(in) :: doWholeDomainO

        integer :: i,j,s,k
        real(rp) :: kT,vm
        logical :: doWholeDomain
            !! Whether we should apply certain BCs (moments2eta mapping) to entire active domain or not
        logical, dimension(Grid%shGrid%isg:Grid%shGrid%ieg,&
                           Grid%shGrid%jsg:Grid%shGrid%jeg) :: doBC
            !! Whether we do moments to eta mapping at given grid cell

        integer :: psphIdx
            !! Index of plasmasphere
        real(rp) :: etaBelow
            !! Amount of eta below lowest lambda bound (every i,j gets a copy)
        real(rp) :: pressBelow
            !! Pressure of eta below lowest lambda bound (every i,j gets a copy)

        if (present(doWholeDomainO)) then
            doWholeDomain = doWholeDomainO
        else
            doWholeDomain = .false.
        endif

        ! Determine where to do BCs
        ! Will definitely be its own function later to do ellipse fitting, restriction of fast flows, etc
        if(doWholeDomain) then
            where (State%active .eq. RAIJUBUFFER .or. State%active .eq. raijuACTIVE)
                doBC = .true.
            elsewhere
                doBC = .false.
            end where
        else
            where (State%active .eq. RAIJUBUFFER)
                doBC = .true.
            elsewhere
                doBC = .false.
            end where
        endif

        psphIdx = spcIdx(Grid, F_PSPH)
        
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,s,vm,kT, etaBelow, pressBelow)
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
                        State%Davg(i,j,s), kT, vm, etaBelow)
                    
                    ! etaBelow has the amount of eta that is below the lowest lambda channel bound
                    !! TODO: Check to see if we are missing too much pressure
                    ! Maybe we want to put it in plasmasphere channel cause its cold H+
                    if (Model%doExcesstoPsph .and. Grid%spc(s)%mapExtraToPsph) then
                        State%eta(i,j,Grid%spc(psphIdx)%kStart) = State%eta(i,j,Grid%spc(psphIdx)%kStart) + etaBelow
                    endif

                enddo  ! s
            enddo  ! j
        enddo  ! i

        ! Do first round of determining active shells for each k
        ! We do this everywhere
        do s=1,Grid%nSpc
            call setActiveShellsByContribution(Grid%shGrid, Grid%spc(s), State%bVol, &
                State%eta(:,:,Grid%spc(s)%kStart:Grid%spc(s)%kEnd), &
                State%activeShells(:,Grid%spc(s)%kStart:Grid%spc(s)%kEnd), &  ! This is what we're writing to
                worthyFracO=Model%worthyFrac)
        enddo  ! s

    end subroutine applyRaijuBCs

!------
! Active I Shell calculations
!------

    subroutine setActiveShellsByContribution(shGrid, spc, bVol, etas, activeShellsOut, nSpacesO, worthyFracO)
        !! For each lambda channel, calculate active shells based on how much they contribute to the total pressure and/or density
        type(ShellGrid_T), intent(in) :: shGrid
        type(raijuSpecies_T), intent(in) :: spc
        real(rp), dimension(shGrid%isg:shGrid%ieg,shGrid%jsg:shGrid%jeg), intent(in) :: bVol
        real(rp), dimension(shGrid%isg:shGrid%ieg,shGrid%jsg:shGrid%jeg,spc%kStart:spc%kEnd), intent(in) :: etas
        logical, dimension(shGrid%isg:shGrid%ieg, spc%kStart:spc%kEnd), intent(inout) :: activeShellsOut
        integer, optional, intent(in) :: nSpacesO
            !! Number of i spaces between last good value and active i for species
        real(rp), optional, intent(in) :: worthyFracO

        integer :: i,j,k, iL, iU
        integer :: nSpaces
        real(rp) :: worthyFrac
        real(rp) :: alamc, kDen, kPress
        real(rp), dimension(shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg) :: spcDen, spcPress
        logical, dimension(shGrid%isg:shGrid%ieg, shGrid%jsg:shGrid%jeg) :: as2D

        if(present(nSpacesO)) then
            nSpaces = nSpacesO
        else
            nSpaces = nSpacesDef
        endif

        if(present(worthyFracO)) then
            worthyFrac = worthyFracO
        else
            worthyFrac = fracWorthyDef
        endif
        
        as2D = .false.
        activeShellsOut = .false.

        ! Start by getting total density and pressure for each i,j point
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j)
        do j=shGrid%jsg,shGrid%jeg
            do i=shGrid%isg,shGrid%ieg
                spcDen  (i,j) = SpcEta2Den  (spc, etas(i,j,:), bVol(i,j))
                spcPress(i,j) = spcEta2Press(spc, etas(i,j,:), bVol(i,j))
            enddo
        enddo

        ! Then calculate active shells for each k
        do k=spc%kStart,spc%kEnd

            ! Setup for this k
            as2D = .false.
            alamc = 0.5*abs(spc%alami(k) + spc%alami(k+1))

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,kDen,kPress)
            do j=shGrid%jsg,shGrid%jeg
                do i=shGrid%isg,shGrid%ieg

                    kDen   = etak2Den  (etas(i,j,k), bVol(i,j))
                    kPress = etak2Press(etas(i,j,k), alamc, bVol(i,j))
            
                    if (kDen/spcDen(i,j) > worthyFrac .or. kPress/spcPress(i,j) > worthyFrac) then
                        as2D(i,j) = .true.
                    endif
                enddo
            enddo

            ! Now collapse j
            !! Re-use as2D, first j element 
            do i=shGrid%isg,shGrid%ieg
                as2D(i,shGrid%jsg) = any(as2D(i,:))
            enddo

            ! Finally, assign to activeShells, including buffer spaces
            do i=shGrid%isg,shGrid%ieg
                iL = max(i-nSpaces, shGrid%isg)
                iU = min(i+nSpaces, shGrid%ieg)
                activeShellsOut(i,k) = any(as2D(iL:iU, shGrid%jsg))
            enddo

        enddo
    end subroutine setActiveShellsByContribution


end module raijuBCs