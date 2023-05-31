module sifBCs

    use planethelper

    use sifdefs
    use siftypes
    use sifetautils

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
        logical, dimension(Grid%shGrid%is:Grid%shGrid%ie,&
                           Grid%shGrid%js:Grid%shGrid%je) :: doBC

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

        
        do s=1,Grid%nSpc

            do i=Grid%shGrid%is,Grid%shGrid%ie
                do j=Grid%shGrid%js,Grid%shGrid%je
                    ! Skip if we should leave point alone
                    if(.not. doBC(i,j)) then
                        write(*,*)i,j
                        cycle
                    endif

                    vm = State%bvol(i,j)**(-2./3.)
                
                    kT = DP2kT(State%Davg(i,j,s), State%Pavg(i,j,s))  ! [keV]
                    call DkT2SpcEta(Model,Grid%spc(s), &
                        State%eta(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd),&
                        State%Davg(i,j,s), kT, vm)
                    !write(*,*)s,i,j,kT,maxval(State%eta(i,j,Grid%spc(s)%kStart:Grid%spc(s)%kEnd))
                enddo  ! j
            enddo  ! i

            ! Do first round of determining active shells for each k
            do k=1,Grid%spc(s)%N
                State%activeShells(:,k) = setLambdaActiveShells(Grid%shGrid, Grid%spc(s), State%bvol(i,j), &
                        State%eta(:,:,Grid%spc(s)%kStart:Grid%spc(s)%kStart), k)
            enddo  ! k

        enddo  ! s

    end subroutine applySifBCs


    function setLambdaActiveShells(shGrid, spc, bVol, etas, k, nSpacesO) result(activeShells)
        !! Determine which i shells should be evolved and which should not, for a given lambda channel
        type(ShellGrid_T), intent(in) :: shGrid
        type(SIFSpecies_T), intent(in) :: spc
        real(rp), intent(in) :: bVol
        integer, intent(in) :: k
            !! index for lambda dimension
        real(rp), dimension(:,:,:), intent(in) :: etas
        integer, optional, intent(in) :: nSpacesO
            !! Number of i spaces between last good value and active i for species

        integer :: i,j, iL, iU
        integer :: nSpaces
        logical, dimension(shGrid%is:shGrid%ie) :: activeShells

        if(present(nSpacesO)) then
            nSpaces = nSpacesO
        else
            nSpaces = nSpacesDef
        endif

        activeShells = .false.

        ! Determine shells that have enough stuff to be evolved
        do i=shGrid%is,shGrid%ie
            do j=shGrid%js,shGrid%je
                if (isWorthy(spc, bVol, etas(i,j,:), k)) then
                    ! While we're here, turn on this and adjacent i shells
                    iL = max(i-nSpaces, 1)
                    iU = min(i+nSpaces, shGrid%ie)
                    activeShells(iL:iU) = .true.
                    exit  ! No need to evaluate the rest of j
                endif
            enddo
        enddo

    end function setLambdaActiveShells


    function isWorthy(spc, bVol, etas, k, fracO) result(isW)
        !! Determine if lambda @ index k is contributing a certain percentage to the total pressure or density
        !! Evaluated at a single spatial point
        type(SIFSpecies_T), intent(in) :: spc
        real(rp), intent(in) :: bVol
        real(rp), dimension(spc%N), intent(in) :: etas
        integer, intent(in) :: k
        real(rp), optional, intent(in) :: fracO
            ! Fraction of total den/press that channel must contribute in order to be 

        real(rp) :: frac
        real(rp) :: alamc
        real(rp) :: kDen, kPress, spcDen, spcPress
        logical :: isW

        if(present(fracO)) then
            frac = fracO
        else
            frac = fracWorthyDef
        endif

        kDen = etak2Den(etas(k), bVol)
        spcDen = SpcEta2Den(spc, etas, bVol)

        alamc = 0.5*(spc%alami(k) + spc%alami(k+1))

        kPress = etak2Press(etas(k), alamc, bVol)
        spcPress = spcEta2Press(spc, etas, bVol)

        if (kDen/spcDen > frac .or. kPress/spcPress > frac) then
            isW = .true.
        else
            isW = .false.
        endif


    end function isWorthy

end module sifBCs