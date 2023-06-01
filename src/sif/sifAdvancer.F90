module sifadvancer

    use planethelper

    use sifdefs
    use siftypes
    use sifetautils

    implicit none

    contains


!------
! Active I Shell calculations
!------
    function setLambdaActiveShells(shGrid, spc, bVol, etas, k, nSpacesO, worthyFracO) result(activeShells)
        !! Determine which i shells should be evolved and which should not, for a given lambda channel
        type(ShellGrid_T), intent(in) :: shGrid
        type(SIFSpecies_T), intent(in) :: spc
        real(rp), dimension(shGrid%isg:shGrid%ieg,shGrid%jsg:shGrid%jeg), intent(in) :: bVol
        integer, intent(in) :: k
            !! index for lambda dimension
        real(rp), dimension(shGrid%isg:shGrid%ieg,shGrid%jsg:shGrid%jeg,spc%kStart:spc%kEnd), intent(in) :: etas
        integer, optional, intent(in) :: nSpacesO
            !! Number of i spaces between last good value and active i for species
        real(rp), optional, intent(in) :: worthyFracO

        integer :: i,j, iL, iU
        integer :: nSpaces
        real(rp) :: worthyFrac
        logical, dimension(shGrid%isg:shGrid%ieg) :: activeShells

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

        activeShells = .false.

        ! Determine shells with cells that have enough stuff to be evolved
        !$OMP PARALLEL DO default(shared) collapse(1) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,iL,iU)
        do i=shGrid%isg,shGrid%ieg
            do j=shGrid%jsg,shGrid%jeg
                if (isWorthy(spc, bVol(i,j), etas(i,j,:), k, worthyFrac)) then
                    ! While we're here, turn on this and adjacent i shells
                    iL = max(i-nSpaces, shGrid%isg)
                    iU = min(i+nSpaces, shGrid%ieg)
                    activeShells(iL:iU) = .true.
                    !exit  ! No need to evaluate the rest of j
                endif
            enddo
        enddo

    end function setLambdaActiveShells


    function isWorthy(spc, bVol, etas, k, fracO) result(isW)
        !! Determine if lambda @ index k is contributing a certain percentage to the total pressure or density
        !! Evaluated at a single spatial point
        type(SIFSpecies_T), intent(in) :: spc
        real(rp), intent(in) :: bVol
        real(rp), dimension(spc%kStart:spc%kEnd), intent(in) :: etas
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

        alamc = 0.5*abs(spc%alami(k) + spc%alami(k+1))

        kPress = etak2Press(etas(k), alamc, bVol)
        spcPress = spcEta2Press(spc, etas, bVol)

        if (kDen/spcDen > frac .or. kPress/spcPress > frac) then
            isW = .true.
        else
            isW = .false.
        endif
    end function isWorthy

!------
! Cell Velocity calculations
!------

! Cell velocity will not change during sub-stepping, so can calculate for everyone at once

    subroutine calcVEffective(Grid, State, vEff, vGCO, vExBO, vCorotO)
        ! Calculates vEffective as the sum of vExB, vCorot, and vGC
        ! Returns vEffectie as vEff, must be provided
        ! Optionally returns individual velocities if requested
        type(sifGrid_T) , intent(in) :: Grid
        type(sifState_T), intent(in) :: State
        real(rp), dimension(:,:,:), allocatable, intent(out) :: vEff
        real(rp), dimension(:,:,:), allocatable, optional, intent(out) :: vGCO
        real(rp), dimension(:,:)  , allocatable, optional, intent(out) :: vExBO
        real(rp), dimension(:,:)  , allocatable, optional, intent(out) :: vCorotO
        

        real(rp), dimension(:,:,:), allocatable :: vGC
        real(rp), dimension(:,:), allocatable :: vExB, vCorot
        

        associate(sh=>Grid%shGrid, spc=>Grid%spc)

            ! Allocate everyone
            allocate( vEff   (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( vGC    (sh%isg:sh%ieg, sh%jsg:sh%jeg, Grid%Nk) )
            allocate( vExB   (sh%isg:sh%ieg, sh%jsg:sh%jeg) )
            allocate( vCorot (sh%isg:sh%ieg, sh%jsg:sh%jeg) )


        end associate

    end subroutine calcVEffective

end module sifadvancer