!! Routines to handle packing coupling objects when one way driving with output files

module sifowdcpl

    use imagtubes

    use siftypes
    use sifcpltypes

    !!Temporary. Eventually we will just use shellGrid stuff
    use calcdbtypes
    use calcdbio
    use mixinterp
    use mixgeom

    implicit none

    type(Map_T) :: m2sMap

    contains

    

    !> This function takes updated model states and does the operations
    !> necessary to update the cplBase%fromV object
    subroutine packFromV(fromV, vApp, rmState, sApp)
        type(sif_fromV_T), intent(inout) :: fromV
        type(voltApp_T), intent(inout) :: vApp
        type(rmState_T) :: rmState
        type(sifApp_T) , intent(in) :: sApp

        ! Update coupling time
        fromV%tLastUpdate = sApp%State%t

        ! Using chimp, populate imagTubes
        call genImagTubes(fromV, vApp, sApp)

        ! Set potential
        call mix_map_grids(m2sMap, rmState%nPot, fromV%pot)

    end subroutine packFromV


    subroutine genImagTubes(fromV, vApp, sApp)
        type(sif_fromV_T), intent(inout) :: fromV
        type(voltApp_T), intent(in   ) :: vApp
        type(sifApp_T ), intent(in) :: sApp

        integer :: i,j
        real(rp) :: seedR
        type(fLine_T) :: fLine
        ! Get field line info and potential from voltron
        ! And put the data into sif's fromV coupling object

        associate(sh=>sApp%Grid%shGrid , &
            planet=>sApp%Model%planet, &
            ebApp =>vApp%ebTrcApp)

            seedR =  planet%ri_m/planet%rp_m + TINY
            ! Do field line tracing, populate fromV%ijTubes
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%is,sh%ie+1
                do j=sh%js,sh%je+1
                    call CleanStream(fromV%fLines(i,j))

                    call MHDTube(ebApp, planet,   & !ebTrcApp, planet
                        sh%th(i), sh%ph(j), seedR, &  ! colat, lon, r
                        fromV%ijTubes(i,j), fromV%fLines(i,j), &  ! IMAGTube_T, fLine_T
                        doShiftO=.true.,doShueO=.false. &
                        )
                enddo
            enddo
        end associate

    end subroutine genImagTubes



!----
!  !! TEMPORARY, replace with shellGrid interp stuff later
!------

    subroutine InitMixMap(shGrid, mixState)
        !! Adapted from rcmXimag.F90:rcmGrid
        !! Take a shell grid and generate map from remix grid to shellGrid
        type(ShellGrid_T), intent(in) :: shGrid
        type(rmState_T), intent(in) :: mixState

        integer :: i
        type(mixGrid_T) :: mixGrid, mixedGrid  ! Mix grid from file, shGrid converted to mixGrid
        real(rp), dimension(:,:), allocatable :: colat2D,lon2D

        allocate(colat2D(shGrid%Nt,shGrid%Np))  ! +1 because we're doing corners
        allocate(lon2D  (shGrid%Nt,shGrid%Np))

        do i=1,shGrid%Np
            colat2D(:,i) = shGrid%thc(shGrid%is:shGrid%ie)
        enddo

        do i=1,shGrid%Nt
            lon2D(i,:) = shGrid%phc(shGrid%js:shGrid%je)
        enddo

        call init_grid_fromXY(mixGrid, mixState%XY(:,:,XDIR),mixState%XY(:,:,YDIR),.false.,.true.)
        call init_grid_fromTP(mixedGrid, colat2D, lon2D,.false.,.true.)
        call mix_set_map(mixGrid, mixedGrid, m2sMap)
        !write(*,*) map%M(:,:,1)

    end subroutine InitMixMap


end module sifowdcpl