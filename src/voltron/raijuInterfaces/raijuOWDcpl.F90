!! Routines to handle packing coupling objects when one way driving with output files

module raijuowdcpl
    
    ! Base
    use ebtypes
    use shellinterp
    use remixReader

    ! Raiju
    use raijutypes
    use raijucpltypes
    use raijuGrids
    
    ! Voltron
    use imagtubes

    implicit none


    contains

    !> This function takes updated model states and does the operations
    !> necessary to update the cplBase%fromV object
    subroutine packFromV(fromV, vApp, rmReader, raiApp)
        type(raiju_fromV_T), intent(inout) :: fromV
        type(voltApp_T), intent(inout) :: vApp
        type(rmReader_T) :: rmReader
        type(raijuApp_T) , intent(in) :: raiApp

        integer :: i,j
        real(rp), dimension(:,:), allocatable :: tmpPot

        ! Update coupling time
        fromV%tLastUpdate = raiApp%State%t

        ! Using chimp, populate imagTubes
        call genImagTubes(fromV, vApp, raiApp)

        ! Set potential
        call InterpShellVar_TSC_SG(rmReader%shGr, rmReader%nsPot(1), fromV%shGr, fromV%pot)

    end subroutine packFromV


    subroutine genImagTubes(fromV, vApp, sApp)
        type(raiju_fromV_T), intent(inout) :: fromV
        type(voltApp_T), intent(in   ) :: vApp
        type(raijuApp_T ), intent(in) :: sApp

        integer :: i,j
        real(rp) :: seedR
        type(magLine_T) :: magLine
        ! Get field line info and potential from voltron
        ! And put the data into RAIJU's fromV coupling object

        associate(sh=>sApp%Grid%shGrid , &
            planet=>sApp%Model%planet, &
            ebApp =>vApp%ebTrcApp)

            seedR =  planet%ri_m/planet%rp_m + TINY
            ! Do field line tracing, populate fromV%ijTubes
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%isg,sh%ieg+1
                do j=sh%jsg,sh%jeg+1
                    call CleanLine(fromV%magLines(i,j))

                    call MHDTube(ebApp, planet,   & !ebTrcApp, planet
                        sh%th(i), sh%ph(j), seedR, &  ! colat, lon, r
                        fromV%ijTubes(i,j), fromV%magLines(i,j), &  ! IMAGTube_T, magLine_T
                        doShiftO=.true.)
                enddo
            enddo
        end associate

    end subroutine genImagTubes


end module raijuowdcpl