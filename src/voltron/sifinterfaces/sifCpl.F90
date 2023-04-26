module sifCpl
    !! Contains functions used to pass data between voltron and SIF

    use imagtubes

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
        type(fLine_T) :: fLine
        ! Get field line info and potential from voltron
        ! And put the data into SIF

        associate(fromV=>cplBase%fromV, &
            sh    =>sApp%Grid%shGrid , &
            planet=>sApp%Model%planet, &
            ebApp =>vApp%ebTrcApp)
            ! Do field line tracing, poopulate fromV%ijTubes
            do i=sh%is,sh%ie+1
                do j=sh%js,sh%je+1
                    call MHDTube(ebApp, planet,   & !ebTrcApp, planet
                        sh%th(i), sh%ph(j), planet%ri_m/planet%rp_m, &  ! colat, lon, r
                        fromV%ijTubes(i,j), fLine, &  ! IMAGTube_T, fLine_T
                        doShiftO=.true. &
                        )
                enddo
            enddo
            call imagTubes2SIF(fromV%ijTubes, sApp%Model, sApp%Grid, sApp%State)

        end associate
    end subroutine sifCpl_Volt2SIF

    subroutine sifCpl_SIF2Volt(cplBase, vApp, sApp)
        type(sif_cplBase_T), intent(inout) :: cplBase
        type(voltApp_T), intent(inout) :: vApp
        type(sifApp_T) , intent(in   ) :: sApp
    end subroutine sifCpl_SIF2Volt


! Helpers

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
            do i=sh%is,sh%ie+1
                do j=sh%js,sh%je+1
                    State%xyzMin(i,j,:)  = ijTubes(i,j)%X_bmin / Model%planet%rp_m  ! xyzMn in Rp
                enddo
            enddo
            ! Cell-centered quantities
            do i=sh%is,sh%ie
                do j=sh%js,sh%je
                    State%Pavg(i,j,1)    = Vcc2D(ijTubes(i:i+1,j:j+1)%Pave) * 1.0e+9  ! Pa -> nPa
                    State%Davg(i,j,1)    = Vcc2D(ijTubes(i:i+1,j:j+1)%Nave) * 1.0e-6  ! #/m^3 -> #/cc
                    State%Bmin(i,j,ZDIR) = Vcc2D(ijTubes(i:i+1,j:j+1)%bmin) * 1.0e+9  ! Tesla -> nT
                    !State%topo(i,j)      = Vcc2D(ijTubes(i:i+1,j:j+1)%topo)
                    !State%active(i,j)    = Vcc2D(ijTubes(i:i+1,j:j+1)%topo) !! TODO: Handle active
                    State%thc(i,j)  = Vcc2D(PI/2-ijTubes(i:i+1,j:j+1)%latc)
                    State%phc(i,j)       = Vcc2D(ijTubes(i:i+1,j:j+1)%lonc)
                    State%bvol(i,j)      = Vcc2D(ijTubes(i:i+1,j:j+1)%Vol)  ! Rp/T
                enddo
            enddo
        end associate

        contains

        function Vcc2D(v) result (vcc)
            !! 2D cell center of a value
            real(rp), dimension(2,2), intent(in) :: v

            real(rp) :: vcc

            vcc = 0.25*(v(1,1) + v(2,1) + v(1,2) + v(2,2))
        end function Vcc2D
        
    end subroutine imagTubes2SIF


end module sifCpl