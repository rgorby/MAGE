module sifCpl
    !! Contains functions used to pass data between voltron and SIF

    use math
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
        real(rp) :: seedR
        type(fLine_T) :: fLine
        ! Get field line info and potential from voltron
        ! And put the data into SIF

        associate(fromV=>cplBase%fromV, &
            sh    =>sApp%Grid%shGrid , &
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
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%is,sh%ie+1
                do j=sh%js,sh%je+1
                    State%xyzMin(i,j,:)  = ijTubes(i,j)%X_bmin / Model%planet%rp_m  ! xyzMin in Rp
                    State%topo(i,j)      = ijTubes(i,j)%topo
                enddo
            enddo
            ! Cell-centered quantities
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%is,sh%ie
                do j=sh%js,sh%je
                    State%Pavg(i,j,1)    = toCenter2D(ijTubes(i:i+1,j:j+1)%Pave) * 1.0e+9  ! Pa -> nPa
                    State%Davg(i,j,1)    = toCenter2D(ijTubes(i:i+1,j:j+1)%Nave) * 1.0e-6  ! #/m^3 -> #/cc
                    State%Bmin(i,j,ZDIR) = toCenter2D(ijTubes(i:i+1,j:j+1)%bmin) * 1.0e+9  ! Tesla -> nT
                    State%thc(i,j)  = toCenter2D(PI/2-ijTubes(i:i+1,j:j+1)%latc)
                    State%phc(i,j)       = toCenter2D(ijTubes(i:i+1,j:j+1)%lonc)
                    State%bvol(i,j)      = toCenter2D(ijTubes(i:i+1,j:j+1)%Vol) * 1.0e-9  ! Rp/T -> Rp/nT
                enddo
            enddo
        end associate

        ! TODO: Save calculating active domain for last, once we have all other info
        !State%active(i,j)   = ??

        
    end subroutine imagTubes2SIF


end module sifCpl