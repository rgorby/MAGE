module voltappHelper

    use kdefs
    use voltTypes
    use imagtubes
    use kai2geo
    use mixinterfaceutils

    implicit none

    contains

    subroutine initVoltState(vApp)
        type(voltApp_T), intent(inout) :: vApp

        associate(sh=>vApp%shGrid, State=>vApp%State)

            ! Allocations
            allocate(State%ijTubes(sh%is:sh%ie+1, sh%js:sh%je+1))
            call init_TubeShell(sh, State%tubeShell)
            call initShellVar(sh, SHGR_CORNER, State%potential_total)
            call initShellVar(sh, SHGR_CORNER, State%potential_corot)

            call initShellVar(sh, SHGR_CC, State%bIonoMag)
            call initShellVar(sh, SHGR_CC, State%bIonoRad)

            ! Default values
            call initMagVars(sh, vApp%planet, State)
        end associate 


        contains

        subroutine initMagVars(sh, planet, State)
            !! Default routine to init magnetic field-related variables on voltron grid. 
            !! TODO: Will make this more formal / someone else's job once we have different mag field models like IGRF
            type(ShellGrid_T), intent(in) :: sh
            type(planet_T   ), intent(in) :: planet
            type(voltState_t), intent(inout) :: State

            integer :: j
            real(rp), dimension(sh%isg:sh%ieg) :: cosThc
                !! cos(th) at cell centers
            real(rp), dimension(sh%isg:sh%ieg) :: BMagTh
                !! [nT] mag field along theta direction

            cosThc = cos(sh%thc)
            do j=sh%jsg,sh%jeg
                State%bIonoMag%data(:,j) = planet%magMoment*G2nT &
                                         /sh%radius**3.0 &
                                         * sqrt(1.0+3.0*cosThc**2.0)  ! [nT]
                State%bIonoRad%data(:,j) = planet%magMoment*G2nT &
                                         /sh%radius**3.0 &
                                         * 2*cosThc  ! [nT]
            enddo


        end subroutine initMagVars
        
    end subroutine initVoltState

    subroutine updateVoltPotential(vApp)
        class(voltApp_T), intent(inout) :: vApp

        type(ShellGridVar_T) :: iono_var
            !! We hand this to mixVarToVoltron to give us back the remix potential on the voltron ShellGrid

        ! Get ionospheric potential from remix, add to coration to get total potential in SM coordinates
        call initShellVar(vApp%shGrid, SHGR_CORNER, iono_var)
        call calcCorotPotential(vApp%planet, vApp%shGrid, vApp%State%potential_corot,doGeoCorotO=vApp%doGeoCorot)
        call mixVarToVoltron(vApp%remixApp, POT, vApp%shGrid, iono_var)
        ! Combine iono and corotation potential to form total potential in SM frame
        vApp%State%potential_total%data = iono_var%data + vApp%State%potential_corot%data
        vApp%State%potential_total%mask = .true.

    end subroutine updateVoltPotential

    subroutine calcCorotPotential(planet, sh, pCorot,doGeoCorotO)
        ! Calculate corotation potential to store on voltron ShellGrid
        type(planet_T), intent(in) :: planet
        type(ShellGrid_T), intent(in) :: sh
        type(ShellGridVar_T), intent(inout) :: pCorot
            !! [kV] corotation potential we calculate
        logical, intent(in), optional :: doGeoCorotO
        
        integer :: i,j

        pCorot%data = 0.0
        pCorot%mask = .true.

        if (pCorot%loc .ne. SHGR_CORNER) then
            write(*,*)"Error in voltappHelper.F90's calcCorotPotential:"
            write(*,*)"Only handling corner-located corot var for now"
            stop
        endif

        if (present(doGeoCorotO)) then
            if (doGeoCorotO) then
                do j=sh%jsg,sh%jeg+1
                    do i=sh%isg,sh%ieg+1
                        ! geopack corotation returns in kV
                        call geocorotation_from_SMTP(sh%th(i), sh%ph(j), pCorot%data(i,j))
                    enddo
                enddo

                return

            endif
        endif

        ! IfdoGeoCorotO was not provided, or it was false, we default to corotation on aligned dipole and rotational axis
        do j=sh%jsg,sh%jeg+1
            pCorot%data(:,j) = -planet%psiCorot*(planet%rp_m/planet%ri_m)*sin(sh%th)**2  ! [kV]
        enddo

    end subroutine calcCorotPotential

end module voltappHelper