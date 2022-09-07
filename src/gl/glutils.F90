module glutils
    use gltypes
    use math
    implicit none

    contains

        !> Ensures close to zero numbers are large enough 
        !> from trigonometric functions
        !> 
        subroutine zero2tiny2d(arr)
            real(rp), dimension(:, :), intent(inout) :: arr
            real(rp), dimension(:, :), allocatable :: tmparr
            integer, dimension(2) :: dims
            !debug
            !write(555,*)"[EP] in zero2tiny2d: printing dim1 and dim2"
            !write(555,*) dim1, dim2
            !write(555,*) tiny
            dims = shape(arr)
            allocate (tmparr(dims(1), dims(2)))
            tmparr = 1.0

            where (arr .lt. 0) tmparr = -1.0
            where (abs(arr) .lt. BIGTINY) arr = BIGTINY*tmparr

            deallocate (tmparr)
        end subroutine zero2tiny2d

        !> Calculate coordinate transformation for "Cap"
        !> model coordinates
        !> 
        subroutine calcCoords(Model, State)
            type(glModel_T), intent(inout) :: Model
            type(glState_T), intent(inout) :: State
            
            !;  the self-similar transformation:
            State%rsquig = State%rpb/Model%phiss

            ! ;  the radial coordinate is now transformed introducing expansion of
            ! ;  the field
            State%rlam = State%rsquig + Model%apar

            ! ;  The force balance equations will actually
            ! ;  be solved in yet another system, this one with a shifted center,
            ! ;  This shifted symmetry introduces asymmetry.
            ! ;  We must first move to Cartesian coordinates.
            ! ;
            ! ;  note I have changed the coordinates to a right-hand system
            ! ;  (finally)
            ! ;  the shift xo is along the x axis, so now one needs to set 
            ! ;  central meridian to -90 to see the bubble at the West limb

            State%xtr = State%rlam*sin(State%thpb)*cos(State%phpb) - Model%xo
            State%ytr = State%rlam*sin(State%thpb)*sin(State%phpb)
            State%ztr = State%rlam*cos(State%thpb)
            
            ! now transform to spherical coords
            ! rcap used for inside bubble stream function 
            State%rcap = sqrt(State%xtr*State%xtr + State%ytr*State%ytr + State%ztr*State%ztr)
            
            ! ; we want to rotate the coordinates in the bubble 
            ! ;  about the x-axis by an angle sigma (parameter)
            ! ;
            ! ; note the counterclockwise rotation
            ! ; this is necessary to move into the bubble coordinates
            ! ; which are rotated sigma clockwise from the physical coordinates
            ! ;
            State%ytilde = State%ytr*cos(Model%sigma) - State%ztr*sin(Model%sigma)
            State%ztilde = State%ztr*cos(Model%sigma) + State%ytr*sin(Model%sigma)

            State%rat = (State%ztilde/State%Rcap)
            where (State%rat .gt. 1.) State%rat = 1.
            where (State%rat .lt. -1.) State%rat = -1.

            State%Thcap = acos(State%rat)
            State%Phcap = atan2(State%ytilde, State%xtr)

            where (State%Thcap /= State%Thcap) State%Thcap = 0.
            where (State%Phcap /= State%Phcap) State%Phcap = 0.
            where (State%rcap /= State%rcap) State%rcap = 0.

        end subroutine calcCoords

        !> emergence times:
        !> velocity at each point depends on the point's position at t=0: 
        !> v=sqrt(eta)*Rsun_km*velmult*r(t=0) km/s
        !> (here sqrt(eta)*Rsun_km~258.55)
        !> critical points are
        !>  x0-apar+[-2, -1, 0, 1, 2]*rbub/2=
        !>   =x0-apar+2*rbub/2+[-4, -3, -2, -1, 0]*rbub/2=
        !>   =frontheight-[4, 3, 2, 1, 0]*rbub/2
        !>
        !> they need to travel the following distance to emerge to frontheight:
        !> frontheight-(frontheight-[4, 3, 2, 1, 0]*rbub/2)=
        !>  =[4, 3, 2, 1, 0]*rbub/2
        !>
        !>  convention: emergence_time=0 means a collapsed critical point
        function calcEmergenceTimes(Model, State) result(etimes)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(in) :: State

            real(rp), dimension(5) :: etimes,  tmp1, tmp2 ! to calculate emergence times
            integer :: i

            tmp1 = [4, 3, 2, 1, 0]*Model%r0/2 ! distance to travel
            tmp2 = sqrt(eta0)*Model%velmult*(Model%frontheight + tmp1)*3600 ! velocity in [solar radii]/hr
    
            do i = 1, 5
                if (tmp1(i) .ge. Model%frontheight) then ! collapsed, r le 0, or tmp1(i) ge frontheight
                    etimes(i) = 0
                else
                    etimes(i) = tmp1(i)/tmp2(i)
                end if
            end do

            etimes(5) = 0 ! in DOICMEM, frontheight is at 0.1AU always -- in the future for other tasks, we might need to change this
        end function

        !> 
        !> xyz must have dimensions of xyz(Model%Ni, Model%Nj, Model%Nk, NDIM)
        function solutionSphereToCartesian(xyz, Model, State, SolutionSphere) result(SolutionCartesian)
            type(glState_T) :: State
            type(glModel_T) :: Model
            real(rp), dimension(:,:,:,:), intent(in) :: xyz
            type(glSolution_T), intent(in) :: SolutionSphere
            type(glSolution_T) :: SolutionCartesian
            integer :: i,j,k

            call allocSolution(State, SolutionCartesian)

            if (SolutionSphere%CoordSystem .eq. SPHERICAL) then
                do i = 1, State%Ni
                    do j = 1, State%Nj
                        do k = 1, State%Nk                        
                            SolutionCartesian%b(i,j,k,:) = rtp2xyz(xyz(i,j,k,:), SolutionSphere%b(i,j,k,:))
                            SolutionCartesian%v(i,j,k,:) = rtp2xyz(xyz(i,j,k,:), SolutionSphere%v(i,j,k,:))
                            SolutionCartesian%j(i,j,k,:) = rtp2xyz(xyz(i,j,k,:), SolutionSphere%j(i,j,k,:))
                        end do
                    end do
                end do
                SolutionCartesian%pres = SolutionSphere%pres
                SolutionCartesian%dens = SolutionSphere%dens
                SolutionCartesian%pres = SolutionSphere%pres
                SolutionCartesian%inside_mask = SolutionSphere%inside_mask
                SolutionCartesian%CoordSystem = CARTESIAN
            else 
                if (Model%isLoud) write(*,*) "Solution is not curently in Spherical Coordinates, please check your coodinate system choice"
            end if
        end function solutionSphereToCartesian

        !> Allocate and Intialize Solution Variables
        !>
        subroutine allocSolution(State, Solution)
            type(glState_T), intent(in) :: State
            type(glSolution_T), intent(inout) :: Solution

            if (.not. allocated(Solution%dens)) allocate(Solution%dens(State%Ni, State%Nj, State%Nk))
            if (.not. allocated(Solution%pres)) allocate(Solution%pres(State%Ni,State%Nj, State%Nk))
            if (.not. allocated(Solution%temp)) allocate(Solution%temp(State%Ni, State%Nj, State%Nk))
            if (.not. allocated(Solution%b)) allocate(Solution%b(State%Ni, State%Nj, State%Nk, NDIM))
            if (.not. allocated(Solution%v)) allocate(Solution%v(State%Ni, State%Nj, State%Nk, NDIM))
            if (.not. allocated(Solution%j)) allocate(Solution%j(State%Ni, State%Nj, State%Nk, NDIM))
            if (.not. allocated(Solution%inside_mask)) allocate(Solution%inside_mask(State%Ni, State%Nj, State%Nk)) 

            Solution%dens = 0.
            Solution%pres = 0.
            Solution%temp = 0.
            Solution%b = 0.
            Solution%v = 0.
            Solution%j = 0.
            Solution%inside_mask = 0
        end subroutine allocSolution

        !> Allocate and Initalize State Variables
        !> 
        subroutine allocState(State)
            type(glState_T),  intent(inout) :: State

            if (.not. allocated(State%xyz)) allocate (State%xyz(State%Ni,State%Nj,State%Nk,NDIM))
            if (.not. allocated(State%r)) allocate (State%r(State%Ni))
            if (.not. allocated(State%rpb)) allocate (State%rpb(State%Nj, State%Nk))
            if (.not. allocated(State%thpb)) allocate (State%thpb(State%Nj, State%Nk))
            if (.not. allocated(State%phpb)) allocate (State%phpb(State%Nj, State%Nk))
            if (.not. allocated(State%rout)) allocate (State%rout(State%Nj, State%Nk))
            if (.not. allocated(State%thout)) allocate (State%thout(State%Nj, State%Nk))
            if (.not. allocated(State%phout)) allocate (State%phout(State%Nj, State%Nk))
            if (.not. allocated(State%rcap)) allocate (State%rcap(State%Nj, State%Nk))
            if (.not. allocated(State%thcap)) allocate (State%thcap(State%Nj, State%Nk))
            if (.not. allocated(State%phcap)) allocate (State%phcap(State%Nj, State%Nk))
            if (.not. allocated(State%rsquig)) allocate (State%rsquig(State%Nj, State%Nk))
            if (.not. allocated(State%rlam)) allocate (State%rlam(State%Nj, State%Nk))
            if (.not. allocated(State%xtr)) allocate (State%xtr(State%Nj, State%Nk))
            if (.not. allocated(State%ytr)) allocate (State%ytr(State%Nj, State%Nk))
            if (.not. allocated(State%ztr)) allocate (State%ztr(State%Nj, State%Nk))
            if (.not. allocated(State%F)) allocate (State%F(State%Nj, State%Nk))
            if (.not. allocated(State%xtilde)) allocate (State%xtilde(State%Nj, State%Nk))
            if (.not. allocated(State%ytilde)) allocate (State%ytilde(State%Nj, State%Nk))
            if (.not. allocated(State%ztilde)) allocate (State%ztilde(State%Nj, State%Nk))
            if (.not. allocated(State%rtilde)) allocate (State%rtilde(State%Nj, State%Nk))
            if (.not. allocated(State%glpi)) allocate (State%glpi(State%Nj, State%Nk))
            if (.not. allocated(State%bpresin)) allocate (State%bpresin(State%Nj, State%Nk))
            if (.not. allocated(State%presin)) allocate (State%presin(State%Nj, State%Nk))
            if (.not. allocated(State%presin_0)) allocate (State%presin_0(State%Nj, State%Nk))
            if (.not. allocated(State%pbackin)) allocate (State%pbackin(State%Nj, State%Nk))
            if (.not. allocated(State%dbackin)) allocate (State%dbackin(State%Nj, State%Nk))
            if (.not. allocated(State%densin)) allocate (State%densin(State%Nj, State%Nk))
            if (.not. allocated(State%densback)) allocate (State%densback(State%Nj, State%Nk))
            if (.not. allocated(State%DensbackHEonly)) allocate (State%DensbackHEonly(State%Nj, State%Nk))
            if (.not. allocated(State%ptot)) allocate (State%ptot(State%Nj, State%Nk))
            if (.not. allocated(State%presback)) allocate (State%presback(State%Nj, State%Nk))
            if (.not. allocated(State%bmag)) allocate (State%bmag(State%Nj, State%Nk))
            if (.not. allocated(State%streamout)) allocate (State%streamout(State%Nj, State%Nk))
            if (.not. allocated(State%densout)) allocate (State%densout(State%Nj, State%Nk))
            if (.not. allocated(State%presout)) allocate (State%presout(State%Nj, State%Nk))
            if (.not. allocated(State%brlambout)) allocate (State%brlambout(State%Nj, State%Nk))
            if (.not. allocated(State%bthlambout)) allocate (State%bthlambout(State%Nj, State%Nk))
            if (.not. allocated(State%bphlambout)) allocate (State%bphlambout(State%Nj, State%Nk))
            if (.not. allocated(State%tderivout)) allocate (State%tderivout(State%Nj, State%Nk))
            if (.not. allocated(State%tderivR)) allocate (State%tderivR(State%Nj, State%Nk))
            if (.not. allocated(State%tderivmu)) allocate (State%tderivmu(State%Nj, State%Nk))
            if (.not. allocated(State%tderiv)) allocate (State%tderiv(State%Nj, State%Nk))
            if (.not. allocated(State%jrlambout)) allocate (State%jrlambout(State%Nj, State%Nk))
            if (.not. allocated(State%jthlambout)) allocate (State%jthlambout(State%Nj, State%Nk))
            if (.not. allocated(State%jphlambout)) allocate (State%jphlambout(State%Nj, State%Nk))
            if (.not. allocated(State%cavinside)) allocate (State%cavinside(State%Nj, State%Nk))
            if (.not. allocated(State%blittlerlamb)) allocate (State%blittlerlamb(State%Nj, State%Nk))
            if (.not. allocated(State%blittlethlamb)) allocate (State%blittlethlamb(State%Nj, State%Nk))
            if (.not. allocated(State%blittlephlamb)) allocate (State%blittlephlamb(State%Nj, State%Nk))
            if (.not. allocated(State%jlittlerlamb)) allocate (State%jlittlerlamb(State%Nj, State%Nk))
            if (.not. allocated(State%jlittlethlamb)) allocate (State%jlittlethlamb(State%Nj, State%Nk))
            if (.not. allocated(State%jlittlephlamb)) allocate (State%jlittlephlamb(State%Nj, State%Nk))
            if (.not. allocated(State%stream)) allocate (State%stream(State%Nj, State%Nk))
            if (.not. allocated(State%bstrengthphys)) allocate (State%bstrengthphys(State%Nj, State%Nk))

            State%xyz = 0.
            State%r = 0.            
            State%rpb = 0.
            State%thpb = 0.
            State%phpb = 0.
            State%rout = 0.
            State%thout = 0.
            State%phout = 0.
            State%rcap = 0.
            State%thcap = 0.
            State%phcap = 0.
            State%rsquig = 0.
            State%rlam = 0.
            State%xtr = 0.
            State%ytr = 0.
            State%ztr = 0.
            State%F = 0.
            State%xtilde = 0.
            State%ytilde = 0.
            State%ztilde = 0.
            State%rtilde = 0.
            State%glpi = 0.
            State%densin = 0.
            State%densback = 0.
            State%DensbackHEonly = 0.
            State%densout = 0.
            State%bpresin = 0.
            State%presin = 0.
            State%presin_0 = 0.
            State%pbackin = 0.
            State%dbackin = 0.
            State%ptot = 0.
            State%presback  = 0.
            State%bmag = 0.
            State%streamout = 0.
            State%presout = 0.
            State%brlambout = 0.
            State%bthlambout = 0.
            State%bphlambout = 0.
            State%tderivout = 0.
            State%tderivR = 0.
            State%tderivmu = 0.
            State%tderiv = 0.
            State%jrlambout = 0.
            State%jthlambout = 0.
            State%jphlambout = 0.
            State%cavinside = 0.
            State%blittlerlamb = 0.
            State%blittlethlamb = 0.
            State%blittlephlamb = 0.
            State%jlittlerlamb = 0.
            State%jlittlethlamb = 0.
            State%jlittlephlamb = 0.
            State%stream = 0.
            State%bstrengthphys = 0.
        end subroutine allocState

        !> Deallocate State Variables
        !> 
        subroutine deallocState(State)
            type(glState_T),  intent(inout) :: State
            if (allocated(State%r)) deallocate (State%r)
            if (allocated(State%rpb)) deallocate (State%rpb)
            if (allocated(State%thpb)) deallocate (State%thpb)
            if (allocated(State%phpb)) deallocate (State%phpb)
            if (allocated(State%rout)) deallocate (State%rout)
            if (allocated(State%thout)) deallocate (State%thout)
            if (allocated(State%phout)) deallocate (State%phout)
            if (allocated(State%rcap)) deallocate (State%rcap)
            if (allocated(State%Thcap)) deallocate (State%Thcap)
            if (allocated(State%Phcap)) deallocate (State%Phcap)
            if (allocated(State%rsquig)) deallocate (State%rsquig)
            if (allocated(State%rlam)) deallocate (State%rlam)
            if (allocated(State%xtr)) deallocate (State%xtr)
            if (allocated(State%ytr)) deallocate (State%ytr)
            if (allocated(State%ztr)) deallocate (State%ztr)
            if (allocated(State%F)) deallocate (State%F)
            if (allocated(State%xtilde)) deallocate (State%xtilde)
            if (allocated(State%ytilde)) deallocate (State%ytilde)
            if (allocated(State%ztilde)) deallocate (State%ztilde)
            if (allocated(State%rtilde)) deallocate (State%rtilde)
            if (allocated(State%thtilde)) deallocate (State%thtilde)
            if (allocated(State%phtilde)) deallocate (State%phtilde)
        end subroutine deallocState
end module glutils