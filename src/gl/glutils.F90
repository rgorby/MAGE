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
        !> etimes is returned in [s], offset by Tstart_transient [s]
        function calcEmergenceTimes(Model, State, gT0) result(etimes)
            type(glModel_T), intent(in) :: Model
            type(glState_T), intent(in) :: State
            real(rp), intent(in) :: gT0
            real(rp), dimension(5) :: etimes,  tmp1, tmp2 ! to calculate emergence times
            integer :: i

            tmp1 = [4, 3, 2, 1, 0]*Model%r0/2 ! distance to travel in [solar radii]
            tmp2 = sqrt(eta0)*Model%velmult*(Model%frontheight + tmp1) ! velocity in [solar radii]/h
    
            do i = 1, 5
                if (tmp1(i) .ge. Model%frontheight) then ! collapsed, r le 0, or tmp1(i) ge frontheight
                    etimes(i) = 0.
                else
                    etimes(i) = ((tmp1(i)/tmp2(i))*3600. +  Model%Tstart_transient)/gT0
                end if
            end do

            etimes(5) = 0 ! in DOICMEM, frontheight is at 0.1AU always -- in the future for other tasks, we might need to change this
            if(Model%isLoud) write(*,"(1x,A40,2x,7F)") "Emergence times, lastP, eta0*velmult: ",  etimes, maxval(etimes), sqrt(eta0)*Model%velmult
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
                do i = State%is, State%ie
                    do j = State%js, State%je
                        do k = State%ks, State%ke                  
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

        !> Utility function to print model parameters
        !>
        !> 
        subroutine printModelParameters(Model, str)
            type(glModel_T), intent(in) :: Model
            character(*), optional, intent(in) :: str
            if(present(str)) write(*,"(1x,A24,A24)"), "Model Parameters ", str
            write(*,"(1x,A40)"), "----------------------------------------"
            write(*,"(1x,A24,2x,1F)") "frontheight: ", Model%frontheight
            write(*,"(1x,A24,2x,1F)") "legsang: ", Model%legsang
            write(*,"(1x,A24,2x,1F)") "topmorph: ", Model%topmorph
            write(*,"(1x,A24,2x,1F)") "latitude: ", Model%latitude
            write(*,"(1x,A24,2x,1F)") "longitude: ", Model%longitude
            write(*,"(1x,A24,2x,1F)") "sigma: ", Model%sigma
            write(*,"(1x,A24,2x,1F)") "cmer: ", Model%cmer
            write(*,"(1x,A24,2x,1F)") "Bmax: ", Model%Bmax
            write(*,"(1x,A24,2x,1F)") "cmeV: ", Model%cmeV
            write(*,"(1x,A24,2x,1F)") "vel_fh: ", Model%vel_fh
            write(*,"(1x,A24,2x,1F)") "apar: ", Model%apar
            write(*,"(1x,A24,2x,1F)") "r0: ", Model%r0
            write(*,"(1x,A24,2x,1F)") "x0: ", Model%x0
            write(*,"(1x,A24,2x,1F)") "ao: ", Model%ao
            write(*,"(1x,A24,2x,1F)") "bmagmax: ", Model%bmagmax
            write(*,"(1x,A24,2x,1F)") "alnot: ", Model%alnot
            write(*,"(1x,A24,2x,1F)") "alnotrbubuse: ", Model%alnotrbubuse
            write(*,"(1x,A24,2x,1F)") "alnotrbub: ", Model%alnotrbub
            write(*,"(1x,A24,2x,1F)") "muse: ", Model%muse
            write(*,"(1x,A24,2x,1F)") "xo: ", Model%xo
            write(*,"(1x,A24,2x,1F)") "rbub: ", Model%rbub
            write(*,"(1x,A24,2x,1F)") "anot: ", Model%anot
            write(*,"(1x,A24,2x,1F)") "eta: ", Model%eta
            write(*,"(1x,A24,2x,1F)") "phiss: ", Model%phiss
            write(*,"(1x,A24,2x,1F)") "aa: ", Model%aa
            write(*,"(1x,A24,2x,1F)") "bb: ", Model%bb
            write(*,"(1x,A24,2x,1F)") "cc: ", Model%cc
            write(*,"(1x,A24,2x,1F)") "dd: ", Model%dd
            write(*,"(1x,A24,2x,1F)") "ee: ", Model%ee
            write(*,"(1x,A24,2x,1F)") "ff: ", Model%ff
            write(*,"(1x,A24,2x,1F)") "velmult: ", Model%velmult
            write(*,"(1x,A24,2x,1F)") "alpha: ", Model%alpha
            write(*,"(1x,A24,2x,1F)") "s_eta0: ", Model%s_eta0
            write(*,"(1x,A24,2x,1F)") "s_eta: ", Model%s_eta
            write(*,"(1x,A40)"), "----------------------------------------"
        end subroutine printModelParameters

        !> Allocate and Intialize Solution Variables
        !>
        subroutine allocSolution(State, Solution)
            type(glState_T), intent(in) :: State
            type(glSolution_T), intent(inout) :: Solution

            if (.not. allocated(Solution%dens)) allocate(Solution%dens(State%is:State%ie, State%js:State%je, State%ks:State%ke))
            if (.not. allocated(Solution%pres)) allocate(Solution%pres(State%is:State%ie, State%js:State%je, State%ks:State%ke))
            if (.not. allocated(Solution%temp)) allocate(Solution%temp(State%is:State%ie, State%js:State%je, State%ks:State%ke))
            if (.not. allocated(Solution%b)) allocate(Solution%b(State%is:State%ie, State%js:State%je, State%ks:State%ke, NDIM))
            if (.not. allocated(Solution%v)) allocate(Solution%v(State%is:State%ie, State%js:State%je, State%ks:State%ke, NDIM))
            if (.not. allocated(Solution%j)) allocate(Solution%j(State%is:State%ie, State%js:State%je, State%ks:State%ke, NDIM))
            if (.not. allocated(Solution%inside_mask)) allocate(Solution%inside_mask(State%is:State%ie, State%js:State%je, State%ks:State%ke)) 

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

            if (.not. allocated(State%xyz)) allocate (State%xyz(State%is:State%ie,State%js:State%je,State%ks:State%ke,NDIM))
            if (.not. allocated(State%r)) allocate (State%r(State%is:State%ie))
            if (.not. allocated(State%rpb)) allocate (State%rpb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%thpb)) allocate (State%thpb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%phpb)) allocate (State%phpb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%rout)) allocate (State%rout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%thout)) allocate (State%thout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%phout)) allocate (State%phout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%rcap)) allocate (State%rcap(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%thcap)) allocate (State%thcap(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%phcap)) allocate (State%phcap(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%rsquig)) allocate (State%rsquig(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%rlam)) allocate (State%rlam(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%xtr)) allocate (State%xtr(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%ytr)) allocate (State%ytr(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%ztr)) allocate (State%ztr(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%F)) allocate (State%F(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%xtilde)) allocate (State%xtilde(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%ytilde)) allocate (State%ytilde(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%ztilde)) allocate (State%ztilde(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%rtilde)) allocate (State%rtilde(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%glpi)) allocate (State%glpi(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%bpresin)) allocate (State%bpresin(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%presin)) allocate (State%presin(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%presin_0)) allocate (State%presin_0(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%pbackin)) allocate (State%pbackin(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%dbackin)) allocate (State%dbackin(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%densin)) allocate (State%densin(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%densback)) allocate (State%densback(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%DensbackHEonly)) allocate (State%DensbackHEonly(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%ptot)) allocate (State%ptot(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%presback)) allocate (State%presback(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%bmag)) allocate (State%bmag(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%streamout)) allocate (State%streamout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%densout)) allocate (State%densout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%presout)) allocate (State%presout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%brlambout)) allocate (State%brlambout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%bthlambout)) allocate (State%bthlambout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%bphlambout)) allocate (State%bphlambout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%tderivout)) allocate (State%tderivout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%tderivR)) allocate (State%tderivR(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%tderivmu)) allocate (State%tderivmu(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%tderiv)) allocate (State%tderiv(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%jrlambout)) allocate (State%jrlambout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%jthlambout)) allocate (State%jthlambout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%jphlambout)) allocate (State%jphlambout(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%cavinside)) allocate (State%cavinside(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%blittlerlamb)) allocate (State%blittlerlamb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%blittlethlamb)) allocate (State%blittlethlamb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%blittlephlamb)) allocate (State%blittlephlamb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%jlittlerlamb)) allocate (State%jlittlerlamb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%jlittlethlamb)) allocate (State%jlittlethlamb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%jlittlephlamb)) allocate (State%jlittlephlamb(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%stream)) allocate (State%stream(State%js:State%je, State%ks:State%ke))
            if (.not. allocated(State%bstrengthphys)) allocate (State%bstrengthphys(State%js:State%je, State%ks:State%ke))

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