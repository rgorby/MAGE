!K: RCM Clawpack code, adapted by Kareem from what Frank sent me
MODULE rcmclaw
    use rcm_precision
    implicit none

    PUBLIC claw2ez95

    INCLUDE 'rcm_include.h'

    integer, parameter, private :: maux   =     3,   &  ! aux array has maux components
                                   mx     = isize,   &  !cells in x direction
                                   my     = jsize-3, &  !cells in y direction
                                   max_allowed_steps = 200  ! max number of time steps we allow per call to claw2

    real(rprec),  parameter, private :: max_allowed_cfl   = 1.d0
    real(rprec),  parameter, private :: desired_cfl       = 0.9d0  ! Courant number (allowed, desired)

    real(rprec), parameter, private :: dx      = (mx - 0) / float(mx),  &  ! ---grid spacing
                                       dy      = (my - 0) / float(my)
    integer, parameter, private :: RCMCLAWLIM = -1 !NOTE: Hard-coded limiter below
    
    CONTAINS

    subroutine claw2ez95(rcm_time_step_s,qIN,aIN,bIN,cIN,final_iterations)

    !     An easy-to-use clawpack driver routine for simple applications
    !     Documentation is available at
    !                 http://www.amath.washington.edu/~claw/doc.html

    !     Authors: Randall J. LeVeque, Grady I. Lemoine
    !     Modified by Colby Lemon to customize for RCM

        implicit none

        real(rprec), intent(in)   :: rcm_time_step_s
        
        real(rprec), intent(in),    dimension(mx, my) :: aIN,bIN,cIN
        real(rprec), intent(inout), dimension(mx, my) :: qIN
        integer, intent(out) :: final_iterations

        real(rprec) :: max_allowed_dt, actual_max_cfl
        real(rprec) :: aux(maux, -1:mx+2, -1:my+2)
        real(rprec) :: q(-1:mx+2, -1:my+2)

    !---set auxiliary arrays
        aux(1, 1:mx, 1:my) = aIN(1:mx, 1:my)   ! loc_didt from RCM-E, colatitude (i) velocity at grid point/cell (i,j)
        aux(2, 1:mx, 1:my) = bIN(1:mx, 1:my)   ! loc_djdt from RCM-E, MLT (j) velocity at grid point/cell (i,j)
        aux(3, 1:mx, 1:my) = cIN(1:mx, 1:my)   ! loc_rate from RCM-E, the loss term

        q(1:mx,1:my) = qIN
        max_allowed_dt = rcm_time_step_s  ! try doing the whole thing in 1 step. What do we have to lose?

        call claw2(q, aux, rcm_time_step_s, max_allowed_dt, actual_max_cfl, RCMCLAWLIM, final_iterations)
        final_iterations = 0

        if (final_iterations > max_allowed_steps) then
            WRITE(*,*) 'ERROR: maximum time steps exceeded in clawpack'
            WRITE(*,*) max_allowed_dt, actual_max_cfl, RCMCLAWLIM, final_iterations
            STOP
        end if

        qIN(:,:) = q(1:mx,1:my)
    end subroutine claw2ez95

!====
!Below here are the clawpack routines
    subroutine claw2(q, aux, t_end, max_allowed_dt, actual_max_cfl,limiter, iterations)

        implicit none
        real(rprec),  intent(inout), dimension(-1:mx+2, -1:my+2) :: q    ! for RCM, this is eeta
        real(rprec),  intent(inout), dimension(maux, -1:mx+2, -1:my+2) :: aux  !inout because BCs change outer grid cell values
        real(rprec),  intent(in)             :: t_end, max_allowed_dt
        real(rprec),  intent(out)            :: actual_max_cfl
        integer, intent(in)             :: limiter
        integer, intent(out)            :: iterations

        real(rprec) :: t, t_old, dt, cfl, q_old(-1:mx+2, -1:my+2)
        integer :: n, ibc, jbc

        t = 0.0; iterations = 0; actual_max_cfl = 0.d0; dt = max_allowed_dt  ! initialize variables

        do n = 1,max_allowed_steps

            DO ibc=1,2                                     ! colatitude BC: zero-order extrapolation
                q  (   1-ibc,:) = q  (   1,:)
                aux(:, 1-ibc,:) = aux(:, 1,:)
                q  (  mx+ibc,:) = q  (  mx,:)
                aux(:,mx+ibc,:) = aux(:,mx,:)
            ENDDO

            DO jbc=1,2                                     ! local time BC: periodic
               q   (  :,1-jbc)  = q   (  :,my+1-jbc)
               aux (:,:,1-jbc)  = aux (:,:,my+1-jbc)
               q   (  :,my+jbc) = q   (  :,jbc)
               aux (:,:,my+jbc) = aux (:,:,jbc)
            END DO

            t_old = t; q_old = q                             ! save initial state in case we need to go back

            dt = min(dt, t_end - t_old)                      ! don't overstep t_end

            call step2(q,aux,dt,limiter,cfl)                 ! take one step on the conservation law

            iterations = iterations + 1; t = t_old + dt      ! update counters

            if (cfl <= max_allowed_cfl) then                 ! was Courant number too large?
                actual_max_cfl = max(cfl, actual_max_cfl)    ! no: accept this step
            else
                t = t_old; q = q_old                         ! yes: reject this step
                dt = min(max_allowed_dt, dt*desired_cfl/cfl) ! take a smaller step
                continue
            endif

            q(:,:) = q(:,:)*exp(-aux(3,:,:)*dt)              ! aux(3,:,:) is the loss term from the RCM

            if (cfl == 0.d0) then
                dt = max_allowed_dt
            else
                dt = min(max_allowed_dt, dt*desired_cfl/cfl) ! choose new time step
            endif

            if (t >= t_end) exit                             ! are we done yet?
        end do

        if (t < t_end .AND. iterations == max_allowed_steps) then
            iterations = max_allowed_steps + 1               ! Signal to caller we've exceeded max_allowed_steps
        end if

    end subroutine claw2

    subroutine step2(q,aux,dt,limiter,cfl)
    ! Take one time step, updating q.
    ! On entry, q is the initial data for this step
    ! On exit,  q is updated for the time step.

        implicit none

        real(rprec), intent(inout) :: q(        -1:mx+2, -1:my+2)
        real(rprec), intent(in)    :: aux(maux, -1:mx+2, -1:my+2), dt
        integer, intent(in)   :: limiter
        real(rprec), intent(out)   :: cfl

        real(rprec) ::    gadd(2,-1:mx+2, -1:my+2)
        real(rprec), dimension(  -1:mx+2, -1:my+2) :: cfl2d, u, v, cqxx, xwave, ywave, amdq, apdq
        real(rprec) :: dtdx, dtdy
        integer :: i, j

    ! The method used is specified by order and transverse_order:
    ! order == 1 if only first order increment waves are to be used;
    ! order == 2 if second order correction terms are to be added, with a flux limiter as specified by mthlim.
        integer, parameter :: order = 2

    ! transverse_order ...
    ! = 0 if no transverse propagation is to be applied. Increment and perhaps correction waves are propagated normal to the interface.
    ! = 1 if transverse propagation of increment waves (but not correction waves, if any) is to be applied.
    ! = 2 if transverse propagation of correction waves is also to be included.
        integer, parameter :: transverse_order = 2

        dtdx = dt/dx; dtdy = dt/dy
        u(:,:) = aux(1,:,:)   ! x-velocity
        v(:,:) = aux(2,:,:)   ! y-velocity

        do i = 0, mx+2
            xwave(i,:) = q(i,:) - q(i-1,:)   ! Note that the i'th Riemann problem has left state q(:,i-1) and right state q(:,i)
        end do
        do j = 0, my+2
            ywave(:,j) = q(:,j) - q(:,j-1)
        end do

    !-------------- X-dimension Flux calculation (previously in subroutine "flux2") ---------------
    ! rpn2 subroutine was here: solve Riemann problem at each interface and compute Godunov updates

        amdq(:,:) = min(u, 0.d0) * xwave  ! The flux difference df = s*wave  all goes in the downwind direction:
        apdq(:,:) = max(u, 0.d0) * xwave

        q(1:mx,:) = q(1:mx,:) - dtdx*(apdq(1:mx,:) + amdq(2:mx+1,:))

        cfl2d = abs(dtdx*u)
        cfl = maxval(cfl2d)

        ! modify F fluxes for second order q_{xx} correction terms:
        if (order > 1) then
            if (limiter /= 0) call flux_limiter(1,limiter,xwave,u,cfl2d)  !   apply limiter to waves

            cqxx(:,:) = abs(u) * (1.d0 - abs(u)*dtdx) * xwave

            q(1:mx,:) = q(1:mx,:) - dtdx*0.5*(cqxx(2:mx+1,:) - cqxx(1:mx,:))
        end if

        ! modify G fluxes for transverse propagation
        ! rpt2 is embedded here, instead of as a separate subroutine
        if (transverse_order > 0) then
            if (transverse_order == 2) then  ! incorporate cqxx into amdq and apdq so that it is split also.
                do j = 0, my+1; do i = 0, mx+1
                    gadd(1,i,j) = -0.5d0 * dtdx * min(v(i,j  ), 0.d0) * (amdq(i+1,j) + apdq(i,j) + cqxx(i+1,j) - cqxx(i,j))
                    gadd(2,i,j) = -0.5d0 * dtdx * max(v(i,j+1), 0.d0) * (amdq(i+1,j) + apdq(i,j) + cqxx(i+1,j) - cqxx(i,j))
                end do; end do
            else
                do j = 0, my+1; do i = 0, mx+1
                    gadd(1,i,j) = -0.5d0 * dtdx * min(v(i,j  ), 0.d0) * (amdq(i+1,j) + apdq(i,j))
                    gadd(2,i,j) = -0.5d0 * dtdx * max(v(i,j+1), 0.d0) * (amdq(i+1,j) + apdq(i,j))
                end do; end do
            end if

            do j=1,my
                q(:,j) = q(:,j) + dtdy*(gadd(1,:,j) - gadd(2,:,j) - gadd(1,:,j+1) + gadd(2,:,j-1))
            end do
        end if

    !-------------- Y-dimension Flux calculation (previously in subroutine "flux2") ---------------
    ! rpn2 subroutine was here: solve Riemann problem at each interface and compute Godunov updates

        amdq(:,:) = min(v, 0.d0) * ywave  ! The flux difference df = s*wave  all goes in the downwind direction:
        apdq(:,:) = max(v, 0.d0) * ywave

        q(:,1:my) = q(:,1:my) - dtdy*(apdq(:,1:my) + amdq(:,2:my+1))

        cfl2d = abs(dtdy*v)
        cfl = max(cfl, maxval(cfl2d))

        ! modify F fluxes for second order q_{xx} correction terms:
        if (order > 1) then
            if (limiter /= 0) call flux_limiter(2,limiter,ywave,v,cfl2d)  !   apply limiter to waves

            cqxx(:,:) = abs(v) * (1.d0 - abs(v)*dtdy) * ywave

            q(:,1:my) = q(:,1:my) - dtdy*0.5*(cqxx(:,2:my+1) - cqxx(:,1:my))
        end if

        ! modify G fluxes for transverse propagation
        ! rpt2 is embedded here, instead of as a separate subroutine
        if (transverse_order > 0) then
            if (transverse_order == 2) then  ! incorporate cqxx into amdq and apdq so that it is split also.
                do j = 0, my+1; do i = 0, mx+1
                    gadd(1,i,j) = -0.5d0 * dtdy * min(u(i  ,j), 0.d0) * (amdq(i,j+1) + apdq(i,j) + cqxx(i,j+1) - cqxx(i,j))
                    gadd(2,i,j) = -0.5d0 * dtdy * max(u(i+1,j), 0.d0) * (amdq(i,j+1) + apdq(i,j) + cqxx(i,j+1) - cqxx(i,j))
                end do; end do
            else
                do j = 0, my+1; do i = 0, mx+1
                    gadd(1,i,j) = -0.5d0 * dtdy * min(u(i  ,j), 0.d0) * (amdq(i,j+1) + apdq(i,j))
                    gadd(2,i,j) = -0.5d0 * dtdy * max(u(i+1,j), 0.d0) * (amdq(i,j+1) + apdq(i,j))
                end do; end do
            end if

            do i = 1, mx
                q(i,:) = q(i,:) + dtdx*(gadd(1,i,:) - gadd(2,i,:) - gadd(1,i+1,:) + gadd(2,i-1,:))
            end do
        end if

    end subroutine step2

    subroutine flux_limiter(ixy,limiter,wave,s,cfl2d)

    ! Apply a limiter to the waves.
    ! The limiter is computed by comparing the 2-norm of each wave with
    ! the projection of the wave from the interface to the left or
    ! right onto the current wave.  For a linear system this would
    ! correspond to comparing the norms of the two waves.  For a
    ! nonlinear problem the eigenvectors are not colinear and so the
    ! projection is needed to provide more limiting in the case where the
    ! neighboring wave has large norm but points in a different direction
    ! in phase space.
    !
    ! The specific limiter used in each family is determined by the
    ! value of the corresponding element of the array mthlim, as used in
    ! the function philim.
    ! Note that a different limiter may be used in each wave family.
    !
    ! dotl and dotr denote the inner product of wave with the wave to
    ! the left or right.  The norm of the projections onto the wave are then
    ! given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm
    ! of wave.

        implicit none

        integer, intent(in) :: ixy, limiter
        real(rprec), intent(inout) :: wave(-1:mx+2,-1:my+2)
        real(rprec), intent(in) :: s(-1:mx+2,-1:my+2), cfl2d(-1:mx+2,-1:my+2)

        real(rprec), dimension(-1:mx+2,-1:my+2) :: dotl, dotr, r, c
        real(rprec), dimension(-1:mx+2,-1:my+2) :: cfmod1, cfmod2, theta, beta, ultra, s1, s2, phimax
        integer :: i,j

    ! The following loops should be conceptually equivalent to (but somewhat faster than):
    !        dotr( 0:mx+1,:) = wave( 0:mx+1,:)*wave(1:mx+2,:)
    !        dotl( 0:mx+1,:) = wave(-1:mx,  :)*wave(0:mx+1,:)
    !        where (s > 0.d0)
    !            r = dotl/wave**2
    !        elsewhere
    !            r = dotr/wave**2
    !        end where

        if (ixy == 1) then
            dotr(-1,:) = 0.0
            do j = 0, my+1; do i = 0, mx+1
                dotl(i,j) = dotr(i-1,j)
                dotr(i,j) = wave(i,j)*wave(i+1,j)

                if (wave(i,j) == 0) cycle

                if (s(i,j) > 0.0) then
                    r(i,j) = dotl(i,j)/wave(i,j)**2
                else
                    r(i,j) = dotr(i,j)/wave(i,j)**2
                end if
    !            wave = wave*max(0.d0, min(1.d0, 2.d0*r), min(2.d0, r))
            end do; end do
        else
            dotr(:,-1) = 0.0
            do i = 0, mx+1; do j = 0, my+1
                dotl(i,j) = dotr(i,j-1)
                dotr(i,j) = wave(i,j)*wave(i,j+1)

                if (wave(i,j) == 0) cycle

                if (s(i,j) > 0.0) then
                    r(i,j) = dotl(i,j)/wave(i,j)**2
                else
                    r(i,j) = dotr(i,j)/wave(i,j)**2
                end if
    !            wave = wave*max(0.d0, min(1.d0, 2.d0*r), min(2.d0, r))
            end do; end do
        end if
        

    ! Hard code limiter above (and comment out the code below) to get slight speedup
        !K: Hard-coding limiter to avoid "illegal instruction" error w/ OMP
        wave = wave*max(0.d0, min(1.d0, r))

        ! select case(limiter)

        !     case(1)  ! Minmod
        !         wave = wave*max(0.d0, min(1.d0, r))

        !     case(2)  ! Superbee
        !         wave = wave*max(0.d0, min(1.d0, 2.d0*r), min(2.d0, r))

        !     case(3)  ! Van Leer
        !         wave = wave*(r + abs(r)) / (1.d0 + abs(r))

        !     case(4)  ! Monotonized - Centered
        !         c = (1.d0 + r)/2.d0
        !         wave = wave*max(0.d0, min(c, 2.d0, 2.d0*r))

        !     case(5)  ! Beam Warming
        !         wave = wave*r

        !     case(6)  ! Fromm
        !         wave = wave*0.5*(1 + r)

        !     case(7)  ! Albada 2
        !         wave = wave*(r**2 + r) / (1 + r**2)

        !     case(8)  ! Albada 3
        !         wave = wave*0.5*(1+r)*(1 - (abs(1-r)**3) / (1+abs(r)**3))

        !     case(9)  ! Roe's Ultrabee (CFL-superbee)
        !         wave = wave*max(0.0, min(1.0, 2.0*r/max(0.001, cfl2d)), &
        !                              min(r,   2.0/max(0.001, 1-cfl2d)))

        !     case(10)  ! Modified CFL-superbee with Beta = 2/3, theta = 1.0
        !         beta = 0.666666667; theta = 1.0

        !         ultra = max(0.0, min(theta*2.0/max(0.001, 1-cfl2d),  &
        !                            r*theta*2.0/max(0.001, cfl2d)))
        !         s2 = (1.0 + cfl2d)/3.0

        !         wave = wave*max(0.0, min(ultra, max(1.0 + (s2 - beta/2.0) * (r - 1.0),  &
        !                                             1.0 + (s2 + beta/2.0) * (r - 1.0)) ))

        !     case(11)  ! Arora Roe
        !         cfmod1 = max(0.001, cfl2d); cfmod2 = min(0.999, cfl2d)
        !         theta = 0.95

        !         s1 = 2.0 / cfmod1
        !         s2 = (1.0 + cfl2d) / 3.0
        !         phimax = 2.0 / (1.0 - cfmod2)

        !         wave = wave*min( max(1.0 + s2*(r-1.0), (1.0 - theta)*s1), &
        !                          max(theta*s1*r, (1.0 - theta)*phimax*r), &
        !                          theta*phimax)
        ! end select

    end subroutine flux_limiter

END MODULE rcmclaw
!  Solves a 2D hyperbolic system of conservation laws of the general form
!
!     capa * q_t + A q_x + B q_y = psi
!
!  The "capacity function" capa(x,y) and source term psi are optional (see below).
!
!  More complete documentation:  http://www.amath.washington.edu/~claw
!   A description of the input parameters is at the bottom of this module file
!
!  The user must supply the following subroutines:
!
!  bc2:        subroutine to specify the boundary conditions
!  rpn2, rpt2: subroutines specifying the Riemann solvers.
!  src2:       subroutine to solve capa * q_t = psi over a single time step.
!
!  b4step2  The routine b4step2 is called each time step and
!           can be supplied by the user in order to perform
!           other operations that are necessary every time
!           step.  For example, if the variables stored in
!           the aux arrays are time-dependent then these
!           values can be set. (b4step2 not used by RCM)


!  Description of parameters...
!  ----------------------------

!    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
!        On input:  initial data at time t_start.
!        On output: final solution at time rcm_time_step_s.
!        q(m,i,j) = value of mth quantity in the (i,j) cell.
!        Values within the physical domain are in q(m,i,j)
!                for i = 1,2,...,mx   and j = 1,2,...,my.
!        mbc extra cells on each end are needed for boundary conditions
!        as specified in the routine bc2.

!    aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
!        Array of auxiliary variables that are used in specifying the problem.
!        If method(7) = 0 then there are no auxiliary variables and aux
!                         can be a dummy variable.
!        If method(7) = maux > 0 then there are maux auxiliary variables
!                         and aux must be dimensioned as above.

!        Capacity functions are one particular form of auxiliary variable.
!        These arise in some applications, e.g. the
!        determinant of the Jacobian if a mapped grid is used, or a density
!        or porosity function in some advection problems.
!        See Clawpack Note # 5 for examples.

!        If method(6) = 0 then there is no capacity function.
!        If method(6) = mcapa > 0  then there is a capacity function and
!            capa(i,j), the "capacity" of the (i,j) cell, is assumed to be
!            stored in aux(mcapa,i,j).
!            In this case we require method(7).ge.mcapa.

!    mwaves is the number of waves that result from the
!           solution of each Riemann problem.  Often mwaves = meqn but
!           for some problems these may be different, e.g. for the Euler
!           equations meqn = 4 but mwaves = 3 since there are only 3
!           distinct wave speeds.

!    mbc is the number of "ghost cells" that must be added on to each
!       side of the domain to handle boundary conditions.  The cells
!       actually in the physical domain are labelled from 1 to mx in x and
!       from 1 to my in y.  The arrays are dimensioned actually indexed
!       from 1-mbc to mx+mbc and from 1-mbc to my+mbc.
!       For the methods currently implemented, mbc = 2 should be used.
!       If the user implements another method that has a larger stencil and
!       hence requires more ghost cells, a larger value of mbc could be used.
!       q is extended from the physical domain to the ghost cells by the
!       user-supplied routine bc2.

!    mx is the number of grid cells in the x-direction, in the
!       physical domain.  In addition there are mbc grid cells
!       along each edge of the grid that are used for boundary
!       conditions.

!    my is the number of grid cells in the y-direction, in the
!       physical domain.  In addition there are mbc grid cells
!       along each edge of the grid that are used for boundary
!       conditions.



!    dx = grid spacing in x.
!         (for a computation in ax <= x <= bx,  set dx = (bx-ax)/mx.)

!    dy = grid spacing in y.
!         (for a computation in ay <= y <= by,  set dy = (by-ay)/my.)

!    t_start = initial time (here we always start at t=0)

!    rcm_time_step_s = Desired final time (on input).
!              If rcm_time_step_s<t_start, then claw2 returns after a single successful
!                 time step has been taken (single-step mode).
!              Otherwise, as many steps are taken as needed to reach rcm_time_step_s,
!                 up to a maximum of nv(1).
!         = Actual time reached (on output).

!    dtv(1:3) = array of values related to the time step:
!         dtv(3) = smallest dt that was used (on output)
!         dtv(4) = largest dt that was used (on output)
!         dtv(5) = dt used in the last step (on output)

!    cflv(1:4) = array of values related to Courant number:
!         cflv(1) = maximum Courant number to be allowed.
!                   With variable time steps the step is retracted and a
!                   smaller step taken if the Courant
!                   number is larger than this value.
!                   With fixed time steps the routine aborts.
!                   Usually cflv(1)=1.0 should work
!                   (or cflv(1)=0.5 if method(3)=0).
!         cflv(2) = unused if method(1) = 0.
!                 = desired Courant number if method(1) = 1.
!                   Should be somewhat less than cflv(1), e.g. 0.9
!         cflv(3) = largest Courant number observed (on output).
!         cflv(4) = Courant number in last step (on output).

!    nv(1:2) = array of values related to the number of time steps:
!         nv(1) = unused if method(1) = 0
!               = maximum number of time steps allowed if method(1) = 1
!         iterations = number of time steps taken (on output).

!    method(1:7) = array of values specifying the numerical method to use
!                  and also indicating whether source terms, capacity
!                  function, auxiliary variables are present in the equation.

!         method(1) = 0 if fixed size time steps are to be taken.
!                       In this case, dt = dtv(1) in all steps.
!                   = 1 if variable time steps are to be used.
!                       In this case, dt = dtv(1) in the first step and
!                       thereafter the value cflv(2) is used to choose the
!                       next time step based on the maximum wave speed seen
!                       in the previous step.  Note that since this value
!                       comes from the previous step, the Courant number will
!                       not in general be exactly equal to the desired value
!                       If the actual Courant number in the next step is
!                       greater than cflv(1), then this step is redone with a
!                       smaller dt.

!         method(2) = 1 if only first order increment waves are to be used.
!                   = 2 if second order correction terms are to be added, with
!                       a flux limiter as specified by mthlim.


!         method(3) = 0 if no transverse propagation is to be applied.
!                       Increment and perhaps correction waves are propagated
!                       normal to the interface.
!                   = 1 if transverse propagation of increment waves
!                       (but not correction waves, if any) is to be applied.
!                   = 2 if transverse propagation of correction waves is also
!                       to be included.

!                   = -1 if dimensional splitting is to be used instead
!                        of the multi-dimensional wave-propagation.  The
!                        Godunov splitting is used which consists of
!                        sweeping first in x and then in y, with a step of
!                        length dt in each.  The routine bc2 is called
!                        before either sweep to set boundary data, and in
!                        the x-sweep goes over the rows of ghost cells too
!                        so that proper boundary conditions should be set
!                        for the y-sweeps by this process.  Dimensional
!                        splitting is somewhat faster than the unsplit
!                        method and works as well for many (though not all)
!                        problems.

!                   = -2 if dimensional splitting is to be used with the
!                        Strang splitting, consisting of
!                           sweep in x over time dt/2
!                           sweep in y over time dt
!                           sweep in x over time dt/2
!                        This is not recommended because it is slower than
!                        the Godunov splitting and does not appear to be
!                        appreciably better.  Moreover, the boundary
!                        conditions will not be properly set for the final
!                        x-sweep.  (The code could be modified to achieve
!                        this by sweeping over more ghost cells.)

!         method(4) = 0 to suppress printing
!                   = 1 to print dt and Courant number every time step

!         method(5) = 0 if there is no source term psi.  In this case
!                       the subroutine src2 is never called so a dummy
!                       parameter can be given.
!                   = 1 if there is a source term.  In this case
!                       the subroutine src2 must be provided and a
!                       fractional step method is used.
!                       In each time step the following sequence is followed:
!                            call bc to extend data to ghost cells
!                            call step2 to advance hyperbolic eqn by dt
!                            call src2 to advance source terms by dt
!                   = 2 if there is a source term and Strang splitting is to
!                       be used instead of the Godunov splitting above.
!                       In each time step the following sequence is followed:
!                            call bc to extend data to ghost cells
!                            call src2 to advance source terms by dt/2
!                            call step2 to advance hyperbolic equation by dt
!                            call src2 to advance source terms by dt/2
!                       For most problems 1 is recommended rather than 2
!                       since it is less expensive and works essentially as
!                       well on most problems.


!         method(6) = 0 if there is no capacity function capa.
!                   = mcapa > 0 if there is a capacity function.  In this case
!                       aux(mcapa,i,j) is the capacity of cell (i,j) and you
!                       must also specify method(7) .ge. mcapa and set aux.

!         method(7) = 0 if there is no aux array used.
!                   = maux > 0  if there are maux auxiliary variables.

!         The recommended choice of methods for most problems is
!            method(1) = 1,  method(2) = 2,  method(3) = 2.


!    mthlim(1:mwaves) = array of values specifying the flux limiter to be used
!                     in each wave family mw.  Often the same value will be used
!                     for each value of mw, but in some cases it may be
!                     desirable to use different limiters.  For example,
!                     for the Euler equations the superbee limiter might be
!                     used for the contact discontinuity (mw=2) while another
!                     limiter is used for the nonlinear waves.  Several limiters
!                     are built in and others can be added by modifying the
!                     subroutine philim.

!        mthlim(mw) = 0 for no limiter
!                   = 1 for minmod
!                   = 2 for superbee
!                   = 3 for van Leer
!                   = 4 for monotonized centered


! =========================================================================

!  Copyright 1994 -- 1999 R. J. LeVeque

!  This software is made available for research and instructional use only.
!  You may copy and use this software without charge for these non-commercial
!  purposes, provided that the copyright notice and associated text is
!  reproduced on all copies.  For all other uses (including distribution of
!  modified versions), please contact the author at the address given below.

!  *** This software is made available "as is" without any assurance that it
!  *** will work for your purposes.  The software may in fact have defects, so
!  *** use the software at your own risk.

!  --------------------------------------
!    CLAWPACK Version 4.1,  August, 2002
!    Webpage: http://www.amath.washington.edu/~claw
!  --------------------------------------

!    Author:  Randall J. LeVeque
!             Applied Mathematics
!             Box 352420
!             University of Washington,
!             Seattle, WA 98195-2420
!             rjl@amath.washington.edu
! =========================================================================

