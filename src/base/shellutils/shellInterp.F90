! Various data structures and routines to do interpolation from (and to?) a spherical shell grid
! Routines are fairly broken up to try and allow for a reasonable balance of flexibility and performance
! TODO: Add routine to take list of scattered points and interpolate to ShellGrid_T
module shellInterp
    use kdefs
    use math
    use shellgrid
    use shellutils

    implicit none

    integer, parameter, private :: NumTSC = 9

    contains


    ! InterpShellVar - TODO
    ! InterpShellVar_ChildToParent - TODO
    ! InterpShellVar_ParentToChild - TODO
    ! InterpShellVar_TSC_SG
    ! InterpShellVar_TSC_pnt

    subroutine InterpShellVar(sgSource, sgVar, sgDest, varOut)
        !! This is meant to be the highest-abstraction option for interpolation
        !! Not expected to be used very frequently, but maybe it makes sense for some simple variables
        !! Basically, the source and destination are well-defined within ShellGrid's intended use, 
        !!   so we can make the best decisions here on how to interpolate this variable
        type(ShellGrid_T)   , intent(in) :: sgSource
            !! Source shellGrid that sgVar lives on
        type(ShellGridVar_T), intent(in) :: sgVar
            !! the variable we are interpolating
        type(ShellGrid_T)   , intent(in) :: sgDest
            !! The destination grid
        type(ShellGridVar_T), intent(inout) :: varOut
            !! sgVar interpolated ontp sgDest


    end subroutine InterpShellVar


    subroutine InterpShellVar_TSC_SG(sgSource, sgVar, sgDest, varOut, dThetaO, dPhiO)
        !! Interpolate a ShellGridVar to another ShellGrid using the Triangle-Schaped Cloud method
        type(ShellGrid_T)   , intent(in) :: sgSource
            !! Source shellGrid that sgVar lives on
        type(ShellGridVar_T), intent(in) :: sgVar
            !! the variable we are interpolating
        type(ShellGrid_T)   , intent(in) :: sgDest
            !! The destination grid
        type(ShellGridVar_T), intent(inout) :: varOut
            !! sgVar interpolated ontp sgDest
        real(rp), dimension(:), optional, intent(in) :: dThetaO
            !! Cell width in theta direction
        real(rp), dimension(:), optional, intent(in) :: dPhiO
            !! Width in theta direction


        integer :: extraPnt
        integer :: i,j
        real(rp), dimension(:), allocatable :: dTheta
        real(rp), dimension(:), allocatable :: dPhi
        logical :: goodInterp

        if ( present(dThetaO) ) then
            ! First, make sure they are the same shape
            if (size(dTheta) .ne. sgVar%Ni) then
                write(*,*) "WARNING: InterpShellVar_TSC_SG got a mismatched dThetaO shape. Dying."
                stop
            endif
            ! Otherwise, we can safely copy one to the other
            allocate(dTheta(sgVar%isv:sgVar%iev))
            dTheta = dThetaO
        else
            ! If dTheta not present, we calculate it ourselves

            if (sgVar%loc == SHCC .or. sgVar%loc == SHFPH) then
                ! Note that the location id is the variable's 1D location w.r.t. the theta axis
                call calcdx_TSC(sgSource%th, sgSource%isg, sgSource%ieg, SHCC, dTheta)
            elseif (sgVar%loc == SHCORNER .or. sgVar%loc == SHFTH) then
                call calcdx_TSC(sgSource%th, sgSource%isg, sgSource%ieg, SHCORNER, dTheta)
            endif
            !! Note: calcdx_TSC sets the dimension of dTheta and dPhi
        endif

        if ( present(dPhiO) ) then
            ! First, make sure they are the same shape
            if (size(dPhiO) .ne. sgVar%Nj) then
                write(*,*) "WARNING: InterpShellVar_TSC_SG got a mismatched dPhiO shape. Dying."
                stop
            endif
            ! Otherwise, we can safely copy one to the other
            allocate(dPhi(sgVar%jsv:sgVar%jev))
            dPhi = dPhiO
        else
            ! If dPhi not present, we calculate it ourselves
            if (sgVar%loc == SHCC .or. sgVar%loc == SHFTH) then
                call calcdx_TSC(sgSource%ph, sgSource%jsg, sgSource%jeg, SHCC, dPhi)
            elseif (sgVar%loc == SHCORNER .or. sgVar%loc == SHFPH) then
                call calcdx_TSC(sgSource%ph, sgSource%jsg, sgSource%jeg, SHCORNER, dPhi)
            endif
        endif

        ! Now that we have our dTheta and dPhi, we can start interpolating

        ! Which destination grid locations we loop over depends on the destination variable location
        select case(varOut%loc)
            case(SHCC)
                !do j=varOut%jsv,varOut%jev
                !    do i=varOut%isv,varOut%iev
                !^^^ This indexing works just fine, but I'm not gonna do it cause its less clear what we're actually looping over
                !$OMP PARALLEL DO default(shared) collapse(1) &
                !$OMP schedule(dynamic) &
                !$OMP private(i,j)
                do j=sgDest%jsg,sgDest%jeg
                    do i=sgDest%isg,sgDest%ieg
                        if (.not. varOut%mask(i,j)) cycle
                        ! NOTE/TODO: This is where we would do transformations of destination grid's theta dn phi to source grid
                        ! in the case where they have different coordinate systems
                        call InterpShellVar_TSC_pnt( \
                                sgSource, sgVar,\
                                sgDest%thc(i), sgDest%phc(j),\
                                varOut%data(i,j), \
                                dTheta, dPhi,\
                                goodInterp)
                        ! TODO: Handle case where goodInterp is false here
                        ! Probably will be model dependent. Maybe we return a 2D goodInterp array if an optional array is provided to us
                    enddo
                enddo
            case(SHCORNER)
                !$OMP PARALLEL DO default(shared) collapse(1) &
                !$OMP schedule(dynamic) &
                !$OMP private(i,j)
                do j=sgDest%jsg,sgDest%jeg+1
                    do i=sgDest%isg,sgDest%ieg+1
                        if (.not. varOut%mask(i,j)) cycle
                        ! NOTE/TODO: This is where we would do transformations of destination grid's theta dn phi to source grid
                        ! in the case where they have different coordinate systems
                        call InterpShellVar_TSC_pnt( \
                                sgSource, sgVar,\
                                sgDest%th(i), sgDest%ph(j),\
                                varOut%data(i,j), \
                                dTheta, dPhi,\
                                goodInterp)
                    enddo
                enddo
            case(SHFTH)
                !$OMP PARALLEL DO default(shared) collapse(1) &
                !$OMP schedule(dynamic) &
                !$OMP private(i,j)
                do j=sgDest%jsg,sgDest%jeg
                    do i=sgDest%isg,sgDest%ieg+1
                        if (.not. varOut%mask(i,j)) cycle
                        ! NOTE/TODO: This is where we would do transformations of destination grid's theta dn phi to source grid
                        ! in the case where they have different coordinate systems
                        call InterpShellVar_TSC_pnt( \
                                sgSource, sgVar,\
                                sgDest%th(i), sgDest%phc(j),\
                                varOut%data(i,j), \
                                dTheta, dPhi,\
                                goodInterp)
                    enddo
                enddo
            case(SHFPH)
                !$OMP PARALLEL DO default(shared) collapse(1) &
                !$OMP schedule(dynamic) &
                !$OMP private(i,j)
                do j=sgDest%jsg,sgDest%jeg+1
                    do i=sgDest%isg,sgDest%ieg
                        if (.not. varOut%mask(i,j)) cycle
                        ! NOTE/TODO: This is where we would do transformations of destination grid's theta dn phi to source grid
                        ! in the case where they have different coordinate systems
                        call InterpShellVar_TSC_pnt( \
                                sgSource, sgVar,\
                                sgDest%thc(i), sgDest%ph(j),\
                                varOut%data(i,j), \
                                dTheta, dPhi,\
                                goodInterp)
                    enddo
                enddo
        end select

        ! Fill j ghosts with active cell data
        call wrapJ_SGV(sgDest, varOut)
        

    end subroutine InterpShellVar_TSC_SG


    subroutine InterpShellVar_TSC_pnt(sgsource, sgVar, th, pin, Qinterp, dThetaO, dPhiO, goodInterpO)
        !! Given the source information, interpolate sgVar to point (t,p) and return as Qout
        type(ShellGrid_T   ), intent(in) :: sgSource
            !! Source ShellGrid
        type(ShellGridVar_T), intent(in) :: sgVar
            !! Variable relative to provided ShellGrid
        real(rp), intent(in) :: th
            !! Theta coordinate with respect to source grid
        real(rp), intent(in) :: pin
            !! Phi coordinate with respect to source grid
        real(rp), intent(out) :: Qinterp
            !! Interpolated value we are returning
        real(rp), dimension(sgVar%isv:sgVar%iev), optional, intent(in) :: dThetaO
            !! Cell width in theta, centered at variable positions on source grid
        real(rp), dimension(sgVar%jsv:sgVar%jev), optional, intent(in) :: dPhiO
            !! Cell width in phi, centered at variable positions on source grid
        logical, optional, intent(inout) :: goodInterpO
            !! True if we are returning a meaningful interpolated value

        
        real(rp) :: ph
            !! Cleaned-up phi location we actually use
        integer :: i0, j0
            !! i and j locations of point t,p
            !! Whether they are respect to corner, center, or a face depends on sgVar%loc
        real(rp) :: dTh, dPh
            !! dTheta and dPhi at location i0,j0 , assuming we are in domain
        real(rp) :: t0, p0
            !! Theta and Phi values at point i0,j0
        real(rp) :: eta, zeta
        real(rp), dimension(-1:+1) :: wE,wZ
        real(rp), dimension(NumTSC) :: Ws,Qs
        logical , dimension(NumTSC) :: isGs
        integer :: ipnt,jpnt,n,di,dj

        ! Default return values
        Qinterp = 0.0
        if (present(goodInterpO)) goodInterpO = .false.

        ! Check bounds
        ph = modulo(pin,2*PI)

        ! Also make sure the requested theta is in bounds
        if ( (th<0.).or.(th>PI) ) then
            write(*,*) "ERROR in InterpShellVar_TSC_pnt: Theta should be in the range [0,PI]. Quitting..."
            stop
        end if

        call getShellILoc(sgSource, sgVar%loc, th, i0, t0)  ! Sets i0 and t0 to closest i/theta values to interp point t
        call getShellJLoc(sgSource, sgVar%loc, ph, j0, p0)  ! Sets j0 and p0 to closest i/phi   values to interp point p

        if (i0 > sgVar%iev .or. i0 < sgVar%isv) then
            return
        endif

        if (j0 > sgVar%jev .or. j0 < sgVar%jsv) then
            write(*,*) "ERROR in InterpShellVar_TSC_pnt: Phi out of bounds. idx=",j0
            write(*,*) "This wasn't supposed to be possible, good job."
            return
        endif

        if (j0 > sgVar%jev .or. j0 < sgVar%jsv) then
            write(*,*) "phi oob",j0
            write(*,*)"How did you manage that?"
            return
        endif

        !if ( (i0 < sgSource%isg) .or. (i0 > sgSource%isg + sgVar%Ni - 1) ) then
        !    write(*,*)"ERROR: InterpShellVar_TSC_pnt can't handle points outside of grid yet"
        !    stop
        !endif

        ! Make sure data is good at this point
        if (.not. sgVar%mask(i0,j0)) then
            return
        endif

        ! If still here we're gonna do something, so we can tell our caller we are returning a valid value
        if (present(goodInterpO)) goodInterpO = .true.

        ! Determine dTheta and dPhi for this ij point
        if (present(dThetaO)) then
            ! First, make sure dTheta and dPhi are defined at sgVar locations     
            if ( size(dThetaO) .ne. sgVar%Ni ) then
                write(*,*)"ERROR in InterpShellVar_TSC_pnt: dTheta != sgVar%Ni"
                write(*,*) size(dThetaO), sgVar%Ni
                stop
            endif
            dTh = dThetaO(i0)
        else
            ! dThetaO array not provided, so we calculate it ourselves
            if (sgVar%loc == SHCC .or. sgVar%loc == SHFPH) then
                dTh = Diff1D_4halfh(sgSource%th, sgSource%isg, sgSource%ieg  , i0)
            else if (sgVar%loc == SHCORNER .or. sgVar%loc == SHFTH) then
                dTh = Diff1D_4h    (sgSource%th, sgSource%isg, sgSource%ieg+1, i0)
            endif
        endif

        if (present(dPhiO)) then
            if ( size(dPhiO) .ne. sgVar%Nj ) then
                write(*,*)"ERROR in InterpShellVar_TSC_pnt: dPhi != sgVar%Nj"
                write(*,*) size(dPhiO), sgVar%Nj
                stop
            endif
            dPh = dPhiO(j0)
        else
            ! dPhiO array not provided, so we calculate it ourselves
            if (sgVar%loc == SHCC .or. sgVar%loc == SHFTH) then
                dPh = Diff1D_4halfh(sgSource%ph, sgSource%jsg, sgSource%jeg  , j0)
            else if (sgVar%loc == SHCORNER .or. sgVar%loc == SHFPH) then
                dPh = Diff1D_4h    (sgSource%ph, sgSource%jsg, sgSource%jeg+1, j0)
            endif
        endif


        if (sgSource%doNP .and. (i0==sgSource%is)) then
            call interpPole(sgSource,sgVar,th,ph,Qinterp)
            ! Handle north pole and return
            write(*,*) "Not implemented!"
            stop
        endif

        ! First, if active grid has poles 
        if (sgSource%doSP .and. (i0==sgSource%ie)) then
            ! Handle south pole and return
            write(*,*) "Not implemented!"
            stop
        endif

        ! Now, if ghost grid has poles
        if (sgSource%ghostSP .and. (i0==sgsource%ieg)) then
            write(*,*) "Not implemented!"
            stop
        endif

        if (sgSource%ghostNP .and. (i0==sgSource%isg)) then
            write(*,*) "Not implemented!"
            stop
        endif

        ! Note: If still here we know i0 isn't on the boundary

        eta  = (th - t0)/dTh
        zeta = (ph - p0)/dPh

        call ClampMapVar(eta)
        call ClampMapVar(zeta)

        ! Calculate weights
        call TSCweight1D(eta ,wE)
        call TSCweight1D(zeta,wZ)

        ! Collect weights/values
        n = 1
        do dj=-1,+1
            do di=-1,+1
                ipnt = i0+di
                jpnt = j0+dj
                
                ! Wrap around boundary
                if (jpnt<1)           jpnt = sgSource%Np
                if (jpnt>sgSource%Np) jpnt = 1

                ! Do zero-grad for theta
                if (ipnt<1)           ipnt = 1
                if (ipnt>sgSource%Nt) ipnt = sgSource%Nt

                Qs(n) = sgVar%data(ipnt,jpnt)
                Ws(n) = wE(di)*wZ(dj)
                isGs(n) = sgVar%mask(ipnt,jpnt)

                if (.not. isGs(n)) Ws(n) = 0.0

                n = n + 1
            enddo
        enddo

        ! Renormalize
        Ws = Ws/sum(Ws)

        ! Get final value
        Qinterp = dot_product(Qs,Ws)

        ! Have some internal functions
        contains
        
        ! Clamps mapping in [-0.5,0.5]
        subroutine ClampMapVar(ez)
            REAL(rp), intent(inout) :: ez
            if (ez<-0.5) ez = -0.5
            if (ez>+0.5) ez = +0.5
        end subroutine ClampMapVar

        ! 1D triangular shaped cloud weights
        ! 1D weights for triangular shaped cloud interpolation
        ! Assuming on -1,1 reference element, dx=1
        ! Check for degenerate cases ( |eta| > 0.5 )
        subroutine TSCweight1D(eta,wE)
            real(rp), intent(in)  :: eta
            real(rp), intent(out) :: wE(-1:1)

            wE(-1) = 0.5*(0.5-eta)**2.0
            wE( 1) = 0.5*(0.5+eta)**2.0
            wE( 0) = 0.75 - eta**2.0

        end subroutine TSCweight1D

    end subroutine InterpShellVar_TSC_pnt


    subroutine getShellILoc(shGr, varLoc, t, iLoc, tLocO)
        type(ShellGrid_T), intent(in) :: shGr
        integer :: varLoc
            !! Location id of the source variable
        real(rp), intent(in) :: t
        integer, intent(out) :: iLoc
        real(rp), optional, intent(out) :: tLocO

        real(rp) :: tLoc

        if (varLoc == SHCC .or. varLoc == SHFPH) then
            !! Variable is defined at center w.r.t. theta direction
            if ( (t>shGr%maxGTheta) ) then                
                iLoc = shGr%ieg+ceiling((t-shGr%maxGTheta)/(shGr%th(shGr%ieg+1)-shGr%th(shGr%ieg)))
                tLoc = shGr%thc(shGr%ieg)  ! Just return the last available theta value
                !write(*,*)"theta going out of bounds",t,shGr%maxGTheta
            else if ( (t<shGr%minGTheta) ) then
                iLoc = shGr%isg-ceiling((shGr%minGTheta-t)/(shGr%th(shGr%isg+1)-shGr%th(shGr%isg)))
                tLoc = shGr%thc(shGr%isg)
                !write(*,*)"theta going out of bounds",t,shGr%minGTheta
            else
                ! If still here then the lat bounds are okay, find closest lat cell center
                iLoc = minloc( abs(shGr%thc-t),dim=1 )
                tLoc = shGr%thc(iLoc)
            endif

            if (present(tLocO)) then
                tLocO = tLoc
            endif

        elseif (varLoc == SHCORNER .or. varLoc == SHFTH) then
            !! Variable is defined at corners w.r.t. theta direction
            if ( (t>shGr%maxTheta) ) then
                iLoc = shGr%ieg+1 + floor( 0.5 + (t-shGr%maxGTheta)/(shGr%th(shGr%ieg+1)-shGr%th(shGr%ieg)) )
                tLoc = shGr%th(shGr%ieg+1)  ! Just return the last available theta value
                !write(*,*)"theta going out of bounds",t,shGr%maxGTheta
            else if ( (t < shGr%minTheta)) then
                iLoc = shGr%isg   - floor( 0.5 + (shGr%minGTheta-t)/(shGr%th(shGr%isg+1)-shGr%th(shGr%isg)) )
                tLoc = shGr%th(shGr%isg)
                !write(*,*)"theta going out of bounds",t,shGr%maxGTheta
            else
                ! If still here then the lat bounds are okay, find closest lat cell corner
                iLoc = minloc( abs(shGr%th-t),dim=1 )
                tLoc = shGr%th(iLoc)
            endif

            if (present(tLocO)) then
                tLocO = tLoc
            endif

        endif

    end subroutine getShellILoc


    subroutine getShellJLoc(shGr, varLoc, pin, jLoc, pLoc)
        type(ShellGrid_T), intent(in) :: shGr
        integer :: varLoc
            !! Location id of the source variable
        real(rp), intent(in) :: pin
        integer, intent(out) :: jLoc
        real(rp), optional, intent(out) :: pLoc

        real(rp) :: p, dp, dJ

        p = modulo(pin,2*PI)

        ! note, shellGrid only implements [0,2pi] grids
        ! but do this check here in case it's needed in the future
        if ( (p>shGr%maxPhi) .or. (p<shGr%minPhi) ) then
            ! Point not on this grid, get outta here
            write(*,*) "ERROR in getShellJLoc, phi outside of bounds"
            write(*,*) p, shGr%minPhi, shGr%maxPhi
            stop
        endif

        if (varLoc == SHCC .or. varLoc == SHFTH) then
            !! Variable is defined at centers w.r.t. phi direction
            if (shGr%isPhiUniform) then
                ! note this is faster, thus preferred
                dp = shGr%phc(2)-shGr%phc(1)
                dJ = p/dp
                jLoc = floor(dJ) + 1
            else
                jLoc = minloc( abs(shGr%phc-p),dim=1 ) ! Find closest lat cell center
            endif

            if (present(pLoc)) then
                pLoc = shGr%phc(jLoc)
            endif

        elseif (varLoc == SHCORNER .or. varLoc == SHFPH) then
            !! Variable is defined at corners w.r.t. phi direction
            if (shGr%isPhiUniform) then
                ! note this is faster, thus preferred
                dp = shGr%ph(2)-shGr%ph(1)
                dJ = p/dp + 0.5
                jLoc = floor(dJ) + 1
            else
                jLoc = minloc( abs(shGr%ph-p),dim=1 ) ! Find closest lat cell center
            endif

            if (present(pLoc)) then
                pLoc = shGr%ph(jLoc)
            endif

        endif

    end subroutine getShellJLoc

    !! Big TODO here
    subroutine interpPole(shGr,Qin,t,pin,Qinterp)
        type(ShellGrid_T), intent(in) :: shGr
        type(ShellGridVar_T), intent(in)  :: Qin
        real(rp), intent(out) :: Qinterp
        real(rp), intent(in)  :: t,pin
        real(rp) :: f0,f1,f2,I1,I2
        integer :: j,pole,iind ! the latter is the index of the pole cell (first/last for NORTH/SOUTH)
        integer :: jpi2,jpi32,jpi  ! which cell do pi/2, 3pi/2 and pi points belong to

        Qinterp = 0.0

        ! first, find out which pole we're at
        ! note, if we're inside this function, we already know we're at one of the poles
        if ( (t.ge.0).and.(t.le.shGr%th(shGr%is+1)) ) then
            pole = NORTH
            iind = shGr%is
        else if ( (t.le.PI).and.(t.ge.shGr%th(shGr%ie)) ) then
            pole = SOUTH
            iind = shGr%ie
        else
            write(*,*) "Inside interpPole. Shouldn't be here. Quitting..."
        end if

        write(*,*) 'which pole ',pole,iind


        ! represent the function near pole to first order in theta as
        ! f(t,p) = f0 + f1*cos(p)*t + f2*sin(p)*t 
        ! (Lewis&Bellan, J. Math. Phys. 31, 2592 (1990); 
        ! https://doi.org/10.1063/1.529009
        ! 
        ! then, 2pi*f0 = \int_0^2pi f(t(i=1),p), where t(i=1) is the cell center of first cell in i-direction
        ! to calculate f1, define
        ! I1 = \int_(-pi/2)^(pi/2) f(t,p)dp = - \int_(pi/2)^(3pi/2) f(t,p)dp = 2*f1*t+f0*pi
        ! compute I1 = 0.5*(\int_(-pi/2)^(pi/2) f(t,p)dp - \int_(pi/2)^(3pi/2))
        ! to take all points around the ring into account
        ! then f1 = (I1-f0*pi)/(2*t)
        !
        ! similarly,
        ! f2 = (I2-f0*pi)/(2*t), where
        ! I2 = 0.5*(\int_0^pi f(t,p)dp - \int_pi^(2pi) f(t,p)dp)

        f0 = 0.
        do j=1,shGr%Np
            f0 = f0 + (shGr%ph(j+1)-shGr%ph(j))*Qin%data(iind,j)/(2.*PI)

            ! find which cells do pi/2, 3pi/2 and pi points belong to
            ! this can be done a priori but it doesn't add to the computation
            if ( (shGr%ph(j).le.0.5*pi).and.(shGr%ph(j+1).gt.0.5*pi) ) jpi2  = j
            if ( (shGr%ph(j).le.    pi).and.(shGr%ph(j+1).gt.    pi) ) jpi   = j
            if ( (shGr%ph(j).le.1.5*pi).and.(shGr%ph(j+1).gt.1.5*pi) ) jpi32 = j
        end do
        
        write(*,*) 'pi indices ',jpi2,jpi,jpi32

    end subroutine interpPole


!------
! Low-level interp helpers
!------

    subroutine calcdx_TSC(x, isg, ieg, loc, dx)
        !! Calculates dx for positional array x using a 4-point stencil
        !! This is desirable for TSC so that we can have better accuracy 
        !! in the case of a non-uniformly space grid
        
        ! Note, even though we are not doing much in this function, we are 
        ! breaking it out here in case we want to do higher-order for TSC or something later
        real(rp), dimension(isg:ieg+1), intent(in) :: x
            !! 1D spatial grid, should always be cell corners
        integer, intent(in) :: isg, ieg
            !! Start/end indices of cell centers
        integer, intent(in) :: loc
            !! Location identifier, MUST BE SHCC OR SHCORNER
            !! This is the location of the points we are calculating dx at relative to x
        real(rp), dimension(:), allocatable, intent(out) :: dx
            !! 'cell width' we return

        integer :: i

        if (allocated(dx)) deallocate(dx)

        if (loc == SHCORNER) then
            
            allocate(dx(isg:ieg+1))
            do i=isg,ieg+1
                dx(i) = Diff1D_4h(x, isg, ieg+1, i)
            enddo

        else if (loc == SHCC) then

            allocate(dx(isg:ieg))
            do i=isg,ieg
                dx(i) = Diff1D_4halfh(x, isg, ieg, i)
            enddo

        else
            write(*,*) "ERROR: Invalid location id in calcdx_TSC. Must be SHCC or SHCORNER"
            stop
        endif

    end subroutine calcdx_TSC


    function Diff1D_4halfh(Q,is,ie,i0) result(Qp)
        !! Use 4-point stencil to calculate first derivative of coordinates
        !! This is for the case where we are using corners to calculate the difference at cell center
        !! adapted from chimp/ebinit.F90
        real(rp), intent(in) :: Q(is:ie+1)
            !! Cell corners
        integer, intent(in) :: is,ie
            !! Start and end indices of bounding grid
        integer, intent(in) :: i0
            !! Index of the point we are evaluating, offset half a cell from its lower bound Q(i0)
        real(rp) :: Qp
        real(rp) :: Qblk(4),c(4)

        ! Note that we have fewer cases than Diff1D_4h
        ! That's because even when i0 is the last cell-centered coordinate,
        ! we still have 1 usable corner value before we reach the void
        if (i0 == is) then
            ! Q coordinates at -0.5,0.5,1.5,2.5 relative to our point
            Qblk = [Q(is), Q(is+1), Q(is+2), Q(is+3)]
            c = [-23.0, 21.0, 3.0, -1.0]/25.0
        else if (i0 == ie) then
            ! Q coordinates at -2.5,-1.5,-0.5,0.5 relative to our point
            Qblk = [Q(ie-2), Q(ie-1), Q(ie), Q(ie+1)]
            c = [1.0, -3.0, -21.0, 23.0]/25.0
        else
            ! Q coordinates at -1.5,-0.5,0.5,1.5 relative to our point
            Qblk = [Q(i0-1), Q(i0), Q(i0+1), Q(i0+2)]
            c = [1.0, -27.0, 27.0, -1.0]/24.0
        endif
        Qp = dot_product(Qblk,c)

    end function Diff1D_4halfh
    

    function Diff1D_4h(Q,is,ie,i0) result(Qp)
        !! Use 4-point stencil to calculate first derivative of coordinates
        !! This is for the case where position we are calculating the difference for is at Q(i0)
        !! (In contrast to e.g. using cell corner values to calculate the difference located at the cell center)
        !! adapted from chimp/ebinit.F90
        real(rp), intent(in) :: Q(is:ie)
            !! Cell corners
        integer, intent(in) :: is,ie,i0
        real(rp) :: Qp
        real(rp) :: Qblk(4),c(4)
        if (i0 == is) then
            !Forward
            Qblk = [Q(is),Q(is+1),Q(is+2),Q(is+3)]
            c = [-11.0,18.0,-9.0,2.0]/6.0
        else if (i0 == is+1) then
            !1 back
            Qblk = [Q(is),Q(is+1),Q(is+2),Q(is+3)]
            c = [-2.0,-3.0,6.0,-1.0]/6.0
        else if (i0 == ie) then
            Qblk = [Q(ie-3),Q(ie-2),Q(ie-1),Q(ie)]
            c = [-2.0,9.0,-18.0,11.0]/6.0
        else if (i0 == ie-1) then
            Qblk = [Q(ie-3),Q(ie-2),Q(ie-1),Q(ie)]
            c = [1.0,-6.0,3.0,2.0]/6.0
        else
            !Centered
            Qblk = [Q(i0-2),Q(i0-1),Q(i0+1),Q(i0+2)]
            c = [1.0,-8.0,8.0,-1.0]/12.0
        endif
        Qp = dot_product(Qblk,c)
        
    end function Diff1D_4h

end module shellInterp
