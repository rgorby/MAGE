! Various data structures and routines to do interpolation from (and to?) a spherical shell grid
! Routines are fairly broken up to try and allow for a reasonable balance of flexibility and performance
! TODO: Add routine to take list of scattered points and interpolate to ShellGrid_T
module shellInterp
    use kdefs
    use math
    use shellgrid

    implicit none

    integer, parameter, private :: NumTSC = 9

    contains


    ! InterpShellVar
    ! InterpShellVar_ChildToParent
    ! InterpShellVar_ParentToChild
    ! InterpShellVar_TSC_SG
    ! InterpShellVar_TSC_pnt

    subroutine InterpShellVar(sgVar, sgSource, sgDest, varOut)
        !! This is meant to be the highest-abstraction option for interpolation
        !! Not expected to be used very frequently, but maybe it makes sense for some simple variables
        !! Basically, the source and destination are well-defined within ShellGrid's intended use, 
        !!   so we can make the best decisions here on how to interpolate this variable
        type(ShellGridVar_T), intent(in) :: sgVar
            !! the variable we are interpolating
        type(ShellGrid_T)   , intent(in) :: sgSource
            !! Source shellGrid that sgVar lives on
        type(ShellGrid_T)   , intent(in) :: sgDest
            !! The destination grid
        type(ShellGridVar_T), intent(inout) :: varOut
            !! sgVar interpolated ontp sgDest


    end subroutine InterpShellVar


    subroutine InterpShellVar_TSC_SG(sgVar, sgSource, sgDest, varOut, dThetaO, dPhiO)
        !! Interpolate a ShellGridVar to another ShellGrid using the Triangle-Schaped Cloud method
        type(ShellGridVar_T), intent(in) :: sgVar
            !! the variable we are interpolating
        type(ShellGrid_T)   , intent(in) :: sgSource
            !! Source shellGrid that sgVar lives on
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
        real(rp), dimension(:), allocatable :: dTheta, dPhi

        if ( present(dThetaO) ) then
            ! First, make sure they are the same shape
            if (size(dTheta) .ne. size(dThetaO)) then
                write(*,*) "WARNING: InterpShellVar_TSC_SG got a mismatched dThetaO shape. Dying."
                stop
            endif
            ! Otherwise, we can safely copy one to the other
            ! First, allocate the right size
            extraPnt = merge(1, 0, sgVar%loc == SHCORNER .or. sgVar%loc == SHFTH)
            allocate(dTheta(sgSource%isg:sgSource%ieg+extraPnt))
            ! Then copy
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
            if (size(dPhi) .ne. size(dPhiO)) then
                write(*,*) "WARNING: InterpShellVar_TSC_SG got a mismatched dPhiO shape. Dying."
                stop
            endif
            ! Otherwise, we can safely copy one to the other
            ! First, allocate the right size
            extraPnt = merge(1, 0, sgVar%loc == SHCORNER .or. sgVar%loc == SHFPH)
            allocate(dPhi(sgSource%jsg:sgSource%jeg+extraPnt))
            ! Then copy
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
                do j=sgDest%jsg,sgDest%jeg
                    do i=sgDest%isg,sgDest%ieg
                        varOut%data(i,j) = InterpShellVar_TSC_pnt( \
                                            sgVar, sgSource,\
                                            dTheta, dPhi,\
                                            sgDest%thc(i), sgDest%phc(j))
                    enddo
                enddo
            case(SHCORNER)
                do j=sgDest%jsg,sgDest%jeg+1
                    do i=sgDest%isg,sgDest%ieg+1
                        varOut%data(i,j) = InterpShellVar_TSC_pnt( \
                                            sgVar, sgSource,\
                                            dTheta, dPhi,\
                                            sgDest%th(i), sgDest%ph(j))
                    enddo
                enddo
            case(SHFTH)
                do j=sgDest%jsg,sgDest%jeg
                    do i=sgDest%isg,sgDest%ieg+1
                        varOut%data(i,j) = InterpShellVar_TSC_pnt( \
                                            sgVar, sgSource,\
                                            dTheta, dPhi,\
                                            sgDest%th(i), sgDest%phc(j))
                    enddo
                enddo
            case(SHFPH)
                do j=sgDest%jsg,sgDest%jeg+1
                    do i=sgDest%isg,sgDest%ieg
                        varOut%data(i,j) = InterpShellVar_TSC_pnt( \
                                            sgVar, sgSource,\
                                            dTheta, dPhi,\
                                            sgDest%thc(i), sgDest%ph(j))
                    enddo
                enddo
        end select
        

    end subroutine InterpShellVar_TSC_SG


    function InterpShellVar_TSC_pnt(sgVar, sgSource, dTheta, dPhi, t, pin) result(Qout)
        !! Given the source information, interpolate sgVar to point (t,p) and return as Qout
        type(ShellGridVar_T), intent(in) :: sgVar
        type(ShellGrid_T   ), intent(in) :: sgSource
        real(rp), dimension(:), intent(in) :: dTheta
        real(rp), dimension(:), intent(in) :: dPhi
        real(rp), intent(in) :: t
        real(rp), intent(in) :: pin

        real(rp) :: Qout
        real(rp) :: p
        integer :: i0, j0
            !! i and j locations of point t,p
            !! Whether they are respect to corner, center, or a face depends on sgVar%loc

        ! First, make sure dTheta and dPhi are defined at sgVar locations        
        if ( size(dTheta) .ne. sgVar%Ni ) then
            write(*,*)"ERROR in InterpShellVar_TSC_pnt: dTheta != sgVar%Ni"
            write(*,*) size(dTheta), sgVar%Ni
            stop
        endif

        if ( size(dPhi) .ne. sgVar%Nj ) then
            write(*,*)"ERROR in InterpShellVar_TSC_pnt: dPhi != sgVar%Nj"
            write(*,*) size(dPhi), sgVar%Nj
            stop
        endif

        ! Check bounds
        p = modulo(pin,2*PI)

        ! Also make sure the requested theta is in bounds
        if ( (t<0.).or.(t>PI) ) then
            write(*,*) "ERROR in InterpShellVar_TSC_pnt: Theta should be in the range [0,PI]. Quitting..."
            stop
        end if

        i0 = getShellILoc(sgSource, sgVar%loc, t)
        j0 = getShellJLoc(sgSource, sgVar%loc, p)

        if ( (i0 < sgSource%isg) .or. (i0 > sgSource%isg + sgVar%Ni - 1) ) then
            write(*,*)"ERROR: InterpShellVar_TSC_pnt can't handle points outside of grid yet"
            stop
        endif

        Qout = 0.0
    end function InterpShellVar_TSC_pnt


    function getShellILoc(shGr, varLoc, t) result(iLoc)
        type(ShellGrid_T), intent(in) :: shGr
        integer :: varLoc
            !! Location id of the source variable
        real(rp), intent(in) :: t

        integer :: iLoc

        !if (varLoc == SHCC     .or. varLoc == SHFPH) thLoc = SHCC
        !if (varLoc == SHCORNER .or. varLoc == SHFTH) thLoc = SHCORNER

        if (varLoc == SHCC .or. varLoc == SHFPH) then
            !! Variable is defined at center w.r.t. theta direction
            if ( (t>shGr%maxGTheta) ) then                
                iLoc = shGr%ieg+ceiling((t-shGr%maxGTheta)/(shGr%th(shGr%ieg+1)-shGr%th(shGr%ieg)))
                return
            endif
    
            if ( (t<shGr%minGTheta) ) then
                iLoc = shGr%isg-ceiling((shGr%minGTheta-t)/(shGr%th(shGr%isg+1)-shGr%th(shGr%isg)))
                return
            endif
    
            ! If still here then the lat bounds are okay, find closest lat cell center
            iLoc = minloc( abs(shGr%thc-t),dim=1 )

        elseif (varLoc == SHCORNER .or. varLoc == SHFTH) then
            !! Variable is defined at corners w.r.t. theta direction
            if ( (t>shGr%maxTheta) ) then
                iLoc = shGr%ieg+1 + floor( 0.5 + (t-shGr%maxGTheta)/(shGr%th(shGr%ieg+1)-shGr%th(shGr%ieg)) )
            endif

            if ( (t < shGr%minTheta)) then
                iLoc = shGr%isg   - floor( 0.5 + (shGr%minGTheta-t)/(shGr%th(shGr%isg+1)-shGr%th(shGr%isg)) )
            endif

            ! If still here then the lat bounds are okay, find closest lat cell corner
            iLoc = minloc( abs(shGr%th-t),dim=1 )
        endif   

    end function getShellILoc


    function getShellJLoc(shGr, varLoc, pin) result(jLoc)
        type(ShellGrid_T), intent(in) :: shGr
        integer :: varLoc
            !! Location id of the source variable
        real(rp), intent(in) :: pin

        real(rp) :: p, dp, dJ
        integer :: jLoc

        p = modulo(pin,2*PI)

        ! note, shellGrid only implements [0,2pi] grids
        ! but do this check here in case it's needed in the future
        if ( (p>shGr%maxPhi) .or. (p<shGr%minPhi) ) then
            ! Point not on this grid, get outta here
            write(*,*) "ERROR in getShellJLoc, phi somehow outside of bounds"
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
            end if
        elseif (varLoc == SHCORNER .or. varLoc == SHFPH) then
            !! Variable is defined at corners w.r.t. phi direction
            if (shGr%isPhiUniform) then
                ! note this is faster, thus preferred
                dp = shGr%ph(2)-shGr%ph(1)
                dJ = p/dp + 0.5
                jLoc = floor(dJ) + 1
            else
                jLoc = minloc( abs(shGr%ph-p),dim=1 ) ! Find closest lat cell center
            end if
        endif

    end function getShellJLoc

    ! Interpolate on grid shGr a cell-centered variable at point t(heta),p(hi)
    ! Qin is the shellVar
    ! Result is Qinterp
    ! Optional : isGood (Nt,Np), a mask for good/bad data
    ! Optional : isGoodP, whether Qinterp is a good value
    subroutine InterpShell(shGr,Qin,t,pin,Qinterp,isGoodP,isGood)
        type(ShellGrid_T), intent(in) :: shGr
            !! ShellGrid of source grid
        type(ShellGridVar_T), intent(in)  :: Qin
            !! Variable stored on source grid
        real(rp), intent(out) :: Qinterp
            !! Value we return
        real(rp), intent(in)  :: t,pin
            !! Theta and Phi location of point we are interpolating to
        logical , intent(out), optional :: isGoodP
        logical , intent(in) , optional :: isGood(shGr%Nt,shGr%Np) ! TODO: consider if this should include ghosts

        integer :: i0,j0,ij0(2),di,dj
        integer :: ip,jp,n
        real(rp) :: p,dt,dp,eta,zeta
        real(rp), dimension(NumTSC) :: Ws,Qs
        logical , dimension(NumTSC) :: isGs
        real(rp), dimension(-1:+1) :: wE,wZ


        if (Qin%loc .ne. SHCC) then
            write(*,*) "InterpShell only interpolates cell-centered data right now, goodbye"
            stop
        endif

        ! if ( (t>shGr%maxTheta).and.(.not.shGr%bcsApplied(SOUTH)) ) then
        !     ! Point inside ghosts but BCs not applied, get outta here
        !     return
        ! endif

        ! if ( (t<shGr%minTheta).and.(.not.shGr%bcsApplied(NORTH)) ) then
        !     ! Point inside ghosts but BCs not applied, get outta here
        !     return
        ! endif


        Qinterp = 0.0
        if (present(isGoodP)) then
            isGoodP = .false.
        endif

        ! make sure phi is inside the [0,2pi] interval
        ! we do this inside GetShellIJ, but do it here also just in case phi is used here
        p = modulo(pin,2*PI)

        ! also make sure the requested theta is in bounds
        if ( (t<0.).or.(t>PI) ) then
            write(*,*) "Inside InterpShell."
            write(*,*) "Theta should be in the range [0,PI]. Quitting..."
            stop
        end if
        
        ! Note, check inside if the point is in our grid (including ghosts)
        ! Note, all points are within the longitude bounds because we only use periodic grids
        ! If the point is outside of latitude bounds, ij0(1) will be set to a nominal cell number extrapolated outside of bounds
        ! (i.e., how many cells outside of the boundary are we)
        call GetShellIJ(shGr,t,p,ij0) !Find the i,j cell this point is in

        i0 = ij0(1)
        j0 = ij0(2)

        if (present(isGood)) then
            ! Check cell is good
            if (.not. isGood(i0,j0)) return
        endif

        ! Have central cell and know that it's good
        if (present(isGoodP)) then
            isGoodP = .true.
        end if

        ! Trap for near-pole cases

        if (shGr%doNP .and. (i0==1)) then
            call interpPole(shGr,Qin,t,pin,Qinterp)

            ! Handle north pole and return
            write(*,*) "Not implemented!"
            stop
        endif

        ! First, if active grid has poles 
        if (shGr%doSP .and. (i0==shGr%Nt)) then
            ! Handle south pole and return
            write(*,*) "Not implemented!"
            stop
        endif

        ! Now, if ghost grid has poles
        if (shGr%ghostSP .and. (i0==shGr%ieg)) then
        endif

        if (shGr%ghostNP .and. (i0==shGr%isg)) then
        endif

        ! Note: If still here we know i0 isn't on the boundary

        ! Calculate local mapping
        dt = shGr%th(i0+1)-shGr%th(i0) 
        dp = shGr%ph(j0+1)-shGr%ph(j0)

        eta  = ( t - shGr%thc(i0) )/dt
        zeta = ( p - shGr%phc(j0) )/dp


        call ClampMapVar(eta)
        call ClampMapVar(zeta)

        ! Calculate weights
        call TSCweight1D(eta ,wE)
        call TSCweight1D(zeta,wZ)

        ! Now loop over surrounding cells and get weights/values
        n = 1
        do dj=-1,+1
            do di=-1,+1
                ip = i0+di
                jp = j0+dj
                ! Wrap around boundary
                if (jp<1)       jp = shGr%Np
                if (jp>shGr%Np) jp = 1

                ! Do zero-grad for theta
                if (ip<1)         ip = 1
                if (ip>shGr%Nt) ip = shGr%Nt

                Qs(n) = Qin%data(ip,jp)
                Ws(n) = wE(di)*wZ(dj)
                
                if (present(isGood)) then
                    isGs(n) = isGood(ip,jp)
                else
                    isGs(n) = .true.
                endif
                if (.not. isGs(n)) Ws(n) = 0.0

                n = n + 1
            enddo
        enddo !dj

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

    end subroutine InterpShell

    ! For a shGr type find the ij cell that the point t(theta)/p(hi) is in
    ! NOTE: Returns 0,0 if the point isn't in the grid
    subroutine GetShellIJ(shGr,t,pin,ij0)
        type(ShellGrid_T), intent(in) :: shGr
        real(rp), intent(in) :: t,pin
        integer, intent(out) :: ij0(2)

        real(rp) :: p,dp,dJ
        integer  :: iX,jX

        ! make sure longitude is inside the [0,2pi] interval
        p = modulo(pin,2*PI)

        ! note, shellGrid only implements [0,2pi] grids
        ! but do this check here in case it's needed in the future
        if ( (p>shGr%maxPhi) .or. (p<shGr%minPhi) ) then
            ! Point not on this grid, get outta here
            return
        endif

        ij0 = 0

        ! First get lon part, b/c we always assume periodic
        ! so always on the grid in the phi direction
        if (shGr%isPhiUniform) then
            ! note this is faster, thus preferred
            dp = shGr%phc(2)-shGr%phc(1)
            dJ = p/dp
            jX = floor(dJ) + 1
        else
            jX = minloc( abs(shGr%phc-p),dim=1 ) ! Find closest lat cell center
        end if

        ! Now do latitude, this is more complicated

        ! Do some short circuiting
        if ( (t>shGr%maxGTheta) ) then
            ! Point outside ghost grid
            ! get an idea of where we are using the last available deltaTheta
            ! and get outta here
            
            iX = shGr%ieg+ceiling((t-shGr%maxGTheta)/(shGr%th(shGr%ieg+1)-shGr%th(shGr%ieg)))
            return
        endif

        if ( (t<shGr%minGTheta) ) then
            ! Point outside ghost grid
            ! get an idea of where we are using the last available deltaTheta
            ! and get outta here
            
            iX = shGr%isg-ceiling((shGr%minGTheta-t)/(shGr%th(shGr%isg+1)-shGr%th(shGr%isg)))
            return
        endif

        ! If still here then the lat bounds are okay, let's do this

        ! Get lat part
        iX = minloc( abs(shGr%thc-t),dim=1 ) ! Find closest lat cell center

        ! ! Impose bounds just in case
        ! iX = max(iX,1)
        ! iX = min(iX,shGr%Nt)
        ! jX = max(jX,1)
        ! jX = min(jX,shGr%Np)

        ij0 = [iX,jX]

    end subroutine GetShellIJ

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
            write(*,*) "Inside interPole. Shouldn't be here. Quitting..."
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
        integer, intent(in) :: isg, ieg
            !! Start/end indices of cell centers
        real(rp), dimension(isg:ieg+1), intent(in) :: x
            !! 1D spatial grid, should always be cell corners
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
        !! This is for the case where position we are calculating the difference for is at Q(i0+1/2)
        !! e.g. we ar using corners to calculate the difference at cell center
        !! adapted from chimp/ebinit.F90
        integer, intent(in) :: is,ie
            !! Start and end indices of bounding grid
        integer, intent(in) :: i0
            !! Index of the pointe we are evaluating, offset half a cell from its lower bound Q(i0)
        real(rp), intent(in) :: Q(is:ie)
        real(rp) :: Qp
        real(rp) :: Qblk(4),c(4)

        ! Note that we have fewer cases than Diff1D_4h
        ! That's because even when i0 is the last cell-centered coordinate,
        ! we still have 1 usable corner value before we reach the void
        if (i0 == is) then
            ! Q coordinates at -0.5,0.5,1.5,2.5 relative to our point
            Qblk = [Q(is), Q(is+1), Q(is+2), Q(is+3)]
            c = [-23.0, 21.0, 3.0, -1.0]/25.0
        else if (i0 == ie-1) then
            ! Q coordinates at -2.5,-1.5,-0.5,0.5 relative to our point
            Qblk = [Q(is-2), Q(is-1), Q(is), Q(is+1)]
            c = [1.0, -3.0, -21.0, 23.0]/25.0
        else
            ! Q coordinates at -1.5,-0.5,0.5,1.5 relative to our point
            Qblk = [Q(is-1), Q(is), Q(is+1), Q(is+2)]
            c = [1.0, -27.0, 27.0, -1.0]/24.0
        endif
        Qp = dot_product(Qblk,c)

    end function Diff1D_4halfh
    

    function Diff1D_4h(Q,is,ie,i0) result(Qp)
        !! Use 4-point stencil to calculate first derivative of coordinates
        !! This is for the case where position we are calculating the difference for is at Q(i0)
        !! (In contrast to e.g. using cell corner values to calculate the difference located at the cell center)
        !! adapted from chimp/ebinit.F90
        integer, intent(in) :: is,ie,i0
        real(rp), intent(in) :: Q(is:ie)
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
