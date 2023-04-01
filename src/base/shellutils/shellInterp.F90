! Various data structures and routines to do interpolation from (and to?) a spherical shell grid
! TODO: Add routine to take list of scattered points and interpolate to ShellGrid_T
module shellInterp
    use kdefs
    use math
    use shellgrid

    implicit none

    integer, parameter, private :: NumTSC = 9

    contains

    ! Interpolate on grid shGr a cell-centered variable at point t(heta),p(hi)
    ! The cell-centered variable is selected via the Qind parameter
    ! Result is Qinterp
    ! Optional : isGood (Nt,Np), a mask for good/bad data
    ! Optional : isGoodP, whether Qinterp is a good value
    subroutine InterpShell(shGr,Qind,t,pin,Qinterp,isGoodP,isGood)
        type(ShellGrid_T), intent(in) :: shGr
        integer, intent(in)  :: Qind
        real(rp), intent(out) :: Qinterp
        real(rp), intent(in)  :: t,pin
        logical , intent(out), optional :: isGoodP
        logical , intent(in) , optional :: isGood(shGr%Nt,shGr%Np) ! TODO: consider if this should include ghosts

        integer :: i0,j0,ij0(2),di,dj
        integer :: ip,jp,n
        real(rp) :: p,dt,dp,eta,zeta
        real(rp), dimension(NumTSC) :: Ws,Qs
        logical , dimension(NumTSC) :: isGs
        real(rp), dimension(-1:+1) :: wE,wZ



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
        if (shGr%doSP .and. (i0==shGr%Nt)) then
            ! Handle south pole and return
            write(*,*) "Not implemented!"
            stop
        endif

        if (shGr%doNP .and. (i0==1)) then
            ! Handle north pole and return
            write(*,*) "Not implemented!"
            stop
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

                Qs(n) = shGr%shellVars(Qind)%cellData(ip,jp)
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
            
            iX = shGr%ieg+1+floor((t-shGr%maxGTheta)/(shGr%th(shGr%ieg+1)-shGr%th(shGr%ieg)))
            return
        endif

        if ( (t<shGr%minGTheta) ) then
            ! Point outside ghost grid
            ! get an idea of where we are using the last available deltaTheta
            ! and get outta here
            
            iX = shGr%isg-1-floor((shGr%minGTheta-t)/(shGr%th(shGr%isg+1)-shGr%th(shGr%isg)))
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

end module shellInterp
