! Various data structures and routines to do interpolation from (and to?) a spherical shell grid
! TODO: Add routine to take list of scattered points and interpolate to ShellGrid_T
module shellInterp
    use kdefs
    use math

    implicit none

    integer, parameter, private :: NumTSC = 9

    !Interpolate on grid shGr a cell-centered variable (Q) at point lat,lon
    !Result is Qp
    !Optional : isGood (NLat,NLon), a mask for good/bad data
    !Optional : isGoodP, whether Qp is a good value
    ! subroutine InterpShell(shGr,Q,lat,lonin,Qp,isGoodP,isGood)
    !     type(ShellGrid_T), intent(in) :: shGr
    !     real(rp), intent(in)  :: Q(shGr%NLat,shGr%NLon)
    !     real(rp), intent(out) :: Qp
    !     real(rp), intent(in)  :: lat,lonin
    !     logical , intent(out), optional :: isGoodP
    !     logical , intent(in) , optional :: isGood(shGr%NLat,shGr%NLon)

    !     integer :: i0,j0,ij0(2),di,dj
    !     integer :: ip,jp,n
    !     real(rp) :: dlat,dphi,eta,zeta,lon
    !     real(rp), dimension(NumTSC) :: Ws,Qs
    !     logical , dimension(NumTSC) :: isGs
    !     real(rp), dimension(-1:+1) :: wE,wZ

    !     Qp = 0.0
    !     if (present(isGoodP)) then
    !         isGoodP = .false.
    !     endif

    !     !Do some short circuiting
    !     if ( (lat>shGr%maxLat) .or. (lat<shGr%minLat) ) then
    !         !Point not on this grid, get outta here
    !         return
    !     endif

    !     if (lonin<0) then
    !         lon = lonin+2*PI
    !     else
    !         lon = lonin
    !     endif

    !     !Now we know this point is in our grid
    !     call GetShellIJ(shGr,lat,lon,ij0) !Find the i,j cell this point is in

    !     i0 = ij0(1)
    !     j0 = ij0(2)

    !     if (present(isGood)) then
    !         !Check cell is good
    !         if (.not. isGood(i0,j0)) return
    !     endif

    ! !Have central cell and know that it's good
    !     isGoodP = .true.
    !     !Trap for near-pole cases
    !     if (shGr%doSP .and. (i0==1)) then
    !         !Handle south pole and return
    !         write(*,*) "Not implemented!"
    !         stop
    !     endif

    !     if (shGr%doNP .and. (i0==shGr%NLat)) then
    !         !Handle north pole and return
    !         write(*,*) "Not implemented!"
    !         stop
    !     endif

    !     !Note: If still here we know i0 isn't on the boundary

    !     !Calculate local mapping
    !     dlat = shGr%LatI(i0+1)-shGr%LatI(i0) 
    !     dphi = shGr%LonI(j0+1)-shGr%LonI(j0)

    !     eta  = ( lat - shGr%LatC(i0) )/dlat
    !     zeta = ( lon - shGr%LonC(j0) )/dphi

    !     call ClampMapVar(eta)
    !     call ClampMapVar(zeta)

    !     !Calculate weights
    !     call TSCweight1D(eta ,wE)
    !     call TSCweight1D(zeta,wZ)

    ! !Now loop over surrounding cells and get weights/values

    !     n = 1
    !     do dj=-1,+1
    !         do di=-1,+1
    !             ip = i0+di
    !             jp = j0+dj
    !             !Wrap around boundary
    !             if (jp<1)         jp = shGr%NLon
    !             if (jp>shGr%NLon) jp = 1

    !             !Do zero-grad for lat
    !             if (ip<1)         ip = 1
    !             if (ip>shGr%NLat) ip = shGr%NLat

    !             Qs(n) = Q(ip,jp)
    !             Ws(n) = wE(di)*wZ(dj)
                
    !             if (present(isGood)) then
    !                 isGs(n) = isGood(ip,jp)
    !             else
    !                 isGs(n) = .true.
    !             endif
    !             if (.not. isGs(n)) Ws(n) = 0.0

    !             n = n + 1
    !         enddo
    !     enddo !dj

    !     !Renormalize
    !     Ws = Ws/sum(Ws)

    ! !Get final value
    !     Qp = dot_product(Qs,Ws)
        
    ! !Have some internal functions
    ! contains
    
    !     !Clamps mapping in [-0.5,0.5]
    !     subroutine ClampMapVar(ez)
    !       REAL(rp), intent(inout) :: ez
    !       if (ez<-0.5) ez = -0.5
    !       if (ez>+0.5) ez = +0.5
    !     end subroutine ClampMapVar

    !     !1D triangular shaped cloud weights
    !     !1D weights for triangular shaped cloud interpolation
    !     !Assuming on -1,1 reference element, dx=1
    !     !Check for degenerate cases ( |eta| > 0.5 )
    !     subroutine TSCweight1D(eta,wE)
    !         real(rp), intent(in)  :: eta
    !         real(rp), intent(out) :: wE(-1:1)

    !         wE(-1) = 0.5*(0.5-eta)**2.0
    !         wE( 1) = 0.5*(0.5+eta)**2.0
    !         wE( 0) = 0.75 - eta**2.0

    !     end subroutine TSCweight1D

    ! end subroutine InterpShell

    ! !For a shGr type find the ij cell that the point lat/lon is in
    ! !NOTE: Returns 0,0 if the point isn't in the grid
    ! subroutine GetShellIJ(shGr,lat,lonin,ij0)
    !     type(ShellGrid_T), intent(in) :: shGr
    !     real(rp), intent(in) :: lat,lonin
    !     integer, intent(out) :: ij0(2)

    !     real(rp) :: lon,dp,dJ
    !     integer  :: iX,jX

    !     if (lonin<0) then
    !         lon = lonin+2*PI
    !     else
    !         lon = lonin
    !     endif

    !     ij0 = 0

    ! !Do some short circuiting
    !     if ( (lat>shGr%maxLat) .or. (lat<shGr%minLat) ) then
    !         !Point not on this grid, get outta here
    !         return
    !     endif

    ! !If still here then the lat bounds are okay, let's do this

    !     !Get lat part
    !     iX = minloc( abs(shGr%LatC-lat),dim=1 ) !Find closest lat cell center

    !     !Now get lon part, assume uniform phi spacing
    !     dp = shGr%LonC(2)-shGr%LonC(1)
    !     dJ = lon/dp
    !     jX = floor(dJ) + 1

    !     !Impose bounds just in case
    !     iX = max(iX,1)
    !     iX = min(iX,shGr%NLat)
    !     jX = max(jX,1)
    !     jX = min(jX,shGr%NLon)

    !     ij0 = [iX,jX]

    ! end subroutine GetShellIJ

end module shellInterp
