!Routines to handle CHIMP 2D<->3D mapping for deep coupling
module ebsquish
    use kdefs
    use gamtypes
    use ebtypes
    use volttypes
    use streamline
    use earthhelper
    
    implicit none

    !Projection type
    abstract interface
        subroutine Projection_T(ebModel,ebState,xyz,t,x1,x2)
            Import :: rp,NDIM,chmpModel_T,ebState_T
            type(chmpModel_T), intent(in) :: ebModel
            type(ebState_T)  , intent(in) :: ebState
            real(rp), intent(in)  :: xyz(NDIM), t
            real(rp), intent(out) :: x1,x2
        end subroutine Projection_T
    end interface

    real(rp), parameter, private :: startEps = 0.05
    real(rp), parameter, private :: rEps = 0.125
    real(rp), private :: Rinner

    contains

    !Find i-index of outer boundary of coupling domain
    function ShellBoundary(gModel,Gr,R) result(iMax)
        type(Model_T), intent(in) :: gModel
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(in) :: R
        integer :: iMax
        integer :: i

        !Just using dayside line
        !Avoiding findloc because it's not well supported by GNU
        do i=Gr%is,Gr%ie
            iMax = i
            if (Gr%xyz(i,Gr%js,Gr%ks,XDIR)>=R) exit
        enddo

    end function ShellBoundary

    !Do 3D->2D mapping from Gamera volume to 2D inner mag model grid
    subroutine Squish(vApp)
        type(voltApp_T), intent(inout) :: vApp
        
        integer :: i,j,k,Nk,nSkp
        real(rp) :: t,x1,x2
        real(rp), dimension(NDIM) :: xyz,xy0
        procedure(Projection_T), pointer :: ProjectXYZ
        
        ProjectXYZ => NULL()

        select case (vApp%prType)
        case(LPPROJ)
            ProjectXYZ => Proj2LP
        case(LLPROJ)
            ProjectXYZ => Proj2LL
        case DEFAULT
            write(*,*) 'Unkown projection type, bailing ...'
            stop
        end select

        associate(ebModel=>vApp%ebTrcApp%ebModel,ebGr=>vApp%ebTrcApp%ebState%ebGr,ebState=>vApp%ebTrcApp%ebState, &                  
                  xyzSquish=>vApp%chmp2mhd%xyzSquish,isGood=>vApp%chmp2mhd%isGood)

        t = ebState%eb1%time

        Rinner = norm2(ebGr%xyz(ebGr%is,ebGr%js,ebGr%ks,XDIR:ZDIR))

        xyzSquish = 0.0
        isGood = .false.

        if (vApp%doQkSquish) then
            nSkp = 2 !Stride through grid for projections
        else
            nSkp = 1
        endif

        !Force iDeep to be even
        vApp%iDeep = ceiling(vApp%iDeep/2.0)

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,xyz,x1,x2)
        do k=ebGr%ks,ebGr%ke+1,nSkp
            do j=ebGr%js,ebGr%je+1,nSkp
                do i=ebGr%is,vApp%iDeep+1,nSkp
                    xyz = ebGr%xyz(i,j,k,XDIR:ZDIR)
                    if (norm2(xyz) <= vApp%rTrc) then
                        !Do projection
                        call ProjectXYZ(ebModel,ebState,xyz,t,x1,x2)
                    else
                        !Set null projection because outside radius
                        x1 = 0.0
                        x2 = 0.0
                    endif
                    
                    xyzSquish(i,j,k,1) = x1
                    xyzSquish(i,j,k,2) = x2

                    if ( (abs(x1)>0.0) .and. (abs(x2)>0.0) ) then
                        isGood(i,j,k) = .true.
                    endif

                enddo
            enddo
        enddo

        if (vApp%doQkSquish) then
            call FillSkips(ebModel,ebGr,vApp%iDeep,xyzSquish,isGood)
        endif

        Nk = ebGr%ke-ebGr%ks+1

        !Ensure same value for degenerate axis points
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i)
        do i=ebGr%is,vApp%iDeep+1
            !x1 (regular average)
            xyzSquish(i,ebGr%js  ,:,1) = sum(xyzSquish(i,ebGr%js  ,ebGr%ks:ebGr%ke,1))/Nk
            xyzSquish(i,ebGr%je+1,:,1) = sum(xyzSquish(i,ebGr%je+1,ebGr%ks:ebGr%ke,1))/Nk
            !x2 (circular average)
            xyzSquish(i,ebGr%js  ,:,2) = CircMean(xyzSquish(i,ebGr%js  ,ebGr%ks:ebGr%ke,2))
            xyzSquish(i,ebGr%je+1,:,2) = CircMean(xyzSquish(i,ebGr%je+1,ebGr%ks:ebGr%ke,2))
        enddo

        vApp%chmp2mhd%iMax = vApp%iDeep

        end associate
    end subroutine Squish

    subroutine FillSkips(ebModel,ebGr,iDeep,xyzSquish,isGood)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebGrid_T)   , intent(in) :: ebGr
        integer          , intent(in) :: iDeep
        real(rp), intent(inout) :: xyzSquish(ebGr%is:ebGr%ie+1,ebGr%js:ebGr%je+1,ebGr%ks:ebGr%ke+1,2)
        logical , intent(inout) :: isGood   (ebGr%is:ebGr%ie+1,ebGr%js:ebGr%je+1,ebGr%ks:ebGr%ke+1)

    end subroutine FillSkips

    !Project XYZ to R,phi at Z=0 plane
    subroutine Proj2LP(ebModel,ebState,xyz,t,x1,x2)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebState_T)  , intent(in) :: ebState
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(in) :: t
        real(rp), intent(out) :: x1,x2

        real(rp) :: L,phi,z
        real(rp), dimension(NDIM) :: xyzSeed,xy0

        ! trap for when we're within epsilon of the inner boundary
        ! (really, it's probably only the first shell of nodes at R=Rinner_boundary that doesn't trace correctly)
        if ( (norm2(xyz)-Rinner)/Rinner < startEps ) then
           ! dipole-shift to startEps
           xyzSeed = DipoleShift(xyz,norm2(xyz)+startEps)
        else
           xyzSeed = xyz
        end if

        ! note, this assumes internally tracing along -B for z>0 and along B for z<0
        ! if that fails, it returns xy0=HUGE
        call getEquatorProjection(ebModel,ebState,xyzSeed,t,xy0)

        !Map projection to L,phi
        L = sqrt(xy0(XDIR)**2.0+xy0(YDIR)**2.0)
        phi = atan2(xy0(YDIR),xy0(XDIR))
        if (phi<0) phi = phi+2*PI
        
        x1 = L
        x2 = phi
    end subroutine Proj2LP

    !Project XYZ to lat-lon on ionosphere
    !TODO: Define cutoff radius within which we just dipole map for speed
    subroutine Proj2LL(ebModel,ebState,xyz,t,x1,x2)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebState_T)  , intent(in) :: ebState
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(in) :: t
        real(rp), intent(out) :: x1,x2

        real(rp), dimension(NDIM) :: xE,xIon,xyz0
        real(rp) :: dX,rC
        logical :: isGood

        x1 = 0.0
        x2 = 0.0

        ! trap for when we're within epsilon of the inner boundary
        ! (really, it's probably only the first shell of nodes at R=Rinner_boundary that doesn't trace correctly)
        if ( (norm2(xyz)-Rinner)/Rinner < startEps ) then
           ! dipole-shift to startEps
           xyz0 = DipoleShift(xyz,norm2(xyz)+startEps)
        else
           xyz0 = xyz
        end if

        !Use one-sided projection routine from chimp
        !Trace along field line (i.e. to northern hemisphere)
        call project(ebModel,ebState,xyz0,t,xE,+1,toEquator=.false.)

        dX = norm2(xyz0-xE)
        rC = Rinner*(1.+rEps)
        isGood = (dX>TINY) .and. (norm2(xE) <= rC) .and. (xE(ZDIR) > 0)

        if (isGood) then
            !Get invariant lat/lon
            x1 = InvLatitude(xE)
            x2 = atan2(xE(YDIR),xE(XDIR))
            if (x2 < 0) x2 = x2 + 2*PI
        else
            !Set 0/0 for projection failure
            x1 = 0.0
            x2 = 0.0
        endif

    end subroutine Proj2LL

end module ebsquish
