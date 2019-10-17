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
        
        integer :: i,j,k
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

        associate(ebModel=>vApp%ebTrcApp%ebModel,ebGr=>vApp%ebTrcApp%ebState%ebGr,ebState=>vApp%ebTrcApp%ebState)
        t = ebState%eb1%time
        
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyz,x1,x2)
        do k=ebGr%ks,ebGr%ke
            do j=ebGr%js,ebGr%je
                do i=ebGr%is,ebGr%is+vApp%iDeep
                    xyz = ebGr%xyzcc(i,j,k,:)
                    call ProjectXYZ(ebModel,ebState,xyz,t,x1,x2)

                    vApp%chmp2mhd%xyzSquish(i,j,k,1) = x1
                    vApp%chmp2mhd%xyzSquish(i,j,k,2) = x2

                enddo
            enddo
        enddo

        vApp%chmp2mhd%iMax = vApp%iDeep

        end associate
    end subroutine Squish

    !Project XYZ to R,phi at Z=0 plane
    subroutine Proj2LP(ebModel,ebState,xyz,t,x1,x2)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebState_T)  , intent(in) :: ebState
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(in) :: t
        real(rp), intent(out) :: x1,x2

        real(rp) :: L,phi,z
        real(rp), dimension(NDIM) :: xy0

        call getProjection(ebModel,ebState,xyz,t,xy0)
        !Map projection to L,phi

        z = abs(xy0(ZDIR))
        L = sqrt(xy0(XDIR)**2.0+xy0(YDIR)**2.0)
        phi = atan2(xy0(YDIR),xy0(XDIR))
        if (phi<0) phi = phi+2*PI

        if (z/L > 1.0e-3) then
            !Probably failed to get to equator, set L=0
            x1 = 0.0
            x2 = phi
        else
            x1 = L
            x2 = phi
        endif
    end subroutine Proj2LP

    !Project XYZ to lat-lon on ionosphere
    subroutine Proj2LL(ebModel,ebState,xyz,t,x1,x2)
        type(chmpModel_T), intent(in) :: ebModel
        type(ebState_T)  , intent(in) :: ebState
        real(rp), dimension(NDIM), intent(in) :: xyz
        real(rp), intent(in) :: t
        real(rp), intent(out) :: x1,x2

        real(rp), dimension(NDIM) :: xE,xIon

        x1 = 0.0
        x2 = 0.0
        !Use one-sided projection routine from chimp
        !Trace along field line (i.e. to northern hemisphere)
        call project(ebModel,ebState,xyz,t,xE,+1,toEquator=.false.)

        if ( norm2(xE) <= 3.0 .and. (xE(ZDIR)>0) ) then
            !Endpoint is close-ish to earth and in northern hemisphere, this is closed
            x1 = InvLatitude(xE) !Invariant latitude 
            x2 = atan2(xE(YDIR),xE(XDIR))
            if (x2 < 0) x2 = x2 + 2*PI
        else
            !Endpoint is far, this is probably open line
            x1 = 0.0
            x2 = 0.0
        endif

    end subroutine Proj2LL

end module ebsquish
