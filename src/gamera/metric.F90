!Various routines to get cell-level metric quantities

module metric
    use gamtypes
    use math

    implicit none

    contains

    !Given (i,j,k) return cell centered coordinates (xc,yc,zc) of cell i,j,k
    subroutine cellCenter(Gr,i,j,k,xc,yc,zc)
        integer, intent(in) :: i,j,k
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(out) :: xc,yc,zc

        xc = 0.125*( Gr%xyz(i,j,k  ,XDIR) + Gr%xyz(i+1,j,k  ,XDIR) + Gr%xyz(i+1,j+1,k  ,XDIR) + Gr%xyz(i,j+1,k  ,XDIR) + &
                     Gr%xyz(i,j,k+1,XDIR) + Gr%xyz(i+1,j,k+1,XDIR) + Gr%xyz(i+1,j+1,k+1,XDIR) + Gr%xyz(i,j+1,k+1,XDIR) )

        yc = 0.125*( Gr%xyz(i,j,k  ,YDIR) + Gr%xyz(i+1,j,k  ,YDIR) + Gr%xyz(i+1,j+1,k  ,YDIR) + Gr%xyz(i,j+1,k  ,YDIR) + &
                     Gr%xyz(i,j,k+1,YDIR) + Gr%xyz(i+1,j,k+1,YDIR) + Gr%xyz(i+1,j+1,k+1,YDIR) + Gr%xyz(i,j+1,k+1,YDIR) )

        zc = 0.125*( Gr%xyz(i,j,k  ,ZDIR) + Gr%xyz(i+1,j,k  ,ZDIR) + Gr%xyz(i+1,j+1,k  ,ZDIR) + Gr%xyz(i,j+1,k  ,ZDIR) + &
                     Gr%xyz(i,j,k+1,ZDIR) + Gr%xyz(i+1,j,k+1,ZDIR) + Gr%xyz(i+1,j+1,k+1,ZDIR) + Gr%xyz(i,j+1,k+1,ZDIR) )

    end subroutine cellCenter

    !Calculate face center (in direction d, shifted down) of cell i,j,k
    subroutine faceCenter(Gr,i,j,k,xfc,yfc,zfc,d)
        integer, intent(in) :: i,j,k,d
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(out) :: xfc,yfc,zfc

        ! select direction
        select case(d)
            case(IDIR)
                xfc = 0.25*( Gr%xyz(i,j,k,XDIR) + Gr%xyz(i,j+1,k,XDIR) + Gr%xyz(i,j+1,k+1,XDIR) + Gr%xyz(i,j,k+1,XDIR) )
                yfc = 0.25*( Gr%xyz(i,j,k,YDIR) + Gr%xyz(i,j+1,k,YDIR) + Gr%xyz(i,j+1,k+1,YDIR) + Gr%xyz(i,j,k+1,YDIR) )
                zfc = 0.25*( Gr%xyz(i,j,k,ZDIR) + Gr%xyz(i,j+1,k,ZDIR) + Gr%xyz(i,j+1,k+1,ZDIR) + Gr%xyz(i,j,k+1,ZDIR) )
            case(JDIR)
                xfc = 0.25*( Gr%xyz(i,j,k,XDIR) + Gr%xyz(i+1,j,k,XDIR) + Gr%xyz(i+1,j,k+1,XDIR) + Gr%xyz(i,j,k+1,XDIR) )
                yfc = 0.25*( Gr%xyz(i,j,k,YDIR) + Gr%xyz(i+1,j,k,YDIR) + Gr%xyz(i+1,j,k+1,YDIR) + Gr%xyz(i,j,k+1,YDIR) )
                zfc = 0.25*( Gr%xyz(i,j,k,ZDIR) + Gr%xyz(i+1,j,k,ZDIR) + Gr%xyz(i+1,j,k+1,ZDIR) + Gr%xyz(i,j,k+1,ZDIR) )
            case(KDIR)
                xfc = 0.25*( Gr%xyz(i,j,k,XDIR) + Gr%xyz(i+1,j,k,XDIR) + Gr%xyz(i+1,j+1,k,XDIR) + Gr%xyz(i,j+1,k,XDIR) )
                yfc = 0.25*( Gr%xyz(i,j,k,YDIR) + Gr%xyz(i+1,j,k,YDIR) + Gr%xyz(i+1,j+1,k,YDIR) + Gr%xyz(i,j+1,k,YDIR) )
                zfc = 0.25*( Gr%xyz(i,j,k,ZDIR) + Gr%xyz(i+1,j,k,ZDIR) + Gr%xyz(i+1,j+1,k,ZDIR) + Gr%xyz(i,j+1,k,ZDIR) )
        end select

    end subroutine faceCenter

    !d-Face coordinates in order f0=SW,f1=SE,f2=NW,f3=NE
    subroutine faceCoords(Model,Grid,i,j,k,d,f0,f1,f2,f3)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        integer, intent(in) :: i,j,k,d
        real(rp), intent(out), dimension(NDIM) :: f0,f1,f2,f3

        f0(XDIR) = Grid%xyz(i,j,k,XDIR)
        f0(YDIR) = Grid%xyz(i,j,k,YDIR)
        f0(ZDIR) = Grid%xyz(i,j,k,ZDIR)
        
        select case(d)
        case(IDIR)
            !I face
            f1(XDIR) = Grid%xyz(i,j+1,k,XDIR)
            f1(YDIR) = Grid%xyz(i,j+1,k,YDIR)
            f1(ZDIR) = Grid%xyz(i,j+1,k,ZDIR)

            f3(XDIR) = Grid%xyz(i,j+1,k+1,XDIR)
            f3(YDIR) = Grid%xyz(i,j+1,k+1,YDIR)
            f3(ZDIR) = Grid%xyz(i,j+1,k+1,ZDIR)

            f2(XDIR) = Grid%xyz(i,j,k+1,XDIR)
            f2(YDIR) = Grid%xyz(i,j,k+1,YDIR)
            f2(ZDIR) = Grid%xyz(i,j,k+1,ZDIR)
        case(JDIR)
            !J face
            f1(XDIR) = Grid%xyz(i,j,k+1,XDIR)
            f1(YDIR) = Grid%xyz(i,j,k+1,YDIR)
            f1(ZDIR) = Grid%xyz(i,j,k+1,ZDIR)

            f3(XDIR) = Grid%xyz(i+1,j,k+1,XDIR)
            f3(YDIR) = Grid%xyz(i+1,j,k+1,YDIR)
            f3(ZDIR) = Grid%xyz(i+1,j,k+1,ZDIR)

            f2(XDIR) = Grid%xyz(i+1,j,k,XDIR)
            f2(YDIR) = Grid%xyz(i+1,j,k,YDIR)
            f2(ZDIR) = Grid%xyz(i+1,j,k,ZDIR)
        case(KDIR)
            !K face
            f1(XDIR) = Grid%xyz(i+1,j,k,XDIR)
            f1(YDIR) = Grid%xyz(i+1,j,k,YDIR)
            f1(ZDIR) = Grid%xyz(i+1,j,k,ZDIR)

            f3(XDIR) = Grid%xyz(i+1,j+1,k,XDIR)
            f3(YDIR) = Grid%xyz(i+1,j+1,k,YDIR)
            f3(ZDIR) = Grid%xyz(i+1,j+1,k,ZDIR)

            f2(XDIR) = Grid%xyz(i,j+1,k,XDIR)
            f2(YDIR) = Grid%xyz(i,j+1,k,YDIR)
            f2(ZDIR) = Grid%xyz(i,j+1,k,ZDIR)
        end select  
                        

    end subroutine faceCoords

    !d-Edge coordinates in order
    subroutine edgeCoords(Model,Grid,i,j,k,d,e1,e2)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        integer, intent(in) :: i,j,k,d
        real(rp), intent(out), dimension(NDIM) :: e1,e2

        e1 = [Grid%xyz(i,j,k,XDIR),Grid%xyz(i,j,k,YDIR),Grid%xyz(i,j,k,ZDIR)]
        select case(d)
        case(IDIR)
            !I edge
            e2 = [Grid%xyz(i+1,j,k,XDIR),Grid%xyz(i+1,j,k,YDIR),Grid%xyz(i+1,j,k,ZDIR)]
        case(JDIR)
            !J edge
            e2 = [Grid%xyz(i,j+1,k,XDIR),Grid%xyz(i,j+1,k,YDIR),Grid%xyz(i,j+1,k,ZDIR)]
        case(KDIR)
            !K edge
            e2 = [Grid%xyz(i,j,k+1,XDIR),Grid%xyz(i,j,k+1,YDIR),Grid%xyz(i,j,k+1,ZDIR)]
        end select
    end subroutine edgeCoords

    !8 xyz corners of cell i,j,k
    subroutine cellCoords(Model,Grid,i,j,k,xyzC)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        integer, intent(in) :: i,j,k
        real(rp), intent(out), dimension(8,NDIM) :: xyzC

        xyzC(1,:) = [Grid%xyz(i  ,j  ,k  ,XDIR),Grid%xyz(i  ,j  ,k  ,YDIR),Grid%xyz(i  ,j  ,k  ,ZDIR)]
        xyzC(2,:) = [Grid%xyz(i+1,j  ,k  ,XDIR),Grid%xyz(i+1,j  ,k  ,YDIR),Grid%xyz(i+1,j  ,k  ,ZDIR)]
        xyzC(3,:) = [Grid%xyz(i+1,j+1,k  ,XDIR),Grid%xyz(i+1,j+1,k  ,YDIR),Grid%xyz(i+1,j+1,k  ,ZDIR)]
        xyzC(4,:) = [Grid%xyz(i  ,j+1,k  ,XDIR),Grid%xyz(i  ,j+1,k  ,YDIR),Grid%xyz(i  ,j+1,k  ,ZDIR)]
        xyzC(5,:) = [Grid%xyz(i  ,j  ,k+1,XDIR),Grid%xyz(i  ,j  ,k+1,YDIR),Grid%xyz(i  ,j  ,k+1,ZDIR)]
        xyzC(6,:) = [Grid%xyz(i+1,j  ,k+1,XDIR),Grid%xyz(i+1,j  ,k+1,YDIR),Grid%xyz(i+1,j  ,k+1,ZDIR)]
        xyzC(7,:) = [Grid%xyz(i+1,j+1,k+1,XDIR),Grid%xyz(i+1,j+1,k+1,YDIR),Grid%xyz(i+1,j+1,k+1,ZDIR)]
        xyzC(8,:) = [Grid%xyz(i  ,j+1,k+1,XDIR),Grid%xyz(i  ,j+1,k+1,YDIR),Grid%xyz(i  ,j+1,k+1,ZDIR)]

    end subroutine cellCoords

    !Gets i,j,k basis vectors for normal plane to dN (e1,e2)
    subroutine getNormalPlane(dN,e1,e2)
        integer, intent(in) :: dN
        integer, intent(out) :: e1(NDIM), e2(NDIM)

        !Select direction
        select case(dN)
            case(IDIR)
                e1 = [0,1,0]
                e2 = [0,0,1]
            case(JDIR)
                e1 = [1,0,0]
                e2 = [0,0,1]
            case(KDIR)
                e1 = [1,0,0]
                e2 = [0,1,0]
        end select
    end subroutine getNormalPlane
    
!Metric quantities for edge ijk -> cell xyz calculation
!Cell-centered i,j,k hat vectors
!----------------------------------------------
    function ijkVec(Model,Grid,i,j,k,d) result(xCC)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        integer, intent(in) :: i,j,k,d
        real(rp), dimension(NDIM) :: xCC

        ! select direction
        select case(d)
          case(IDIR)
            ! i-vector at cell center
            xCC(XDIR) = 0.25*( (Grid%xyz(i+1,j,k,XDIR) + Grid%xyz(i+1,j+1,k,XDIR) + Grid%xyz(i+1,j,k+1,XDIR) + Grid%xyz(i+1,j+1,k+1,XDIR)) - &
                               (Grid%xyz(i,  j,k,XDIR) + Grid%xyz(i,  j+1,k,XDIR) + Grid%xyz(i,  j,k+1,XDIR) + Grid%xyz(i,  j+1,k+1,XDIR)) )
            xCC(YDIR) = 0.25*( (Grid%xyz(i+1,j,k,YDIR) + Grid%xyz(i+1,j+1,k,YDIR) + Grid%xyz(i+1,j,k+1,YDIR) + Grid%xyz(i+1,j+1,k+1,YDIR)) - &
                               (Grid%xyz(i,  j,k,YDIR) + Grid%xyz(i,  j+1,k,YDIR) + Grid%xyz(i,  j,k+1,YDIR) + Grid%xyz(i,  j+1,k+1,YDIR)) )
            xCC(ZDIR) = 0.25*( (Grid%xyz(i+1,j,k,ZDIR) + Grid%xyz(i+1,j+1,k,ZDIR) + Grid%xyz(i+1,j,k+1,ZDIR) + Grid%xyz(i+1,j+1,k+1,ZDIR)) - &
                               (Grid%xyz(i,  j,k,ZDIR) + Grid%xyz(i,  j+1,k,ZDIR) + Grid%xyz(i,  j,k+1,ZDIR) + Grid%xyz(i,  j+1,k+1,ZDIR)) )
          case(JDIR)
            ! j-vector at cell center
            xCC(XDIR) = 0.25*( (Grid%xyz(i,j+1,k,XDIR) + Grid%xyz(i+1,j+1,k,XDIR) + Grid%xyz(i,j+1,k+1,XDIR) + Grid%xyz(i+1,j+1,k+1,XDIR)) - &
                               (Grid%xyz(i,j,  k,XDIR) + Grid%xyz(i+1,j,  k,XDIR) + Grid%xyz(i,j,  k+1,XDIR) + Grid%xyz(i+1,j,  k+1,XDIR)) )
            xCC(YDIR) = 0.25*( (Grid%xyz(i,j+1,k,YDIR) + Grid%xyz(i+1,j+1,k,YDIR) + Grid%xyz(i,j+1,k+1,YDIR) + Grid%xyz(i+1,j+1,k+1,YDIR)) - &
                               (Grid%xyz(i,j,  k,YDIR) + Grid%xyz(i+1,j,  k,YDIR) + Grid%xyz(i,j,  k+1,YDIR) + Grid%xyz(i+1,j,  k+1,YDIR)) )  
            xCC(ZDIR) = 0.25*( (Grid%xyz(i,j+1,k,ZDIR) + Grid%xyz(i+1,j+1,k,ZDIR) + Grid%xyz(i,j+1,k+1,ZDIR) + Grid%xyz(i+1,j+1,k+1,ZDIR)) - &
                               (Grid%xyz(i,j,  k,ZDIR) + Grid%xyz(i+1,j,  k,ZDIR) + Grid%xyz(i,j,  k+1,ZDIR) + Grid%xyz(i+1,j,  k+1,ZDIR)) )  
          case(KDIR)
            ! k-vector at cell center
            xCC(XDIR) = 0.25*( (Grid%xyz(i,j,k+1,XDIR) + Grid%xyz(i+1,j,k+1,XDIR) + Grid%xyz(i,j+1,k+1,XDIR) + Grid%xyz(i+1,j+1,k+1,XDIR)) - &
                               (Grid%xyz(i,j,k  ,XDIR) + Grid%xyz(i+1,j,k  ,XDIR) + Grid%xyz(i,J+1,k  ,XDIR) + Grid%xyz(i+1,j+1,k  ,XDIR)) )
            xCC(YDIR) = 0.25*( (Grid%xyz(i,j,k+1,YDIR) + Grid%xyz(i+1,j,k+1,YDIR) + Grid%xyz(i,j+1,k+1,YDIR) + Grid%xyz(i+1,j+1,k+1,YDIR)) - &
                               (Grid%xyz(i,j,k  ,YDIR) + Grid%xyz(i+1,j,k  ,YDIR) + Grid%xyz(i,J+1,k  ,YDIR) + Grid%xyz(i+1,j+1,k  ,YDIR)) )
            xCC(ZDIR) = 0.25*( (Grid%xyz(i,j,k+1,ZDIR) + Grid%xyz(i+1,j,k+1,ZDIR) + Grid%xyz(i,j+1,k+1,ZDIR) + Grid%xyz(i+1,j+1,k+1,ZDIR)) - &
                               (Grid%xyz(i,j,k  ,ZDIR) + Grid%xyz(i+1,j,k  ,ZDIR) + Grid%xyz(i,J+1,k  ,ZDIR) + Grid%xyz(i+1,j+1,k  ,ZDIR)) )
        end select
        xCC = xCC/norm2(xCC)

    end function ijkVec
    
    subroutine ijkMatrix(Model,Grid,iCC,jCC,kCC,Mijk)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), dimension(NDIM), intent(in) :: iCC,jCC,kCC
        real(rp), dimension(NDIM,NDIM), intent(out) :: Mijk

        Mijk(1,1) = kCC(ZDIR)*jCC(YDIR) - kCC(YDIR)*jCC(ZDIR);
        Mijk(2,1) = kCC(XDIR)*jCC(ZDIR) - kCC(ZDIR)*jCC(XDIR);
        Mijk(3,1) = kCC(YDIR)*jCC(XDIR) - kCC(XDIR)*jCC(YDIR);

        Mijk(1,2) = kCC(YDIR)*iCC(ZDIR) - kCC(ZDIR)*iCC(YDIR);
        Mijk(2,2) = kCC(ZDIR)*iCC(XDIR) - kCC(XDIR)*iCC(ZDIR);
        Mijk(3,2) = kCC(XDIR)*iCC(YDIR) - kCC(YDIR)*iCC(XDIR);

        Mijk(1,3) = jCC(ZDIR)*iCC(YDIR) - jCC(YDIR)*iCC(ZDIR);
        Mijk(2,3) = jCC(XDIR)*iCC(ZDIR) - jCC(ZDIR)*iCC(XDIR);
        Mijk(3,3) = jCC(YDIR)*iCC(XDIR) - jCC(XDIR)*iCC(YDIR);

    end subroutine ijkMatrix

!Ring singularity helper functions
!--------------------------------------
    !Take inner-most two rings of vector quantity and reinterpolate inner ring to match outer
    !Important for IJK->XYZ w/o ghosts
    subroutine FixRAVec_S(Model,Grid,Qxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp), intent(inout) :: Qxyz(2,Grid%ks:Grid%ke,NDIM)

        real(rp) :: Q0(NDIM)
        integer :: k
        real(rp) :: w1,w2

        !Check if we're doing ring avg and have either pole
        if ( Model%doRing .and. Model%Ring%doS ) then

            select case (Model%Ring%GridID)
                !------------------
                case ("lfm")
                    !Calculate Qxyz at pole (assuming we have full ring)
                    Q0 = sum ( Qxyz(2,Grid%ks:Grid%ke,:), dim=1 )/(Grid%ke-Grid%ks+1)

                    !Now loop around ring and interpolate between pole and j=2
                    !Using lazy weighting, 0 (pole) -> 0.5 (1cc) -> 1.5 (2cc)
                    w1 = 2.0/3.0
                    w2 = 1.0/3.0

                    do k=Grid%ks,Grid%ke
                        Qxyz(1,k,:) = w1*Q0 + w2*Qxyz(2,k,:)
                    enddo

            end select

        endif

    end subroutine FixRAVec_S

    !Same as above, but do for -X side pole (ie 1->2 instead of 2->1)
    !Important for IJK->XYZ w/o ghosts
    subroutine FixRAVec_E(Model,Grid,Qxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp), intent(inout) :: Qxyz(2,Grid%ks:Grid%ke,NDIM)

        real(rp) :: Q0(NDIM)
        integer :: k
        real(rp) :: w1,w2

        !Check if we're doing ring avg and have either pole
        if ( Model%doRing .and. Model%Ring%doE ) then

            select case (Model%Ring%GridID)
                !------------------
                case ("lfm")
                    !Calculate Qxyz at pole (assuming we have full ring)
                    Q0 = sum ( Qxyz(1,Grid%ks:Grid%ke,:), dim=1 )/(Grid%ke-Grid%ks+1)

                    !Now loop around ring and interpolate between pole and j=2
                    !Using lazy weighting, 0 (pole) -> 0.5 (1cc) -> 1.5 (2cc)
                    w1 = 2.0/3.0
                    w2 = 1.0/3.0

                    do k=Grid%ks,Grid%ke
                        Qxyz(2,k,:) = w1*Q0 + w2*Qxyz(1,k,:)
                    enddo

            end select

        endif

    end subroutine FixRAVec_E

!Lazy quadrature functions
!--------------------------------------
    !Lazy function for volume/barycenter quadrature
    subroutine CellDV(x,y,z,D,Vx,Vy,Vz,P)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: D,Vx,Vy,Vz,P

        D = 1
        Vx = x
        Vy = y
        Vz = z
        P = 0
    end subroutine CellDV

    !Lazy functions for surface area/face centers
    subroutine rVec(x,y,z,Ax,Ay,Az)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Ax,Ay,Az

        Ax = x
        Ay = y
        Az = z
        
    end subroutine rVec

    subroutine IdVec(x,y,z,Ax,Ay,Az)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Ax,Ay,Az

        Ax = 1.0
        Ay = 1.0
        Az = 1.0
        
    end subroutine IdVec

end module metric
