!Routines to map native source grids to Bios-Savart (ground system) grids
module calcdbremap
    use ebtabutils
    use kdefs
    use chmpdefs
    use ebtypes
    use calcdbtypes
    use ebinterp
    use chmpfields
    use geopack

    implicit none

    integer, parameter, private :: i0 = 2
    logical, parameter, private :: doTest = .false.

    contains

    !Map ground coordinates to SM
    subroutine remapGR(Model,t,ebState,gGr)
        type(chmpModel_T), intent(in)    :: Model
        real(rp)         , intent(in)    :: t
        type(ebState_T)  , intent(in)    :: ebState
        type(grGrid_T)   , intent(inout) :: gGr

        real(rp) :: mjd
        integer :: i,j,k
        real(rp), dimension(NDIM) :: grXYZ,smXYZ

        if (.not. gGr%doGEO) then
            gGr%SMxyzC = gGr%GxyzC
            return
        endif

        mjd = MJDAt(ebState%ebTab,Model%t)
        call MJDRecalc(mjd) !Setup geopack for this time

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,grXYZ,smXYZ)
        do k=1,gGr%Nz
            do j=1,gGr%NLon
                do i=1,gGr%NLat
                    grXYZ = gGr%GxyzC(i,j,k,XDIR:ZDIR)
                    call GEO2SM(grXYZ(XDIR),grXYZ(YDIR),grXYZ(ZDIR),smXYZ(XDIR),smXYZ(YDIR),smXYZ(ZDIR))
                    gGr%SMxyzC(i,j,k,XDIR:ZDIR) = smXYZ
                enddo
            enddo
        enddo

    end subroutine remapGR

    !Load Bios Savart grids
    subroutine packBS(Model,t,ebState,ionGrid,facGrid,magBS,ionBS,facBS)
        type(chmpModel_T), intent(in) :: Model
        real(rp)         , intent(in) :: t
        type(ebState_T)  , intent(in) :: ebState
        type(ionGrid_T)  , intent(in) :: ionGrid
        type(facGrid_T)  , intent(in) :: facGrid
        type(BSGrid_T), intent(inout) :: magBS,ionBS,facBS

        real(rp) :: w1,w2,dV
        integer :: i,j,k,l,n
        real(rp), dimension(:,:,:,:), allocatable :: Jxyz

    !----
    !Magnetospheric part
        
    !Start by getting Jxyz at current time
        !Jxyz are SM currents at time t
        associate(ebGr=>ebState%ebGr)
        allocate(Jxyz(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        call GetTWgts(Model,ebState,t,w1,w2)
        !$OMP PARALLEL WORKSHARE
        Jxyz = w1*ebState%eb1%Jxyz + w2*ebState%eb2%Jxyz
        !$OMP END PARALLEL WORKSHARE
        
    !Now do remap and store into BSGrid
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,n,dV)
        do k=ebGr%ks,ebGr%ke
            do j=ebGr%js,ebGr%je
                do i=ebGr%is,ebGr%ie
                    !Get n-d => 1D index
                    n = ijk2n(i,j,k,ebGr%is,ebGr%ie,ebGr%js,ebGr%je,ebGr%ks,ebGr%ke)
                    dV  = ebGr%dV(i,j,k) !Volume element, Re^3

                    if (i <= ebGr%is+i0-1) then
                        dV = 0.0
                        !Zap this contribution
                    endif

                    magBS%XYZcc(n,XDIR:ZDIR) = ebGr%xyz(i,j,k,:) !Cell-centered MHD grid coordinates (SM)
                    magBS%Jxyz (n,XDIR:ZDIR) = Jxyz(i,j,k,:) !SM current at cell-center
                    magBS%dV   (n)           = dV
                    
                    if (doTest) then
                        call MagJTest( magBS%XYZcc(n,XDIR:ZDIR),magBS%Jxyz (n,XDIR:ZDIR) )
                    endif
                enddo
            enddo
        enddo

        end associate

    !----
    !Ionospheric part
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,n)
        do k=1,2 !Hemisphere
            do j=1,ionGrid%Nth !Theta
                do i=1,ionGrid%Np !phi
                    !Get n-d => 1D index
                    n = ijk2n(i,j,k,1,ionGrid%Np,1,ionGrid%Nth,1,2)

                    ionBS%XYZcc(n,XDIR:ZDIR) = ionGrid%XYZcc(i,j,k,:)
                    ionBS%Jxyz (n,XDIR:ZDIR) = ionGrid%Jxyz (i,j,k,:)
                    ionBS%dV   (n)           = ionGrid%dS   (i,j,k)

                enddo
            enddo
        enddo

    !----
    !FAC part
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,n,l)
        do l=1,2 !Hemisphere
            do k=1,facGrid%rSegs
                do j=1,facGrid%Nth
                    do i=1,facGrid%Np
                        !Get n-d => 1D index
                        n = ijkl2n(i,j,k,l,1,ionGrid%Np,1,ionGrid%Nth,1,facGrid%rSegs,1,2)
                        
                        facBS%XYZcc(n,XDIR:ZDIR) = facGrid%XYZcc(i,j,k,l,:)
                        facBS%Jxyz (n,XDIR:ZDIR) = facGrid%Jxyz (i,j,k,l,:)
                        facBS%   dV(n)           = facGrid%dV   (i,j,k,l)
                    enddo !i
                enddo
            enddo
        enddo !l

    end subroutine packBS

    !Some test routines
    subroutine MagJTest(xyz,Jxyz)
        real(rp), intent(in)  :: xyz(NDIM)
        real(rp), intent(out) :: Jxyz(NDIM)
        real(rp) :: x,absY,absZ
        Jxyz = 0.0
        x = xyz(XDIR)
        absY = abs(xyz(YDIR))
        absZ = abs(xyz(ZDIR))

        if ( (absY<128.0) .and. (absZ<5.0) .and. (x>-20) .and. (x<-10) ) then
            Jxyz(YDIR) = 1.0*1.0e+3
        endif

        if ( (absY<128.0) .and. (absZ<5.0) .and. (x<20) .and. (x>10) ) then
            Jxyz(YDIR) = 1.0*1.0e+3
        endif

    end subroutine MagJTest
!=========
!Index squashing routines

    !TODO: Rewrite these routines to be smarter
    function ijk2n(i,j,k,is,ie,js,je,ks,ke) result(n)
        integer, intent(in) :: i,j,k,is,ie,js,je,ks,ke
        integer :: n

        integer :: ip,jp,kp
        integer :: Ni,Nj

        ip = i-is+1
        jp = j-js+1
        kp = k-ks+1
        Ni = ie-is+1
        Nj = je-js+1

        !Convert to n
        n = ip + (jp-1)*Ni + (kp-1)*Ni*Nj
    end function ijk2n

    function ijkl2n(i,j,k,l,is,ie,js,je,ks,ke,ls,le) result(n)
        integer, intent(in) :: i,j,k,l,is,ie,js,je,ks,ke,ls,le
        integer :: n

        integer :: ip,jp,kp,lp
        integer :: Ni,Nj,Nk

        ip = i-is+1
        jp = j-js+1
        kp = k-ks+1
        lp = l-ls+1
        Ni = ie-is+1
        Nj = je-js+1
        Nk = ke-ks+1

        !Convert to n
        n = ip + (jp-1)*Ni + (kp-1)*Ni*Nj + (lp-1)*Ni*Nj*Nk
    end function ijkl2n

end module calcdbremap