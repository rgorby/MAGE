  ! various help functions for earth

module earthhelper
  use kdefs

  implicit none
  
  contains
    !Gallagher plasmasphere density model
    !Yoinked from Slava's code
    function psphD(L) result(D)
        ! approx based on gallagher et al, figure 1
        ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 105, NO. A8, PAGES 18,819-18,833, AUGUST 1, 2000
        ! returns values in ples/cc

        implicit none
        real(rp),intent(in) :: L
        real(rp)  :: D
        real(rp), parameter :: L0 = 4.5
        real(rp), parameter :: alpha = 10.
        real(rp), parameter :: a1 = -0.25
        real(rp), parameter :: b1 = 2.4
        real(rp), parameter :: a2 = -0.5
        real(rp), parameter :: b2 = 4.5
        real ::f,q

        f = 0.5*(1.0+tanh(alpha*(L-L0)))
        q = f*(a1*L + b1) + (1.0-f)*(a2*L + b2)
        D = 10.**q
    end function psphD

    !Calculate invariant latitude (RADIANS) for x,y,z vector (in Rx)
    function InvLatitude(r) result(invlat)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: invlat

        real(rp) :: z,rad,lat,Leq

        z = r(ZDIR)
        rad = norm2(r)

        lat = abs( asin(z/rad) )
        Leq = rad/( cos(lat)*cos(lat) )
        invlat = abs(acos(sqrt(1.0/Leq)))

    end function InvLatitude

    !Calculate dipole L shell for point
    function DipoleL(r) result(Leq)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: Leq

        real(rp) :: z,rad,lat
        z = r(ZDIR)

        rad = norm2(r)
        lat = abs( asin(z/rad) )
        Leq = rad/( cos(lat)*cos(lat) )
    end function DipoleL
    
    !Take point xyz0 and push along dipole to point with radius r
    function DipoleShift(xyz0,r) result(xyz)
        real(rp), intent(in) :: xyz0(NDIM), r
        real(rp), dimension(NDIM) :: xyz

        real(rp) :: L,mlat,mlon
        !Find L of this point
        L = DipoleL(xyz0)

        !Avoid bad values if L<r, push as far as possible
        L = max(L,r)

        !Use r = L*cos^2(latitude)
        mlat = abs(acos(sqrt(r/L)))
        mlon = atan2(xyz0(YDIR),xyz0(XDIR)) !No change in longitude
        if (mlon<0) mlon = mlon+2*PI

        !Get cartesian coordinates
        xyz(XDIR) = r*cos(mlat)*cos(mlon)
        xyz(YDIR) = r*cos(mlat)*sin(mlon)
        xyz(ZDIR) = r*sin(mlat)

    end function DipoleShift

end module earthhelper
