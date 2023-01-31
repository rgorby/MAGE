!Various kaiju wrappers around geopack things
module kai2geo
    use kdefs
    use geopack
    use planethelper
    use math

    implicit none

    logical, parameter, private :: doSpyro = .false. !Whether to sneak spyro

    contains
    
    !Returns SM-theta,phi of rotation axis
    subroutine GeoAxisTP(axTheta,axPhi)
        real(rp), intent(out) :: axTheta,axPhi

        real(rp), dimension(NDIM) :: xyzGEO,xyzSM

        if (doSpyro) then
            !Sneak the spyro assumption in by setting geo axis to Z-hat
            axPhi = 0.0
            axTheta = 0.0
            return
        endif

        xyzGEO = [0.0,0.0,1.0] !Coordinates of axis in geo
        call GEO2SM(xyzGEO,xyzSM)
        axPhi = atan2(xyzSM(YDIR),xyzSM(XDIR))
        axTheta = acos( xyzSM(ZDIR)/norm2(xyzSM) )

    end subroutine GeoAxisTP

    !Return geocorotation potential [kV] from ionospheric SM colat/lon [radians]
    subroutine geocorotation_from_SMTP(smclatIn,smlon,crpot)
        real(rp),intent(in) :: smclatIn,smlon
        real(rp),intent(out) :: crpot
        
        real(rp) :: r,smlat,xyz(NDIM)

        r =  (RionE*1.0e+6)/REarth !Radius of ionosphere in Re
        smlat = pi/2-smclatIn

        ! convert to XYZ
        xyz(XDIR) = r*cos(smlat)*cos(smlon)
        xyz(YDIR) = r*cos(smlat)*sin(smlon)
        xyz(ZDIR) = r*sin(smlat)

        call geocorotation_from_SMXYZ(xyz,crpot)

    end subroutine geocorotation_from_SMTP

    subroutine geocorotation_from_GEOTP(geoclat,geolon,crpot)
        real(rp), intent(in)  :: geoclat,geolon
        real(rp), intent(out) :: crpot

        real(rp) :: smclat,smlon

        call geoTP2smTP(geoclat,geolon,smclat,smlon)
        call geocorotation_from_SMTP(smclat,smlon,crpot)
        
    end subroutine geocorotation_from_GEOTP

    subroutine geocorotation_from_SMXYZ(xyz,crpot)
        real(rp), intent(in)  :: xyz(NDIM)
        real(rp), intent(out) :: crpot

        real(rp) :: r,smlon,smclat,axT,axP,A0
        !Using Eqn 14 from Hones & Bergeson 1965
        r = norm2(xyz)
        smlon  = atan2(xyz(YDIR),xyz(XDIR))
        smclat = acos (xyz(ZDIR)/norm2(xyz))
        !Get theta/phi of axis
        call GeoAxisTP(axT,axP)

        !NOTE: Negative sign gets absorbed into order of A0 x (...)
        A0 = EarthPsi0/r
        crpot = A0*( sin(smclat)*cos(smclat)*cos(smlon-axP)*sin(axT) - cos(axT)*sin(smclat)*sin(smclat) )

    end subroutine geocorotation_from_SMXYZ

    subroutine smTP2geoTP(incolat,inlon,outcolat,outlon)
        real(rp), intent(in)  :: incolat,inlon
        real(rp), intent(out) :: outcolat,outlon

        real(rp) :: r,tempt,tempp
        real(rp),dimension(NDIM) :: indim,outdim

        ! Calculate radius of ionosphere in km
        r =  (RionE*1.0e+6)/REarth

        ! temporarily convert to lat
        tempt = pi/2-incolat
      
        ! use temp values everywhere for consistency
        tempp = inlon

        ! convert to XYZ    
        indim(XDIR) = r*cos(tempt)*cos(tempp)
        indim(YDIR) = r*cos(tempt)*sin(tempp)
        indim(ZDIR) = r*sin(tempt)

        ! Do transform
        call SM2GEO(indim,outdim)
        ! Convert back to colatitude
        ! I'm going to be lazy here and not even care to convert to funky colat output
        outcolat = acos(outdim(ZDIR)/r) ! convert back to standard colatitude
        outlon = modulo((atan2(outdim(YDIR),outdim(XDIR))+2*pi),(2*pi)) 
    end subroutine smTP2geoTP

    subroutine geoTP2smTP(geocolat,geolon,smcolat,smlon)
        real(rp), intent(in)  :: geocolat,geolon
        real(rp), intent(out) :: smcolat,smlon

        real(rp) :: r,geolat
        real(rp), dimension(NDIM) :: geoXYZ,smXYZ
        r =  (RionE*1.0e+6)/REarth !Radius of ionosphere in Re

        geolat = PI/2 - geocolat
        !Convert to xyz
        geoXYZ(XDIR) = r*cos(geolat)*cos(geolon)
        geoXYZ(YDIR) = r*cos(geolat)*sin(geolon)
        geoXYZ(ZDIR) = r*sin(geolat)

        !Do transform
        call GEO2SM(geoXYZ,smXYZ)

        !Convert back to TP
        smcolat = acos(smXYZ(ZDIR)/r)
        smlon   = katan2(smXYZ(YDIR),smXYZ(XDIR))

    end subroutine geoTP2smTP

    !Take the signed invariant latitude and calculated the mirror ratio between R=1 and R=rad
    !Using dipole mapping across the gap
    function IGRFMirrorRatio(mlat,mlon,rad) result(RM)
        real(rp), intent(in) :: mlat,mlon,rad
        real(rp) :: RM

        real(rp) :: rIon
        real(rp), dimension(NDIM) :: xyzIon,xyzMir

        rIon = 1.0
        xyzIon(XDIR) = rIon*cos(mlat)*cos(mlon)
        xyzIon(YDIR) = rIon*cos(mlat)*sin(mlon)
        xyzIon(ZDIR) = rIon*sin(mlat)

        xyzMir = DipoleShift(xyzIon,rad)
        RM = norm2( IGRF_SM(xyzIon) )/norm2( IGRF_SM (xyzMir) )
    end function IGRFMirrorRatio
end module kai2geo
