
!Various routines to handle flux-tube calculations to feed inner magnetosphere models
module imagtubes
    use volttypes
    use streamline
    use chmpdbz, ONLY : DipoleB0
    use imaghelper
    use planethelper

	implicit none

    type IMAGTube_T
        real(rp) :: Vol,bmin,beta_average,Pave,Nave
            !! Flux tube volume, minimum b along FL
            !! Average plasma beta, average pressure, average density
        real(rp) :: X_bmin(NDIM)
        integer(ip) :: topo
        real(rp) :: latc,lonc !Conjugate lat/lon
        real(rp) :: Lb, Tb !Arc length/bounce time
        real(rp) :: losscone,rCurv,wIMAG,TioTe0=4.0_rp
    end type IMAGTube_T

    contains

    subroutine MHDTube(ebTrcApp,planet,colat,lon,r,ijTube,bTrc,nTrcO,doShiftO)
        type(ebTrcApp_T), intent(in) :: ebTrcApp
        type(planet_T), intent(in) :: planet
        real(rp), intent(in) :: colat, lon, r
        type(IMAGTube_T), intent(out) :: ijTube
        type(magLine_T), intent(inout) :: bTrc
        integer, intent(in), optional :: nTrcO
        logical, intent(in), optional :: doShiftO

        real(rp) :: t, bMin,bIon
        real(rp), dimension(NDIM) :: x0, bEq, xyzIon
        real(rp), dimension(NDIM) :: xyzC,xyzIonC
        integer :: OCb
        real(rp) :: bD,bP,dvB,bBeta,rCurv
        real(rp) :: VaMKS,CsMKS,VebMKS !Speeds in km/s
        real(rp) :: TiEV !Temperature in ev

        logical :: doShift,doShue

        if (present(doShiftO)) then
            doShift = doShiftO
        else
            doShift = .false.
        endif

        ! Take seed point in spherical coordinates
        xyzIon(XDIR) = r*sin(colat)*cos(lon)
        xyzIon(YDIR) = r*sin(colat)*sin(lon)
        xyzIon(ZDIR) = r*cos(colat)
        if (doShift) then
            !! TODO: Decide if this is how we want to do it
            x0 = DipoleShift(xyzIon, planet%ri_m/planet%rp_m + TINY)
        else
            x0 = xyzIon
        endif
        
        bIon = norm2(DipoleB0(xyzIon))*oBScl*1.0e-9 !EB=>T, ionospheric field strength


        associate(ebModel=>ebTrcApp%ebModel,ebGr=>ebTrcApp%ebState%ebGr,ebState=>ebTrcApp%ebState)

    !Now do field line trace
        t = ebState%eb1%time !Time in CHIMP units
            !!TODO: What do we do when we want times in between steps? Somethign similar to what slice/chop do to output
        
        if (present(nTrcO)) then
            call genLine(ebModel,ebState,x0,t,bTrc,nTrcO,doShueO=.true.,doNHO=.true.)
        else
            call genLine(ebModel,ebState,x0,t,bTrc,      doShueO=.true.,doNHO=.true.)
        endif


        !Topology
        !OCB =  0 (solar wind), 1 (half-closed), 2 (both ends closed)
        OCb = FLTop(ebModel,ebGr,bTrc)
        ijTube%topo = OCb
        
        if (OCb /= 2) then
            ! If the field line hasn't closed after 15 minutes we are legally allowed to leave
            ijTube%X_bmin = 0.0
            ijTube%bmin = 0.0
            ijTube%Vol = -1.0
            ijTube%Pave = 0.0
            ijTube%Nave = 0.0
            ijTube%beta_average = 0.0
            ijTube%latc = 0.0
            ijTube%lonc = 0.0
            ijTube%Lb   = 0.0
            ijTube%Tb   = 0.0
            ijTube%losscone = 0.0
            ijTube%rCurv = 0.0
            ijTube%wIMAG = 0.0
            return
        endif

    !Get diagnostics from closed field line
        !Minimal surface (bEq in Rp, bMin in EB)
        call FLEq(ebModel,bTrc,bEq,bMin)
        bMin = bMin*oBScl*1.0e-9 !EB=>Tesla
        bEq = bEq*planet%rp_m ! Rp -> meters

        !Plasma quantities
        !dvB = Flux-tube volume (Re/EB)
        call FLThermo(ebModel,ebGr,bTrc,bD,bP,dvB,bBeta)
        !Converts Rp/EB => Rp/T
        dvB = dvB/(oBScl*1.0e-9)
        bP = bP*1.0e-9 !nPa=>Pa
        bD = bD*1.0e+6 !#/cc => #/m3

        ijTube%X_bmin = bEq
        ijTube%bmin = bMin
        ijTube%Vol = dvB
        ijTube%Pave = bP
        ijTube%Nave = bD
        ijTube%beta_average = bBeta

        call FLConj(ebModel,ebGr,bTrc,xyzC)
        if (doShift) then
            xyzIonC = DipoleShift(xyzC, planet%ri_m)
        else
            xyzIonC = xyzC
        endif
        ijTube%latc = asin(xyzIonC(ZDIR)/norm2(xyzIonC))
        ijTube%lonc = modulo( atan2(xyzIonC(YDIR),xyzIonC(XDIR)),2*PI )
        ijTube%Lb = FLArc(ebModel,ebGr,bTrc)
        !NOTE: Bounce timescale may be altered to use IMAG hot density
        ijTube%Tb = FLAlfvenX(ebModel,ebGr,bTrc)
        !ijTube%Tb = FLFastX(ebModel,ebGr,bTrc)

        !write(*,*) 'Bounce compare: = ', FLFastX(ebModel,ebGr,bTrc)/FLAlfvenX(ebModel,ebGr,bTrc)
        
        ijTube%losscone = asin(sqrt(bMin/bIon))

        !Get curvature radius and ExB velocity [km/s]
        call FLCurvRadius(ebModel,ebGr,ebState,bTrc,rCurv,VebMKS)
        ijTube%rCurv = rCurv

    !Get confidence interval
        !VaMKS = flux tube arc length [km] / Alfven crossing time [s]
        VaMKS = (ijTube%Lb*planet%rp_m*1.0e-3)/ijTube%Tb 
        !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
        TiEV = (1.0e+3)*DP2kT(bD*1.0e-6,bP*1.0e+9) !Temp in eV
        CsMKS = 9.79*sqrt((5.0/3)*TiEV)

        ijTube%wIMAG = VaMKS/( sqrt(VaMKS**2.0 + CsMKS**2.0) + VebMKS)

        end associate

    end subroutine MHDTube

end module imagtubes
