
!Various routines to handle flux-tube calculations to feed inner magnetosphere models
module imagtubes
    use volttypes
    use rcm_mhd_interfaces
    use streamline
    use chmpdbz, ONLY : DipoleB0

	implicit none

    !TODO: Get rid of these ugly globals
    real(rp) :: RIonRCM !Units of Rp
    real(rp) :: Rp_m
    real(rp) :: planetM0g

    logical, parameter :: doNewBeta = .false.

    !Some threshold values for poisoning tubes
    !TODO: Make these XML parameters
    real(rp), private :: wImag_C = 0.20 ![0,1]
    real(rp), private :: bMin_C  = 5.0 !nT

!Information taken from MHD flux tubes
    !Pave = Average pressure [Pa]
    !Nave = Average density [#/m3]
    !Vol  = Flux-tube volume [Re/T]
    !bmin = Min field strength [T]
    !X_bmin = Location of Bmin [m]
    !beta_average = Average plasma beta
    !Potential = MIX potential [Volts]
    !iopen = Field line topology (-1: Closed, 1: Open)
    !Lb = Field line length [Re]
    !rCurv = radius of curvature [Re]
    type RCMTube_T
        real(rp) :: Vol,bmin,beta_average,Pave,Nave,pot
        real(rp) :: X_bmin(NDIM)
        integer(ip) :: iopen
        real(rp) :: latc,lonc !Conjugate lat/lon
        real(rp) :: Lb, Tb !Arc length/bounce time
        real(rp) :: losscone,rCurv,wIMAG
    end type RCMTube_T

    !Parameters used to set RCM ICs to feed back into MHD
    type RCMIC_T
        logical :: doIC = .false.
        real(rp) :: dst0,ktRC !Values for inner magnetosphere
        real(rp) :: vSW,dSW  !Solar wind values used to set plamsa sheet
        real(rp) :: dPS,ktPS !Plasmasheet values
        
    end type RCMIC_T
    
    type(RCMIC_T) :: RCMICs

    contains

!--------------
!MHD=>RCM routines
    !MHD flux-tube
    subroutine MHDTube(vApp,lat,lon,ijTube)
        type(voltApp_T), intent(in) :: vApp
        real(rp), intent(in) :: lat,lon
        type(RCMTube_T), intent(out) :: ijTube

        type(fLine_T) :: bTrc
        real(rp) :: t, bMin,bIon
        real(rp), dimension(NDIM) :: x0, bEq, xyzIon
        real(rp), dimension(NDIM) :: xyzC,xyzIonC
        integer :: OCb
        real(rp) :: bD,bP,dvB,bBeta,rCurv

        real(rp) :: VaMKS,CsMKS,VebMKS !Speeds in km/s
        real(rp) :: TiEV !Temperature in ev
    !First get seed for trace
        !Assume lat/lon @ Earth, dipole push to first cell + epsilon
        xyzIon(XDIR) = RIonRCM*cos(lat)*cos(lon)
        xyzIon(YDIR) = RIonRCM*cos(lat)*sin(lon)
        xyzIon(ZDIR) = RIonRCM*sin(lat)
        x0 = DipoleShift(xyzIon,vApp%mhd2chmp%Rin+TINY)
        bIon = norm2(DipoleB0(xyzIon))*oBScl*1.0e-9 !EB=>T, ionospheric field strength

    !Now do field line trace
        associate(ebModel=>vApp%ebTrcApp%ebModel,ebGr=>vApp%ebTrcApp%ebState%ebGr,ebState=>vApp%ebTrcApp%ebState)

        t = ebState%eb1%time !Time in CHIMP units
        
        call genStream(ebModel,ebState,x0,t,bTrc)
        !Topology
        !OCB =  0 (solar wind), 1 (half-closed), 2 (both ends closed)
        OCb = FLTop(ebModel,ebGr,bTrc)
        
        if (OCb /= 2) then
            !Not closed line, set some values and get out
            ijTube%X_bmin = 0.0
            ijTube%bmin = 0.0
            ijTube%iopen = RCMTOPOPEN
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
        bEq = bEq*Rp_m !Re=>meters

        !Plasma quantities
        !dvB = Flux-tube volume (Re/EB)
        call FLThermo(ebModel,ebGr,bTrc,bD,bP,dvB,bBeta)
        !Converts Re/EB => Re/T
        dvB = dvB/(oBScl*1.0e-9)
        bP = bP*1.0e-9 !nPa=>Pa
        bD = bD*1.0e+6 !#/cc => #/m3

        ijTube%X_bmin = bEq
        ijTube%bmin = bMin
        ijTube%iopen = RCMTOPCLOSED
        ijTube%Vol = dvB
        ijTube%Pave = bP
        ijTube%Nave = bD
        ijTube%beta_average = bBeta

        !Find conjugate lat/lon @ RIonRCM
        call FLConj(ebModel,ebGr,bTrc,xyzC)
        xyzIonC = DipoleShift(xyzC,RIonRCM)
        ijTube%latc = asin(xyzIonC(ZDIR)/norm2(xyzIonC))
        ijTube%lonc = modulo( atan2(xyzIonC(YDIR),xyzIonC(XDIR)),2*PI )
        ijTube%Lb = FLArc(ebModel,ebGr,bTrc)
        ijTube%Tb = FLAlfvenX(ebModel,ebGr,bTrc)
        ijTube%losscone = asin(sqrt(bMin/bIon))

        !Get curvature radius and ExB velocity [km/s]
        call FLCurvRadius(ebModel,ebGr,ebState,bTrc,rCurv,VebMKS)
        ijTube%rCurv = rCurv

    !Get confidence interval
        !VaMKS = flux tube arc length [km] / Alfven crossing time [s]
        VaMKS = (ijTube%Lb*Rp_m*1.0e-3)/ijTube%Tb 
        !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
        TiEV = (1.0e+3)*DP2kT(bD*1.0e-6,bP*1.0e+9) !Temp in eV
        CsMKS = 9.79*sqrt((5.0/3)*TiEV)

        ijTube%wIMAG = VaMKS/( sqrt(VaMKS**2.0 + CsMKS**2.0) + VebMKS)

    !Do different calculation of beta if desired
        if (doNewBeta) then
            bBeta = 2.0*(CsMKS/VaMKS)**2.0 !alternative beta definition
            ijTube%beta_average = bBeta
        endif

        end associate

    end subroutine MHDTube

    !Dipole flux tube info
    subroutine DipoleTube(vApp,lat,lon,ijTube)
        type(voltApp_T), intent(in) :: vApp
        real(rp), intent(in) :: lat,lon
        type(RCMTube_T), intent(out) :: ijTube

        real(rp) :: L,colat
        real(rp) :: mdipole

        mdipole = ABS(planetM0g)*G2T ! dipole moment in T
        colat = PI/2 - lat
        L = 1.0/(sin(colat)**2.0)
        ijTube%Vol = 32./35.*L**4.0/mdipole
        ijTube%X_bmin(XDIR) = L*cos(lon)*Rp_m !Rp=>meters
        ijTube%X_bmin(YDIR) = L*sin(lon)*Rp_m !Rp=>meters
        ijTube%X_bmin(ZDIR) = 0.0
        ijTube%bmin = mdipole/L**3.0

        ijTube%iopen = RCMTOPCLOSED
        
        ijTube%pot = 0.0

        ijTube%beta_average = 0.0
        ijTube%Pave = 0.0
        ijTube%Nave = 0.0

        ijTube%latc = -lat
        ijTube%lonc = lon
        ijTube%Lb   = L !Just lazily using L shell
        ijTube%Tb   = 0.0
        ijTube%losscone = 0.0
        ijTube%rCurv = L/3.0

        ijTube%wIMAG = 1.0 !Much imag
    end subroutine DipoleTube

    !Do some trickkery on the tubes if they seem too weird for RCM
    subroutine TrickyTubes(RCMApp)
        type(rcm_mhd_T), intent(inout) :: RCMApp

        integer :: Ni,Nj,i,j
        logical :: isBadTube

        Ni = RCMApp%nLat_ion
        Nj = RCMApp%nLon_ion
        
        !Loop over grid and poison cells we wouldn't trust RCM with
        do j=1,Nj
            do i=1,Ni
                if (RCMApp%iopen(i,j) /= RCMTOPOPEN) then
                    !Check this cell
                    isBadTube = (RCMApp%wImag(i,j)<=wImag_C) .or. (RCMApp%Bmin(i,j)*1.0e+9 <= bMin_C)
                    
                    if (isBadTube) then
                        !Poison this cell
                        RCMApp%iopen(i,j) = RCMTOPOPEN
                        RCMApp%Vol(i,j) = 0.0 
                    endif
                endif
            enddo !i loop
        enddo !j loop
        
    end subroutine TrickyTubes
    
    !Smooth RCM tube data as needed
    subroutine SmoothTubes(RCMApp)
        type(rcm_mhd_T), intent(inout) :: RCMApp

        integer :: Ni,Nj
        logical, dimension(:,:), allocatable :: isG

        Ni = RCMApp%nLat_ion
        Nj = RCMApp%nLon_ion
        allocate(isG(Ni,Nj))
        isG = .not. (RCMApp%iopen == RCMTOPOPEN)

        
        !Smooth some tubes
        call Smooth2D(RCMApp%Vol) !Flux-tube volume
        call Smooth2D(RCMApp%pot) !Electrostatic potential
        
        contains
        subroutine Smooth2D(Q)
            real(rp), dimension(Ni,Nj), intent(inout) :: Q
            
            real(rp), dimension(Ni,Nj) :: Qs
            integer :: i,j,di,dj,ip,jp
            real(rp), dimension(-2:+2,-2:+2) :: Q55
            logical , dimension(-2:+2,-2:+2) :: G55

            Qs(:,:) = Q

            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,di,dj,ip,jp,Q55,G55)
            do j=1,Nj
                do i=3,Ni-3
                    if (isG(i,j)) then
                        !Pack local 5x5 stencil
                        do dj=-2,+2
                            do di=-2,+2
                                ip = i+di
                                jp = WrapJ(j+dj)
                                Q55(di,dj) = Q  (ip,jp)
                                G55(di,dj) = isG(ip,jp)
                            enddo !di
                        enddo !dj
                        !Calc smoothed value and store
                        Qs(i,j) = SmoothOperator55(Q55,G55)
                    endif !isG(i,j)
                enddo !i
            enddo !j
            !Store values
            Q = Qs
        end subroutine Smooth2D

        !Wrap j index around periodic boundary
        function WrapJ(j) result(jp)
            integer, intent(in) :: j
            integer :: jp
            jp = j
            if (jp>Nj) jp = jp-Nj+1
            if (jp<1 ) jp = jp+Nj-1
        end function WrapJ
    end subroutine SmoothTubes

    !Rewire MHD=>RCM info to set RCM's cold start ICs
    subroutine HackTubes(RCMApp,vApp)
        type(rcm_mhd_T), intent(inout) :: RCMApp
        type(voltApp_T), intent(in)    :: vApp

        integer :: i,j
        real(rp) :: llBC,lat,colat,lon,LPk
        logical :: isLL

        real(rp) :: L,Pmhd,Dmhd,P0_rc,N0_rc,N0_ps,P0_ps,N,P
        !Loop through active region and reset things
        llBC = vApp%mhd2chmp%lowlatBC

        LPk = LPk_QTRC()

        do j=1,RCMApp%nLon_ion
            do i=1,RCMApp%nLat_ion
                !Grab coordinates
                colat = RCMApp%gcolat(i)
                lat = PI/2 - colat
                lon = RCMApp%glong(j)
                
                !Decide if we're below low-lat BC or not
                isLL = (lat <= llBC)

                if (isLL) cycle
                
                !Get L,Pmhd,Dmhd (convert back to our units; nPa,#/cc,Re)
                Pmhd = rcmPScl*RCMApp%Pave(i,j)
                Dmhd = rcmNScl*RCMApp%Nave(i,j)
                L    = norm2( RCMApp%X_bmin(i,j,:) )/Rp_m

            !Quiet-time ring current
                P0_rc = P_QTRC(L)
                !Get density from pressure and target temperature
                N0_rc = PkT2Den(P0_rc,RCMICs%ktRC)
            !Statistical plasma sheet numbers
                N0_ps = RCMICs%dPS
                P0_ps = DkT2P(N0_ps,RCMICs%kTPS)

                if ( P0_ps>P0_rc ) then
                    !Use PS values
                    P = P0_ps
                    N = N0_ps
                else
                    !Use RC values
                    P = P0_rc
                    N = N0_rc
                endif

                !Now test against MHD
                P = max(P,Pmhd)
                N = max(N,Dmhd)

                !Now store them
                RCMApp%Pave(i,j) = P/rcmPScl
                RCMApp%Nave(i,j) = N/rcmNScl

            enddo
        enddo

    end subroutine HackTubes

end module imagtubes
