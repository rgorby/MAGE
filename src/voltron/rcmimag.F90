!Routines to handle RCM inner magnetosphere model
!NOTES: 
!-Figure out flux-tube volume units
!-Work on upating legacy Fortran
!-Work on OMP bindings
!-Streamline console noise

module rcmimag
    use volttypes
    use files
    use earthhelper
    use rcm_mhd_interfaces
    use rcm_mix_interface
    use streamline
    use clocks
    use rcm_mhd_mod, ONLY : rcm_mhd
    use rcm_mhd_io
    use msphutils, only : MagMoment

    use cmiutils, only : SquishCorners
    
    implicit none

    real(rp) :: RIonRCM !Units of Rp
    real(rp), private :: rEqMin = 0.0
    real(rp), private :: PPDen = 50.0 !Plasmapause density
    character(len=strLen), private :: h5File

    real(rp), private :: Rp_m
    real(rp), private :: planetM0g

    logical, parameter, private :: doKillRCMDir = .true. !Whether to always kill RCMdir before starting

    !Information taken from MHD flux tubes
    !TODO: Figure out RCM boundaries

    !Pave = Average pressure [Pa]
    !Nave = Average density [#/m3]
    !Vol  = Flux-tube volume [Re/T]
    !bmin = Min field strength [T]
    !X_bmin = Location of Bmin [m]
    !beta_average = Average plasma beta
    !Potential = MIX potential [Volts]
    !iopen = Field line topology (-1: Closed, 1: Open)
    !Lb = Field line length [Re]
    type RCMTube_T
        real(rp) :: Vol,bmin,beta_average,Pave,Nave,pot
        real(rp) :: X_bmin(NDIM)
        integer(ip) :: iopen
        real(rp) :: latc,lonc !Conjugate lat/lon
        real(rp) :: Lb
    end type RCMTube_T

    real(rp), dimension(:,:), allocatable, private :: mixPot

    !Parameters for smoothing toMHD boundary
    type SmoothOperator_T
        integer :: nIter,nRad
    end type SmoothOperator_T

    type(SmoothOperator_T) :: SmoothOp

    type, extends(innerMagBase_T) :: rcmIMAG_T

        ! rcm coupling variable
        type(rcm_mhd_T) :: rcmCpl

        contains

        ! over-ride the base functions with RCM versions
        procedure :: doInit => initRCM
        procedure :: doAdvance => advanceRCM
        procedure :: doEval => EvalRCM
        procedure :: doIO => doRCMIO
        procedure :: doRestart => doRCMRestart

    end type

    contains

    !Initialize RCM inner magnetosphere model
    subroutine initRCM(imag,iXML,isRestart,rad_planet_m,rad_iono_m,M0g,vApp)
        class(rcmIMAG_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        real(rp), intent(in) :: rad_planet_m,rad_iono_m, M0g ! Specific planet aprameters
        type(voltApp_T), intent(inout) :: vApp

        character(len=strLen) :: RunID
        real(rp) :: t0

        associate(RCMApp => imag%rcmCpl, & !type rcm_mhd_T
                  imag2mix => vApp%imag2mix, &
                  t0 => vApp%time, &
                  dtCpl => vApp%DeepDT)
        !Set radii in RCMApp
        RCMApp%planet_radius = rad_planet_m
        RCMApp%iono_radius = rad_iono_m
        Rp_m = rad_planet_m ! for local use
        RIonRCM = rad_iono_m/rad_planet_m

        planetM0g = M0g
        
        call iXML%Set_Val(RunID,"/gamera/sim/runid","sim")
        RCMApp%rcm_runid = trim(RunID)

        if (isRestart) then
            if (doKillRCMDir) then
                !Kill RCMFiles directory even on restart
                call KillRCMDir()
            endif

            !Get t0 and nRes necessary for RCM restart
            call RCMRestartInfo(RCMApp,iXML,t0)

            write(*,*) 'Restarting RCM @ t = ', t0
            vApp%time = t0 !Set vApp's time to correct value from restart
            call rcm_mhd(t0,dtCpl,RCMApp,RCMRESTART,iXML=iXML)
            coldstart = .false. ! set to false is it is a restart
        else
            t0 = vApp%time
            call KillRCMDir()
            write(*,*) 'Initializing RCM ...'
            call rcm_mhd(t0,dtCpl,RCMApp,RCMINIT,iXML=iXML)
        endif

        call init_rcm_mix(RCMApp,imag2mix)

        !Start up IO
        call initRCMIO(RCMApp,isRestart)

        call iXML%Set_Val(SmoothOp%nIter,"imag/nIter",4)
        call iXML%Set_Val(SmoothOp%nRad ,"imag/nRad" ,8)

        end associate

        contains
            subroutine KillRCMDir()
                write(*,*) 'Nuking RCMfiles/ ...'
                CALL SYSTEM("rm -rf RCMfiles > /dev/null 2>&1")
            end subroutine KillRCMDir
    end subroutine initRCM

    !Advance RCM from Voltron data
    subroutine AdvanceRCM(imag,vApp,tAdv)
        class(rcmIMAG_T), intent(inout) :: imag
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        integer :: i,j,n,nStp
        real(rp) :: colat,lat,lon
        real(rp) :: dtAdv
        type(RCMTube_T) :: ijTube

        real(rp) :: llBC,maxRad
        logical :: isLL

        associate(RCMApp => imag%rcmCpl)

        !Lazily grabbing rDeep here, convert to RCM units
        !rEqMin = vApp%rDeep*Re_cgs*1.0e-2 !Re=>meters
            rEqMin = vApp%rDeep*Rp_m !Re=>meters

        llBC = vApp%mhd2chmp%lowlatBC

        call Tic("MAP_RCMMIX")
    !Get potential from mix
        call map_rcm_mix(vApp,mixPot)
        call Toc("MAP_RCMMIX")

        call Tic("RCM_TUBES")
    !Load RCM tubes
       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP schedule(guided) &
       !$OMP private(i,j,colat,lat,lon,isLL,ijTube)
        do j=1,RCMApp%nLon_ion
            do i=1,RCMApp%nLat_ion
                colat = RCMApp%gcolat(i)
                lat = PI/2 - colat
                lon = RCMApp%glong(j)
                
                !Decide if we're below low-lat BC or not
                isLL = (lat <= llBC)

                if (isLL) then
                    !Use mocked up values
                    call DipoleTube(vApp,lat,lon,ijTube)
                else
                    !Trace through MHD
                    call MHDTube(vApp,lat,lon,ijTube)
                endif

                !Stuff data into RCM
                RCMApp%Vol(i,j)          = ijTube%Vol
                RCMApp%bmin(i,j)         = ijTube%bmin
                RCMApp%iopen(i,j)        = ijTube%iopen
                RCMApp%beta_average(i,j) = ijTube%beta_average
                RCMApp%Pave(i,j)         = ijTube%Pave
                RCMApp%Nave(i,j)         = ijTube%Nave
                RCMApp%X_bmin(i,j,:)     = ijTube%X_bmin

                RCMApp%latc(i,j)         = ijTube%latc
                RCMApp%lonc(i,j)         = ijTube%lonc
                RCMApp%Lb(i,j)           = ijTube%Lb
                RCMApp%Tb(i,j)           = AlfvenBounce(ijTube%Nave,ijTube%bmin,ijTube%Lb)
                !mix variables are stored in this order (longitude,colatitude), hence the index flip
                RCMApp%pot(i,j)          = mixPot(j,i)
            enddo
        enddo
        call Toc("RCM_TUBES")

        call Tic("AdvRCM")
    !Advance from vApp%time to tAdv
        dtAdv = tAdv-vApp%time !RCM-DT
        if(coldstart)then
         call rcm_mhd(vApp%time,dtAdv,RCMApp,RCMCOLDSTART)
         coldstart=.false.
        else
         call rcm_mhd(vApp%time,dtAdv,RCMApp,RCMADVANCE)
        end if
        !Update timming data - not needed, frt
!        call rcm_mhd(vApp%time,0.0_rp,RCMApp,RCMWRITETIMING)
        call Toc("AdvRCM")

        !Set ingestion region
        call SetIngestion(RCMApp)

    !Pull data from RCM state for conductance calculations
        vApp%imag2mix%isClosed = (RCMApp%iopen == RCMTOPCLOSED)
        vApp%imag2mix%latc = RCMApp%latc
        vApp%imag2mix%lonc = RCMApp%lonc
        vApp%imag2mix%eflux = RCMApp%flux
        vApp%imag2mix%eavg  = RCMApp%eng_avg

        vApp%imag2mix%iflux = 0.0
        vApp%imag2mix%iavg  = 0.0

        vApp%imag2mix%isFresh = .true.

    !Find maximum extent of closed field region
        maxRad = maxval(norm2(RCMApp%X_bmin,dim=3),mask=vApp%imag2mix%isClosed)
        !maxRad = maxRad/(Re_cgs*1.0e-2)
        maxRad = maxRad/Rp_m
        vApp%rTrc = 1.25*maxRad

        end associate        

        contains
            !Calculate Alfven bounce timescale
            !D = #/m3, B = T, L = Rp
            function AlfvenBounce(D,B,L) result(dTb)
                real(rp), intent(in) :: D,B,L
                real(rp) :: dTb

                real(rp) :: Va,nCC,bNT

                if ( (D<TINY) .or. (L<TINY) ) then
                    dTb = 0.0
                    return
                endif
                nCC = D*rcmNScl !Get n in #/cc
                bNT = B*1.0e+9 !Convert B to nT
                Va = 22.0*bNT/sqrt(nCC) !km/s, from NRL plasma formulary
                !dTb = (L*Re_km)/Va
                dTb = (L*Rp_m*1.0e-3)/Va
            end function AlfvenBounce
    end subroutine AdvanceRCM

    !Set region of RCM grid that's "good" for MHD ingestion
    subroutine SetIngestion(RCMApp)
        type(rcm_mhd_T), intent(inout) :: RCMApp
        integer :: n,i,j,iC

        integer , dimension(:), allocatable :: jBnd
        real, dimension(:), allocatable :: jRad,jRadG
        integer :: NSmth,NRad
        real(rp) :: RadC,rIJ

        NSmth = SmoothOp%nIter
        NRad  = SmoothOp%nRad

        allocate(jBnd (  RCMApp%nLon_ion  ))
        
        allocate(jRad (  RCMApp%nLon_ion  ))
        allocate(jRadG(1-NRad:RCMApp%nLon_ion+NRad))

        do j=1,RCMApp%nLon_ion
            do i = 1,RCMApp%nLat_ion
                if (RCMApp%toMHD(i,j) .and. (RCMApp%iopen(i,j) == RCMTOPCLOSED)) then
                    exit
                endif
            enddo
            jBnd(j) = min(i+1,RCMApp%nLat_ion)
            jRad(j) = norm2( RCMApp%X_bmin(jBnd(j),j,XDIR:YDIR) )
        enddo

        do n=1,NSmth
            jRadG(1:RCMApp%nLon_ion) = jRad
            jRadG(1-NRad:0) = jRad(RCMApp%nLon_ion-NRad+1:RCMApp%nLon_ion)
            jRadG(RCMApp%nLon_ion+1:RCMApp%nLon_ion+NRad) = jRad(1:NRad)
            do j=1,RCMApp%nLon_ion
                !Take mean over range
                jRad(j) = sum(jRadG(j-NRad:j+NRad))/(2.0*NRad+1)
                !jRad(j) = (product(jRadG(j-NRad:j+NRad)))**(1.0/(2.0*NRad+1))
                !jRad(j) = minval(jRadG(j-NRad:j+NRad))
            enddo
        enddo

        RCMApp%toMHD = .false.
        do j=1,RCMApp%nLon_ion
            RadC = jRad(j)
            do i = jBnd(j),RCMApp%nLat_ion
                rIJ = norm2(RCMApp%X_bmin(i,j,XDIR:YDIR))
                if ( (rIJ <= RadC) .and. (RCMApp%iopen(i,j) == RCMTOPCLOSED) ) then
                    RCMApp%toMHD(i,j) = .true.
                else
                    RCMApp%toMHD(i,j) = .false.
                endif
            enddo
        enddo

    end subroutine SetIngestion

    !Evaluate eq map at a given point
    !Returns density (#/cc) and pressure (nPa)
    subroutine EvalRCM(imag,x1,x2,t,imW,isEdible)
        class(rcmIMAG_T), intent(inout) :: imag
        real(rp), intent(in) :: x1,x2,t
        real(rp), intent(out) :: imW(NVARIMAG)
        logical, intent(out) :: isEdible

        real(rp) :: colat,nrcm,prcm,npp,ntot
        integer, dimension(2) :: ij0

        associate(RCMApp => imag%rcmCpl, lat => x1, lon => x2)

        !Set defaults
        imW(:) = 0.0
        imW(IMDEN ) = 0.0
        imW(IMPR  ) = 0.0
        imW(IMTSCL) = 0.0
        isEdible = .false.

        colat = PI/2 - lat

        !Do 1st short cut tests
        isEdible =  (colat >= RCMApp%gcolat(1)) .and. (colat <= RCMApp%gcolat(RCMApp%nLat_ion)) &
                    .and. (lat > TINY)

        if (.not. isEdible) return

        !If still here, find mapping (i,j) on RCM grid of point
        call GetRCMLoc(lat,lon,ij0)

        !Do second short cut tests
        isEdible = RCMApp%toMHD(ij0(1),ij0(2))
        if (.not. isEdible) return

        prcm = rcmPScl*RCMApp%Prcm (ij0(1),ij0(2))
        nrcm = rcmNScl*RCMApp%Nrcm (ij0(1),ij0(2))
        npp  = rcmNScl*RCMApp%Npsph(ij0(1),ij0(2))

        ntot = 0.0
        !Decide which densities to include
        if (npp >= PPDen) then
            ntot = ntot + npp
        endif
        if ( (nrcm>TINY) .and. (prcm>TINY) ) then
            ntot = ntot + nrcm
        endif

        !Store data
        imW(IMDEN)  = ntot
        imW(IMPR)   = prcm
        imW(IMTSCL) = RCMApp%Tb(ij0(1),ij0(2))
        imW(IMX1)   = (180.0/PI)*lat
        imW(IMX2)   = (180.0/PI)*lon

        end associate

        contains

        subroutine GetRCMLoc(lat,lon,ij0)
            real(rp), intent(in) :: lat,lon
            integer, intent(out) :: ij0(2)

            integer :: iX,jX,iC,n
            real(rp) :: colat,dp,dcol,dI,dJ

            associate(gcolat=>imag%rcmCpl%gcolat,glong=>imag%rcmCpl%glong, &
                      nLat=>imag%rcmCpl%nLat_ion,nLon=>imag%rcmCpl%nLon_ion)

            !Assuming constant lon spacing
            dp = glong(2) - glong(1)

            !Get colat point
            colat = PI/2 - lat
!            iC = findloc(gcolat >= colat,.true.,dim=1) - 1
! bypass as findloc does not work for gfortran<9
           !  Work-around code        
           do n=1,nLat
            if (gcolat(n) > colat) exit
            enddo
            iC = n-1
            dcol = gcolat(iC+1)-gcolat(iC)
            dI = (colat-gcolat(iC))/dcol
            if (dI <= 0.5) then
                iX = iC
            else
                iX = iC+1
            endif

            !Get lon point
            dJ = lon/dp
            if ( (dJ-floor(dJ)) <= 0.5 ) then
                jX = floor(dJ)+1
            else
                jX = floor(dJ)+2
            endif

            !Impose bounds just in case
            iX = max(iX,1)
            iX = min(iX,nLat)
            jX = max(jX,1)
            jX = min(jX,nLon)

            ij0 = [iX,jX]

            end associate
        end subroutine GetRCMLoc

    end subroutine EvalRCM

!--------------
!MHD=>RCM routines
    !MHD flux-tube
    subroutine MHDTube(vApp,lat,lon,ijTube)
        type(voltApp_T), intent(in) :: vApp
        real(rp), intent(in) :: lat,lon
        type(RCMTube_T), intent(out) :: ijTube

        type(fLine_T) :: bTrc
        real(rp) :: t, bMin
        real(rp), dimension(NDIM) :: x0, bEq, xyzIon
        real(rp), dimension(NDIM) :: xyzC,xyzIonC
        integer :: OCb
        real(rp) :: bD,bP,dvB,bBeta

    !First get seed for trace
        !Assume lat/lon @ Earth, dipole push to first cell
        xyzIon(XDIR) = RIonRCM*cos(lat)*cos(lon)
        xyzIon(YDIR) = RIonRCM*cos(lat)*sin(lon)
        xyzIon(ZDIR) = RIonRCM*sin(lat)
        x0 = DipoleShift(xyzIon,vApp%mhd2chmp%Rin)
        
    !Now do field line trace
        associate(ebModel=>vApp%ebTrcApp%ebModel,ebGr=>vApp%ebTrcApp%ebState%ebGr,ebState=>vApp%ebTrcApp%ebState)

        t = ebState%eb1%time !Time in CHIMP units
        call genStream(ebModel,ebState,x0,t,bTrc)

    !Get diagnostics from field line
        !Minimal surface (bEq in Rp, bMin in EB)
        call FLEq(ebModel,bTrc,bEq,bMin)
        bMin = bMin*oBScl*1.0e-9 !EB=>Tesla
        !bEq = bEq*Re_cgs*1.0e-2 !Re=>meters
        bEq = bEq*Rp_m !Re=>meters

        !Plasma quantities
        !dvB = Flux-tube volume (Re/EB)
        call FLThermo(ebModel,ebGr,bTrc,bD,bP,dvB,bBeta)
        !Converts Re/EB => Re/T
        dvB = dvB/(oBScl*1.0e-9)
        bP = bP*1.0e-9 !nPa=>Pa
        bD = bD*1.0e+6 !#/cc => #/m3
        !Topology
        !OCB =  0 (solar wind), 1 (half-closed), 2 (both ends closed)
        OCb = FLTop(ebModel,ebGr,bTrc)
  
    !Scale and store information
        if (OCb == 2) then
            !Closed field line
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
        else
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

        !mdipole = EarthM0g*G2T ! dipole moment in T
        mdipole = ABS(planetM0g)*G2T ! dipole moment in T
        colat = PI/2 - lat
        L = 1.0/(sin(colat)**2.0)
        ijTube%Vol = 32./35.*L**4.0/mdipole
        !ijTube%X_bmin(XDIR) = L*cos(lon)*Re_cgs*1.0e-2 !Re=>meters
        !ijTube%X_bmin(YDIR) = L*sin(lon)*Re_cgs*1.0e-2 !Re=>meters
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

        !ijTube%Nave = psphD(L)*1.0e+6 !#/cc => #/m3

    end subroutine DipoleTube

!IO wrappers
    subroutine doRCMIO(imag,nOut,MJD,time)
        class(rcmIMAG_T), intent(inout) :: imag
        integer, intent(in) :: nOut
        real(rp), intent(in) :: MJD,time

        call WriteRCM(imag%rcmCpl,nOut,MJD,time)
    end subroutine doRCMIO

    subroutine doRCMRestart(imag,nRes,MJD,time)
        class(rcmIMAG_T), intent(inout) :: imag
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time

        call WriteRCMRestart(imag%rcmCpl,nRes,MJD,time)
    end subroutine doRCMRestart
    
end module rcmimag
