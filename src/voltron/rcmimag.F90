!Routines to handle RCM inner magnetosphere model
!NOTES: 
!-Figure out flux-tube volume units
!-Work on upating legacy Fortran
!-Work on OMP bindings
!-Streamline console noise

module rcmimag
    use volttypes
    use ioh5
    use files
    use earthhelper
    use rcm_mhd_interfaces
    use rcm_mix_interface
    use streamline
    use clocks

    implicit none

    type(rcm_mhd_T), private :: RCMApp

    !Scaling parameters
    real(rp), parameter, private :: rcmPScl = 1.0e+9 !Convert Pa->nPa
    real(rp), parameter, private :: rcmNScl = 1.0e-6 !Convert #/m3 => #/cc
    real(rp), parameter, private :: IMGAMMA = 5.0/3.0

    real(rp), parameter :: RIonRCM = (RionE/REarth)*1.0e+6 !Units of Re
    real(rp), private :: rEqMin = 0.0
    real(rp), private :: PPDen = 50.0 !Plasmapause density
    integer, parameter :: MAXRCMIOVAR = 30
    character(len=strLen), private :: h5File


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
    type RCMTube_T
        real(rp) :: Vol,bmin,beta_average,Pave,Nave,pot
        real(rp) :: X_bmin(NDIM)
        integer(ip) :: iopen
        real(rp) :: latc,lonc !Conjugate lat/lon
    end type RCMTube_T

    type RCMEllipse_T
        !Ellipse parameters (in m)
        real(rp) :: xSun,xTail,yDD

        !Safety parameters
        real(rp) :: xSScl,xTScl,yScl
        logical  :: isDynamic !Whether to update parameters
    end type RCMEllipse_T

    real(rp), dimension(:,:), allocatable, private :: mixPot
    type(RCMEllipse_T) :: RCMEll

    contains

    !Initialize RCM inner magnetosphere model
    subroutine initRCM(iXML,isRestart,imag2mix,t0,dtCpl,nRes)
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        type(imag2Mix_T), intent(inout) :: imag2mix
        real(rp), intent(in) :: t0,dtCpl
        integer, intent(in), optional :: nRes

        character(len=strLen) :: RunID,RCMH5
        logical :: fExist

        call iXML%Set_Val(RunID,"/gamera/sim/runid","sim")
        RCMApp%rcm_runid = trim(RunID)

        if (isRestart) then
            RCMApp%rcm_nRes = nRes
            write(*,*) 'Restarting RCM @ t = ', t0
            call rcm_mhd(t0,dtCpl,RCMApp,RCMRESTART)
        else
            write(*,*) 'Initializing RCM ...'
            call rcm_mhd(t0,dtCpl,RCMApp,RCMINIT)
        endif
        call init_rcm_mix(RCMApp,imag2mix)

        h5File = trim(RunID) // ".mhdrcm.h5" !MHD-RCM coupling data
        RCMH5  = trim(RunID) // ".rcm.h5" !RCM data

        fExist = CheckFile(h5File)
        write(*,*) 'RCM outputting to ',trim(h5File)
        if ( (.not. isRestart) .or. (isRestart .and. (.not.fExist)) ) then
            !Not a restart or it is a restart and no file
            call CheckAndKill(h5File) !For non-restart but file exists
            call CheckAndKill(RCMH5)

            !Create base file
            call initRCMIO()
        endif
        
    !Get ellipse quantities
        !Get values in Re, then convert to meters
        call iXML%Set_Val(RCMEll%xSun ,"/RCM/ellipse/xSun" ,8.0)
        call iXML%Set_Val(RCMEll%xTail,"/RCM/ellipse/xTail",8.0)
        call iXML%Set_Val(RCMEll%yDD  ,"/RCM/ellipse/yDD"  ,8.0)
        !Scale
        RCMEll%xSun  = REarth*RCMEll%xSun
        RCMEll%xTail = REarth*RCMEll%xTail
        RCMEll%yDD   = REarth*RCMEll%yDD
        call iXML%Set_Val(RCMEll%isDynamic,"/RCM/ellipse/isDynamic"  ,.true.)
        !Get safety parameters (only for dynamic ellipse)
        call iXML%Set_Val(RCMEll%xSScl ,"/RCM/ellipse/xSScl" ,0.70)
        call iXML%Set_Val(RCMEll%xTScl ,"/RCM/ellipse/xTScl" ,0.90)
        call iXML%Set_Val(RCMEll%yScl  ,"/RCM/ellipse/yScl"  ,0.70)

    end subroutine initRCM

    !Advance RCM from Voltron data
    subroutine AdvanceRCM(vApp,tAdv)
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        integer :: i,j,n,nStp
        real(rp) :: colat,lat,lon
        real(rp) :: dtAdv
        type(RCMTube_T) :: ijTube

        real(rp) :: llBC,maxRad
        logical :: isLL

        !Lazily grabbing rDeep here, convert to RCM units
        rEqMin = vApp%rDeep*Re_cgs*1.0e-2 !Re=>meters

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
        do i=1,RCMApp%nLat_ion
            do j=1,RCMApp%nLon_ion
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

                !mix variables are stored in this order (longitude,colatitude), hence the index flip
                RCMApp%pot(i,j)          = mixPot(j,i)
                
            enddo
        enddo

        !Set ingestion region
        call SetIngestion()

        call Toc("RCM_TUBES")

        call Tic("AdvRCM")
    !Advance from vApp%time to tAdv
        dtAdv = tAdv-vApp%time !RCM-DT
        call rcm_mhd(vApp%time,dtAdv,RCMApp,RCMADVANCE)
        !Update timming data
        call rcm_mhd(vApp%time,0.0_rp,RCMApp,RCMWRITETIMING)
        call Toc("AdvRCM")

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
        maxRad = maxRad/(Re_cgs*1.0e-2)
        vApp%rTrc = 1.05*maxRad
        
    end subroutine AdvanceRCM

    !Set region of RCM grid that's "good" for MHD ingestion
    subroutine SetIngestion()
        integer :: i,j,iC
        
        real(rp) :: x0,ell,a,b
        logical :: jClosed

        !Start by looping from high-lat downwards until we find full ring of closed lines
        do i = 1,RCMApp%nLat_ion
            jClosed = all(RCMApp%iopen(i,:) == RCMTOPCLOSED)
            if (jClosed) exit
        enddo
        iC = i

        !Construct new ellipse if we're doing dynamic
        if (RCMEll%isDynamic) then
            !Now find maximum sunward point on this ring
            RCMEll%xSun  = maxval(RCMApp%X_bmin(iC,:,XDIR))
            RCMEll%xTail = minval(RCMApp%X_bmin(iC,:,XDIR))
            RCMEll%yDD   = maxval(abs(RCMApp%X_bmin(iC,:,YDIR)))

            ! RCMEll%xSun  = maxval(    RCMApp%X_bmin(:,:,XDIR) ,mask=(RCMApp%iopen(:,:) == RCMTOPCLOSED) )
            ! RCMEll%xTail = minval(    RCMApp%X_bmin(:,:,XDIR) ,mask=(RCMApp%iopen(:,:) == RCMTOPCLOSED) )
            ! RCMEll%yDD   = maxval(abs(RCMApp%X_bmin(:,:,YDIR)),mask=(RCMApp%iopen(:,:) == RCMTOPCLOSED) )

            !Rescale to give some breathing room
            RCMEll%xSun  = RCMEll%xSScl*RCMEll%xSun 
            RCMEll%xTail = RCMEll%xTScl*RCMEll%xTail
            RCMEll%yDD   = RCMEll%yScl *RCMEll%yDD  

            !Constrain ellipse parameters by reqmin
            !TODO: Test unconstrained tail
            if (    RCMEll%xSun   >= rEqMin) RCMEll%xSun  =  rEqMin
            if (abs(RCMEll%xTail) >= rEqMin) RCMEll%xTail = -rEqMin
            if (    RCMEll%yDD    >= rEqMin) RCMEll%yDD   =  rEqMin
            ! write(*,*) 'iC = ', iC
            ! write(*,*) 'xSun  = ', RCMEll%xSun/REarth
            ! write(*,*) 'xTail = ', RCMEll%xTail/REarth
            ! write(*,*) 'yDD   = ', RCMEll%yDD/REarth

        endif

        !Set derived quantities
        x0 = (RCMEll%xSun + RCMEll%xTail)/2
        a  = (RCMEll%xSun - RCMEll%xTail)/2
        b  =  RCMEll%yDD


       !$OMP PARALLEL DO default(shared) &
       !$OMP private(ell)
        !do i=iC,RCMApp%nLat_ion
        do i=1,RCMApp%nLat_ion
            do j=1,RCMApp%nLon_ion
                ell = ((RCMApp%X_bmin(i,j,XDIR)-x0)/a)**2.0 + (RCMApp%X_bmin(i,j,YDIR)/b)**2.0
                if ( (ell <= 1) .and. (RCMApp%iopen(i,j) == RCMTOPCLOSED) ) RCMApp%toMHD(i,j) = .true.
            enddo
        enddo


    end subroutine SetIngestion



    !Evaluate eq map at a given point
    !Returns density (#/cc) and pressure (nPa)
    subroutine EvalRCM(lat,lon,llC,t,imW)
        real(rp), intent(in) :: lat,lon,t
        real(rp), intent(in) :: llC(2,2,2,2)
        real(rp), intent(out) :: imW(NVARIMAG)

        real(rp) :: nrcm,prcm,npp,ntot
        integer  :: n
        logical  :: isGood,isGoods(8)
        real(rp) :: lls(8,2),colats(8)
        integer  :: ijs(8,2)

        !Set defaults
        imW(:) = 0.0
        imW(IMDEN ) = 0.0
        imW(IMPR  ) = 0.0
        imW(IMTSCL) = 0.0

        colat = PI/2 - lat

        !Repack
        lls(:,1) = reshape(llC(:,:,:,1),[8])
        lls(:,2) = reshape(llC(:,:,:,2),[8])
        colats = PI/2 - lls(:,1)

        !Do 1st short cut tests
        isGood = all(colats >= RCMApp%gcolat(1)) .and. all(colats <= RCMApp%gcolat(RCMApp%nLat_ion)) &
                 .and. all(lls(:,1) > TINY)
        if (.not. isGood) return

        !If still here, find mapping (i,j) on RCM grid of each corner
        call CornerLocs(lls,ijs)

        !Do second short cut tests
        do n=1,8
            isGoods(n) = RCMApp%toMHD(ijs(n,1),ijs(n,2))
        enddo
        isGood = all(isGoods)

        if (.not. isGood) return
        
        prcm = CornerAvg(ijs,RCMApp%Prcm )*rcmPScl
        npp  = CornerAvg(ijs,RCMApp%Npsph)*rcmNScl
        nrcm = CornerAvg(ijs,RCMApp%Nrcm )*rcmNScl

        if ( (npp >= PPDen) .and. (prcm > TINY) ) then
            ntot = npp + nrcm
        else
            ntot = nrcm
        endif

        !Store data
        imW(IMDEN)  = ntot
        imW(IMPR)   = prcm
        imW(IMTSCL) = 1.0
        imW(IMX1)   = (180.0/PI)*lat
        imW(IMX2)   = (180.0/PI)*lon

        contains

            !Get RCM cells for corner lat/lons
            subroutine CornerLocs(lls,ijs)
                real(rp), intent(in) :: lls(8,2)
                integer, intent(out) :: ijs(8,2)

                integer :: n,i0,j0
                real(rp) :: colat,lon
                do n=1,8
                    colat = PI/2 - lls(n,1)
                    lon   = lls(n,2)

                    i0 = minloc( abs(colat-RCMApp%gcolat),dim=1 )
                    j0 = minloc( abs(lon  -RCMApp%glong ),dim=1 )
                    ijs(n,:) = [i0,j0]
                enddo

            end subroutine CornerLocs

            !Average a quantity over the corners
            function CornerAvg(ijs,Q) result(Qavg)
                integer, intent(in) :: ijs(8,2)
                real(rp), intent(in) :: Q(RCMApp%nLat_ion,RCMApp%nLon_ion)
                real(rp) :: Qavg

                integer :: n,i0,j0

                Qavg = 0.0
                do n=1,8
                    i0 = ijs(n,1)
                    j0 = ijs(n,2)
                    Qavg = Qavg + Q(i0,j0)
                enddo
                Qavg = Qavg/8.0
            end function CornerAvg

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
        !Minimal surface (bEq in Re, bMin in EB)
        call FLEq(ebModel,bTrc,bEq,bMin)
        bMin = bMin*oBScl*1.0e-9 !EB=>Tesla
        bEq = bEq*Re_cgs*1.0e-2 !Re=>meters

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
            !xyzIonC(ZDIR) uses abs(mlat), so make negative
            ijTube%latc = asin(-xyzIonC(ZDIR)/norm2(xyzIonC))
            ijTube%lonc = modulo( atan2(xyzIonC(YDIR),xyzIonC(XDIR)),2*PI )
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

        mdipole = EarthM0g*G2T ! dipole moment in T

        colat = PI/2 - lat
        L = 1.0/(sin(colat)**2.0)
        ijTube%Vol = 32./35.*L**4.0/mdipole
        ijTube%X_bmin(XDIR) = L*cos(lon)*Re_cgs*1.0e-2 !Re=>meters
        ijTube%X_bmin(YDIR) = L*sin(lon)*Re_cgs*1.0e-2 !Re=>meters
        ijTube%X_bmin(ZDIR) = 0.0
        ijTube%bmin = mdipole/L**3.0

        ijTube%iopen = RCMTOPCLOSED
        
        ijTube%pot = 0.0

        ijTube%beta_average = 0.0
        ijTube%Pave = 0.0
        ijTube%Nave = 0.0

        ijTube%latc = -lat
        ijTube%lonc = lon
        !ijTube%Nave = psphD(L)*1.0e+6 !#/cc => #/m3

    end subroutine DipoleTube


!--------------
!Kaiju RCM IO Routines
    subroutine initRCMIO()
        type(IOVAR_T), dimension(MAXRCMIOVAR) :: IOVars

        real(rp), dimension(:,:), allocatable :: iLat,iLon

        integer :: i,j,NLat,NLon
        real(rp) :: dLat,dLon,clMin,clMax

        NLat = RCMApp%nLat_ion
        NLon = RCMApp%nLon_ion

        clMin = RCMApp%gcolat(1)
        clMax = RCMApp%gcolat(NLat)
        dLat = (clMax-clMin)/NLat
        dLon = (2*PI-0.0)/NLon

        allocate(iLat(NLat+1,NLon+1))
        allocate(iLon(NLat+1,NLon+1))

        do j=1,NLon+1
            iLon(:,j) = 0.0 + (j-1)*dLon
        enddo
        dLat = (RCMApp%gcolat(2)-RCMApp%gcolat(1))
        iLat(1,:) = clMin-0.5*dLat
        do i=2,NLat
            dLat = (RCMApp%gcolat(i)-RCMApp%gcolat(i-1))
            iLat(i,:) = iLat(i-1,:) + dLat
        enddo
        !Replicate last dlat
        iLat(NLat+1,:) = iLat(NLat,:) + dLat

        iLat = 90.0-iLat*180.0/PI !Turn colat into lat
        iLon = iLon*180.0/PI

        !Reset IO chain
        call ClearIO(IOVars)
        
        !Flipping lat/lon
        call AddOutVar(IOVars,"X",iLon)
        call AddOutVar(IOVars,"Y",iLat)
                        

        call WriteVars(IOVars,.true.,h5File)

    end subroutine initRCMIO

    subroutine WriteRCM(nOut,MJD,time)
        integer, intent(in) :: nOut
        real(rp), intent(in) :: MJD,time

        type(IOVAR_T), dimension(MAXRCMIOVAR) :: IOVars
        character(len=strLen) :: gStr

        real(rp) :: rcm2Wolf
        integer, dimension(2) :: DimLL
        integer :: Ni,Nj
        
        rcm2Wolf = (1.0e-9)**(IMGAMMA-1.0) !Convert to Wolf units, RCM: Pa (Re/T)^gam => nPa (Re/nT)^gam
        
        !Reset IO chain
        call ClearIO(IOVars)

        call AddOutVar(IOVars,"N",RCMApp%Nrcm*rcmNScl)
        call AddOutVar(IOVars,"Npsph",RCMApp%Npsph*rcmNScl)
        call AddOutVar(IOVars,"P",RCMApp%Prcm*rcmPScl)
        call AddOutVar(IOVars,"IOpen",RCMApp%iopen*1.0_rp)
        call AddOutVar(IOVars,"bVol",RCMApp%Vol)
        call AddOutVar(IOVars,"pot",RCMApp%pot)
        call AddOutVar(IOVars,"xMin",RCMApp%X_bmin(:,:,XDIR)/REarth)
        call AddOutVar(IOVars,"yMin",RCMApp%X_bmin(:,:,YDIR)/REarth)
        call AddOutVar(IOVars,"zMin",RCMApp%X_bmin(:,:,ZDIR)/REarth)
        call AddOutVar(IOVars,"bMin",RCMApp%Bmin)
        
        call AddOutVar(IOVars,"S",rcm2Wolf*RCMApp%Prcm*(RCMApp%Vol**IMGAMMA) )
        call AddOutVar(IOVars,"beta",RCMApp%beta_average)
        call AddOutVar(IOVars,"Pmhd",RCMApp%Pave*rcmPScl)
        call AddOutVar(IOVars,"Nmhd",RCMApp%Nave*rcmNScl)
        call AddOutVar(IOVars,"latc",RCMApp%latc*180.0/PI)
        call AddOutVar(IOVars,"lonc",RCMApp%lonc*180.0/PI)

        call AddOutVar(IOVars,"eavg",RCMApp%eng_avg*1.0e-3) !ev->keV
        call AddOutVar(IOVars,"eflux",RCMApp%flux)

        call AddOutVar(IOVars,"toMHD",merge(1.0_rp,0.0_rp,RCMApp%toMHD))

        !Trim output for colat/aloct to remove wrapping
        DimLL = shape(colat)
        Ni = DimLL(1)
        Nj = DimLL(2)

        call AddOutVar(IOVars,"colat",colat(:,3:Nj))
        call AddOutVar(IOVars,"aloct",aloct(:,3:Nj))
        !Add attributes
        call AddOutVar(IOVars,"time",time)
        call AddOutVar(IOVars,"MJD",MJD)

        write(gStr,'(A,I0)') "Step#", nOut
        call WriteVars(IOVars,.true.,h5File,gStr)

        !Call RCM output
        RCMApp%rcm_nOut = nOut
        call rcm_mhd(time,TINY,RCMApp,RCMWRITEOUTPUT)
        
    end subroutine WriteRCM

    subroutine WriteRCMRestart(nRes,MJD,time)
        
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time

        RCMApp%rcm_nRes = nRes
        call rcm_mhd(time,TINY,RCMApp,RCMWRITERESTART)
        
    end subroutine WriteRCMRestart
end module rcmimag
