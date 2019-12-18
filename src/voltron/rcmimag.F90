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

    implicit none

    integer(ip), parameter,private :: RCMINIT=0,RCMADVANCE=1,RCMRESTART=2,RCMWRITERESTART=-2,RCMWRITEOUTPUT=-3
    type(rcm_mhd_T), private :: RCMApp

    !Scaling parameters
    real(rp), parameter, private :: rcmPScl = 1.0e+9 !Convert Pa->nPa
    real(rp), parameter, private :: rcmNScl = 1.0e-6 !Convert #/m3 => #/cc
    real(rp), parameter, private :: IMGAMMA = 5.0/3.0

    real(rp), parameter :: RIonRCM = (RionE/REarth)*1.0e+6
    

    integer, parameter :: MAXRCMIOVAR = 30
    character(len=strLen), private :: h5File


    !Information taken from MHD flux tubes
    !TODO: Figure out RCM boundaries
    !TODO: Figure out iopen values

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
    end type RCMTube_T

    real(rp), dimension(:,:), allocatable, private :: mixPot
    contains

    !Initialize RCM inner magnetosphere model
    subroutine initRCM(iXML,isRestart,t0,dtCpl,nRes)
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        real(rp), intent(in) :: t0,dtCpl
        integer, intent(in), optional :: nRes

        character(len=strLen) :: RunID,RCMH5
        logical :: fExist

        if (isRestart) then
            RCMApp%rcm_nRes = nRes
            write(*,*) 'Restarting RCM ...'
            call rcm_mhd(t0,dtCpl,RCMApp,RCMRESTART)
            call init_rcm_mix(RCMApp)
        else
            write(*,*) 'Initializing RCM ...'
            call rcm_mhd(t0,dtCpl,RCMApp,RCMINIT)
            call init_rcm_mix(RCMApp)
        endif

        call iXML%Set_Val(RunID,"/gamera/sim/runid","sim")
        RCMApp%rcm_runid = trim(RunID)
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
        
    end subroutine initRCM

    !Advance RCM from Voltron data
    subroutine AdvanceRCM(vApp,tAdv)
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        integer :: i,j,n,nStp
        real(rp) :: colat,lat,lon
        real(rp) :: dtAdv
        type(RCMTube_T) :: ijTube

        real(rp) :: llBC
        logical :: isLL

        llBC = vApp%mhd2chmp%lowlatBC

    !Get potential from mix
        call map_rcm_mix(vApp,mixPot)

    !Load RCM tubes
       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP schedule(dynamic) &
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

                !mix variables are stored in this order (longitude,colatitude), hence the index flip
                RCMApp%pot(i,j)          = mixPot(j,i)
                
            enddo
        enddo

    !Advance from vApp%time to tAdv
        dtAdv = tAdv-vApp%time !RCM-DT
        call rcm_mhd(vApp%time,dtAdv,RCMApp,RCMADVANCE)

    end subroutine AdvanceRCM

    !Evaluate eq map at a given point
    !Returns density (#/cc) and pressure (nPa)
    subroutine EvalRCM(lat,lon,t,imW)
        real(rp), intent(in) :: lat,lon,t
        real(rp), intent(out) :: imW(NVARIMAG)

        real(rp) :: colat, alpha, beta, LScl
        integer :: i0,j0
        logical :: isGood

        colat = PI/2 - lat
        !Just find closest cell
        i0 = minloc( abs(colat-RCMApp%gcolat),dim=1 )
        j0 = minloc( abs(lon  -RCMApp%glong ),dim=1 )

        !Test for good cell
        isGood = (colat >= minval(RCMApp%gcolat)) .and. (colat <= maxval(RCMApp%gcolat)) &
                 & .and. (RCMApp%iopen(i0,j0) == -1) .and. (lat>TINY)

        imW = 0.0

        if (isGood) then
            beta = RCMApp%beta_average(i0,j0)
            alpha = 1.0/(1.0 + beta*IMGAMMA/2.0)
            LScl = RCMApp%Vol(i0,j0)*RCMApp%bmin(i0,j0) !Lengthscale

            imW(IMDEN ) = RCMApp%Nrcm(i0,j0)*rcmNScl
            imW(IMPR  ) = RCMApp%Prcm(i0,j0)*rcmPScl*alpha
            imW(IMLSCL) = LScl
            imW(IMTSCL) = 1.0
        else
            imW(IMDEN ) = 0.0
            imW(IMPR  ) = 0.0
            imW(IMLSCL) = 0.0
            imW(IMTSCL) = 1.0
        endif
        
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

        end associate


    !Scale and store information
        ijTube%X_bmin = bEq
        ijTube%bmin = bMin
        select case(OCb)
        case(0)
            !Solar wind (is this right?)
            ijTube%iopen = 1
            ijTube%Vol = -dvB
        case(1)
            !Open field
            ijTube%iopen = 1
            ijTube%Vol = -dvB
        case(2)
            !Closed field
            ijTube%iopen = -1
            ijTube%Vol = dvB
        case default
            !WTF? (timeout)
            ijTube%iopen = 1
            ijTube%Vol = -1
        end select

        ijTube%Pave = bP
        ijTube%Nave = bD
        ijTube%beta_average = bBeta
        
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
        ijTube%iopen = -1
        ijTube%pot = 0.0

        ijTube%beta_average = 0.0
        ijTube%Pave = 0.0
        ijTube%Nave = 0.0

        !ijTube%Nave = psphD(L)*1.0e+6 !#/cc => #/m3

    end subroutine DipoleTube

    !Lazy test flux tube
    subroutine oDipoleTube(vApp,lat,lon,ijTube)
        type(voltApp_T), intent(in) :: vApp
        real(rp), intent(in) :: lat,lon
        type(RCMTube_T), intent(out) :: ijTube

        real(rp) :: L,colat

        real(rp) :: mdipole = 3.0e-5 ! dipole moment in T
        real(rp) :: Lmax = 4.0 ! location of Pressure max
        real(rp) :: pmax = 5.0e-8 ! pressure max in Pa
        real(rp) :: pmin = 1.0e-11 ! min BG pressure in Pa
        real(rp) :: nmax = 1.0e7 ! dens in ple/m^3
        real(rp) :: nmin = 1.0e4 ! min dens in ple/m^3
        real(rp) :: potmax = 5.0e4 ! potential max
        real(rp) :: re = 6380.e3
        real(rp) :: colat_boundary

        colat = PI/2 - lat
        L = 1.0/(sin(colat)**2.0)
        ijTube%Vol =32./35.*L**4.0/mdipole
        ijTube%X_bmin(XDIR) = L*cos(lon)*re
        ijTube%X_bmin(YDIR) = L*sin(lon)*re
        ijTube%X_bmin(ZDIR) = 0.0
        ijTube%bmin = mdipole/L**3.0
        ijTube%iopen = -1
        ijTube%beta_average = 0.1
        ijTube%Pave = pmax*exp(-(L-Lmax)**2.0) + pmin
        ijTube%Nave = nmax*exp(-(L-Lmax)**2.0) + nmin
        colat_boundary = PI/4.0
        if (colat < colat_boundary) then
            ijTube%pot = -potmax/2.0*sin(lon)*sin(colat)
        else
            ijTube%pot = -potmax/2.0*sin(lon)*sin(colat_boundary)/sin(colat)
        endif

    end subroutine oDipoleTube

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

        do i=1,NLat+1
            do j=1,NLon+1
                iLat(i,j) = clMin + (i-1)*dLat
                iLon(i,j) = 0.0 + (j-1)*dLon
            enddo
        enddo

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
