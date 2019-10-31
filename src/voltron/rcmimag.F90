!Routines to handle RCM inner magnetosphere model
!NOTES: 
!-Figure out flux-tube volume units
!-add ReMIX potential to MHD=>RCM tubes
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

    integer(ip), parameter,private :: RCMINIT=0,RCMADVANCE=1,RCMWRITEREC=-2,RCMWRITEINDEX=-1
    type(rcm_mhd_T), private :: RCMApp

    !Scaling parameters
    real(rp), private :: rcmPScl = 1.0e+9 !Convert Pa->nPa
    real(rp), private :: rcmNScl = 1.0e-6 !Convert #/m3 => #/cc
    real(rp), parameter :: RIonRCM = (RionE/REarth)*1.0e+6
    integer, parameter :: MAXRCMIOVAR = 20
    character(len=strLen), private :: h5File


    !Information taken from MHD flux tubes
    !TODO: Figure out -volume for open flux tubes?
    !TODO: Figure out RCM boundaries
    !TODO: Figure out iopen values
    !TODO: Figure out units for potential

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
    subroutine initRCM(iXML,isRestart,t0,dtCpl)
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        real(rp), intent(in) :: t0,dtCpl

        character(len=strLen) :: RunID
        logical :: fExist

        if (isRestart) then
            write(*,*) 'This should be a restart'
            write(*,*) 'Restarting RCM ...'
            call rcm_mhd(t0,dtCpl,RCMApp,RCMINIT)
            call init_rcm_mix(RCMApp)
        else
            write(*,*) 'Initializing RCM ...'
            call rcm_mhd(t0,dtCpl,RCMApp,RCMINIT)
            call init_rcm_mix(RCMApp)
        endif

        call iXML%Set_Val(RunID,"/gamera/sim/runid","sim")

        h5File = trim(RunID) // ".rcm.h5"

        fExist = CheckFile(h5File)
        write(*,*) 'RCM outputting to ',trim(h5File)
        if ( (.not. isRestart) .or. (isRestart .and. (.not.fExist)) ) then
            !Not a restart or it is a restart and no file
            call CheckAndKill(h5File) !For non-restart but file exists

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

    !Get potential from mix
        call map_rcm_mix(vApp,mixPot)

    !Load RCM tubes
       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP private(i,j,colat,lat,lon,ijTube)
        do i=1,RCMApp%nLat_ion
            do j=1,RCMApp%nLon_ion
                colat = RCMApp%gcolat(i)
                lat = PI/2 - colat
                lon = RCMApp%glong(j)
                
                !call DipoleTube(vApp,lat,lon,ijTube)
                call MHDTube(vApp,lat,lon,ijTube)

                !Pull data into RCMApp
                RCMApp%Vol(i,j)          = ijTube%Vol
                RCMApp%bmin(i,j)         = ijTube%bmin
                RCMApp%iopen(i,j)        = ijTube%iopen
                RCMApp%beta_average(i,j) = ijTube%beta_average
                RCMApp%Pave(i,j)         = ijTube%Pave
                RCMApp%Nave(i,j)         = ijTube%Nave
                !RCMApp%pot(i,j)          = ijTube%pot
                ! mix variables are stored in this order (longitude,colatitude), hence the index flip
                RCMApp%pot(i,j)          = mixPot(j,i)   
                RCMApp%X_bmin(i,j,:)     = ijTube%X_bmin
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

        real(rp) :: colat
        integer :: i0,j0

        colat = PI/2 - lat

        !Just find closest cell
        i0 = minloc( abs(colat-RCMApp%gcolat),dim=1 )
        j0 = minloc( abs(lon  -RCMApp%glong ),dim=1 )

        imW(IMDEN) = RCMApp%Nrcm(i0,j0)*rcmNScl
        imW(IMPR ) = RCMApp%Prcm(i0,j0)*rcmPScl

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
        !type(RCMTube_T) :: dpTube
        integer :: OCb
        real(rp) :: bD,bP,dvB,bBeta
    !First get seed for trace
        !Assume lat/lon @ Earth, dipole push
        xyzIon(XDIR) = RIonRCM*cos(lat)*cos(lon)
        xyzIon(YDIR) = RIonRCM*cos(lat)*sin(lon)
        xyzIon(ZDIR) = RIonRCM*sin(lat)
        x0 = DipoleShift(xyzIon,2.05_rp)
        
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

        ! !Get dipole tube to test against
        ! call DipoleTube(vApp,lat,lon,dpTube)

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
            !WTF?
            ijTube%iopen = -999
            ijTube%Vol = -999
        end select

        ijTube%Pave = bP
        ijTube%Nave = bD
        ijTube%beta_average = bBeta
        

        !ijTube%pot = dpTube%pot

        ! write(*,*) '---'
        ! write(*,*) 'Lat/Lon = ', lat*180.0/PI,lon*180.0/PI
        ! write(*,*) 'x0 = ', x0
        ! write(*,'(a,2es9.2)') 'Vol = ', ijTube%Vol,dpTube%Vol
        ! write(*,*) 'Den = ', ijTube%Nave,dpTube%Nave
        ! write(*,*) 'Pre = ', ijTube%Pave,dpTube%Pave
        ! write(*,*) 'iop = ', ijTube%iopen,dpTube%iopen
        ! write(*,*) 'bmin = ', ijTube%bmin,dpTube%bmin
        ! write(*,*) 'xEq = ', ijTube%X_bmin,dpTube%X_bmin
        ! write(*,*) 'beta = ', ijTube%beta_average,dpTube%beta_average
        ! write(*,*) 'pot = ', ijTube%pot, dpTube%pot
        ! write(*,*) '---'

    end subroutine MHDTube

    !Lazy test flux tube
    subroutine DipoleTube(vApp,lat,lon,ijTube)
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
        call AddOutVar(IOVars,"S",RCMApp%Prcm*(RCMApp%Vol**(5.0/3.0)) )
        call AddOutVar(IOVars,"beta",RCMApp%beta_average)

        !Add attributes
        call AddOutVar(IOVars,"time",time)
        call AddOutVar(IOVars,"MJD",MJD)

        write(gStr,'(A,I0)') "Step#", nOut
        call WriteVars(IOVars,.true.,h5File,gStr)

    end subroutine WriteRCM

    subroutine WriteRCMRestart(nRes,MJD,time)
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time

        !Do two things, force a record output and output time index
        !call rcm_mhd(time,TINY,RCMApp,RCMWRITEREC)
        !call rcm_mhd(time,TINY,RCMApp,RCMWRITEINDEX)
        write(*,*) 'I should do an RCM restart here ...'
    end subroutine WriteRCMRestart
end module rcmimag
