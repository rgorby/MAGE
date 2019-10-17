!Routines to handle RCM inner magnetosphere model
module rcmimag
    use volttypes
    use ioh5
    use files
    use earthhelper
    use rcm_mhd_interfaces

    implicit none

    integer(ip), parameter,private :: RCMINIT=0,RCMADVANCE=1,RCMFINISH=-1
    type(rcm_mhd_T), private :: RCMApp

    !Scaling parameters
    real(rp), private :: rcmPScl = 1.0e+9 !Convert Pa->nPa
    real(rp), private :: rcmNScl = 1.0 !Convert xxx => #/cc

    !Do I need this stuff?
    real(rp), private :: rcm_boundary_s =35,rcm_boundary_e =2
    real(rp), private :: colat_boundary
    real(rp), private :: ddt


    type RCMTube_T
        real(rp) :: Vol,bmin,beta_average,Pave,Nave,pot
        real(rp) :: X_bmin(NDIM)
        integer(ip) :: iopen
    end type RCMTube_T

    contains

    !Initialize RCM inner magnetosphere model
    subroutine initRCM(iXML,isRestart,t0,dtCpl)
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        real(rp), intent(in) :: t0,dtCpl

        if (isRestart) then
            write(*,*) 'What do I do here?'
            stop
        else
            write(*,*) 'Initializing RCM ...'
            call rcm_mhd(t0,dtCpl,RCMApp,RCMINIT)
        endif

        call iXML%Set_Val(ddt,"rcm/ddt",15.0) !RCM substep [s]

    end subroutine initRCM

    !Advance RCM from Voltron data
    subroutine AdvanceRCM(vApp,tAdv)
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        integer :: i,j,n,rcmbndy,nStp
        real(rp) :: colat,lat,lon
        real(rp) :: dtCum

        type(RCMTube_T) :: ijTube

    !Load RCM tubes
        rcmbndy = 30 !I don't know where this is coming from
        colat_boundary = sin(RCMApp%gcolat(rcmbndy))
        write(*,*) RCMApp%nLat_ion,RCMApp%nLon_ion

        do i=1,RCMApp%nLat_ion
            do j=1,RCMApp%nLon_ion
                colat = RCMApp%gcolat(i)
                lat = PI/2 - colat
                lon = RCMApp%glong(j)

                !Load a dipole-ish tube
                call DipoleTube(lat,lon,ijTube)

                !Pull data into RCMApp
                RCMApp%Vol(i,j)          = ijTube%Vol
                RCMApp%bmin(i,j)         = ijTube%bmin
                RCMApp%iopen(i,j)        = ijTube%iopen
                RCMApp%beta_average(i,j) = ijTube%beta_average
                RCMApp%Pave(i,j)         = ijTube%Pave
                RCMApp%Nave(i,j)         = ijTube%Nave
                RCMApp%pot(i,j)          = ijTube%pot
                RCMApp%X_bmin(i,j,:)     = ijTube%X_bmin
            enddo
        enddo

        !Set RCM boundary
        RCMApp%Vol(1:rcmbndy,:) = -1.0
        RCMApp%iopen(1:rcmbndy,:) = 1 !Open

    !Advance from vApp%time to tAdv
        !Substep until done
        !NOTE: Weird use of real(iprec) on RCM side
        dtCum = 0.0
        nStp = int( (tAdv-vApp%time)/ddt )
        do n=1,nStp
            call rcm_mhd(vApp%time+dtCum,ddt,RCMApp,RCMADVANCE)
            dtCum = dtCum+ddt
        enddo

        !Add some diagnostic stuff here?

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

    !Lazy test flux tube
    subroutine DipoleTube(lat,lon,ijTube)
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

        if (colat < colat_boundary) then
            ijTube%pot = -potmax/2.0*sin(lon)*sin(colat)
        else
            ijTube%pot = -potmax/2.0*sin(lon)*sin(colat_boundary)/sin(colat)
        endif

    end subroutine DipoleTube
end module rcmimag
