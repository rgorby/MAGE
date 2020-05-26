module rcm_mhd_io
    use rcm_mhd_interfaces
    use ioh5
    use xml_input
    use rcm_mhd_mod,  ONLY : rcm_mhd
    use rcm_mod_subs, ONLY : colat, aloct

    integer, parameter   , private :: MAXRCMIOVAR = 30
    character(len=strLen), private :: h5File,RCMH5
    real(rp), parameter  , private :: IMGAMMA = 5.0/3.0

    contains
!--------------
!Kaiju RCM IO Routines
    subroutine initRCMIO(RCMApp,isResO)
        type(rcm_mhd_t), intent(inout) :: RCMApp
        logical, intent(in), optional :: isResO

        type(IOVAR_T), dimension(MAXRCMIOVAR) :: IOVars
        real(rp), dimension(:,:), allocatable :: iLat,iLon

        integer :: i,j,NLat,NLon
        real(rp) :: dLat,dLon,clMin,clMax
        logical :: isRestart,fExist

        !Set isRestart
        if (present(isResO)) then
            isRestart = isResO
        else
            isRestart = .false.
        endif

        !Create file names and nuke old stuff
        h5File = trim(RCMApp%rcm_runid) // ".mhdrcm.h5" !MHD-RCM coupling data
        RCMH5  = trim(RCMApp%rcm_runid) // ".rcm.h5" !RCM data

        fExist = CheckFile(h5File)
        write(*,*) 'RCM outputting to ',trim(h5File)
        if (.not. isRestart) then
            !Kill it all
            call CheckAndKill(h5File) !For non-restart but file exists
            call CheckAndKill(RCMH5)
        endif

        if (isRestart .and. fExist) then
            !File already exists, don't need to init
            return
        endif

        !If we're still here then we need to do work
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
        call AddOutVar(IOVars,"UnitsID","RCMMHD")

        call WriteVars(IOVars,.true.,h5File)

    end subroutine initRCMIO

    subroutine WriteRCM(RCMApp,nOut,MJD,time)
        USE constants, ONLY: nt
        USE rcm_mod_subs, ONLY:isize,jsize,jwrap
        type(rcm_mhd_t), intent(inout) :: RCMApp
        integer, intent(in) :: nOut
        real(rp), intent(in) :: MJD,time

        type(IOVAR_T), dimension(MAXRCMIOVAR) :: IOVars
        character(len=strLen) :: gStr

        real(rp) :: rcm2Wolf
!        integer, dimension(2) :: DimLL
!        integer :: Ni,Nj
        
        rcm2Wolf = nt**(IMGAMMA-1.0) !Convert to Wolf units, RCM: Pa (Re/T)^gam => nPa (Re/nT)^gam
        
        !Reset IO chain
        call ClearIO(IOVars)

        call AddOutVar(IOVars,"N",RCMApp%Nrcm*rcmNScl,uStr="#/cc")
        call AddOutVar(IOVars,"Npsph",RCMApp%Npsph*rcmNScl,uStr="#/cc")
        call AddOutVar(IOVars,"P",RCMApp%Prcm*rcmPScl,uStr="nPa")
        call AddOutVar(IOVars,"IOpen",RCMApp%iopen*1.0_rp)
        call AddOutVar(IOVars,"bVol",RCMApp%Vol*nt,uStr="Re/nT")
        call AddOutVar(IOVars,"pot",RCMApp%pot,uStr="V")
        call AddOutVar(IOVars,"xMin",RCMApp%X_bmin(:,:,XDIR)/REarth,uStr="Re")
        call AddOutVar(IOVars,"yMin",RCMApp%X_bmin(:,:,YDIR)/REarth,uStr="Re")
        call AddOutVar(IOVars,"zMin",RCMApp%X_bmin(:,:,ZDIR)/REarth,uStr="Re")
        call AddOutVar(IOVars,"bMin",RCMApp%Bmin,uStr="T")
        
        call AddOutVar(IOVars,"S",rcm2Wolf*RCMApp%Prcm*(RCMApp%Vol**IMGAMMA),uStr="Wolf")
        call AddOutVar(IOVars,"beta",RCMApp%beta_average)
        call AddOutVar(IOVars,"Pmhd",RCMApp%Pave*rcmPScl,uStr="nPa")
        call AddOutVar(IOVars,"Nmhd",RCMApp%Nave*rcmNScl,uStr="#/cc")
        call AddOutVar(IOVars,"latc",RCMApp%latc*180.0/PI,uStr="deg")
        call AddOutVar(IOVars,"lonc",RCMApp%lonc*180.0/PI,uStr="deg")
        call AddOutVar(IOVars,"Lb"  ,RCMApp%Lb,uStr="Re")
        call AddOutVar(IOVars,"Tb"  ,RCMApp%Tb,uStr="s")

        call AddOutVar(IOVars,"eavg",RCMApp%eng_avg*1.0e-3,uStr="keV") !ev->keV
        call AddOutVar(IOVars,"eflux",RCMApp%flux,uStr="ergs/cm2")
        call AddOutVar(IOVars,"birk",RCMApp%fac,uStr="uA/m2")

        call AddOutVar(IOVars,"toMHD",merge(1.0_rp,0.0_rp,RCMApp%toMHD))

        !Trim output for colat/aloct to remove wrapping
!        DimLL = shape(colat)
!        Ni = DimLL(1)
!        Nj = DimLL(2)

        call AddOutVar(IOVars,"colat",colat(:,jwrap:jsize))
        call AddOutVar(IOVars,"aloct",aloct(:,jwrap:jsize))
        !Add attributes
        call AddOutVar(IOVars,"time",time)
        call AddOutVar(IOVars,"MJD",MJD)

        write(gStr,'(A,I0)') "Step#", nOut
        call WriteVars(IOVars,.true.,h5File,gStr)

        !Call RCM output
        RCMApp%rcm_nOut = nOut
        call rcm_mhd(time,TINY,RCMApp,RCMWRITEOUTPUT)
        
    end subroutine WriteRCM

    subroutine WriteRCMRestart(RCMApp,nRes,MJD,time)
        type(rcm_mhd_t), intent(inout) :: RCMApp
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time

        RCMApp%rcm_nRes = nRes
        call rcm_mhd(time,TINY,RCMApp,RCMWRITERESTART)
        
    end subroutine WriteRCMRestart

    subroutine RCMRestartInfo(RCMApp,xmlInp,t0,isRCMopt)
        type(rcm_mhd_t)  , intent(inout) :: RCMApp
        type(XML_Input_T), intent(in)    :: xmlInp
        real(rp), intent(out) :: t0
        logical, intent(in), optional :: isRCMopt

        integer :: nRes
        character(len=strLen) :: resID,nStr,inH5
        type(IOVAR_T), dimension(MAXRCMIOVAR) :: IOVars
        logical :: doSP,isRCM

        if (present(isRCMopt)) then
            isRCM = isRCMopt
        else
            isRCM = .false.
        endif
        if (isRCM) then
            call xmlInp%Set_Val(resID,"/rcm/restart/resID","msphere")
            call xmlInp%Set_Val(nRes ,"/rcm/restart/nRes" ,-1)
        else
            call xmlInp%Set_Val(resID,"/gamera/restart/resID","msphere")
            call xmlInp%Set_Val(nRes ,"/gamera/restart/nRes" ,-1)
        endif            
        !Get number string
        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') nRes
        endif
        
        inH5 = trim(resID) // ".RCM.Res." // trim(nStr) // ".h5"

        call CheckFileOrDie(inH5,"Restart file not found ...")

        !Get time data out of restart
        doSP = .false. !Restarts are always double precision

        call ClearIO(IOVars) !Reset IO chain
        call AddInVar(IOVars,"time",vTypeO=IOREAL )
        call ReadVars(IOVars,doSP,inH5)
        t0   = GetIOReal(IOVars,"time")

        if (ioExist(inH5,"nRes")) then
            call ClearIO(IOVars) !Reset IO chain
            call AddInVar(IOVars,"nRes",vTypeO=IOINT  )
            call ReadVars(IOVars,doSP,inH5)
            nRes = GetIOInt(IOVars,"nRes")
            RCMApp%rcm_nRes = nRes 
        endif

        RCMApp%rcm_nRes = nRes + 1 !Holds step for *NEXT* restart
    end subroutine RCMRestartInfo

end module rcm_mhd_io
