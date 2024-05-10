module rcm_mhd_io
    use rcm_mhd_interfaces
    use ioh5
    use xml_input
    use rcm_mod_subs, ONLY : colat, aloct
    use rice_housekeeping_module, ONLY : nSkipFL,doFLOut
    use rcmdefs
    
    implicit none

    integer, parameter   , private :: MAXRCMIOVAR = 80
    character(len=strLen), private :: h5File,RCMH5,FLH5
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
        FLH5   = trim(RCMApp%rcm_runid) // ".rcmfl.h5" !RCM field lines
        RCMH5  = trim(RCMApp%rcm_runid) // ".rcm.h5" !RCM data
        
        fExist = CheckFile(h5File)
        write(*,*) 'RCM outputting to ',trim(h5File)

        if (.not. isRestart) then
            !Kill it all
            call CheckAndKill(h5File) !For non-restart but file exists
            call CheckAndKill(FLH5)
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
        
        integer :: NLat,NLon

        NLat = RCMApp%nLat_ion
        NLon = RCMApp%nLon_ion
        
        rcm2Wolf = nt**(IMGAMMA-1.0) !Convert to Wolf units, RCM: Pa (Re/T)^gam => nPa (Re/nT)^gam
        

        !Reset IO chain
        call ClearIO(IOVars)

        call AddOutVar(IOVars,"N",RCMApp%Nrcm*rcmNScl,uStr="#/cc")
        call AddOutVar(IOVars,"Npsph",RCMApp%Npsph*rcmNScl,uStr="#/cc")
        call AddOutVar(IOVars,"P" ,RCMApp%Prcm *rcmPScl,uStr="nPa")
        call AddOutVar(IOVars,"Pe",RCMApp%Percm*rcmPScl,uStr="nPa")
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
        call AddOutVar(IOVars,"Nmhd0",RCMApp%N0*rcmNScl,uStr="#/cc")
        call AddOutVar(IOVars,"oxyfrac",RCMApp%oxyfrac,uStr="fraction")

        call AddOutVar(IOVars,"latc",RCMApp%latc*180.0/PI,uStr="deg")
        call AddOutVar(IOVars,"lonc",RCMApp%lonc*180.0/PI,uStr="deg")
        call AddOutVar(IOVars,"lossc",RCMApp%losscone*180.0/PI,uStr="deg")
        call AddOutVar(IOVars,"Lb"  ,RCMApp%Lb,uStr="Re")
        call AddOutVar(IOVars,"Tb"  ,RCMApp%Tb,uStr="s")
        call AddOutVar(IOVars,"radcurv"  ,RCMApp%radcurv,uStr="Re")
        call AddOutVar(IOVars,"wIMAG"  ,RCMApp%wIMAG,uStr="weight")

        call AddOutVar(IOVars,"eeavg" ,RCMApp%eng_avg(:,:,RCMELECTRON)*1.0e-3,uStr="keV") !ev->keV electrons
        call AddOutVar(IOVars,"eeflux",RCMApp%flux   (:,:,RCMELECTRON),uStr="ergs/cm2")
        call AddOutVar(IOVars,"enflux",RCMApp%nflx   (:,:,RCMELECTRON),uStr="#/cm2/s")
        call AddOutVar(IOVars,"ieavg" ,RCMApp%eng_avg(:,:,RCMPROTON)*1.0e-3,uStr="keV") !ev->keV ions
        call AddOutVar(IOVars,"ieflux",RCMApp%flux   (:,:,RCMPROTON),uStr="ergs/cm2")
        call AddOutVar(IOVars,"influx",RCMApp%nflx   (:,:,RCMPROTON),uStr="#/cm2/s")

        call AddOutVar(IOVars,"birk",RCMApp%fac,uStr="uA/m2",dStr="RCM Vasyliunas FACs")
        call AddOutVar(IOVars,"nTrc",RCMApp%nTrc*1.0_rp,uStr="steps")

        call AddOutVar(IOVars,"TioTe0",RCMApp%TioTe0,dStr="Empirical Ti/Te")

        call AddOutVar(IOVars,"toMHD",merge(1.0_rp,0.0_rp,RCMApp%toMHD))
        call AddOutVar(IOVars,"errD",RCMApp%errD,uStr="X'/X")
        call AddOutVar(IOVars,"errP",RCMApp%errP,uStr="X'/X")

        call AddOutVar(IOVars,"colat",colat(:,jwrap:jsize))
        call AddOutVar(IOVars,"aloct",aloct(:,jwrap:jsize))
        !Add attributes
        call AddOutVar(IOVars,"time",time)
        call AddOutVar(IOVars,"MJD",MJD)

        write(gStr,'(A,I0)') "Step#", nOut
        call WriteVars(IOVars,.true.,h5File,gStr)
        
    end subroutine WriteRCM

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
            call xmlInp%Set_Val(resID,"/Kaiju/rcm/restart/resID","msphere")
            call xmlInp%Set_Val(nRes ,"/Kaiju/rcm/restart/nRes" ,-1)
        else
            call xmlInp%Set_Val(resID,"/Kaiju/gamera/restart/resID","msphere")
            call xmlInp%Set_Val(nRes ,"/Kaiju/gamera/restart/nRes" ,-1)
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

    !Write out field lines
    subroutine WriteRCMFLs(RCMFLs,nOut,MJD,time,Ni,Nj)
        USE ebtypes
        integer, intent(in) :: nOut,Ni,Nj
        real(rp), intent(in) :: MJD,time
        type(magLine_T), intent(in), dimension(Ni,Nj) :: RCMFLs

        type(IOVAR_T), dimension(MAXRCMIOVAR) :: IOVars
        character(len=strLen) :: gStr,lnStr
        integer :: i,j,n
        
        !Bail out if we're not doing this
        if (.not. doFLOut) return

    !Create group and write base data
        write(gStr,'(A,I0)') "Step#", nOut
        call AddOutVar(IOVars,"time",time)
        call AddOutVar(IOVars,"MJD",MJD)

        
        call WriteVars(IOVars,.true.,FLH5,gStr)
        call ClearIO(IOVars)

        !Now loop through and create subgroup for each line (w/ striding)
        !TODO: Avoid the individual write for every line
        n = 0
        do i=1,Ni,nSkipFL
            do j=1,Nj-1,nSkipFL
                write(lnStr,'(A,I0)') "Line#", n
                if (RCMFLs(i,j)%isGood) then
                    call OutLine(RCMFLs(i,j),gStr,lnStr,IOVars)
                    n = n + 1
                endif
            enddo
        enddo

    end subroutine WriteRCMFLs

    !Write out individual line
    subroutine OutLine(fL,gStr,lnStr,IOVars)
        USE ebtypes
        USE gdefs
        type(magLine_T), intent(in) :: fL
        character(len=strLen), intent(in) :: gStr,lnStr
        type(IOVAR_T), intent(inout), dimension(MAXRCMIOVAR) :: IOVars
        integer :: i,Np,Npp,n0
        
        call ClearIO(IOVars)
        Np = fL%Nm + fL%Np + 1
        if (Np<=nSkipFL) return
        n0 = fL%Nm

        !Add scalar stuff
        !Record seed point
        call AddOutVar(IOVars,"x0",fL%x0(XDIR))
        call AddOutVar(IOVars,"y0",fL%x0(YDIR))
        call AddOutVar(IOVars,"z0",fL%x0(ZDIR))

        !Do striding through field line points
        Npp = size(fL%xyz(0:-n0:-nSkipFL,XDIR))

        call AddOutVar(IOVars,"xyz",transpose(fL%xyz(0:-n0:-nSkipFL,XDIR:ZDIR)))
        call AddOutVar(IOVars,"Np",Npp)
        call AddOutVar(IOVars,"n0",1) !Seed point is now the first point

        !Only output some of the variables
        call AddOutVar(IOVars,"B",fL%magB(0:-n0:-nSkipFL),uStr="nT")
        call AddOutVar(IOVars,"D",fL%Gas (0:-n0:-nSkipFL,DEN     ,BLK),uStr="#/cc")
        call AddOutVar(IOVars,"P",fL%Gas (0:-n0:-nSkipFL,PRESSURE,BLK),uStr="nPa" )

        !Write output chain
        call WriteVars(IOVars,.true.,FLH5,gStr,lnStr)
        call ClearIO(IOVars)

    end subroutine OutLine

    subroutine WriteMHD2IMagRestart(RCMApp,nRes,MJD,time)
        type(rcm_mhd_t)  , intent(inout) :: RCMApp
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time

        type(IOVAR_T), dimension(MAXRCMIOVAR) :: IOVars
        character(len=strLen) :: ResF,lnResF

        write (ResF, '(A,A,I0.5,A)') trim(RCMApp%rcm_runid), ".mhd2imag.Res.", nRes, ".h5"
        call CheckAndKill(ResF)
        call ClearIO(IOVars)

        !Main attributes
        call AddOutVar(IOVars,"nRes",nRes)
        call AddOutVar(IOVars,"MJD" ,MJD)
        call AddOutVar(IOVars,"time",time)
        call AddOutVar(IOVars,"llBC",RCMApp%llBC)
        call AddOutVar(IOVars,"dtCpl",RCMApp%dtCpl)
        call AddOutVar(IOVars,"planet_radius",RCMApp%planet_radius)
        call AddOutVar(IOVars,"iono_radius"  ,RCMApp%iono_radius)

        !Variables
        call AddOutVar(IOVars,"gcolat"      ,RCMApp%gcolat      )
        call AddOutVar(IOVars,"glong"       ,RCMApp%glong       )  
        call AddOutVar(IOVars,"pot"         ,RCMApp%pot         )
        call AddOutVar(IOVars,"eng_avg"     ,RCMApp%eng_avg     )    
        call AddOutVar(IOVars,"flux"        ,RCMApp%flux        )  
        call AddOutVar(IOVars,"nflx"        ,RCMApp%nflx        )  
        call AddOutVar(IOVars,"fac"         ,RCMApp%fac         )  
        call AddOutVar(IOVars,"Pave"        ,RCMApp%Pave        )   
        call AddOutVar(IOVars,"Nave"        ,RCMApp%Nave        )
        call AddOutVar(IOVars,"N0"          ,RCMApp%N0          ) 
        call AddOutVar(IOVars,"Vol"         ,RCMApp%Vol         )
        call AddOutVar(IOVars,"X_bmin"      ,RCMApp%X_bmin      )
        call AddOutVar(IOVars,"Bmin"        ,RCMApp%Bmin        )
        call AddOutVar(IOVars,"beta_average",RCMApp%beta_average)
        call AddOutVar(IOVars,"iopen"       ,RCMApp%iopen*1.0_rp)

        call AddOutVar(IOVars,"Prcm"    ,RCMApp%Prcm    )
        call AddOutVar(IOVars,"Nrcm"    ,RCMApp%Nrcm    )
        call AddOutVar(IOVars,"Npsph"   ,RCMApp%Npsph   )
        call AddOutVar(IOVars,"sigmap"  ,RCMApp%sigmap  )
        call AddOutVar(IOVars,"sigmah"  ,RCMApp%sigmah  )
        call AddOutVar(IOVars,"oxyfrac" ,RCMApp%oxyfrac )
        call AddOutVar(IOVars,"Percm"   ,RCMApp%Percm   )
        call AddOutVar(IOVars,"latc"    ,RCMApp%latc    )
        call AddOutVar(IOVars,"lonc"    ,RCMApp%lonc    )
        call AddOutVar(IOVars,"TioTe0"  ,RCMApp%TioTe0  )
        call AddOutVar(IOVars,"Lb"      ,RCMApp%Lb      )
        call AddOutVar(IOVars,"Tb"      ,RCMApp%Tb      )
        call AddOutVar(IOVars,"losscone",RCMApp%losscone)
        call AddOutVar(IOVars,"radcurv" ,RCMApp%radcurv )
        call AddOutVar(IOVars,"wIMAG"   ,RCMApp%wIMAG   )

        call AddOutVar(IOVars,"toMHD",merge(1.0_rp,0.0_rp,RCMApp%toMHD))
        
        !Let er rip
        call WriteVars(IOVars,.false.,ResF)
        !Create link to latest restart
        write (lnResF, '(A,A,A,A)') trim(RCMApp%rcm_runid), ".mhd2imag.Res.", "XXXXX", ".h5"
        call MapSymLink(ResF,lnResF)

    end subroutine WriteMHD2IMagRestart

    subroutine ReadMHD2IMagRestart(RCMApp,nRes)
        type(rcm_mhd_t)  , intent(inout) :: RCMApp
        integer, intent(in) :: nRes

        type(IOVAR_T), dimension(MAXRCMIOVAR) :: IOVars
        character(len=strLen) :: ResF,nStr
        logical :: fExist
        real(rp), dimension(:,:), allocatable :: iopenX,toMHDX
        integer :: NLat,NLon

        NLat = RCMApp%nLat_ion
        NLon = RCMApp%nLon_ion

        !Get number string
        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') nRes
        endif

        write (ResF, '(A,A,A,A)') trim(RCMApp%rcm_runid), ".mhd2imag.Res.", trim(nStr), ".h5"
        write(*,*) 'Trying to read MHD2Imag restart from ', trim(ResF)
        inquire(file=ResF,exist=fExist)

        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open MHD2Imag restart file, skipping ...'
            return
        else
            write(*,*) 'Found MHD2Imag restart, reading ...'
        endif

    !Read data if still here
        call ClearIO(IOVars)
        call AddInVar(IOVars,"gcolat"      )
        call AddInVar(IOVars,"glong"       )
        call AddInVar(IOVars,"pot"         )
        call AddInVar(IOVars,"eng_avg"     )
        call AddInVar(IOVars,"flux"        )
        call AddInVar(IOVars,"nflx"        )
        call AddInVar(IOVars,"fac"         )
        call AddInVar(IOVars,"Pave"        )
        call AddInVar(IOVars,"Nave"        )
        call AddInVar(IOVars,"Vol"         )
        call AddInVar(IOVars,"X_bmin"      )
        call AddInVar(IOVars,"Bmin"        )
        call AddInVar(IOVars,"beta_average")
        call AddInVar(IOVars,"iopen"       )
        call AddInVar(IOVars,"Prcm"       )
        call AddInVar(IOVars,"Nrcm"       )
        call AddInVar(IOVars,"Npsph"      )
        call AddInVar(IOVars,"sigmap"     )
        call AddInVar(IOVars,"sigmah"     )
        call AddInVar(IOVars,"oxyfrac"    )
        call AddInVar(IOVars,"Percm"      )
        call AddInVar(IOVars,"latc"       )
        call AddInVar(IOVars,"lonc"       )
        call AddInVar(IOVars,"Lb"         )
        call AddInVar(IOVars,"Tb"         )
        call AddInVar(IOVars,"losscone"   )
        call AddInVar(IOVars,"radcurv"    )
        call AddInVar(IOVars,"wIMAG"      )
        call AddInVar(IOVars,"toMHD"      )
        !Get data
        call ReadVars(IOVars,.false.,ResF)
  
    !Unpack data
        !1D
        !Disabling reading geometric stuff (use recalculated value)
        !call IOArray1DFill(IOVars,"gcolat"      ,RCMApp%gcolat      )
        !call IOArray1DFill(IOVars,"glong"       ,RCMApp%glong       )

        !2D
        call IOArray2DFill(IOVars,"pot"         ,RCMApp%pot         )
        call IOArray2DFill(IOVars,"fac"         ,RCMApp%fac         )  
        call IOArray2DFill(IOVars,"Pave"        ,RCMApp%Pave        )   
        call IOArray2DFill(IOVars,"Nave"        ,RCMApp%Nave        ) 
        call IOArray2DFill(IOVars,"Vol"         ,RCMApp%Vol         )
        call IOArray2DFill(IOVars,"Bmin"        ,RCMApp%Bmin        )
        call IOArray2DFill(IOVars,"beta_average",RCMApp%beta_average)
        call IOArray2DFill(IOVars,"Prcm"        ,RCMApp%Prcm        )
        call IOArray2DFill(IOVars,"Nrcm"        ,RCMApp%Nrcm        )
        call IOArray2DFill(IOVars,"Npsph"       ,RCMApp%Npsph       )
        call IOArray2DFill(IOVars,"sigmap"      ,RCMApp%sigmap      )
        call IOArray2DFill(IOVars,"sigmah"      ,RCMApp%sigmah      )
        call IOArray2DFill(IOVars,"oxyfrac"     ,RCMApp%oxyfrac     )
        call IOArray2DFill(IOVars,"Percm"       ,RCMApp%Percm       )
        call IOArray2DFill(IOVars,"latc"        ,RCMApp%latc        )
        call IOArray2DFill(IOVars,"lonc"        ,RCMApp%lonc        )
        call IOArray2DFill(IOVars,"Lb"          ,RCMApp%Lb          )
        call IOArray2DFill(IOVars,"Tb"          ,RCMApp%Tb          )
        call IOArray2DFill(IOVars,"losscone"    ,RCMApp%losscone    )
        call IOArray2DFill(IOVars,"radcurv"     ,RCMApp%radcurv     )
        call IOArray2DFill(IOVars,"wIMAG"       ,RCMApp%wIMAG       )

        !3D
        call IOArray3DFill(IOVars,"eng_avg"     ,RCMApp%eng_avg     )    
        call IOArray3DFill(IOVars,"flux"        ,RCMApp%flux        )  
        if(ioExist(ResF,"nflx")) then
            call IOArray3DFill(IOVars,"nflx"        ,RCMApp%nflx        )
        endif
        call IOArray3DFill(IOVars,"X_bmin"      ,RCMApp%X_bmin      )

        !Weird data
        allocate(iopenX(NLat,NLon))
        allocate(toMHDX(NLat,NLon))
        call IOArray2DFill(IOVars,"iopen",iopenX)
        RCMApp%iopen = nint(iopenX)
        call IOArray2DFill(IOVars,"toMHD",toMHDX)
        RCMApp%toMHD = (toMHDX>0.5)

        !Handle some optional values
        if (ioExist(ResF,"N0")) then
            call ClearIO(IOVars)
            call AddInVar(IOVars,"N0")
            call ReadVars(IOVars,.false.,ResF)
            call IOArray2DFill(IOVars,"N0",RCMApp%N0)
        else
            RCMApp%N0 = 0.0
        endif

        if (ioExist(ResF,"TioTe0")) then
            call ClearIO(IOVars)
            call AddInVar(IOVars,"TioTe0")
            call ReadVars(IOVars,.false.,ResF)
            call IOArray2DFill(IOVars,"TioTe0",RCMApp%TioTe0)
        else
            RCMApp%TioTe0 = tiote_RCM
        endif
        write(*,*) 'Finished reading MHD2Imag restart ...'
    end subroutine ReadMHD2IMagRestart
end module rcm_mhd_io
