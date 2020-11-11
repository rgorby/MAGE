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
    use imagtubes
    use rcm_mhd_interfaces
    use rcm_mix_interface
    use clocks
    use kronos
    use rcm_mhd_mod, ONLY : rcm_mhd
    use rcm_mhd_io
    use rcmdefs, only : DenPP0
    
    
    implicit none

    real(rp), private :: rTrc0 = 2.0 !Padding factor for RCM domain to ebsquish radius
    logical , private, parameter :: doKillRCMDir = .true. !Whether to always kill RCMdir before starting
    integer, parameter, private :: MHDPad = 0 !Number of padding cells between RCM domain and MHD ingestion
    logical , private :: doWolfLim  = .false. !Whether to do wolf-limiting
    logical , private :: doIsoWolf  = .true. !Preserve RCM temp. when doing wolf-limiting
    logical , private :: doBounceDT = .true. !Whether to use Alfven bounce in dt-ingest
    logical , private :: doWIMTScl = .false. !Whether to modulate ingestion timescale by wIM
    real(rp), private :: nBounce = 1.0 !Scaling factor for Alfven transit

    real(rp), private :: wIM_C = 0.0 !Critical wIM for MHD ingestion inclusion

    real(rp), dimension(:,:), allocatable, private :: mixPot


    type, extends(innerMagBase_T) :: rcmIMAG_T

        ! rcm coupling variable
        type(rcm_mhd_T) :: rcmCpl

        contains

        ! over-ride the base functions with RCM versions
        procedure :: doInit => initRCM
        procedure :: doAdvance => advanceRCM
        procedure :: doEval => EvalRCM
        procedure :: doIO => doRCMIO
        procedure :: doConIO => doRCMConIO
        procedure :: doRestart => doRCMRestart

    end type

    contains

    !Initialize RCM inner magnetosphere model
    subroutine initRCM(imag,iXML,isRestart,rad_planet_m,rad_iono_m,M0g,vApp)
        class(rcmIMAG_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        real(rp), intent(in) :: rad_planet_m,rad_iono_m, M0g ! Specific planet parameters
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

        call iXML%Set_Val(doWolfLim ,"/gamera/source/doWolfLim" ,doWolfLim )
        call iXML%Set_Val(doBounceDT,"/gamera/source/doBounceDT",doBounceDT)
        call iXML%Set_Val(nBounce   ,"/gamera/source/nBounce"   ,nBounce   )

        if (isRestart) then
            if (doKillRCMDir) then
                !Kill RCMFiles directory even on restart
                call ResetRCMDir()
            endif

            !Get t0 and nRes necessary for RCM restart
            call RCMRestartInfo(RCMApp,iXML,t0)

            write(*,*) 'Restarting RCM @ t = ', t0
            vApp%time = t0 !Set vApp's time to correct value from restart
            call rcm_mhd(t0,dtCpl,RCMApp,RCMRESTART,iXML=iXML)
            doColdstart = .false. ! set to false if it is a restart
        else
            t0 = vApp%time
            call ResetRCMDir()
            write(*,*) 'Initializing RCM ...'
            call InitRCMICs(imag,vApp,iXML)
            call rcm_mhd(t0,dtCpl,RCMApp,RCMINIT,iXML=iXML)
            doColdStart = .true.
        endif

        call init_rcm_mix(RCMApp,imag2mix)

        !Start up IO
        call initRCMIO(RCMApp,isRestart)

        end associate

        contains
            subroutine ResetRCMDir()
                write(*,*) 'Reset RCMfiles/ ...'
                CALL SYSTEM("rm -rf RCMfiles > /dev/null 2>&1")
                CALL SYSTEM("mkdir RCMfiles > /dev/null 2>&1")
                CALL SYSTEM("touch RCMfiles/rcm.printout > /dev/null 2>&1")
                CALL SYSTEM("touch RCMfiles/rcm.index > /dev/null 2>&1")
            end subroutine ResetRCMDir
            
    end subroutine initRCM

    !Setup ICs to pass to RCM if asked to
    subroutine InitRCMICs(imag,vApp,iXML)
        class(rcmIMAG_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        type(voltApp_T), intent(in) :: vApp

        real(rp) :: t0
        
        call iXML%Set_Val(RCMICs%doIC,"imag/doInit",.false.)
        t0 = 0.0
        if (RCMICs%doIC) then
            !Want initial dst0
            RCMICs%dst0 = GetSWVal("symh",vApp%tilt%wID,t0)
            call iXML%Set_Val(RCMICs%ktRC,"imag/ktRC",30.0)
            RCMICs%vSW = abs(GetSWVal("Vx",vApp%tilt%wID,t0))
            RCMICs%dSW = GetSWVal("D",vApp%tilt%wID,t0)

            !Set PS values (see Borovsky paper)
            RCMICs%dPS  = 0.292*(RCMICs%dSW**0.49)
            RCMICs%kTPS = -3.65 + 0.0190*RCMICs%vSW*1.0e-3 !m/s=>km/s
            RCMICs%kTPS = max(RCMICs%kTPS,TINY)

            !Tune RC pressure profile, using just dst(T=0) (will try to do better later)
            call SetQTRC(RCMICs%dst0)
        else
            !Zero out any additional ring current
            call SetQTRC(0.0_rp)
        endif

        contains

        function GetSWVal(vID,fID,t0) result(qSW)
            character(len=*), intent(in) :: vID,fID
            real(rp), intent(in) :: t0
            real(rp) :: qSW

            type(TimeSeries_T) :: tsQ
            tsQ%wID = trim(fID)
            call tsQ%initTS(trim(vID),doLoudO=.false.)
            qSW = tsQ%evalAt(t0)
        end function GetSWVal

    end subroutine InitRCMICs

    !Advance RCM from Voltron data
    subroutine AdvanceRCM(imag,vApp,tAdv)
        class(rcmIMAG_T), intent(inout) :: imag
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        integer :: i,j,n,nStp
        real(rp) :: colat,lat,lon
        real(rp) :: dtAdv
        type(RCMTube_T) :: ijTube

        real(rp) :: maxRad
        logical :: isLL,doHackIC
        
        associate(RCMApp => imag%rcmCpl)

        RCMApp%llBC  = vApp%mhd2chmp%lowlatBC
        RCMApp%dtCpl = vApp%DeepDT
        call Tic("MAP_RCMMIX")
    !Get potential from mix
        call map_rcm_mix(vApp,mixPot)
        call Toc("MAP_RCMMIX")

        call Tic("RCM_TUBES")
    !Load RCM tubes
       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP schedule(dynamic) &
       !$OMP private(i,j,colat,lat,lon,isLL,ijTube)
        do j=1,RCMApp%nLon_ion
            do i=1,RCMApp%nLat_ion
                colat = RCMApp%gcolat(i)
                lat = PI/2 - colat
                lon = RCMApp%glong(j)
                
                !Decide if we're below low-lat BC or not
                isLL = (lat <= RCMApp%llBC)
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
                RCMApp%losscone(i,j)     = ijTube%losscone
                RCMApp%Lb(i,j)           = ijTube%Lb
                RCMApp%radcurv(i,j)      = ijTube%rCurv
                !Get some kind of bounce timscale, either real integrated value or lazy average
                !RCMApp%Tb(i,j)           = AlfvenBounce(ijTube%Nave,ijTube%bmin,ijTube%Lb)
                RCMApp%Tb(i,j)           = ijTube%Tb
                RCMApp%wIMAG(i,j)        = ijTube%wIMAG

                !mix variables are stored in this order (longitude,colatitude), hence the index flip
                RCMApp%pot(i,j)          = mixPot(j,i)

                !Set composition
                RCMApp%oxyfrac(i,j)      = 0.0
                
            enddo
        enddo
        doHackIC = (vApp%time <= vApp%DeepDT) .and. RCMICs%doIC !Whether to hack MHD/RCM coupling for ICs

        call Toc("RCM_TUBES")
    !Do any tube hacking needed before sending tubes to RCM
        if (doHackIC) then
        !Tune values to send to RCM for its cold start
            !Setup quiet time ring current to hit target using both current BSDst and target dst
            !Replacing first RC estimate w/ Dst at end of blow-in period
            call SetQTRC(RCMICs%dst0-vApp%BSDst)
            call HackTubes(RCMApp,vApp)
        endif

        !Coverup some bad tubes
        call TrickyTubes(RCMApp)
        
        !Smooth out FTV/potential on tubes b/c RCM will take gradient
        call SmoothTubes(RCMApp,vApp)

    !Advance from vApp%time to tAdv
        call Tic("AdvRCM")
        dtAdv = tAdv-vApp%time !RCM-DT
        if (doColdstart) then
            write(*,*) 'Cold-starting RCM @ t = ', vApp%time
            call rcm_mhd(vApp%time,dtAdv,RCMApp,RCMCOLDSTART)
            doColdstart = .false.
        else
            call rcm_mhd(vApp%time,dtAdv,RCMApp,RCMADVANCE)
        end if

        call Toc("AdvRCM")

    !Set ingestion region
        if (doHackIC) then
            !For ICs ingest from entire closed field region just this once
            RCMApp%toMHD = .not. (RCMApp%iopen == RCMTOPOPEN)
        else
            call SetIngestion(RCMApp)
            !Find maximum extent of RCM domain (RCMTOPCLOSED but not RCMTOPNULL)
            maxRad = maxval(norm2(RCMApp%X_bmin,dim=3),mask=(RCMApp%iopen == RCMTOPCLOSED))
            maxRad = maxRad/Rp_m
            vApp%rTrc = rTrc0*maxRad
        endif

    !Pull data from RCM state for conductance calculations
        !NOTE: this is not the closed field region, this is actually the RCM domain
        vApp%imag2mix%inIMag = (RCMApp%iopen == RCMTOPCLOSED)
        vApp%imag2mix%latc = RCMApp%latc
        vApp%imag2mix%lonc = RCMApp%lonc

    ! electrons precipitation
        vApp%imag2mix%eflux = RCMApp%flux(:,:,1)
        vApp%imag2mix%eavg  = RCMApp%eng_avg(:,:,1)
    ! ion precipitation
        vApp%imag2mix%iflux = RCMApp%flux(:,:,2)
        vApp%imag2mix%iavg  = RCMApp%eng_avg(:,:,2)

        vApp%imag2mix%isFresh = .true.

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
                dTb = (L*Rp_m*1.0e-3)/Va
            end function AlfvenBounce
    end subroutine AdvanceRCM

    !Set region of RCM grid that's "good" for MHD ingestion
    subroutine SetIngestion(RCMApp)
        type(rcm_mhd_T), intent(inout) :: RCMApp

        integer , dimension(:), allocatable :: jBnd
        integer :: i,j
        logical :: inMHD,isClosed
        
        RCMApp%toMHD(:,:) = .false.
        !Testing lazy quick boundary
        allocate(jBnd (  RCMApp%nLon_ion  ))
        !Now find nominal current boundary
        jBnd(:) = RCMApp%nLat_ion-1
        do j=1,RCMApp%nLon_ion
            do i = RCMApp%nLat_ion,1,-1
                inMHD = RCMApp%toMHD(i,j)
                isClosed = (RCMApp%iopen(i,j) == RCMTOPCLOSED)
                if ( .not. isClosed ) then
                    jBnd(j) = min(i+1+MHDPad,RCMApp%nLat_ion)
                    exit
                endif

            enddo !i loop
            RCMApp%toMHD(:,j) = .false.
            RCMApp%toMHD(jBnd(j):,j) = .true.

        enddo
        
        !Finally, upscale ingestion timescales based on wImag
        if (doWIMTScl) then
            where (RCMApp%wImag>TINY)
                RCMApp%Tb = RCMApp%Tb/RCMApp%wImag
            endwhere
        endif

    end subroutine SetIngestion

    !Enforce Wolf-limiting on an MHD/RCM thermodynamic state
    subroutine WolfLimit(nrc,prc,npp,nmhd,pmhd,beta,nlim,plim)
        real(rp), intent(in)  :: nrc,prc,npp,nmhd,pmhd,beta
        real(rp), intent(out) :: nlim,plim

        real(rp) :: nrcm,prcm
        real(rp) :: alpha

        nlim = 0.0
        plim = 0.0

        !Test RC contribution
        if ( (nrc>TINY) .and. (prc>TINY) ) then
            nrcm = nrc !Ring current contribution
            prcm = prc !Only pressure contribution from RC
        else
            nrcm = 0.0
            prcm = 0.0
        endif

        !Test plasmasphere density
        if (npp >= DenPP0) then
            !Plasmasphere good, add it to total rcm density
            nrcm = nrcm + npp
        endif

        !Now have total density/pressure from RC+PP
        if (.not. doWolfLim) then
            nlim = nrcm
            plim = prcm
            return
        endif

        !Get scaling term
        alpha = 1.0 + beta*5.0/6.0
        if ( (nrcm>TINY) .and. (prcm>TINY) ) then

            !Apply wolf-limiting
            plim = ( pmhd*(alpha-1) + prcm )/alpha
            nlim = nrcm - (0.5*beta/alpha)*(nmhd/pmhd)*(prcm-pmhd)
            if (nlim <= TINY) then
                !Something went bad, nuke both
                nlim = 0.0
                plim = 0.0
            endif
            return
        endif

        !If still here at least one of density/pressure is degenerate
        !Test for den only, i.e. plasmasphere
        if (nrcm>TINY) then
            !Take raw density
            nlim = nrcm
            plim = 0.0
            return
        endif
        !Return null values
        nlim = 0.0
        plim = 0.0

    end subroutine WolfLimit

    !Evaluate eq map at a given point
    !Returns density (#/cc) and pressure (nPa)
    subroutine EvalRCM(imag,x1,x2,t,imW,isEdible)
        class(rcmIMAG_T), intent(inout) :: imag
        real(rp), intent(in) :: x1,x2,t
        real(rp), intent(out) :: imW(NVARIMAG)
        logical, intent(out) :: isEdible

        real(rp) :: colat,nrcm,prcm,npp,pScl,beta,pmhd,nmhd,wIM
        real(rp) :: plim,nlim,Tb
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
        beta =  RCMApp%beta_average(ij0(1),ij0(2))
        pmhd = rcmPScl*RCMApp%Pave (ij0(1),ij0(2))
        nmhd = rcmNScl*RCMApp%Nave (ij0(1),ij0(2))
        wIM  =         RCMApp%wImag(ij0(1),ij0(2))
        Tb   =         RCMApp%Tb   (ij0(1),ij0(2))

        nlim = 0.0
        plim = 0.0

        if (doWolfLim) then
            call WolfLimit(nrcm,prcm,npp,nmhd,pmhd,beta,nlim,plim)
        else
            !Just lazyily use same function w/ beta=0
            call WolfLimit(nrcm,prcm,npp,nmhd,pmhd,0.0_rp,nlim,plim)
        endif

        !Store values
        imW(IMDEN) = nlim
        imW(IMPR)  = plim

        if (doBounceDT) then
            !Use Alfven bounce timescale
            imW(IMTSCL) = nBounce*RCMApp%Tb(ij0(1),ij0(2))
        else if (doWIMTScl) then
            !Use current coupling timescale, modulated by wImag
            imW(IMTSCL) = RCMApp%dtCpl/max( RCMApp%wImag(ij0(1),ij0(2)), TINY )
        endif

        imW(IMX1)   = rad2deg*lat
        imW(IMX2)   = rad2deg*lon

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
!Use findloc w/ intel for speed
#ifdef __INTEL_COMPILER
            iC = findloc(gcolat >= colat,.true.,dim=1) - 1
#else 
!Bypass as findloc does not work for gfortran<9    
           !Work-around code        
            do n=1,nLat
                if (gcolat(n) >= colat) exit
            enddo
            iC = n-1
#endif
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
    

    !Some quick console diagnostics for RCM info
    subroutine doRCMConIO(imag,MJD,time)
        class(rcmIMAG_T), intent(inout) :: imag
        real(rp), intent(in) :: MJD, time

        integer :: i0,j0,maxIJ(2)

        real(rp) :: maxPRCM,maxD,maxDP,maxPMHD,maxDMHD,maxL,maxMLT,maxBeta
        real(rp) :: limP,limD,wTrust,wTMin,maxT,maxWT,maxLam

        associate(RCMApp => imag%rcmCpl)
    !Start by getting some data
        !Pressure peak info
        maxIJ = maxloc(RCMApp%Prcm,mask=RCMApp%toMHD)
        i0 = maxIJ(1); j0 = maxIJ(2)

        maxPRCM  = RCMApp%Prcm (i0,j0)*rcmPScl
        maxPMHD  = RCMApp%Pave (i0,j0)*rcmPScl
        maxBeta  = RCMApp%beta_average(i0,j0)
        maxD     = RCMApp%Nrcm (i0,j0)*rcmNScl
        maxDMHD  = RCMApp%Nave (i0,j0)*rcmNScl
        maxDP = RCMApp%Npsph(i0,j0)*rcmNScl

        maxL = norm2(RCMApp%X_bmin(i0,j0,XDIR:YDIR))/Rp_m
        maxMLT = atan2(RCMApp%X_bmin(i0,j0,YDIR),RCMApp%X_bmin(i0,j0,XDIR))*180.0/PI
        if (maxMLT<0) maxMLT = maxMLT+360.0
        maxWT = 100*RCMApp%wImag(i0,j0)

        maxLam = (1.0e-3)*RCMApp%MaxAlam/( (RCMApp%vol(i0,j0)*1.0e-9)**(2.0/3.0) ) !Max RCM energy channel [keV]

        !Get pressure weighted confidence
        wTrust = sum(RCMApp%Prcm*RCMApp%wIMAG,mask=RCMApp%toMHD)/sum(RCMApp%Prcm,mask=RCMApp%toMHD)
        wTrust = 100.0*wTrust
        !Get min confidence in MHD domain
        wTMin = 100.0*minval(RCMApp%wIMAG,mask=RCMApp%toMHD)

    !Do some output
        if (maxPRCM<TINY) return

        write(*,*) ANSIYELLOW
        write(*,*) 'RCM'
        !write (*, '(a, f8.2,a,f6.2,a,f6.2,a)')      '  Trust    = ' , wTrust, '% (P-AVG) / ', wTMin, '% (MIN) / ', maxWT, '% (@ MAX)'
        if (doWolfLim) then
            call WolfLimit(maxD,maxPRCM,maxDP,maxDMHD,maxPMHD,maxBeta,limD,limP)
            write (*, '(a, f8.3,a,f8.3,a)')      '  Max RC-P = ' , maxPRCM, ' (RCM) / ', limP, ' (LIM) [nPa]'
            !Get temperature from RCM raw (unlimited)
            maxT = DP2kT(maxD,maxPRCM)
            !maxT = DP2kT(limD,limP)
        else
            write (*, '(a,1f8.3,a)')             '  Max RC-P = ' , maxPRCM, ' [nPa]'
            maxT = DP2kT(maxD,maxPRCM)
        endif

        write (*, '(a,2f8.3,a)')             '   @ L/MLT = ' , maxL, maxMLT, ' [deg]'
        if (doWolfLim) then
            !Add limited density
            write (*, '(a, f8.3,a,f8.3,a,f8.3,a)')      '      w/ D = ' , maxD, ' (RC) / ', maxDP, ' (PSPH) / ', limD, ' (LIM) [#/cc]' 
        else
            write (*, '(a, f8.3,a,f8.3,a)')      '      w/ D = ' , maxD, ' (RC) / ', maxDP, ' (PSPH) [#/cc]'
        endif
        write (*, '(a,1f8.3,a)')             '      w/ T = ' , maxT, ' [keV]'

        write (*, '(a,1f8.3,a)')             '  Max RC-D = ' , maxval(RCMApp%Nrcm,mask=RCMApp%toMHD)*rcmNScl,' [#/cc]'
        write(*,'(a)',advance="no") ANSIRESET!, ''

        end associate

    end subroutine doRCMConIO
end module rcmimag
