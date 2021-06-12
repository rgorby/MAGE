! Collection of data and objects for the voltron middle man

module voltapp
    use mixtypes
    use ebtypes
    use chmpdefs
    use starter
    use mhd2mix_interface
    use mix2mhd_interface
    use mhd2chmp_interface
    use chmp2mhd_interface
    use ebsquish
    use innermagsphere
    use dates
    use kronos
    use voltio
    use msphutils, only : RadIonosphere
    use gcminterp
    use gcmtypes
    
    implicit none

    contains

    !Initialize Voltron (after Gamera has already been initialized)
    subroutine initVoltron(vApp,gApp,optFilename)
        type(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: inpXML, kaijuRoot
        type(XML_Input_T) :: xmlInp
        type(TimeSeries_T) :: tsMJD
        real(rp) :: gTScl,tSpin,tIO
        logical :: doSpin

        if(present(optFilename)) then
            ! read from the prescribed file
            inpXML = optFilename
            call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")
        else
            !Find input deck
            call getIDeckStr(inpXML)

        endif
        !Start by shutting up extra ranks
        if (.not. vApp%isLoud) call xmlInp%BeQuiet()

        !Create XML reader
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Voltron',.true.)

        ! try to verify that the XML file has "Kaiju" as a root element
        kaijuRoot = ""
        call xmlInp%Get_Key_Val("/Gamera/sim/H5Grid",kaijuRoot)
        if(len(trim(kaijuRoot)) /= 0) then
            write(*,*) "The input XML appears to be of an old style."
            write(*,*) "As of June 12th, 2021 it needs a root element of <Kaiju>."
            write(*,*) "Please modify your XML config file by adding this line at the top:"
            write(*,*) "<Kaiju>"
            write(*,*) "and this line at the bottom:"
            write(*,*) "</Kaiju>"
            write(*,*) "OR (preferred) convert your configuration to an INI file and use"
            write(*,*) " the XMLGenerator.py script to create conforming XML files."
            write(*,*) "Please refer to the python script or"
            write(*,*) " the [Generating XML Files] wiki page for additional info."
            stop
        endif

        !Setup OMP if on separate node (otherwise can rely on gamera)
        if (vApp%isSeparate) then
            call SetOMP(xmlInp)
        endif

        ! read number of squish blocks
        call xmlInp%Set_Val(vApp%ebTrcApp%ebSquish%numSquishBlocks,"coupling/numSquishBlocks",4)

    !Initialize state information
        !Set file to read from and pass desired variable name to initTS
        call xmlInp%Set_Val(vApp%tilt%wID,"/Kaiju/Gamera/wind/tsfile","NONE")
        call vApp%tilt%initTS("tilt",doLoudO=.false.)
        vApp%symh%wID = vApp%tilt%wID
        call vApp%symh%initTS("symh",doLoudO=.false.)

        gTScl = gApp%Model%Units%gT0

        !Use MJD from time series
        tsMJD%wID = vApp%tilt%wID
        call tsMJD%initTS("MJD",doLoudO=.false.)
        gApp%Model%MJD0 = tsMJD%evalAt(0.0_rp) !Evaluate at T=0
        
    !Time options
        call xmlInp%Set_Val(vApp%tFin,'time/tFin',1.0_rp)
        !Sync Gamera to Voltron endtime
        gApp%Model%tFin = vApp%tFin/gTScl
        
        call vApp%IO%init(xmlInp,vApp%time,vApp%ts)
        
        !Pull numbering from Gamera
        vApp%IO%tsNext = gApp%Model%IO%tsNext
        
        !Force Gamera IO times to match Voltron IO
        call IOSync(vApp%IO,gApp%Model%IO,1.0/gTScl)

    !Shallow coupling
        call xmlInp%Set_Val(vApp%ShallowDT ,"coupling/dt" , 0.1_rp)
        vApp%TargetShallowDT = vApp%ShallowDT
        call xmlInp%Set_Val(vApp%doGCM, "coupling/doGCM",.false.)

        if (vApp%doGCM) then
            call init_gcm(vApp%gcm,gApp%Model%isRestart)
        end if

    !Deep coupling
        call xmlInp%Set_Val(vApp%DeepDT, "coupling/dtDeep", -1.0_rp)
        vApp%TargetDeepDT = vApp%DeepDT
        call xmlInp%Set_Val(vApp%rTrc,   "coupling/rTrc"  , 40.0)

        if (vApp%DeepDT>0) then
            vApp%doDeep = .true.
        else
            vApp%doDeep = .false.
        endif

        if(gApp%Model%isRestart) then
            call readVoltronRestart(vApp, xmlInp)
            vApp%IO%tOut = floor(vApp%time/vApp%IO%dtOut)*vApp%IO%dtOut
            vApp%IO%tRes = vApp%time + vApp%IO%dtRes
            vApp%IO%tsNext = vApp%ts
            if(vApp%isSeparate) then
                gApp%Model%ts = vApp%ts
                gApp%Model%t  = vApp%time/gTScl
                gApp% State%time  = gApp%Model%t
                gApp%oState%time  = gApp%Model%t-gApp%Model%dt
            endif
        else
            ! non-restart initialization
            !Check for spinup info
            call xmlInp%Set_Val(doSpin,"spinup/doSpin",.true.)
            ! Deep enabled, not restart, not spinup is an error. Restart or spinup is required
            if (vApp%doDeep .and. (.not. doSpin) ) then
                write(*,*) 'Spinup is required with deep coupling. Please enable the spinup/doSpin option. At least 1 minute of spinup is recommended.'
                stop
            endif
            if (doSpin) then
                call xmlInp%Set_Val(tSpin,"spinup/tSpin",7200.0) !Default two hours
                !Rewind Gamera time to negative tSpin (seconds)
                gApp%Model%t = -tSpin/gTScl
                !Reset State/oState
                gApp% State%time  = gApp%Model%t
                gApp%oState%time  = gApp%Model%t-gApp%Model%dt
                call xmlInp%Set_Val(tIO,"spinup/tIO",0.0) !Time of first restart and output
                gApp%Model%IO%tRes = tIO/gTScl
                gApp%Model%IO%tOut = tIO/gTScl
                vApp%IO%tRes = tIO
                vApp%IO%tOut = tIO
            endif
            vApp%time = gApp%Model%t*gTScl !Time in seconds
            vApp%ts   = gApp%Model%ts !Timestep
            vApp%MJD = T2MJD(vApp%time,gApp%Model%MJD0)
            vApp%ShallowT = vApp%time ! shallow coupling immediately
            !Set first deep coupling (defaulting to 0)
            call xmlInp%Set_Val(vApp%DeepT, "coupling/tDeep", 0.0_rp)
        endif

        if (vApp%doDeep) then
            !Whether to do fast eb-squishing
            call xmlInp%Set_Val(vApp%doQkSquish,"coupling/doQkSquish",.false.)
            call xmlInp%Set_Val(vApp%qkSquishStride,"coupling/qkSquishStride", 2)
            if(vApp%doQkSquish .and. popcnt(vApp%qkSquishStride) /= 1) then
                write(*,*) 'Quick Squish Stride must be a power of 2'
                stop
            endif
            call xmlInp%Set_Val(vApp%chmp2mhd%epsSquish,"ebsquish/epsSquish",0.05)

            !Verify that Gamera has location to hold source info
            if (.not. gApp%Model%doSource) then
                write(*,*) 'Must have GAMERA/source/doSource="T" when running inner magnetosphere model'
                stop
            endif

            !Verify CHIMP data has been set
            if (      (.not. xmlInp%Exists("chimp/units/uid"))     &
               & .or. (.not. xmlInp%Exists("chimp/fields/grType")) &
               & .or. (.not. xmlInp%Exists("chimp/domain/dtype")) ) then   
                write(*,*) 'Necessary CHIMP XML paramters not found, sort that out ...'
                stop
            endif
            
            !Initialize deep coupling type/inner magnetosphere model
            call InitInnerMag(vApp,gApp,xmlInp)
        endif

        !Check for dynamic coupling cadence
        call xmlInp%Set_Val(vApp%doDynCplDT ,"coupling/doDynDT" , .false.)
        
        if(present(optFilename)) then
            ! read from the prescribed file
            call initializeFromGamera(vApp, gApp, optFilename)
        else
            call initializeFromGamera(vApp, gApp)
        endif

        if(.not. vApp%isSeparate) then
            !Do first couplings if the gamera data is local and therefore uptodate
            call Tic("IonCoupling")
            call ShallowUpdate(vApp,gApp,vApp%time)
            call Toc("IonCoupling")
        
            if (vApp%doDeep .and. (vApp%time>=vApp%DeepT)) then
                call Tic("DeepCoupling")
                call DeepUpdate(vApp,gApp)
                call Toc("DeepCoupling")
            endif
        endif

        !Recalculate timestep
        gApp%Model%dt = CalcDT(gApp%Model,gApp%Grid,gApp%State)
        if (gApp%Model%dt0<TINY) gApp%Model%dt0 = gApp%Model%dt
        
        !Bring overview info
        call printConfigStamp()

        !Finally do first output stuff
        !console output
        if (vApp%isSeparate) then
            call consoleOutputVOnly(vApp,gApp,gApp%Model%MJD0)
        else
            call consoleOutputV(vApp,gApp)
        endif
        !file output
        if (.not. gApp%Model%isRestart) then
            if(vApp%isSeparate) then
                call fOutputVOnly(vApp,gApp)
            else
                call fOutputV(vApp, gApp)
            endif
        endif
    end subroutine initVoltron

    !Step Voltron if necessary (currently just updating state variables)
    subroutine stepVoltron(vApp, gApp)
        class(voltApp_T), intent(inout) :: vApp
        class(gamApp_T) , intent(in)    :: gApp

        vApp%time = gApp%Model%t*gApp%Model%Units%gT0 !Time in seconds
        vApp%MJD = T2MJD(vApp%time,gApp%Model%MJD0)
        vApp%ts = gApp%Model%ts

    end subroutine stepVoltron
    
    !Initialize Voltron app based on Gamera data
    subroutine initializeFromGamera(vApp, gApp, optFilename)
        type(voltApp_T), intent(inout) :: vApp
        type(gamApp_T), intent(inout) :: gApp
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: RunID
        type(TimeSeries_T) :: f107

        logical :: isRestart
        real(rp) :: maxF107,Rin
        integer :: n

        isRestart = gApp%Model%isRestart
        RunID = trim(gApp%Model%RunID)
        
        if(vApp%writeFiles) call InitVoltIO(vApp,gApp)
        
    !Remix from Gamera
        !Set mix default grid before initializing
        Rin = norm2(gApp%Grid%xyz(1,1,1,:)) !Inner radius
        call SetMixGrid0(Rin,gApp%Grid%Nkp)

        if(present(optFilename)) then
            ! read from the prescribed file
            call init_mix(vApp%remixApp%ion,[NORTH, SOUTH],optFilename=optFilename,RunID=RunID,isRestart=isRestart,nRes=vApp%IO%nRes,optIO=vApp%writeFiles)
        else
            call init_mix(vApp%remixApp%ion,[NORTH, SOUTH],RunID=RunID,isRestart=isRestart,nRes=vApp%IO%nRes,optIO=vApp%writeFiles)
        endif
        vApp%remixApp%ion%rad_iono_m = RadIonosphere() * gApp%Model%units%gx0 ! [Rp] * [m/Rp]

        !Set F10.7 from time series (using max)
        f107%wID = vApp%tilt%wID
        call f107%initTS("f10.7",doLoudO=.false.)
        maxF107 = f107%getMax()
        

        do n=1,2
            vApp%remixApp%ion(n)%P%f107 = maxF107
        enddo
        write(*,*) 'Using F10.7 = ', maxF107
        write(*,*) 'Using MJD0  = ', gApp%Model%MJD0

        call init_mhd2Mix(vApp%mhd2mix, gApp, vApp%remixApp)
        call init_mix2Mhd(vApp%mix2mhd, vApp%remixApp, gApp)
        !vApp%mix2mhd%mixOutput = 0.0
        
    !CHIMP (TRC) from Gamera
        if (vApp%doDeep) then
            !Verify that there's some place to put deep coupling info
            if (.not. gApp%Model%doSource) then
                write(*,*) 'For deep coupling, GAMERA/source/doSouce="T" must be set!'
                write(*,*) 'Exiting ...'
                stop
            endif

            ! initialize chimp
            associate(ebTrcApp=>vApp%ebTrcApp)
            if (present(optFilename)) then
                call init_volt2Chmp(ebTrcApp,gApp,optFilename=optFilename)
            else
                call init_volt2Chmp(ebTrcApp,gApp)
            endif

            call init_mhd2Chmp(vApp%mhd2chmp, gApp, ebTrcApp)
            call init_chmp2Mhd(vApp%chmp2mhd, ebTrcApp, gApp)

            end associate
            vApp%iDeep = ShellBoundary(gApp%Model,gApp%Grid,vApp%rTrc)
        endif !doDeep

    end subroutine initializeFromGamera

!----------
!Shallow coupling stuff
    subroutine ShallowUpdate(vApp, gApp, time)
        type(gamApp_T), intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp
        real(rp) :: time

        ! convert gamera data to mixInput
        call Tic("G2R")
        call convertGameraToRemix(vApp%mhd2mix, gApp, vApp%remixApp)
        call Toc("G2R")

        ! run remix
        call Tic("ReMIX")
        call runRemix(vApp, time)
        call Toc("ReMIX")

        ! convert mixOutput to gamera data
        call Tic("R2G")
        call convertRemixToGamera(vApp%mix2mhd, vApp%remixApp, gApp)
        call Toc("R2G")

        vApp%ShallowT = vApp%ShallowT + vApp%ShallowDT

    end subroutine ShallowUpdate

    subroutine runRemix(vApp, time)
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: time
        real(rp) :: curTilt

        ! convert gamera inputs to remix
        if (vApp%doDeep) then
            call mapIMagToRemix(vApp%imag2mix,vApp%remixApp)
        endif
        call mapGameraToRemix(vApp%mhd2mix, vApp%remixApp)

        ! determining the current dipole tilt
        call vApp%tilt%getValue(vApp%time,curTilt)

        if (vApp%doGCM .and. time >=0 .and. .not.(vApp%gcm%isRestart)) then
            call coupleGCM2MIX(vApp%gcm,vApp%remixApp%ion,vApp%doGCM,mjd=vApp%MJD,time=vApp%time)
        end if

        ! solve for remix output
        if (time<=0) then
            call run_mix(vApp%remixApp%ion,curTilt,doModelOpt=.false.)
        else if (vApp%doGCM) then
            call run_mix(vApp%remixApp%ion,curTilt,gcm=vApp%gcm)
        else
            call run_mix(vApp%remixApp%ion,curTilt,doModelOpt=.true.)
        endif

        ! get stuff from mix to gamera
        call mapRemixToGamera(vApp%mix2mhd, vApp%remixApp)

    end subroutine runRemix

!----------
!Deep coupling stuff (time coming from vApp%time, so in seconds)
    subroutine DeepUpdate(vApp, gApp)
        type(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        if (.not. vApp%doDeep) then
            !Why are you even here?
            return
        endif

        call PreDeep(vApp, gApp)
          call DoImag(vApp)
          call SquishStart(vApp)
            call Squish(vApp) ! do all squish blocks here
          call SquishEnd(vApp)
        call PostDeep(vApp, gApp)

    end subroutine DeepUpdate

    subroutine PreDeep(vApp, gApp)
        type(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Update coupling time now so that voltron knows what to expect
        vApp%DeepT = vApp%DeepT + vApp%DeepDT

        !Update i-shell to trace within in case rTrc has changed
        vApp%iDeep = ShellBoundary(gApp%Model,gApp%Grid,vApp%rTrc)
        
        !Pull in updated fields to CHIMP
        call Tic("G2C")
        call convertGameraToChimp(vApp%mhd2chmp,gApp,vApp%ebTrcApp)
        call Toc("G2C")

    end subroutine

    subroutine DoImag(vApp)
        class(voltApp_T), intent(inout) :: vApp

        !Advance inner magnetosphere model to tAdv
        call Tic("InnerMag")
        call vApp%imagApp%doAdvance(vApp,vApp%DeepT)
        call Toc("InnerMag")

    end subroutine

    subroutine PostDeep(vApp, gApp)
        type(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Now use imag model and squished coordinates to fill Gamera source terms
        call Tic("IM2G")
        call InnerMag2Gamera(vApp,gApp)
        call Toc("IM2G")

    end subroutine

    subroutine CheckQuickSquishError(vApp, gApp, x2Err, x4Err)
        type(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp
        real(rp), intent(out) :: x2Err, x4Err

        integer :: i,j,k,baseQkSqStr
        logical :: baseDoQkSq
        real(rp), dimension(:,:,:,:), allocatable :: baseXyzSquish
        real(rp), dimension(2) :: posErr

        associate(ebGr=>vApp%ebTrcApp%ebState%ebGr)

        allocate(baseXyzSquish,  MOLD=vApp%chmp2mhd%xyzSquish)
        baseQkSqStr = vApp%qkSquishStride
        baseDoQkSq = vApp%doQkSquish

        ! squish with no quick squish stride
        vApp%qkSquishStride = 1
        vApp%doQkSquish = .false.
        call DeepUpdate(vApp, gApp)
        baseXyzSquish = vApp%chmp2mhd%xyzSquish

        ! squish with 2x quick squish stride
        vApp%qkSquishStride = 2
        vApp%doQkSquish = .true.
        call DeepUpdate(vApp, gApp)
        x2Err = 0
        do i=ebGr%is,vApp%iDeep+1
            do j=ebGr%js,ebGr%je+1
                do k=ebGr%ks,ebGr%ke+1
                    posErr = abs(baseXyzSquish(i,j,k,:) - vApp%chmp2mhd%xyzSquish(i,j,k,:))
                    if(posErr(2) > PI) posErr(2) = 2*PI - posErr(2)
                    x2Err = x2Err + NORM2(posErr)
                enddo
            enddo
        enddo

        ! squish with 4x quick squish stride
        vApp%qkSquishStride = 4
        call DeepUpdate(vApp, gApp)
        x4Err = 0
        do i=ebGr%is,vApp%iDeep+1
            do j=ebGr%js,ebGr%je+1
                do k=ebGr%ks,ebGr%ke+1
                    posErr = abs(baseXyzSquish(i,j,k,:) - vApp%chmp2mhd%xyzSquish(i,j,k,:))
                    if(posErr(2) > PI) posErr(2) = 2*PI - posErr(2)
                    x4Err = x4Err + NORM2(posErr)
                enddo
            enddo
        enddo

        vApp%qkSquishStride = baseQkSqStr
        vApp%doQkSquish = baseDoQkSq
        deallocate(baseXyzSquish)

        end associate

    end subroutine

    !Initialize CHIMP data structure
    subroutine init_volt2Chmp(ebTrcApp,gApp,optFilename)
        type(ebTrcApp_T), intent(inout) :: ebTrcApp
        type(gamApp_T), intent(in) :: gApp
        character(len=*), intent(in), optional     :: optFilename

        character(len=strLen) :: xmlStr
        type(XML_Input_T) :: inpXML
        real(rp) :: xyz0(NDIM)

    !Create input XML object
        if (present(optFilename)) then
            xmlStr = trim(optFilename)
        else
            call getIDeckStr(xmlStr)
        endif
        inpXML = New_XML_Input(trim(xmlStr),"Kaiju/Chimp",.true.)

    !Initialize model
        associate(Model=>ebTrcApp%ebModel,ebState=>ebTrcApp%ebState,ebGr=>ebTrcApp%ebState%ebGr,Gr=>gApp%Grid)
        call setUnits (Model,inpXML)
        Model%T0   = 0.0
        Model%tFin = 0.0
        Model%dt   = 0.0
        Model%t    = 0.0
        call setInterpolation(Model,inpXML)
        Model%doMHD = .true.
        call inpXML%Set_Val(Model%epsds,'tracer/epsds',1.0e-2)    
        call setBackground(Model,inpXML)

    !Initialize ebState
        !CHIMP grid is initialized from Gamera's active corners
        call ebInit_fromMHDGrid(Model,ebState,inpXML,Gr%xyz(Gr%is:Gr%ie+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1,1:NDIM))
        !Replace CHIMP 8-point average centers w/ more accurate Gamera quadrature centers        
        ebGr%xyzcc(ebGr%is:ebGr%ie,ebGr%js:ebGr%je,ebGr%ks:ebGr%ke,:) = Gr%xyzcc(Gr%is:Gr%ie,Gr%js:Gr%je,Gr%ks:Gr%ke,:)

        call InitLoc(Model,ebState%ebGr,inpXML)

        !Do simple test to make sure locator is reasonable
        xyz0 = Gr%xyz(Gr%is+1,Gr%js,Gr%ks,:)
        if (.not. inDomain(xyz0,Model,ebState%ebGr) ) then
            write(*,*) 'Configuration error: CHIMP Domain incorrect'
            stop
        endif

        end associate

    end subroutine init_volt2Chmp

end module voltapp

