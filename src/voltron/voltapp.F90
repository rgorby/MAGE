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
    
    implicit none

    contains

    !Initialize Voltron (after Gamera has already been initialized)
    subroutine initVoltron(vApp,gApp,optFilename)
        type(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp
        character(len=*), optional, intent(in) :: optFilename

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        type(TimeSeries_T) :: tsMJD
        real(rp) :: gTScl,tSpin,tIO
        logical :: doSpin,doDelayIO

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

#ifdef _OPENMP
        if (vApp%isLoud) then
            write(*,*) 'Voltron running threaded'
            write(*,*) '   # Threads = ', omp_get_max_threads()
            write(*,*) '   # Cores   = ', omp_get_num_procs()
        endif
#else
        if (vApp%isLoud) then
            write (*,*) 'Voltron running without threading'
        endif
#endif

    !Create XML reader
        xmlInp = New_XML_Input(trim(inpXML),'Voltron',.true.)

    !Initialize state information
        !Set file to read from and pass desired variable name to initTS
        call xmlInp%Set_Val(vApp%tilt%wID,"/Gamera/wind/tsfile","NONE")
        call vApp%tilt%initTS("tilt")

        gTScl = gApp%Model%Units%gT0

        !Check for spinup info
        call xmlInp%Set_Val(doSpin,"spinup/doSpin",.true.)
        tIO = 0.0
        doDelayIO = .false.
        if (doSpin .and. (.not. gApp%Model%isRestart)) then
            !Doing spinup and not a restart
            call xmlInp%Set_Val(tSpin,"spinup/tSpin",3600.0) !Default two hours
            !Rewind Gamera time to negative tSpin (seconds)
            gApp%Model%t = -tSpin/gTScl 
            !Reset State/oState
            gApp% State%time  = gApp%Model%t
            gApp%oState%time  = gApp%Model%t-gApp%Model%dt

            doDelayIO = .true.
            call xmlInp%Set_Val(tIO,"spinup/tIO",0.0) !Time of first restart
        endif

        vApp%time = gApp%Model%t*gTScl !Time in seconds
        vApp%ts   = gApp%Model%ts !Timestep

        !Use MJD from time series
        tsMJD%wID = vApp%tilt%wID
        call tsMJD%initTS("MJD")
        gApp%Model%MJD0 = tsMJD%evalAt(0.0_rp) !Evaluate at T=0
        
        vApp%MJD = T2MJD(vApp%time,gApp%Model%MJD0)

    !Time options
        call xmlInp%Set_Val(vApp%tFin,'time/tFin',1.0_rp)
        !Sync Gamera to Voltron endtime
        gApp%Model%tFin = vApp%tFin/gTScl
        
    !IO/Restart options
        if (doDelayIO) then
            call vApp%IO%init(xmlInp,tIO)
        else
            call vApp%IO%init(xmlInp,vApp%time)
        endif

        !Pull numbering from Gamera
        vApp%IO%nRes = gApp%Model%IO%nRes
        vApp%IO%nOut = gApp%Model%IO%nOut
        !Force Gamera IO times to match Voltron IO
        call IOSync(vApp%IO,gApp%Model%IO,1.0/gTScl)

    !Shallow coupling
        !Start shallow coupling immediately
        vApp%ShallowT = vApp%time
        call xmlInp%Set_Val(vApp%ShallowDT ,"coupling/dt" , 0.1_rp)

    !Deep coupling
        vApp%DeepT = 0.0_rp
        call xmlInp%Set_Val(vApp%DeepDT, "coupling/dtDeep", -1.0_rp)
        call xmlInp%Set_Val(vApp%rDeep,  "coupling/rDeep" , 10.0_rp)
        call xmlInp%Set_Val(vApp%rTrc,   "coupling/rTrc"  , 2.0*vApp%rDeep)

        if (vApp%DeepDT>0) then
            vApp%doDeep = .true.
        else
            vApp%doDeep = .false.
        endif

        if (vApp%doDeep) then
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
             
            !Set first deep coupling (defaulting to 0)
            call xmlInp%Set_Val(vApp%DeepT, "coupling/tDeep", 0.0_rp)
            !Initialize deep coupling type/inner magnetosphere model
            !call InitInnerMag(vApp,gApp%Model%isRestart,xmlInp)
            call InitInnerMag(vApp,gApp,xmlInp)
        endif

        if(present(optFilename)) then
            ! read from the prescribed file
            call initializeFromGamera(vApp, gApp, optFilename)
        else
            call initializeFromGamera(vApp, gApp)
        endif

        !Do first couplings
        call Tic("IonCoupling")
        call ShallowUpdate(vApp,gApp,vApp%time)
        call Toc("IonCoupling")
        
        if (vApp%doDeep .and. (vApp%time>=vApp%DeepT)) then
            call Tic("DeepCoupling")
            call DeepUpdate(vApp,gApp,vApp%time)
            call Toc("DeepCoupling")
        endif

        !Recalculate timestep
        gApp%Model%dt = CalcDT(gApp%Model,gApp%Grid,gApp%State)
        if (gApp%Model%dt0<TINY) gApp%Model%dt0 = gApp%Model%dt
        
        !Finally do first output stuff
        !console output
        if(vApp%isSeparate) then
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
        real(rp) :: maxF107
        integer :: n

        isRestart = gApp%Model%isRestart
        RunID = trim(gApp%Model%RunID)
        
        call InitVoltIO(vApp,gApp)
        
    !Remix from Gamera
        if(present(optFilename)) then
            ! read from the prescribed file
            call init_mix(vApp%remixApp%ion,[NORTH, SOUTH],optFilename=optFilename,RunID=RunID,isRestart=isRestart)
        else
            call init_mix(vApp%remixApp%ion,[NORTH, SOUTH],RunID=RunID,isRestart=isRestart)
        endif
        vApp%remixApp%ion%rad_iono_m = RadIonosphere() ! Returns in planetary radii
        vApp%remixApp%ion%rad_iono_m = vApp%remixApp%ion%rad_iono_m * gApp%Model%units%gx0 
        !Set F10.7 from time series (using max)
        f107%wID = vApp%tilt%wID
        call f107%initTS("f10.7")
        maxF107 = f107%getMax()
        
        do n=1,2
            vApp%remixApp%ion(n)%P%f107 = maxF107
        enddo
        write(*,*) 'Using F10.7 = ', maxF107
        write(*,*) 'Using MJD0  = ', gApp%Model%MJD0

        call init_mhd2Mix(vApp%mhd2mix, gApp, vApp%remixApp)
        call init_mix2Mhd(vApp%mix2mhd, vApp%remixApp, gApp)
        vApp%mix2mhd%mixOutput = 0.0
        
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

        vApp%ShallowT = time + vApp%ShallowDT

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

        ! solve for remix output
        if (time<=0) then
            call run_mix(vApp%remixApp%ion,curTilt,doModelOpt=.false.)
        else
            call run_mix(vApp%remixApp%ion,curTilt,doModelOpt=.true.)
        endif
        ! get stuff from mix to gamera
        call mapRemixToGamera(vApp%mix2mhd, vApp%remixApp)

    end subroutine runRemix

!----------
!Deep coupling stuff (time coming from vApp%time, so in seconds)
    subroutine DeepUpdate(vApp, gApp, time)
        type(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: time

        real(rp) :: tAdv

        if (.not. vApp%doDeep) then
            !Why are you even here?
            return
        endif

        tAdv = time + vApp%DeepDT !Advance inner magnetosphere through full coupling time 
    
    !Pull in updated fields to CHIMP
        call Tic("G2C")
        call convertGameraToChimp(vApp%mhd2chmp,gApp,vApp%ebTrcApp)
        call Toc("G2C")

    !Advance inner magnetosphere model to tAdv
        call Tic("InnerMag")
        call AdvanceInnerMag(vApp,tAdv)
        call Toc("InnerMag")

    !Squish 3D data to 2D IMAG grid (either RP or lat-lon)
        !Doing field projection at current time
        call Tic("Squish")
        call Squish(vApp)
        call Toc("Squish")

    !Now use imag model and squished coordinates to fill Gamera source terms
        call Tic("IM2G")
        call InnerMag2Gamera(vApp,gApp)
        call Toc("IM2G")
    

    !Setup next coupling
        vApp%DeepT = time + vApp%DeepDT

    end subroutine DeepUpdate


    !Initialize CHIMP data structure
    subroutine init_volt2Chmp(ebTrcApp,gApp,optFilename)
        type(ebTrcApp_T), intent(inout) :: ebTrcApp
        type(gamApp_T), intent(in) :: gApp
        character(len=*), intent(in), optional     :: optFilename

        character(len=strLen) :: xmlStr
        type(XML_Input_T) :: inpXML
        
    !Create input XML object
        if (present(optFilename)) then
            xmlStr = trim(optFilename)
        else
            call getIDeckStr(xmlStr)
        endif
        inpXML = New_XML_Input(trim(xmlStr),"Chimp",.true.)

    !Initialize model
        associate(Model=>ebTrcApp%ebModel,ebState=>ebTrcApp%ebState,Gr=>gApp%Grid)
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
        call InitLoc(Model,ebState%ebGr,inpXML)
        end associate

    end subroutine init_volt2Chmp

end module voltapp

