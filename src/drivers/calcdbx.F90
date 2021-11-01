!Driver for calculating magnetospheric db on ground grid
program calcdbx
    use clocks
    use chmpdefs
    use chmpio
    use ebtypes
    use starter
    use chmpfields
    use math
    use calcdbio
    use calcdbremap

#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    !Main CHIMP data structures
    type(chmpModel_T) :: Model
    type(ebState_T)   :: ebState
    type(XML_Input_T) :: inpXML

    !Holder for remix data
    type(rmState_T) :: rmState

    !Native (SM) source grids
    type(facGrid_T) :: facGrid !FAC grid
    type(ionGrid_T) :: ionGrid !Ionospheric grid

    !Bios-Savart data, source terms in GEO
    type(BSGrid_T) :: facBS,ionBS,magBS

    !Destination (GEO) data
    type(grGrid_T) :: gGr !Ground grid
    
    integer :: NumP
    real(rp) :: wT,cMJD,rSec
    character(len=strLen) :: gStr,utStr
    integer :: iYr,iDoY,iMon,iDay,iHr,iMin

    !write(*,*) "Num threads=",omp_get_max_threads()
    !------------
    !Setup timers
    call initClocks()

    !----------------------------
    !Initialize model and fields
    call goApe(Model,ebState,iXML=inpXML)
    if (.not. Model%doJ) then
        write(*,*) "Must use fields/doJ=T, bailing ..."
        stop
    endif

    !----------------------------
    !Initialize data structures
    call initDBio(Model,ebState,gGr,inpXML,NumP)
    call initRM(Model,ebState,rmState,inpXML)

    call ionGridInit(Model,ebState,rmState,ionGrid)
    call facGridInit(Model,ebState,rmState,ionGrid,facGrid)

    call BSGridInit(Model,ebState,rmState,magBS,ionBS,facBS)

    write(*,*) "Initialization complete, let's do this thing!"
    !Loop from T0 -> tFin
    Model%t = Model%T0
    
    !Do main loop
    do while (Model%t<=Model%tFin)
        cMJD = MJDAt(ebState%ebTab,Model%t)
        call MJDRecalc(cMJD) !Setup geopack for this time

        call Tic("Omega")

    !Read in data and fill native (SM) grids
        call Tic("Step")
        !Update fields to current time
        call updateFields(Model,ebState,Model%t)
        
        !Update remix data to current time
        call updateRemix(Model,ebState,Model%t,rmState)

        !Update FAC/ION grids
        call facGridUpdate(Model,ebState,rmState,facGrid)
        call ionGridUpdate(Model,ebState,rmState,ionGrid)

        call Toc("Step")

    !Pack data into 1D bricks
        call Tic("Pack")
        call packBS(Model,Model%t,ebState,ionGrid,facGrid,magBS,ionBS,facBS)
        call Toc("Pack")

    !Remap ground grid to SM coordinates
        call Tic("Remap")
        call remapGR(Model,Model%t,ebState,gGr)
        call Toc("Remap")

    !Compute BS integrals (SM) on ground and remap dB to GEO
        call Tic("ComputeBS")
        call BS2Gr(Model,Model%t,ebState,magBS,ionBS,facBS,gGr)
        call Toc("ComputeBS")

        !Calc/write DB on grid
        call Tic("Output")
        !Write output grid
        write(gStr,'(A,I0)') "Step#", Model%nOut
        
        call writeDB(Model,ebState,gGr,gStr)
        
        !Setup for next output
        Model%tOut = Model%tOut + Model%dtOut
        Model%nOut = Model%nOut + 1

        if (modulo(Model%ts,Model%tsOut) ==0) then
            cMJD = MJDAt(ebState%ebTab,Model%t)
            call mjd2ut(cMJD,iYr,iDoY,iMon,iDay,iHr,iMin,rSec)
            write(utStr,'(I0.4,a,I0.2,a,I0.2,a,I0.2,a,I0.2,a,I0.2)') iYr,'-',iMon,'-',iDay,' ',iHr,':',iMin,':',nint(rSec)

            wT = readClock("MagDB")
            write(*,'(a,a)')       'UT = ', trim(utStr)
            write(*,'(a,f12.3,a)') '     SMU = ', gGr%SMU, ' [nT]'
            write(*,'(a,2f12.3,a)') '         @ Lat/Lon: ', gGr%SMU_MLat,gGr%SMU_MLon, ' [deg]'
            write(*,'(a,f12.3,a)') '     SML = ', gGr%SML, ' [nT]'
            write(*,'(a,2f12.3,a)') '         @ Lat/Lon: ', gGr%SML_MLat,gGr%SML_MLon, ' [deg]'
            write(*,'(a,f12.3,a)') '     SME = ', gGr%SME, ' [nT]'
            write(*,'(a,f12.3,a)') '     SMR = ', gGr%SMR, ' [nT]'
            write(*,'(a,4f12.3,a)') '         SMR-00/06/12/18 = ', gGr%SMR_00,gGr%SMR_06,gGr%SMR_12,gGr%SMR_18, ' [nT]'
            write(*,'(a,f12.3,a)') '       T = ', Model%t*oTScl, ' ' // trim(tStr)
            write(*,'(a,f8.3)')    '   kDBps = ', 1.0e-3*NumP*Model%tsOut/wT
        endif

        call Toc("Output")
        call Toc("Omega")

        !Timing book keeping
        if (modulo(Model%ts,Model%tsOut) ==0) then
            if (Model%doTimer) call printClocks()
            call cleanClocks()
        endif

        !Update time
        Model%t = Model%t + Model%dt
        Model%ts = Model%ts+1

    enddo


end program calcdbx