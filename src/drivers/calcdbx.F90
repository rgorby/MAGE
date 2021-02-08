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

#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    !Main CHIMP data structures
    type(chmpModel_T) :: Model
    type(ebState_T)   :: ebState
    type(XML_Input_T) :: inpXML

    !New calc-DB structures
    type(rmState_T) :: rmState
    type(sphGrid_T) :: gGr !Ground grid
    type(facGrid_T) :: facGrid !FAC grid
    type(ionGrid_T) :: ionGrid !Ionospheric grid

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
    
    call initDBio(Model,ebState,gGr,inpXML,NumP)
    call initRM(Model,ebState,rmState)

    !Loop from T0 -> tFin
    Model%t = Model%T0
    
    !Do main loop
    do while (Model%t<=Model%tFin)
        call Tic("Omega")

        call Tic("Step")
        !Update fields to current time
        call updateFields(Model,ebState,Model%t)
        
        !Update remix data to current time
        call updateRemix(Model,ebState,Model%t,rmState)

        !Update FAC/ION grids
        call facGridUpdate(Model,ebState,rmState,facGrid)
        
        call Toc("Step")

        call Tic("Compute")
    !Mag DB
        call Tic("MagDB")
        call CalcMagDB(Model,ebState,gGr)
        call Toc("MagDB")
    !Ion DB

    !FAC DB
        call Toc("Compute")

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