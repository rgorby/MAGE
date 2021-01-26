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

    !Main data structures
    type(chmpModel_T) :: Model
    type(ebState_T)   :: ebState
    type(XML_Input_T) :: inpXML


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
    
    call initDBio(Model,ebState,inpXML,NumP)
    
    !Loop from T0 -> tFin
    Model%t = Model%T0

    ! write(*,*) ''
    ! write(*,'(a,I0,a)')  'Tracing ', NumP, ' points'
    ! write(*,'(a,f12.3,a,f12.3)') 'Time Interval = ', Model%T0*oTScl,' / ', Model%tFin*oTScl
    ! write(*,'(a,f8.3)') 'dt = ', Model%dt*oTScl
    
    !Do main loop
    do while (Model%t<=Model%tFin)
        call Tic("Omega")

        call Tic("Step")
        !Update fields to current time
        call updateFields(Model,ebState,Model%t)
        call Toc("Step")

        !Calc/write DB on grid
        call Tic("Output")
        !Write output grid
        write(gStr,'(A,I0)') "Step#", Model%nOut
        call Tic("MagDB")

        call writeDB(Model,ebState,gStr)
        call Toc("MagDB")
        
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