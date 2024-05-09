!Driver for Gamera coupled with Voltron (remix only)

program voltronx
    use clocks
    use gamapp
    use voltapp
    use voltio
    use uservoltic

    implicit none

    type(voltApp_T) :: vApp

    call initClocks()

    vApp%gApp%Model%isLoud = .true.
    vApp%vOptions%gamUserInitFunc => initUser

    call initVoltron(vApp)

    do while (vApp%time < vApp%tFin)
        !Start root timer
        call Tic("Omega")
        
        !Advance voltron models one coupling step
        call Tic("StepVoltron")
        call stepVoltron(vApp)
        call Toc("StepVoltron")
        
        !IO checks
        call Tic("IO")
        !Console output
        if (vApp%IO%doConsole(vApp%ts)) then
            !Using console output from Gamera
            call consoleOutputV(vApp,vApp%gApp)
            !Timing info
            if (vApp%IO%doTimerOut) call printClocks()
            call cleanClocks()
        elseif (vApp%IO%doTimer(vApp%ts)) then
            if (vApp%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif
        
        !Data output
        if (vApp%IO%doOutput(vApp%time)) then
            call fOutputV(vApp,vApp%gApp)
        endif
        !Restart output
        if (vApp%IO%doRestart(vApp%time)) then
            call resOutputV(vApp,vApp%gApp)
        endif

        call Toc("IO")

        call Toc("Omega")

    end do
    write(*,*) "Fin"
    
end program voltronx

