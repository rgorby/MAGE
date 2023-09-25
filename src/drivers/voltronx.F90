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
        
        !Advance Gamera MHD
        !Do any updates to Voltron
        call stepVoltron(vApp)
        
    !Coupling
        call Tic("DeepCoupling")
        if ( (vApp%time >= vApp%DeepT) .and. vApp%doDeep ) then
            call DeepUpdate(vApp, vApp%gApp)
        endif
        call Toc("DeepCoupling")

    !Update coupling DTs
    vApp%DeepDT = vApp%TargetDeepDT
       
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
        
        !Restart output
        if (vApp%IO%doRestart(vApp%time)) then
            call resOutputV(vApp,vApp%gApp)
        endif
        !Data output
        if (vApp%IO%doOutput(vApp%time)) then
            call fOutputV(vApp,vApp%gApp)
        endif

        call Toc("IO")

        call Toc("Omega")

    end do
    write(*,*) "Fin"
    
end program voltronx

