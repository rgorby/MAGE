!Driver for Gamera coupled with Voltron (remix only)

program voltronx
    use clocks
    use gamapp
    use voltapp
    use voltio
    use uservoltic

    implicit none

    type(gamApp_T) :: gApp
    type(voltApp_T) :: vApp
    procedure(StateIC_T), pointer :: userInitFunc => initUser

    call initClocks()

    gApp%Model%isLoud = .true.

    call initGamera(gApp,userInitFunc,doIO=.false.)
    call initVoltron(vApp, gApp)

    do while (vApp%time < vApp%tFin)
        !Start root timer
        call Tic("Omega")
        
    !Advance Gamera MHD
        call stepGamera(gApp)

    !Do any updates to Voltron
        call stepVoltron(vApp,gApp)
        
    !Coupling
        call Tic("DeepCoupling")
        if ( vApp%time >= vApp%DeepT ) then
            call DeepUpdate(vApp, gApp)
        endif
        call Toc("DeepCoupling")

    !Update coupling DTs
    vApp%DeepDT = vApp%TargetDeepDT
       
    !IO checks
        call Tic("IO")
        !Console output
        if (vApp%IO%doConsole(vApp%ts)) then
            !Using console output from Gamera
            call consoleOutputV(vApp,gApp)
            !Timing info
            if (vApp%IO%doTimerOut) call printClocks()
            call cleanClocks()
        elseif (vApp%IO%doTimer(vApp%ts)) then
            if (vApp%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif
        
        !Data output
        if (vApp%IO%doOutput(vApp%time)) then
            call fOutputV(vApp,gApp)
        endif
        !Restart output
        if (vApp%IO%doRestart(vApp%time)) then
            call resOutputV(vApp,gApp)
        endif

        call Toc("IO")

        call Toc("Omega")

    end do
    write(*,*) "Fin"
    
end program voltronx

