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

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = 1
    
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
        if ( (vApp%time >= vApp%DeepT) .and. vApp%doDeep ) then
            call DeepUpdate(vApp, gApp, vApp%time)
        endif
        call Toc("DeepCoupling")

        call Tic("IonCoupling")
        if (vApp%time >= vApp%ShallowT) then
            call ShallowUpdate(vApp, gApp, vApp%time)
        endif
        call Toc("IonCoupling")
        
    !IO checks
        call Tic("IO")
        !Console output
        if (vApp%IO%doConsole(vApp%ts)) then
            !Using console output from Gamera
            call consoleOutputV(vApp,gApp)
        endif
        !Restart output
        if (vApp%IO%doRestart(vApp%time)) then
            call resOutputV(vApp,gApp)
        endif
        !Data output
        if (vApp%IO%doOutput(vApp%time)) then
            call fOutputV(vApp,gApp)
        endif

        call Toc("IO")

    !Do timing info
        if (vApp%IO%doTimer(vApp%ts)) then
            if (vApp%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif

        call Toc("Omega")

    end do

end program voltronx

