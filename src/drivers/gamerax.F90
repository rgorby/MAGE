!Driver for Gamera (uncoupled)

program gamerax
    use clocks
    use gamapp
    use usergamic

    implicit none

    type(gamApp_T) :: gApp
    procedure(StateIC_T), pointer :: userInitFunc => initUser

    !call printConfigStamp()
    call initClocks()

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = .true.
    
    call initGamera(gApp,userInitFunc)

    do while (gApp%Model%t < gApp%Model%tFin)
        call Tic("Omega") !Start root timer
    
    !Step model/s    
        call stepGamera(gApp)

    !Output if necessary
        call Tic("IO")
        
        if (gApp%Model%IO%doConsole(gApp%Model%ts)) then
            call consoleOutput(gApp%Model,gApp%Grid,gApp%State)
        endif

        if (gApp%Model%IO%doOutput(gApp%Model%t)) then
            call fOutput(gApp%Model,gApp%Grid,gApp%State)
        endif

        if (gApp%Model%IO%doRestart(gApp%Model%t)) then
            call resOutput(gApp%Model,gApp%Grid,gApp%State)
        endif

        call Toc("IO")

    !Do timing info
        if (gApp%Model%IO%doTimer(gApp%Model%ts)) then
            if (gApp%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif

        call Toc("Omega")
    end do

end program gamerax
