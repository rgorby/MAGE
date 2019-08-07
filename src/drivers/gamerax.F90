!Driver for Gamera (uncoupled)

program gamerax
    use clocks
    use gamapp

    implicit none

    type(gamApp_T) :: gameraApp

    !call printConfigStamp()
    call initClocks()

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = 1
    
    call initGamera(gameraApp)

    do while (gameraApp%Model%t < gameraApp%Model%tFin)
        call Tic("Omega")
        !Start root timer

        call stepGamera(gameraApp)

        !Do timing info
        if (modulo(gameraApp%Model%ts,gameraApp%Model%tsOut) == 0) then
            if (gameraApp%Model%doTimer) call printClocks()
            call cleanClocks()
        endif
        call Toc("Omega")
    end do

end program gamerax
