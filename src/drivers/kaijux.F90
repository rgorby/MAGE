!Driver for Gamera coupled with Voltron (remix only)

program kaijux
    use clocks
    use gamapp
    use voltapp

    implicit none

    type(gamApp_T) :: gameraApp
    type(voltApp_T) :: voltronApp

    call initClocks()

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = 1
    
    call initGamera(gameraApp)
    call initVoltron(voltronApp, gameraApp)

    do while (gameraApp%Model%t < gameraApp%Model%tFin)
        call Tic("Omega")
        !Start root timer

        call stepGamera(gameraApp)

        if (gameraApp%Model%t >= voltronApp%fastShallowTime) then
            call fastShallowUpdate(voltronApp, gameraApp, gameraApp%Model%t)
        endif

        !Do timing info
        if (modulo(gameraApp%Model%ts,gameraApp%Model%tsOut) == 0) then
            if (gameraApp%Model%doTimer) call printClocks()
            call cleanClocks()
        endif
        call Toc("Omega")
    end do

end program kaijux

