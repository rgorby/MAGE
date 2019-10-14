!Driver for Gamera with MPI decomposition

program gamera_mpix
    use clocks
    use gamapp_mpi
    use usergamic
    use mpidefs

    implicit none

    type(gamAppMpi_T) :: gameraAppMpi
    procedure(StateIC_T), pointer :: userInitFunc => initUser

    ! initialize mpi data type
    call setMpiReal()

    !call printConfigStamp()
    call initClocks()

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = 1
    
    call initGamera_mpi(gameraAppMpi,userInitFunc)

    do while (gameraAppMpi%Model%t < gameraAppMpi%Model%tFin)
        call Tic("Omega")
        !Start root timer

        call stepGamera_mpi(gameraAppMpi)

        !Output if necessary
        call Tic("IO")
        
        if (gameraAppMpi%Model%IO%doConsole(gameraAppMpi%Model%ts)) then
            call consoleOutput(gameraAppMpi%Model,gameraAppMpi%Grid,gameraAppMpi%State)
        endif

        if (gameraAppMpi%Model%IO%doOutput(gameraAppMpi%Model%t)) then
            call fOutput(gameraAppMpi%Model,gameraAppMpi%Grid,gameraAppMpi%State)
        endif

        if (gameraAppMpi%Model%IO%doRestart(gameraAppMpi%Model%t)) then
            call resOutput(gameraAppMpi%Model,gameraAppMpi%Grid,gameraAppMpi%State)
        endif

        call Toc("IO")

    !Do timing info
        if (gameraAppMpi%Model%IO%doTimer(gameraAppMpi%Model%ts)) then
            if (gameraAppMpi%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif

        call Toc("Omega")
    end do

end program gamera_mpix

