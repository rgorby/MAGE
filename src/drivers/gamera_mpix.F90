!Driver for Gamera with MPI decomposition

program gamera_mpix
    use clocks
    use gamapp
    use gamapp_mpi
    use usergamic
    use mpidefs

    implicit none

    type(gamApp_T) :: gameraApp
    procedure(StateIC_T), pointer :: userInitFunc => initUser

    ! initialize mpi data type
    call setMpiReal()

    !call printConfigStamp()
    call initClocks()

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = 1
    
    call initGamera_mpi(gameraApp,userInitFunc)

    do while (gameraApp%Model%t < gameraApp%Model%tFin)
        call Tic("Omega")
        !Start root timer

        call stepGamera_mpi(gameraApp)

        !Do timing info
        if (modulo(gameraApp%Model%ts,gameraApp%Model%tsOut) == 0) then
            if (gameraApp%Model%doTimer) call printClocks()
            call cleanClocks()
        endif
        call Toc("Omega")
    end do

end program gamera_mpix

