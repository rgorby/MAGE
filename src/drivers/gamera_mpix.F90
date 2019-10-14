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

        !Do timing info
        if (modulo(gameraAppMpi%Model%ts,gameraAppMpi%Model%tsOut) == 0) then
            if (gameraAppMpi%Model%doTimer) call printClocks()
            call cleanClocks()
        endif
        call Toc("Omega")
    end do

end program gamera_mpix

