!Driver for Gamera coupled with Voltron (remix only)

program voltronx
    use clocks
    use gamapp
    use voltapp
    use uservoltic

    implicit none

    type(gamApp_T) :: gameraApp
    type(voltApp_T) :: voltronApp
    procedure(StateIC_T), pointer :: userInitFunc => initUser

    call initClocks()

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = 1
    
    call initGamera(gameraApp,userInitFunc)
    call initVoltron(voltronApp, gameraApp)

    do while (gameraApp%Model%t < gameraApp%Model%tFin)
        !Start root timer
        call Tic("Omega")
        
        !Advance Gamera MHD
        call stepGamera(gameraApp)

        !Do any updates to Voltron
        call stepVoltron(voltronApp,gameraApp)
        
        call Tic("DeepCoupling")
        if ( (gameraApp%Model%t >= voltronApp%DeepT) .and. voltronApp%doDeep ) then
            call DeepUpdate(voltronApp, gameraApp, gameraApp%Model%t)
        endif
        call Toc("DeepCoupling")

        call Tic("IonCoupling")
        if (gameraApp%Model%t >= voltronApp%ShallowT) then
            call ShallowUpdate(voltronApp, gameraApp, gameraApp%Model%t)
        endif
        call Toc("IonCoupling")
        
        !Do timing info
        if (modulo(gameraApp%Model%ts,gameraApp%Model%tsOut) == 0) then
            if (gameraApp%Model%doTimer) call printClocks()
            call cleanClocks()
        endif
        call Toc("Omega")
    end do

end program voltronx

