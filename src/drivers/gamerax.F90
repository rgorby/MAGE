!Driver for Gamera (uncoupled)

program gamerax
    use types
    use clocks
    use output
    use init
    use gamutils
    use mhdgroup
    use step

    implicit none

    type(Model_T) :: Model
    type(Grid_T) :: Grid
    type(State_T) :: State, oState

    call initClocks()

    !TODO: Fix this to reset after MPI config to only output from root rank
    verbose = 1
    
    !Initialize Grid/State/Model (Hatch Gamera)
    !Will enforce 1st BCs, caculate 1st timestep, set oldState
    call Hatch(Model,Grid,State,oState)
    call cleanClocks()

    if (.not. Model%isRestart) call fOutput(Model,Grid,State)
    call consoleOutput(Model,Grid,State)


    do while (Model%t < Model%tFin)
        call Tic("Omega")
        !Start root timer

        call Tic("Gamera")
        !Advance system
        call AdvanceMHD(Model,Grid,State,oState,Model%dt)
        call Toc("Gamera")

        !Enforce floors if necessary
        if (Model%doArmor) then
            call Tic("Armor")
            call Armor(Model,Grid,State)
            call Toc("Armor")
        endif

        !Update info
        Model%ts = Model%ts+1
        Model%t = Model%t+Model%dt

        !Call user-defined per-step function
        if (associated(Model%HackStep)) then
            call Tic("HackStep")
            call Model%HackStep(Model,Grid,State)
            call Toc("HackStep")
        endif

        !Calculate new timestep
        call Tic("DT")
        Model%dt = CalcDT(Model,Grid,State)
        call Toc("DT")

        !Enforce BCs
        call Tic("Halos")
        call EnforceBCs(Model,Grid,State)
        call Toc("Halos")

        
        !Output if necessary
        call Tic("IO")
        if (modulo(Model%ts,Model%tsOut) ==0) then
            call consoleOutput(Model,Grid,State)
        endif
        if (Model%t >= Model%tOut) then
            call fOutput(Model,Grid,State)
        endif
        if (Model%doResOut .and. (Model%t >= Model%tRes)) then
            !print *,"RESTART :: ", Model%doResOut, Model%t >= Model%tRes
            call resOutput(Model,Grid,State)
        endif
        call Toc("IO")

        !Do timing info
        if (modulo(Model%ts,Model%tsOut) == 0) then
            if (Model%doTimer) call printClocks()
            call cleanClocks()
        endif
        call Toc("Omega")
    end do

end program gamerax
