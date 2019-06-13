! Main data objects and functions to perform a gamera simulation

module gamapp
    use types
    use step
    use init
    use mhdgroup

    implicit none

    type gamApp_T
        type(Model_T) :: Model
        type(Grid_T) :: Grid
        type(State_T) :: State, oState
    end type gamApp_T

    contains

    subroutine initGamera(gameraApp)
        type(gamApp_T), intent(inout) :: gameraApp

        !Initialize Grid/State/Model (Hatch Gamera)
        !Will enforce 1st BCs, caculate 1st timestep, set oldState
        call Hatch(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState)
        call cleanClocks()

        if (.not. gameraApp%Model%isRestart) call fOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call consoleOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)

    end subroutine initGamera

    subroutine stepGamera(gameraApp)
        type(gamApp_T), intent(inout) :: gameraApp

        call Tic("Gamera")
        !Advance system
        call AdvanceMHD(gameraApp%Model,gameraApp%Grid,gameraApp%State,gameraApp%oState,gameraApp%Model%dt)
        call Toc("Gamera")

        !Enforce floors if necessary
        if (gameraApp%Model%doArmor) then
            call Tic("Armor")
            call Armor(gameraApp%Model,gameraApp%Grid,gameraApp%State)
            call Toc("Armor")
        endif

        !Update info
        gameraApp%Model%ts = gameraApp%Model%ts+1
        gameraApp%Model%t = gameraApp%Model%t+gameraApp%Model%dt

        !Call user-defined per-step function
        if (associated(gameraApp%Model%HackStep)) then
            call Tic("HackStep")
            call gameraApp%Model%HackStep(gameraApp%Model,gameraApp%Grid,gameraApp%State)
            call Toc("HackStep")
        endif

        !Calculate new timestep
        call Tic("DT")
        gameraApp%Model%dt = CalcDT(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call Toc("DT")

        !Enforce BCs
        call Tic("Halos")
        call EnforceBCs(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        call Toc("Halos")


        !Output if necessary
        call Tic("IO")
        if (modulo(gameraApp%Model%ts,gameraApp%Model%tsOut) ==0) then
            call consoleOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        if (gameraApp%Model%t >= gameraApp%Model%tOut) then
            call fOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        if (gameraApp%Model%doResOut .and. (gameraApp%Model%t >= gameraApp%Model%tRes)) then
            !print *,"RESTART :: ", Model%doResOut, Model%t >= Model%tRes
            call resOutput(gameraApp%Model,gameraApp%Grid,gameraApp%State)
        endif
        call Toc("IO")

    end subroutine stepGamera

end module gamapp

