!Driver for CHIMP code
!Trace test-particles through precomputed Gamera EM fields

program pushx
    use clocks
    use chmpdefs
    use tptypes
    use ebtypes
    use starter
    use chmpio
    use pusher
    use chmpfields

    implicit none

    !Main data structures
    type(chmpModel_T) :: Model
    type(tpState_T)   :: tpState
    type(ebState_T)   :: ebState

    !Variables for main loop
    logical :: doPush = .true.
    integer :: n

    !Setup timers
    call initClocks()

    !Initialize main data structures
    call goApe(Model,ebState,tpState)
    Model%doTPOut = .true.
    
    !Prepare for main loop, do first output
    call cleanClocks()

    call Tic("Output")
    call cOutput(Model,ebState,tpState)
    if (Model%t >= Model%tOut) call fOutput(Model,ebState,tpState)
    call Toc("Output")
    
    !Use association for syntax brevity
    associate( TPs=>tpState%TPs )

    do while (doPush)
        !Start root timer
        call Tic("Omega")
        
        !------------------
        !Advance system (could be array of tpStates)
        call Tic("Chimp")
        !$OMP PARALLEL DO default(shared) &
        !$OMP& schedule(dynamic)
        do n=1,tpState%Np
            !Integrate individual particles, advance by dt
            if (TPs(n)%isIn) then
                call PushTP(TPs(n),Model%t,Model%dt,Model,ebState)
            endif
        enddo
        call Toc("Chimp")

        !------------------
        !Update heartbeat/do per-step diagnostics
        call Tic("Step")
        Model%t = Model%t+Model%dt
        Model%ts = Model%ts+1
        tpState%NpT = count(tpState%TPs(:)%isIn)
        
        call updateFields(Model,ebState,Model%t)
        if (Model%doStream) call addIncoming(Model,ebState,tpState)

        call Toc("Step")

        !------------------
        !Do IO if necessary
        call Tic("Output")
        if (modulo(Model%ts,Model%tsOut) ==0) then
            call cOutput(Model,ebState,tpState)
        endif
        if (Model%t >= Model%tOut) then
            call fOutput(Model,ebState,tpState)
        endif
        call Toc("Output")

        !------------------
        !Check for exit conditions
        if (tpState%NpT == 0 .and. (.not. Model%doStream) ) then
            doPush = .false.
            write(*,*) 'No particles active, finishing.'
        endif
        if (Model%t>=Model%tFin) then
            doPush = .false.
            write(*,*) 'End time reached, finishing.'
        endif

        call Toc("Omega")

        !------------------
        !Do timing stuff
        if (modulo(Model%ts,Model%tsOut) == 0) then
            if (Model%doTimer) call printClocks()
            call cleanClocks()
        endif

    enddo !While push loop

    end associate
    !------------------
    !Do finalization housekeeping
    write(*,*) "If you lived here you'd be home by now"

end program pushx