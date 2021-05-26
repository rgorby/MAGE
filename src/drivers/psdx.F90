!Driver for PSD calculations

program psdx
    use clocks
    use chmpdefs
    use ebtypes
    use starter
    use chmpfields
    use psdtypes
    use psdinit
    use psdutils
    use psdio
    use psdcalc
    use xml_input

    implicit none

    !Main data structures
    type(chmpModel_T) :: Model
    type(ebState_T)   :: ebState
    type(XML_Input_T) :: inpXML
    type(PSEq_T)      :: psGr
    type(psdPop_T)    :: psPop

    !------------
    !Setup timers
    call initClocks()

    !----------------------------
    !Initialize
    !Model & fields
    call goApe(Model,ebState,iXML=inpXML)
    !Phase space
    call genPS(Model,psGr,inpXML)
    !PSD population
    call genPop(Model,psPop,inpXML)

    !----------------------------
    !Prep loop and go
    Model%t = Model%T0

    !Do main loop
    do while (Model%t<=Model%tFin)
    !------------------
        !Main work
        call Tic("Omega")

    !Update states
        !Do updates for current time (eb/mhd, psd dV, tp states)
        call Tic("Step") 
        
        call Tic("Fields")
        call updateFields(Model      ,ebState,Model%t)
        call Toc("Fields")
        call Tic("dG")
        call updatePS    (Model,psGr ,ebState,Model%t)
        call Toc("dG")
        call Tic("TPs")
        call updatePop   (Model,psPop,ebState,Model%t)
        call Toc("TPs")
        
        call Toc("Step")

        !Calculate weights for unweighted particles
        call Tic("Weights")
        call CalcWeights(Model,psGr,ebState,psPop)
        call Toc("Weights")

        !Synthesize TPs & weights to PSD
        call Tic("PSD")
        call CalcPSD(Model,psGr,ebState,psPop)
        call Toc("PSD")

        call Tic("Output")
        if (Model%t >= Model%tOut) then
            call fOutPSD(Model,ebState,psGr,psPop)
        endif
        
        call Toc("Output")

        call Toc("Omega")
    !------------------
        !Book keeping

        !Timing info
        if (modulo(Model%ts,Model%tsOut) ==0) then
            if (Model%doTimer) call printClocks()
            call cleanClocks()
        endif

        !Update times
        Model%t = Model%t + Model%dt
        Model%ts = Model%ts+1        
    enddo    

    write(*,*) 'Saving weights ...'
    call fOutWgts(Model,psGr,psPop)
    
    write(*,*) 'Finished PSD!'
end program psdx