!Driver for Gamera (uncoupled)

program gamerax
    use clocks
    use gamapp
    use usergamic

    implicit none

    type(GameraApp_T) :: gApp
    character(len=strLen) :: inpXML
    type(XML_Input_T) :: xmlInp


    !call printConfigStamp()
    call initClocks()

    gApp%Model%isLoud = .true.
    
    gApp%Options%userInitFunc => initUser

    !Find input deck
    call getIDeckStr(inpXML)
    !Create XML reader
    write(*,*) 'Reading input deck from ', trim(inpXML)
    xmlInp = New_XML_Input(trim(inpXML),'Kaiju/Gamera',.true.)

    call gApp%InitModel(xmlInp)
    call gApp%InitIO(xmlInp)

    do while (gApp%Model%t < gApp%Model%tFin)
        call Tic("Omega") !Start root timer
    
    !Step model/s    
        call gApp%AdvanceModel(0.0_rp)

    !Output if necessary
        call Tic("IO")
        
        if (gApp%Model%IO%doConsole(gApp%Model%ts)) then
            call gApp%WriteConsoleOutput()
            !Timing info
            if (gApp%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()
        elseif (gApp%Model%IO%doTimer(gApp%Model%ts)) then
            if (gApp%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif

        if (gApp%Model%IO%doOutput(gApp%Model%t)) then
            call gApp%WriteFileOutput('', 0)
        endif

        if (gApp%Model%IO%doRestart(gApp%Model%t)) then
            call gApp%WriteRestart('')
        endif

        call Toc("IO")

        call Toc("Omega")
    end do
    write(*,*) "Fin"
end program gamerax
