!Driver for Gamera (uncoupled)

program glx
    use clocks
    use giblowapp
    use gloutput

    implicit none

    type(glApp_T) :: glApp

    !call printConfigStamp()
    call initClocks()

    glApp%Model%isLoud = .true.
    
    call initGL(glApp)
    
    if(glApp%Model%isLoud) write(*,*) "Model tFin = ", glApp%Model%tFin

    do while (glApp%Model%time <= glApp%Model%tFin)
        call Tic("Omega") !Start root timer
    
        !Step model/s    
        call step(glApp)

        !Output if necessary
        call Tic("IO")
        
        !TODO: Output
        ! if (glApp%Model%IO%doConsole(glApp%Model%ts)) then
        !     call consoleOutput(glApp%Model,glApp%Grid,glApp%State)
        ! endif

        if (glApp%Model%IO%doOutput(glApp%Model%time)) then
             call fOutput(glApp%Model,glApp%State,glApp%Solution)
        endif

        call Toc("IO")
        glApp%Model%ts = glApp%Model%ts + 1
        glApp%Model%time = glApp%Model%time + glApp%Model%dt
        !Do timing info
        if (glApp%Model%IO%doTimer(glApp%Model%ts)) then
            if (glApp%Model%IO%doTimerOut) call printClocks()
            call cleanClocks()
        endif
        
        ! if (glApp%Model%isDebug) then    
        !     write(*,*) "Radii: ", glApp%State%r
        !     write(*,*) "Inside Solution Values: "
        !     write(*,*) "    dens: ",    pack(glApp%Solution%dens, glApp%Solution%inside_mask .gt. 0)  
        !     write(*,*) "    temp: ",    pack(glApp%Solution%temp, glApp%Solution%inside_mask .gt. 0) 
        !     write(*,*) "    pres: ",    pack(glApp%Solution%pres, glApp%Solution%inside_mask .gt. 0) 
        !     write(*,*) "    br: ",      pack(glApp%Solution%b(:,:,:,XDIR), glApp%Solution%inside_mask .gt. 0)
        !     write(*,*) "    btheta: ",  pack(glApp%Solution%b(:,:,:,YDIR), glApp%Solution%inside_mask .gt. 0)
        !     write(*,*) "    bphi: ",    pack(glApp%Solution%b(:,:,:,ZDIR), glApp%Solution%inside_mask .gt. 0)
        !     write(*,*) "    jr: ",      pack(glApp%Solution%j(:,:,:,XDIR), glApp%Solution%inside_mask .gt. 0)
        !     write(*,*) "    jtheta: ",  pack(glApp%Solution%j(:,:,:,YDIR), glApp%Solution%inside_mask .gt. 0)
        !     write(*,*) "    jphi: ",    pack(glApp%Solution%j(:,:,:,ZDIR), glApp%Solution%inside_mask .gt. 0)
        !     write(*,*) "    vr: ",      pack(glApp%Solution%v(:,:,:,XDIR), glApp%Solution%inside_mask .gt. 0) 
        !     write(*,*) "    vtheta: ",  pack(glApp%Solution%v(:,:,:,YDIR), glApp%Solution%inside_mask .gt. 0)
        !     write(*,*) "    vphi: ",    pack(glApp%Solution%v(:,:,:,ZDIR), glApp%Solution%inside_mask .gt. 0)
        ! end if
        call Toc("Omega")
    end do
    write(*,*) "Fin"
end program glx
