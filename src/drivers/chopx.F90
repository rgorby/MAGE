!Driver for stand-alone EB field slicing

program chopx
    use clocks
    use chmpdefs
    use chopio
    use ebtypes
    use iotable
    use starter
    use chmpfields

    implicit none

    !Main data structures
    type(chmpModel_T) :: Model
    type(ebState_T)   :: ebState
    type(XML_Input_T) :: inpXML

    character(len=strLen) :: gStr
    character(len=strLen) :: utStr

    real(rp) :: mjd
    
    !------------
    !Setup timers
    call initClocks()

    !----------------------------
    !Initialize model and fields
    call goApe(Model,ebState,iXML=inpXML)
    Model%doEBOut = .true.
    call initEB3Dio(Model,ebState,inpXML)

    !Loop from T0 -> tFin
    Model%t = Model%T0

    !Do main loop
    do while (Model%t<=Model%tFin)
        call Tic("Omega")

        call Tic("Step")
        !Update fields to current time
        call updateFields(Model,ebState,Model%t)
        call Toc("Step")


        call Tic("Output")        
        !Write slices
        write(gStr,'(A,I0)') "Step#", Model%nOut
        call writeEB3D(Model,ebState,gStr)
        !Setup for next output
        Model%tOut = Model%tOut + Model%dtOut
        Model%nOut = Model%nOut + 1

        if (modulo(Model%ts,Model%tsOut) ==0) then
            write(*,'(a,f12.3,a)') 'T = ', Model%t*oTScl, ' ' // trim(tStr)
            if (ebState%ebTab%hasMJD) then
                mjd = iotabMJD(ebState%ebTab,Model%t)
                call mjd2utstr(mjd,utStr)
                write (*,'(a,a)')                    '      UT   = ', trim(utStr)
            endif
        endif
        
        call Toc("Output")

        call Toc("Omega")

        !Timing book keeping
        if (modulo(Model%ts,Model%tsOut) ==0) then
            if (Model%doTimer) call printClocks()
            call cleanClocks()
        endif

        !Update time
        Model%t = Model%t + Model%dt
        Model%ts = Model%ts+1

    enddo
    !------------------
    !Do finalization housekeeping
    write(*,*) "If you lived here you'd be home by now"

end program chopx