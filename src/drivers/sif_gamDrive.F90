program sif_gamdrive

    use kdefs
    use planethelper
    use xml_input

    ! Sif stuff
    use sifdefs
    use sifapp

    ! Chimp stuff
    use ebtypes
    use starter
    use chmpfields
    use chopio
    use chmpunits

    ! Voltron stuff
    use volttypes
    use sifCplTypes
    use sifCpl

    implicit none

    type(voltApp_T)   :: vApp
        !! Kinda hacky, just create a voltApp object and initialize only the parts we're gonna use ourselves
    type(sifApp_T   ) :: sApp
    type(sif_cplBase_T) :: sifCplBase

    character(len=strLen) :: XMLStr, gStr
    type(XML_Input_T) :: inpXML    
    
    
    logical :: doChmpOut

    associate(ebModel=>vApp%ebTrcApp%ebModel, ebState=>vApp%ebTrcApp%ebState)
            
        ! Init xml
        call getIDeckStr(XMLStr)

        ! Personal xml flags
        inpXML = New_XML_Input(trim(XMLStr),"Kaiju/SIF",.true.)
        call inpXML%Set_Val(doChmpOut,'driver/doChmpOut',.false.)

        ! Init SIF
        call initSifVolt(vApp, sApp, inpXML)
        call sifCpl_init(vApp, sApp, sifCplBase)

        ! Init CHIMP
        inpXML = New_XML_Input(trim(XMLStr),"Kaiju/Chimp",.true.)
        call goApe(ebModel,ebState,iXML=inpXML)
        
        ebModel%doEBOut = doChmpOut
        if (ebModel%doEBOut) then
            call initEB3Dio(ebModel,ebState,inpXML)
        endif



        ! Ready to loop
        do while ( sApp%State%t < (sApp%Model%tFin + 0.5) )
        ! Output if ready
            if (sApp%State%IO%doOutput(sApp%State%t)) then
                call sifOutput(sApp%Model,sApp%Grid,sApp%State)

                ! Write eb at the same output cadence as imag
                write(gStr,'(A,I0)') "Step#", ebModel%nOut
                call writeEB3D(ebModel,ebState,gStr)
                ebModel%nOut = ebModel%nOut + 1
                ebModel%tOut = inTscl*sApp%State%IO%tOut  ! Idk if we need to set this since chimp isn't controlling its own output
            endif

        ! Update models
            call updateFields(ebModel, ebState, ebModel%t)
            call sifCpl_Volt2SIF(sifCplBase, vApp, sApp)


            ! Advance model times
            sApp%State%t  = sApp%State%t  + sApp%Model%dt
            sApp%State%ts = sApp%State%ts + 1
            ebModel%t  = ebModel%t  + inTscl*sApp%Model%dt
            ebModel%ts = ebModel%ts + 1


        enddo

    end associate

end program sif_gamdrive