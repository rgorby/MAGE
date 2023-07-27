program sifOWDx

    use kdefs
    use planethelper
    use xml_input
    use clocks

    ! Sif stuff
    use sifdefs
    use sifstarter
    use sifBCs
    use sifadvancer

    ! Chimp stuff
    use ebtypes
    use starter
    use chmpfields
    use chopio
    use chmpunits

    ! Remix, actually
    use calcdbtypes
    use calcdbio

    ! Voltron stuff
    use volttypes
    use sifuseric
    use sifCplTypes
    use sifCpl
    use sifowdcpl

    implicit none

    type(voltApp_T)   :: vApp
        !! Kinda hacky, just create a voltApp object and initialize only the parts we're gonna use ourselves
    !Holder for remix data
    type(rmState_T) :: rmState

    type(sifApp_T   ) :: sApp
    type(sif_cplBase_T) :: sifCplBase

    character(len=strLen) :: XMLStr, gStr
    type(XML_Input_T) :: inpXML    
    
    character(len=strLen) :: FLH5
    
    logical :: doChmpOut,doFLOut
    logical :: isFirstCpl = .true.

    real(rp) :: mjd0

    call initClocks()
    call Tic("Omega")

    associate(ebModel=>vApp%ebTrcApp%ebModel, ebState=>vApp%ebTrcApp%ebState)
            
        ! Init xml
        call getIDeckStr(XMLStr)

        ! Personal xml flags
        inpXML = New_XML_Input(trim(XMLStr),"Kaiju/SIF",.true.)
        call inpXML%Set_Val(doChmpOut,'driver/doChmpOut',.false.)
        call inpXML%Set_Val(doFLOut,'driver/doFLOut',.false.)
        call inpXML%Set_Val(sApp%Model%doClockConsoleOut,'driver/doClockOut',.false.)

        ! Init SIF
        call sifInit(sApp, inpXML)
        call sifCpl_init(vApp, sApp, sifCplBase)

        ! Init CHIMP
        inpXML = New_XML_Input(trim(XMLStr),"Kaiju/Chimp",.true.)
        call goApe(ebModel,ebState,iXML=inpXML)
        ebModel%t = inTscl*sApp%Model%t0

        ! Use ebTabs to get MJD
        !> MJD at the start of the simulation (corresponds to sim t=0)
        mjd0 = T2MJD(-1.0*ebState%ebTab%times(1)/inTscl, ebState%ebTab%MJDs(1))
        sApp%State%mjd = mjd0

        ! Init Remix reader
        call initRM  (ebModel,ebState,rmState,inpXML)
        rmState%time = inTscl*sApp%Model%t0

        ! Set mix->SIF map
        call InitMixMap(sApp%Grid%shGrid, rmState)
        
        ! Init outputs
        ebModel%doEBOut = doChmpOut
        if (ebModel%doEBOut) then
            call initEB3Dio(ebModel,ebState,inpXML)
        endif

        FLH5   = trim(sApp%Model%RunID) // ".fl.h5" !RCM field lines
        call CheckAndKill(FLH5)

        ! Ready to loop
        do while ( sApp%State%t < (sApp%Model%tFin + 0.5) )
            
            call Tic("Output")
        ! Output if ready
            if (sApp%State%IO%doOutput(sApp%State%t)) then
                call sifOutput(sApp%Model,sApp%Grid,sApp%State)

                if (ebModel%doEBOut) then
                    ! Write eb at the same output cadence as imag
                    write(gStr,'(A,I0)') "Step#", ebModel%nOut
                    call writeEB3D(ebModel,ebState,gStr)
                    ebModel%nOut = ebModel%nOut + 1
                    ebModel%tOut = inTscl*sApp%State%IO%tOut  ! Idk if we need to set this since chimp isn't controlling its own output
                endif

                call WriteRCMFLs(sifCplBase%fromV%fLines,sApp%State%IO%nOut, &
                        sApp%State%mjd,sApp%State%t, &
                        sApp%Grid%shGrid%Nt,sApp%Grid%shGrid%Np)
            endif
            call Toc("Output")

        ! Update other models
            call Tic("CHIMP update")
            call updateFields(ebModel, ebState, ebModel%t)
            call Toc("CHIMP update")

            call Tic("REMIX update")
            call updateRemix(ebModel,ebState,ebModel%t,rmState)
            call Toc("REMIX update")

            ! Populate sif's fromV object with updated model info
            call Tic("fromV packing")
            call packFromV(sifCplBase%fromV, vApp, rmState, sApp)
            call Toc("fromV packing")

            call Tic("fromV to State")
            ! Now put fomV info into sif's State
            call sifCpl_Volt2SIF(sifCplBase, vApp, sApp)
            call Toc("fromV to State")

        ! Step SIF
            call Tic("SIF Advance")
            call sifAdvance(sApp%Model,sApp%Grid,sApp%State, sApp%Model%dt, isFirstCpl)
            call Toc("SIF Advance")
            !isFirstCpl = .false.

            write(*,*)sApp%State%t
            write(*,*)ebModel%t/inTscl
            write(*,*)rmState%time/inTscl
            write(*,*)"MJDs: "
            write(*,*)sApp%State%mjd
            write(*,*)T2MJD((ebModel%t - ebState%ebTab%times(1))/inTscl, ebState%ebTab%MJDs(1))

            ! Advance model times
            sApp%State%t  = sApp%State%t  + sApp%Model%dt
            sApp%State%ts = sApp%State%ts + 1
            sApp%State%mjd = T2MJD(sApp%State%t,mjd0)

            ebModel%t     = ebModel%t  + inTscl*sApp%Model%dt
            ebModel%ts    = ebModel%ts + 1

            rmState%time  = rmState%time + inTScl*sApp%Model%dt

            if (sApp%Model%doClockConsoleOut) then
                call printClocks()
                call cleanClocks()
            endif
        enddo

    end associate

    call Toc("Omega")
    write(*,*)"Done :)"

    contains

    subroutine sifAdvance(Model, Grid, State, dtCpl, isFirstCplO)
        type(sifModel_T), intent(in) :: Model
        type(sifGrid_T) , intent(in) :: Grid
        type(sifState_T), intent(inout) :: State
        real(rp), intent(in) :: dtCpl
        logical, optional, intent(in) :: isFirstCplO

        logical :: isFirstCpl
        integer :: k

        if (present(isFirstCplO)) then
            isFirstCpl = isFirstCplO
        else
            isFirstCpl = .false.
        endif

        State%dt = dtCpl

        call raijuPreAdvance(Model, Grid, State, isfirstCpl)

        ! Step
        call Tic("AdvanceState")
        call AdvanceState(Model, Grid, State)
        call Toc("AdvanceState")
            ! Push
            ! Losses
        ! etas to moments
        call Tic("Moments Eval")
        call EvalMoments(Grid, State)
        call Toc("Moments Eval")

        

    end subroutine sifAdvance

    subroutine WriteRCMFLs(RCMFLs,nOut,MJD,time,Ni,Nj)
        use ebtypes
        use rice_housekeeping_module, only : nSkipFL
        integer, intent(in) :: nOut,Ni,Nj
        real(rp), intent(in) :: MJD,time
        type(fLine_T), intent(in), dimension(Ni,Nj) :: RCMFLs

        type(IOVAR_T), dimension(40) :: IOVars
        character(len=strLen) :: gStr,lnStr
        integer :: i,j,n
        
        !Bail out if we're not doing this
        if (.not. doFLOut) return

        !Create group and write base data
        write(gStr,'(A,I0)') "Step#", nOut
        call AddOutVar(IOVars,"time",time)
        call AddOutVar(IOVars,"MJD",MJD)

        
        call WriteVars(IOVars,.true.,FLH5,gStr)
        call ClearIO(IOVars)

        !Now loop through and create subgroup for each line (w/ striding)
        !TODO: Avoid the individual write for every line
        n = 0
        do i=1,Ni,nSkipFL
            do j=1,Nj-1,nSkipFL
                write(lnStr,'(A,I0)') "Line#", n
                if (RCMFLs(i,j)%isGood) then
                    call OutLine(RCMFLs(i,j),gStr,lnStr,IOVars)
                    n = n + 1
                endif
            enddo
        enddo
    end subroutine WriteRCMFLs

    ! Do our own fline output cause original is using globals
    !Write out individual line
    subroutine OutLine(fL,gStr,lnStr,IOVars)
        USE ebtypes
        use rice_housekeeping_module, ONLY : nSkipFL
        type(fLine_T), intent(in) :: fL
        character(len=strLen), intent(in) :: gStr,lnStr
        type(IOVAR_T), intent(inout), dimension(40) :: IOVars
        integer :: i,Np,Npp,n0
        
        call ClearIO(IOVars)
        Np = fL%Nm + fL%Np + 1
        if (Np<=nSkipFL) return
        n0 = fL%Nm

        !Add scalar stuff
        !Record seed point
        call AddOutVar(IOVars,"x0",fL%x0(XDIR))
        call AddOutVar(IOVars,"y0",fL%x0(YDIR))
        call AddOutVar(IOVars,"z0",fL%x0(ZDIR))

        !Do striding through field line points
        Npp = size(fL%xyz(0:-n0:-nSkipFL,XDIR))

        call AddOutVar(IOVars,"xyz",transpose(fL%xyz(0:-n0:-nSkipFL,XDIR:ZDIR)))
        call AddOutVar(IOVars,"Np",Npp)
        call AddOutVar(IOVars,"n0",1) !Seed point is now the first point

        !Only output some of the variables
        call AddOutVar(IOVars,"B",fL%lnVars(0)       %V(0:-n0:-nSkipFL),uStr="nT")
        call AddOutVar(IOVars,"D",fL%lnVars(DEN)     %V(0:-n0:-nSkipFL),uStr="#/cc")
        call AddOutVar(IOVars,"P",fL%lnVars(PRESSURE)%V(0:-n0:-nSkipFL),uStr="nPa")

        !Write output chain
        call WriteVars(IOVars,.true.,FLH5,gStr,lnStr)
        call ClearIO(IOVars)

    end subroutine OutLine

end program sifOWDx