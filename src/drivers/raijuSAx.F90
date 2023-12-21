program raijuSAx
    !! Stand-alone RAIJU driver

    use kdefs
    use xml_input
    use clocks
    use geopack, only : mjdRECALC

    ! RAIJU stuff
    use raijudefs
    use raijustarter
    use raijuBCs
    use raijuPreAdvancer
    use raijuAdvancer


    implicit none

    enum, bind(C)
        enumerator :: RaiTest_WM=1
    end enum
    integer :: testCode

    type(raijuApp_T   ) :: raiApp

    character(len=strLen) :: XMLStr, gStr
    type(XML_Input_T) :: inpXML    
    
    character(len=strLen) :: FLH5
    
    logical :: doChmpOut,doFLOut

    real(rp) :: mjd0



    call initClocks()
    call Tic("Omega")
                    
    ! Init xml
    call getIDeckStr(XMLStr)

    ! Personal xml flags
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/RAIJU",.true.)

    !! Break-out testing of specific pieces
    call inpXML%Set_Val(testCode,'testing/code',-1)
    select case(testCode)
        case (-1)  ! Nothing, just pass
            write(*,*) "No test code, doing standard SA driving"
        case(RaiTest_WM)
            write(*,*) "Doing test of Wave Model"
            call test_WM(inpXML)
            stop
        case default
            write(*,*) "Idk what this test code is, try again"
            stop
    end select

    call inpXML%Set_Val(raiApp%Model%doClockConsoleOut,'driver/doClockOut',.false.)

    ! Init RAIJU
    call raijuInit(raiApp, inpXML)

    !> MJD at the start of the simulation (corresponds to sim t=0)
    call inpXML%Set_Val(mjd0,'prob/MJD0',51544.0)  ! default to 2000-01-01T00:00:00
    raiApp%State%mjd = mjd0

    ! If we need geopack, make sure we init
    if (raiApp%Model%doGeoCorot) then
        call mjdRECALC(raiApp%State%mjd)
    endif

    ! Ready to loop
    do while ( raiApp%State%t < (raiApp%Model%tFin + 0.5) )
        
        call Tic("Output")
    ! Output if ready
        if (raiApp%State%IO%doOutput(raiApp%State%t)) then
            call raijuOutput(raiApp%Model,raiApp%Grid,raiApp%State)
        endif
        call Toc("Output")

    ! Update potentials and magnetic field somehow if you want
    ! This used to be where one way driving voltron -> RAIJU happened

    ! Step RAIJU
        call Tic("RAIJU Advance")
        call raijuAdvance(raiApp%Model,raiApp%Grid,raiApp%State, raiApp%Model%dt)
        call Toc("RAIJU Advance")

        !write(*,*)raiApp%State%t
        !write(*,*)"MJDs: "
        !write(*,*)raiApp%State%mjd

        ! Advance model times
        raiApp%State%t  = raiApp%State%t  + raiApp%Model%dt
        raiApp%State%ts = raiApp%State%ts + 1
        raiApp%State%mjd = T2MJD(raiApp%State%t,mjd0)

        ! Update geopack if we need to
        if (raiApp%Model%doGeoCorot) then
            call mjdRECALC(raiApp%State%mjd)
        endif

        if (raiApp%Model%doClockConsoleOut) then
            call printClocks()
            call cleanClocks()
        endif
    enddo


    call Toc("Omega")
    write(*,*)"Main Objective complete!"

    contains

    subroutine raijuAdvance(Model, Grid, State, dtCpl, isFirstCplO)
        type(raijuModel_T), intent(inout) :: Model
        type(raijuGrid_T) , intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
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

    end subroutine raijuAdvance


    subroutine test_WM(iXML)
        type(XML_Input_T), intent(in) :: iXML 
        
        character(len=strLen) :: confName
        type(eLossWM_T) :: eWM
        
        call iXML%Set_Val(confName, "config/fname","raijuconfig.h5")
        call initEWM(eWM, confName, iXML, raiApp%Grid%shGrid)

        write(*,*) "Kp1D (",size(eWM%Kp1D),"):"
        write(*,*) eWM%Kp1D
        write(*,*) "MLT1D (",size(eWM%MLT1D),"):"
        write(*,*) eWM%MLT1D
        write(*,*) "L1D (",size(eWM%L1D),"):"
        write(*,*) eWM%L1D
        write(*,*) "Energy1D (",size(eWM%Energy1D),"):"
        write(*,*) eWM%Energy1D
        write(*,*) "dims(Tau4D):"
        write(*,*) size(eWM%Tau4D)
    end subroutine

end program raijuSAx