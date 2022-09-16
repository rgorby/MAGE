!Driver for mpi decomposed Gamera coupled with Voltron (remix only)

program voltron_mpix
    use clocks
    use gamapp_mpi
    use voltapp_mpi
    use gam2VoltComm_mpi
    use output
    use voltio
    use uservoltic
    use mpi_f08
    use xml_input

    implicit none

    ! this driver requires that each rank be either Gamera OR Voltron, never both
    logical :: isGamera = .true.

    ! gamera data
    type(gamAppMpi_T) :: gApp
    type(gam2VoltCommMpi_T) :: g2vComm

    ! voltron data
    type(voltAppMpi_T) :: vApp

    procedure(StateIC_T), pointer :: userInitFunc => initUser

    integer :: ierror, length, provided, worldSize, worldRank, numHelpers
    type(MPI_Comm) :: gamComm, voltComm
    integer :: required=MPI_THREAD_MULTIPLE
    character( len = MPI_MAX_ERROR_STRING) :: message
    character(len=strLen) :: inpXML, helpersBuf
    logical :: useHelpers, helperQuit
    integer(KIND=MPI_AN_MYADDR) :: tagMax
    logical :: tagSet

    ! initialize MPI
    !Set up MPI with or without thread support
#ifdef _OPENMP
    call MPI_INIT_THREAD(required,provided,ierror)
    if(provided < required) then
        print *,"Not support for MPI_THREAD_MULTIPLE, aborting!"
        call abort()
    end if
    !print *,"MPI + OpenMP !!!"
#else
    !print *," MPI without threading"
    call MPI_INIT(ierror)
#endif

    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    ! initialize mpi data type
    call setMpiReal()

    call initClocks()

    call mpi_comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, tagMax, tagSet, ierror)
    print *, 'Tag Upper-Bound = ', tagMax

    gApp%Model%isLoud = .true.    

    ! need to know how many voltron helpers there are
    call getIDeckStr(inpXML)
    call ReadXmlImmediate(trim(inpXML),'/Kaiju/Voltron/Helpers/useHelpers',helpersBuf,'F',.false.)
    read(helpersBuf,*) useHelpers
    call ReadXmlImmediate(trim(inpXML),'/Kaiju/Voltron/Helpers/numHelpers',helpersBuf,'0',.false.)
    read(helpersBuf,*) numHelpers
    if(.not. useHelpers) numHelpers = 0

    ! create a new MPI communicator for just Gamera
    !    for now this is always all ranks excep the last one (which is reserved for voltron)
    call MPI_Comm_Size(MPI_COMM_WORLD, worldSize, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    call MPI_Comm_Rank(MPI_COMM_WORLD, worldRank, ierror)
    if(ierror /= MPI_Success) then
        call MPI_Error_string( ierror, message, length, ierror)
        print *,message(1:length)
        call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    if(worldRank .ge. (worldSize-1-numHelpers)) then
        ! voltron rank
        isGamera = .false.
        call MPI_Comm_Split(MPI_COMM_WORLD, 1, worldRank, voltComm, ierror)
        if(ierror /= MPI_Success) then
            call MPI_Error_string( ierror, message, length, ierror)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
        end if
    else
        ! gamera rank
        isGamera = .true.
        call MPI_Comm_Split(MPI_COMM_WORLD, 0, worldRank, gamComm, ierror)
        if(ierror /= MPI_Success) then
            call MPI_Error_string( ierror, message, length, ierror)
            print *,message(1:length)
            call mpi_Abort(MPI_COMM_WORLD, 1, ierror)
        end if
    endif

    if(isGamera) then
        call Tic("Omega")
        call initGamera_mpi(gApp,userInitFunc,gamComm,doIO=.false.)
        call initGam2Volt(g2vComm,gApp,MPI_COMM_WORLD)
        call Toc("Omega")

        do while (g2vComm%time < g2vComm%tFin)
            !Start root timer
            call Tic("Omega")
        
            !Advance Gamera MHD
            call stepGamera_mpi(gApp)

            !Update local gamera time units
            call localStepVoltronTime(g2vComm, gApp)

            !Tell gamera to synchronize with voltron (and wait for it if necessary)
            ! If it needs to do shallow or deep coupling, file io, restart io, or
            ! Terminate the simulation
            ! Don't synchronize for console or timing output, not worth it (?)
            if( (g2vComm%time >= g2vComm%DeepT .and. g2vComm%doDeep) .or. &
                (g2vComm%time >= g2vComm%ShallowT) .or. &
                gApp%Model%IO%doRestart(gApp%Model%t) .or. &
                gApp%Model%IO%doOutput(gApp%Model%t) .or. &
                (g2vComm%time >= g2vComm%tFin)) then
                !Do any updates to Voltron
                call Tic("StepVoltron")
                call performStepVoltron(g2vComm,gApp)
                call Toc("StepVoltron")
        
                !Coupling
                if(g2vComm%doDeep .and. g2vComm%time >= g2vComm%DeepT .and. g2vComm%time >= g2vComm%ShallowT) then ! both shallow and deep coupling
                    call Tic("Coupling")
                    call performShallowAndDeepUpdate(g2vComm, gApp)
                    call Toc("Coupling")
                elseif ( g2vComm%time >= g2vComm%DeepT .and. g2vComm%doDeep ) then
                    call Tic("Coupling")
                    call performDeepUpdate(g2vComm, gApp)
                    call Toc("Coupling")
                elseif (g2vComm%time >= g2vComm%ShallowT) then
                    call Tic("Coupling")
                    call performShallowUpdate(g2vComm, gApp)
                    call Toc("Coupling")
                endif
            endif

            !IO checks
            call Tic("IO")
            !NOTE: Does this need to be duplicated from gamera_mpix or can both be done w/ a single subroutine?
            
            !Console output
            if (gApp%Model%IO%doConsole(g2vComm%ts)) then
                !Using console output from Gamera
                call consoleOutput_mpi(gApp)
                if (gApp%Model%IO%doTimerOut .and. debugPrintingRank(gApp)) then
                    write (*,*) "Rank (I,J,K) (", &
                        gApp%Grid%Ri,",", &
                        gApp%Grid%Rj,",", &
                        gApp%Grid%Rk,") is printing debug clock info"
                    call printClocks()
                endif
                call cleanClocks()
            elseif (gApp%Model%IO%doTimer(g2vComm%ts)) then
                if (gApp%Model%IO%doTimerOut .and. debugPrintingRank(gApp)) then
                    write (*,*) "Rank (I,J,K) (", &
                        gApp%Grid%Ri,",", &
                        gApp%Grid%Rj,",", &
                        gApp%Grid%Rk,") is printing debug clock info"
                    call printClocks()
                endif
                call cleanClocks()
            endif

            !Restart output
            if (gApp%Model%IO%doRestart(gApp%Model%t)) then
                call resOutput(gApp%Model, gApp%Grid, gApp%oState, gApp%State)
            endif
            !Data output
            if (gApp%Model%IO%doOutput(gApp%Model%t)) then
                call fOutput(gApp%Model, gApp%Grid, gApp%State)
            endif

            call Toc("IO")

            call Toc("Omega")

        end do

        call endGam2VoltWaits(g2vComm, gApp)

    else 
        ! voltron
        call Tic("Omega")
        call initVoltron_mpi(vApp, userInitFunc, voltComm, MPI_COMM_WORLD)
        call Toc("Omega")

        if(vApp%amHelper) then
            ! do helper loop
            helperQuit = .false.
            do while(.not. helperQuit)
                call Tic("Omega")
                call helpVoltron(vApp, helperQuit)
                call Toc("Omega")
            end do
        else
            ! voltron run loop
            do while (vApp%time < vApp%tFin)
                !Start root timer
                call Tic("Omega")
                !If coupling from Gamera is ready
                if(gameraStepReady(vApp)) then
                    !Do any updates to Voltron
                    call Tic("StepVoltronAndWait")
                    call stepVoltron_mpi(vApp)
                    call Toc("StepVoltronAndWait")

                    !Coupling
                    call Tic("Coupling")
                    if(vApp%doDeep .and. vApp%time >= vApp%DeepT .and. vApp%time >= vApp%ShallowT) then ! both
                        call shallowAndDeepUpdate_Mpi(vApp, vApp%time)
                    elseif (vApp%time >= vApp%DeepT .and. vApp%doDeep ) then
                        call DeepUpdate_mpi(vApp, vApp%time)
                    elseif (vApp%time >= vApp%ShallowT) then
                        call ShallowUpdate_mpi(vApp, vApp%time)
                    endif
                    call Toc("Coupling")

                    !IO checks
                    call Tic("IO")
                    !Console output
                    if (vApp%IO%doConsole(vApp%ts)) then
                        call consoleOutputVOnly(vApp,vApp%gAppLocal,vApp%gAppLocal%Model%MJD0)
                        if (vApp%IO%doTimerOut) then
                            call printClocks()
                        endif
                        call cleanClocks()
                    elseif (vApp%IO%doTimer(vApp%ts)) then
                        if (vApp%IO%doTimerOut) then
                            call printClocks()
                        endif
                        call cleanClocks()
                    endif
                    
                    !Restart output
                    if (vApp%IO%doRestart(vApp%time)) then
                        call resOutputVOnly(vApp,vApp%gAppLocal)
                    endif
                    !Data output
                    if (vApp%IO%doOutput(vApp%time)) then
                        call fOutputVOnly(vApp,vApp%gAppLocal)
                    endif

                    call Toc("IO")

                elseif(deepInProgress(vApp)) then
                    ! If we did not couple, check for deep work to be done
                    call Tic("Coupling")
                    call doDeepBlock(vApp)
                    call Toc("Coupling")
                else
                    ! Gamera wasn't ready and we don't have deep work to do, wait for gamera
                    call waitForGameraStep(vApp)
                endif
                call Toc("Omega")
            enddo
            ! if using helpers, tell them to quit
            if(vApp%useHelpers) call vhReqHelperQuit(vApp)
        endif ! end all voltron loops
        call endVoltronWaits(vApp)
    endif

    call MPI_FINALIZE(ierror)
    write(*,*) "Fin"
    

end program voltron_mpix

