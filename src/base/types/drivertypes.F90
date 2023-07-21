! Module containing types that contain and run apps

module drivertypes

    use basetypes
    use clocks

    implicit none

    type :: AppDriver_T

        ! this array must be allocated and filled before calling RunApps
        type(AppHolder), dimension(:), allocatable :: appPointers

        type (IOClock_T) :: IO
        real(rp) :: t, tFin
        integer :: ts
        character(len=strLen) :: runID

        contains

        procedure :: RunApps
    end type AppDriver_T


    contains

    subroutine RunApps(driver, optionalXml)
        class(AppDriver_T), intent(inout) :: driver
        character(len=*), optional, intent(in) :: optionalXml

        character(len=strLen) :: inpXML
        type(XML_Input_T) :: xmlInp
        integer :: numApps, app

        if(.not. allocated(driver%appPointers)) then
            write (*,*) "You must allocate at least one app within the driver type"
            return
        endif

        numApps = SIZE(driver%appPointers)

        !call printConfigStamp()
        call initClocks()

        if(present(optionalXml)) then
            ! read from the prescribed file
            inpXML = optionalXml
        else
            !Find input deck
            call getIDeckStr(inpXML)
        endif
        call CheckFileOrDie(inpXML,"Error opening input deck, exiting ...")

        !Create XML reader
        write(*,*) 'Reading input deck from ', trim(inpXML)
        xmlInp = New_XML_Input(trim(inpXML),'Kaiju',.true.)
        ! Each model is expected to change the root if they want to

        do app = 1,numApps
            call driver%appPointers(app)%p%InitModel(xmlInp)
            call driver%appPointers(app)%p%InitIO(xmlInp)
        enddo

        do while (driver%t < driver%tFin)
            call Tic("Omega") !Start root timer

            !Step model/s
            do app = 1,numApps
                call driver%appPointers(app)%p%AdvanceModel(0.0_rp)
            enddo

            !Output if necessary
            call Tic("IO")

            if (driver%IO%doConsole(driver%ts)) then
                do app = 1,numApps
                    call driver%appPointers(app)%p%WriteConsoleOutput()
                enddo

                !Timing info
                if (driver%IO%doTimerOut) call printClocks()
                call cleanClocks()

            elseif (driver%IO%doTimer(driver%ts)) then
                if (driver%IO%doTimerOut) call printClocks()
                call cleanClocks()
            endif

            if (driver%IO%doOutput(driver%t)) then
                do app = 1,numApps
                    call driver%appPointers(app)%p%WriteFileOutput(driver%runID, driver%ts)
                enddo
            endif

            if (driver%IO%doRestart(driver%t)) then
                do app = 1,numApps
                   call driver%appPointers(app)%p%WriteRestart(driver%runID)
                enddo
            endif

            call Toc("IO")
            call Toc("Omega")
        end do

        write(*,*) "Fin"

    end subroutine RunApps

end module

