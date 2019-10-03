!Routines to write restart/outputs from Voltron
module voltio
    use gamapp
    use volttypes
    use mixio
    
    implicit none

    integer, parameter, private :: MAXVOLTIOVAR = 10

    contains

    subroutine consoleOutputV(vApp,gApp)
        type(gamApp_T) , intent(in) :: gApp
        type(voltApp_T), intent(in) :: vApp
        real(rp) :: cpcp(2) = 0.0

        !Using console output from Gamera
        call consoleOutput(gApp%Model,gApp%Grid,gApp%State)

        !Augment Gamera console output w/ Voltron stuff
        call getCPCP(vApp%mix2mhd%mixOutput,cpcp)


        write(*,'(a)',advance="no") ANSIBLUE
        !write (*, '(a,f8.3,a)')       '    dt/dt0 = ', 100*Model%dt/dt0, '%'
        write (*, '(a,2f8.3,a)')      '    CPCP  = ' , cpcp(NORTH), cpcp(SOUTH), ' [kV] (N/S)'
        write (*, *) ANSIRESET, ''

    end subroutine consoleOutputV

    subroutine resOutputV(vApp,gApp)
        type(gamApp_T) , intent(inout) :: gApp
        type(voltApp_T), intent(inout) :: vApp

        !Write Gamera restart
        !NOTE: This shouldn't be done in MPI!
        if (.not. gApp%Model%isMPI) then
            call resOutput(gApp%Model,gApp%Grid,gApp%State)
        else
            write(*,*) 'Need to handle Voltron restart case for MPI!'
        endif

        vApp%IO%tRes = vApp%IO%tRes + vApp%IO%dtRes
        vApp%IO%nRes = vApp%IO%nRes + 1            

    end subroutine resOutputV

    subroutine fOutputV(vApp,gApp)
        type(gamApp_T) , intent(inout) :: gApp
        type(voltApp_T), intent(inout) :: vApp

        !Write gamera data
        !NOTE: This shouldn't be done in MPI!
        if (.not. gApp%Model%isMPI) then
            call fOutput(gApp%Model,gApp%Grid,gApp%State) !Gamera
        else
            write(*,*) 'Need to handle Voltron output case for MPI!'
        endif

        !Write ReMIX data
        call writeMix(vApp%remixApp%ion,vApp%IO%nOut,mjd=vApp%MJD,time=vApp%time)

        vApp%IO%tOut = vApp%IO%tOut + vApp%IO%dtOut
        vApp%IO%nOut = vApp%IO%nOut + 1

    end subroutine fOutputV

    subroutine getCPCP(mhdvarsin,cpcp)
        real(rp), dimension(:,:,:,:,:),intent(in) :: mhdvarsin
        real(rp), intent(out) :: cpcp(2)
        cpcp(NORTH) = maxval(mhdvarsin(1,:,:,MHDPSI,NORTH))-minval(mhdvarsin(1,:,:,MHDPSI,NORTH))
        cpcp(SOUTH) = maxval(mhdvarsin(1,:,:,MHDPSI,SOUTH))-minval(mhdvarsin(1,:,:,MHDPSI,SOUTH))

    end subroutine getCPCP
end module voltio

