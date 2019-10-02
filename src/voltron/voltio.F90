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
        call qkMIXOut(vApp%remixApp%ion,vApp%mix2mhd%mixOutput,vApp%time,vApp%IO%nOut)

        vApp%IO%tOut = vApp%IO%tOut + vApp%IO%dtOut
        vApp%IO%nOut = vApp%IO%nOut + 1

    end subroutine fOutputV

    !NOTE: This is a temp remix output which will get replaced by Slava's revamped version
    subroutine qkMIXOut(ion,mhdvarsin,time,step)
        type(mixIon_T),dimension(:),intent(inout) :: ion ! I for ionosphere (is an array of 1 or 2 elements for north and south)
        real(rp), dimension(:,:,:,:,:),intent(in) :: mhdvarsin
        real(rp), intent(in) :: time
        integer, intent(in) :: step

        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars
        character(len=strLen) :: fnstr,fname,vID
        logical :: isThere
        real(rp) :: cpcp(2) = 0.0

        !Save CPCP for diagnostics
        call getCPCP(mhdvarsin,cpcp)

        !Create output name
        write(fnstr,'(I0.6)') step
        fname = 'mixtest'//trim(fnstr)//'.h5'
        call CheckAndKill(fname)

        !Do standard mix output
        call writeMIX(fname,ion)

        !Add extra attribute information to output
        vID = "t"
        isThere = ioExist(fname,trim(vID))
        if (.not. isThere) then
          call ClearIO(IOVars)
          call AddOutVar(IOVars,"t"   ,time)
          call AddOutVar(IOVars,"ts"  ,step)
          call AddOutVar(IOVars,"nCPCP",cpcp(NORTH))
          call AddOutVar(IOVars,"sCPCP",cpcp(SOUTH))
          call WriteVars(IOVars,.true.,fname)
        endif !isThere

    end subroutine qkMIXOut

    subroutine getCPCP(mhdvarsin,cpcp)
        real(rp), dimension(:,:,:,:,:),intent(in) :: mhdvarsin
        real(rp), intent(out) :: cpcp(2)
        cpcp(NORTH) = maxval(mhdvarsin(1,:,:,MHDPSI,NORTH))-minval(mhdvarsin(1,:,:,MHDPSI,NORTH))
        cpcp(SOUTH) = maxval(mhdvarsin(1,:,:,MHDPSI,SOUTH))-minval(mhdvarsin(1,:,:,MHDPSI,SOUTH))

    end subroutine getCPCP
end module voltio

