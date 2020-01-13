!Routines to write restart/outputs from Voltron
module voltio
    use gamapp
    use volttypes
    use mixio
    use clocks
    use innermagsphere
    use wind

    implicit none

    integer, parameter, private :: MAXVOLTIOVAR = 10
    logical, private :: isConInit = .false.
    real(rp), private ::  oMJD = 0.0

    contains

    subroutine consoleOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(in) :: vApp

        !Using console output from Gamera
        call consoleOutput(gApp%Model,gApp%Grid,gApp%State)

        !Using console output from Voltron
        call consoleOutputVOnly(vApp,gApp%Model%MJD0)

    end subroutine consoleOutputV

    subroutine consoleOutputVOnly(vApp,MJD0)
        class(voltApp_T), intent(in) :: vApp
        real(rp), intent(in) :: MJD0

        real(rp) :: cpcp(2) = 0.0

        real(rp) :: dpT,dtWall,cMJD,dMJD,simRate

        integer :: iYr,iDoY,iMon,iDay,iHr,iMin
        real(rp) :: rSec
        character(len=strLen) :: utStr

        !Augment Gamera console output w/ Voltron stuff
        call getCPCP(vApp%mix2mhd%mixOutput,cpcp)

        dpT = vApp%tilt%evalAt(vApp%time)*180.0/PI

        !Figure out some perfromance info
        cMJD = T2MJD(vApp%time,MJD0) !Current MJD


        if (isConInit) then
            !Console output has been initialized
            dMJD = cMJD - oMJD !Elapsed MJD since last console output
            dtWall = kClocks(1)%tElap

            simRate = dMJD*24.0*60.0*60.0/dtWall !Model seconds per wall second
            oMJD = cMJD
        else
            simRate = 0.0
            oMJD = cMJD
            isConInit = .true.
        endif

        !Get MJD info
        call mjd2ut(cMJD,iYr,iDoY,iMon,iDay,iHr,iMin,rSec)
        write(utStr,'(I0.4,a,I0.2,a,I0.2,a,I0.2,a,I0.2,a,I0.2)') iYr,'-',iMon,'-',iDay,' ',iHr,':',iMin,':',nint(rSec)

        write(*,'(a)',advance="no") ANSIBLUE
        write (*, '(a,2f8.3,a)')             '      CPCP = ' , cpcp(NORTH), cpcp(SOUTH), ' [kV, N/S]'
        write (*, '(a,1f8.3,a)')             '      tilt = ' , dpT, ' [deg]'
        write (*,'(a,a)')                    '      UT   = ', trim(utStr)
        write (*, '(a,1f7.3,a)')             '      Running @ ', simRate*100.0, '% of real-time'

        write (*, *) ANSIRESET, ''

    end subroutine consoleOutputVOnly

    !Given vector, get clock/cone angle and magnitude
    function ClockConeMag(V) result(aVec)
        real(rp), dimension(NDIM), intent(in) :: V
        real(rp), dimension(NDIM) :: aVec

        real(rp) :: MagV
        MagV = norm2(V)
        aVec(2) = atan2(V(YDIR),V(ZDIR) )*180.0/PI !Clock angle
        if (aVec(2) < 0) then
            aVec(2) = aVec(2) + 360.0
        endif
        
        if (MagV>TINY) then
            aVec(3) = acos (V(XDIR)/MagV)*180.0/PI
        else
            aVec(3) = 0.0
        endif
        aVec(1) = MagV
    end function ClockConeMag

    subroutine resOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Write Gamera restart
        call resOutput(gApp%Model,gApp%Grid,gApp%State)

        !Write Voltron restart data
        call resOutputVOnly(vApp)

    end subroutine resOutputV

    subroutine resOutputVOnly(vApp)
        class(voltApp_T), intent(inout) :: vApp

        !Write inner mag restart
        if (vApp%doDeep) then
            call InnerMagRestart(vApp,vApp%IO%nRes)
        endif
        if (vApp%time>vApp%IO%tRes) then
            vApp%IO%tRes = vApp%IO%tRes + vApp%IO%dtRes
        endif
        vApp%IO%nRes = vApp%IO%nRes + 1

    end subroutine resOutputVOnly

    subroutine fOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Write gamera data
        call fOutput(gApp%Model,gApp%Grid,gApp%State) !Gamera

        !Write voltron data
        call fOutputVOnly(vApp)

    end subroutine fOutputV

    subroutine fOutputVOnly(vApp)
        class(voltApp_T), intent(inout) :: vApp

        !Write ReMIX data
        call writeMix(vApp%remixApp%ion,vApp%IO%nOut,mjd=vApp%MJD,time=vApp%time)

        !Write inner mag IO if needed
        if (vApp%doDeep) then
            call InnerMagIO(vApp,vApp%IO%nOut)
        endif

        if (vApp%time>vApp%IO%tOut) then
            vApp%IO%tOut = vApp%IO%tOut + vApp%IO%dtOut
        endif
        vApp%IO%nOut = vApp%IO%nOut + 1

    end subroutine fOutputVOnly

    subroutine getCPCP(mhdvarsin,cpcp)
        real(rp), dimension(:,:,:,:,:),intent(in) :: mhdvarsin
        real(rp), intent(out) :: cpcp(2)
        cpcp(NORTH) = maxval(mhdvarsin(1,:,:,MHDPSI,NORTH))-minval(mhdvarsin(1,:,:,MHDPSI,NORTH))
        cpcp(SOUTH) = maxval(mhdvarsin(1,:,:,MHDPSI,SOUTH))-minval(mhdvarsin(1,:,:,MHDPSI,SOUTH))

    end subroutine getCPCP
end module voltio

