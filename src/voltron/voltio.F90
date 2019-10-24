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
        type(gamApp_T) , intent(inout) :: gApp
        type(voltApp_T), intent(in) :: vApp
        real(rp) :: cpcp(2) = 0.0

        real(rp) :: dpT,dtWall,cMJD,dMJD,simRate
        real(rp) :: dSW,pSW
        real(rp), dimension(NDIM) :: xW,bSW,vSW,bAng

        !Using console output from Gamera
        call consoleOutput(gApp%Model,gApp%Grid,gApp%State)

        !Augment Gamera console output w/ Voltron stuff
        call getCPCP(vApp%mix2mhd%mixOutput,cpcp)
        dpT = vApp%tilt%evalAt(vApp%time)*180.0/PI

        !Figure out some perfromance info
        cMJD = T2MJD(vApp%time,gApp%Model%MJD0) !Current MJD

        
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
        
        !Pull solar wind info @ Earth
        xW = 0.0
        select type(pWind=>gApp%Grid%externalBCs(OUTI)%p)
            type is (WindBC_T)
                if (gApp%Model%Ri == gApp%Model%NumRi) then
                    call GetWindAt(pWind,gApp%Model,xW,gApp%Model%t,dSW,pSW,vSW,bSW)
                endif
            class default
                write(*,*) 'WTF?'
        end select
        !dSW = dSW*gApp%Model%Units%gD0
        !pSW = pSW*gApp%Model%Units%gP0
        vSW = vSW*1.0e+2
        bSW = bSW*gApp%Model%Units%gB0
        bAng = ClockConeMag(bSW)

        
        write(*,'(a)',advance="no") ANSIBLUE
        !write (*, '(a,f8.3,a)')       '    dt/dt0 = ', 100*Model%dt/dt0, '%'
        write (*, '(a,f7.2,a,3f8.2,a)')      '     Wind = ' , dSW,     ' [#/cc] / ',vSW,' [km/s, XYZ]'
        write (*, '(a,f7.2,a,2f7.2,a)')      '       IMF = ' , bAng(1), '   [nT] / ',bAng(2),bAng(3),' [deg, Clock/Cone]'
        write (*, '(a,2f8.3,a)')      '     CPCP = ' , cpcp(NORTH), cpcp(SOUTH), ' [kV, N/S]'
        write (*, '(a,1f8.3,a)')      '     tilt = ' , dpT, ' [deg]'
        write (*, '(a,1f7.3,a)')      '     Running @ ', simRate*100.0, '% of real-time'
        
        
        write (*, *) ANSIRESET, ''

    end subroutine consoleOutputV

    !Given vector, get clock/cone angle and magnitude
    function ClockConeMag(V) result(aVec)
        real(rp), dimension(NDIM), intent(in) :: V
        real(rp), dimension(NDIM) :: aVec

        real(rp) :: MagV
        MagV = norm2(V)
        aVec(2) = atan2(V(YDIR),V(ZDIR) )*180.0/PI !Clock angle
        if (MagV>TINY) then
            aVec(3) = acos (V(XDIR)/MagV)*180.0/PI
        else
            aVec(3) = 0.0
        endif
        aVec(1) = MagV
    end function ClockConeMag

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

        !Write inner mag IO if needed
        if (vApp%doDeep) then
            call InnerMagIO(vApp,vApp%IO%nOut)
        endif

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

