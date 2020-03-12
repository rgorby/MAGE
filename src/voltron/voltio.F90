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
        real(rp) :: cpcp(2) = 0.0

        real(rp) :: dpT,dtWall,cMJD,dMJD,simRate
        real(rp) :: dSW,pSW,Dst
        real(rp), dimension(NDIM) :: xW,bSW,vSW,bAng

        integer :: iYr,iDoY,iMon,iDay,iHr,iMin
        real(rp) :: rSec
        character(len=strLen) :: utStr

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
                if (gApp%Grid%hasupperBC(IDIR)) then
                    call GetWindAt(pWind,gApp%Model,xW,gApp%Model%t,dSW,pSW,vSW,bSW)
                endif
            class default
                write(*,*) 'WTF?'
        end select
        vSW = vSW*1.0e+2
        bSW = bSW*gApp%Model%Units%gB0
        bAng = ClockConeMag(bSW)

    !Get MJD info
        call mjd2ut(cMJD,iYr,iDoY,iMon,iDay,iHr,iMin,rSec)
        write(utStr,'(I0.4,a,I0.2,a,I0.2,a,I0.2,a,I0.2,a,I0.2)') iYr,'-',iMon,'-',iDay,' ',iHr,':',iMin,':',nint(rSec)

    !Get Dst estimate
        call EstDST(gApp%Model,gApp%Grid,gApp%State,Dst)

        write(*,'(a)',advance="no") ANSIBLUE
        write (*, '(a,f7.2,a,3f8.2,a)')      '     Wind = ' , dSW,     ' [#/cc] / ',vSW,' [km/s, XYZ]'
        write (*, '(a,f7.2,a,2f7.2,a)')      '       IMF = ' , bAng(1), '   [nT] / ',bAng(2),bAng(3),' [deg, Clock/Cone]'
        write (*, '(a,1f8.3,a)')             '      tilt = ' , dpT, ' [deg]'
        write (*, '(a,2f8.3,a)')             '      CPCP = ' , cpcp(NORTH), cpcp(SOUTH), ' [kV, N/S]'
        write (*, '(a, f8.3,a)')             '    BSDst  ~ ' , Dst, ' [nT]'
        write (*,'(a,a)')                    '      UT   = ', trim(utStr)
        write (*, '(a,1f7.3,a)')             '      Running @ ', simRate*100.0, '% of real-time'
        
        write (*, *) ANSIRESET, ''

    end subroutine consoleOutputV

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

        !Write inner mag restart
        if (vApp%doDeep) then
            call InnerMagRestart(vApp,vApp%IO%nRes)
        endif
        if (vApp%time>vApp%IO%tRes) then
            vApp%IO%tRes = vApp%IO%tRes + vApp%IO%dtRes
        endif
        vApp%IO%nRes = vApp%IO%nRes + 1            

    end subroutine resOutputV

    subroutine fOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Write gamera data
        call fOutput(gApp%Model,gApp%Grid,gApp%State) !Gamera

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

    end subroutine fOutputV

    subroutine getCPCP(mhdvarsin,cpcp)
        real(rp), dimension(:,:,:,:,:),intent(in) :: mhdvarsin
        real(rp), intent(out) :: cpcp(2)
        cpcp(NORTH) = maxval(mhdvarsin(1,:,:,MHDPSI,NORTH))-minval(mhdvarsin(1,:,:,MHDPSI,NORTH))
        cpcp(SOUTH) = maxval(mhdvarsin(1,:,:,MHDPSI,SOUTH))-minval(mhdvarsin(1,:,:,MHDPSI,SOUTH))

    end subroutine getCPCP

    !Use Gamera data to estimate DST
    !(Move this to msphutils?)
    subroutine EstDST(Model,Gr,State,Dst)
        type(Model_T), intent(in)  :: Model
        type(Grid_T) , intent(in)  :: Gr
        type(State_T), intent(in)  :: State
        real(rp)     , intent(out) :: Dst

        integer :: i,j,k
        real (rp), dimension(:,:,:,:), allocatable :: dB,Jxyz !Full-sized arrays

        real(rp), dimension(NDIM) :: xyz,xyz0
        integer :: iMax,iMin

        real(rp) :: dV,r,bs1,bs2,bScl,dBz
        real(rp) :: mu0,d0,u0,B0

        !Very lazy scaling
        mu0 = 4*PI*1.0e-7
        d0 = (1.67e-27)*1.0e+6
        u0 = 1.0e+5
        B0 = sqrt(mu0*d0*u0*u0)*1.0e+9 !nT

        call allocGridVec(Model,Gr,dB  )
        call allocGridVec(Model,Gr,Jxyz)
        
        !Subtract dipole before calculating current
        !$OMP PARALLEL DO default(shared) collapse(2)
        do k=Gr%ksg,Gr%keg
            do j=Gr%jsg,Gr%jeg
                do i=Gr%isg,Gr%ieg
                    dB(i,j,k,:) = State%Bxyz(i,j,k,:) + Gr%B0(i,j,k,:) - MagsphereDipole(Gr%xyzcc(i,j,k,:),Model%MagM0)
                enddo
            enddo
        enddo
        !Calculate current
        call bFld2Jxyz(Model,Gr,dB,Jxyz)

        !Set some lazy config
        xyz0 = 0.0 !Measure at center of Earth
        iMin = Gr%is+4
        iMax = Gr%ie

        !Now do accumulation
        Dst = 0.0
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyz,dV,r,bs1,bs2,bScl,dBz) &
        !$OMP reduction(+:Dst)
        do k=Gr%ksg,Gr%keg
            do j=Gr%jsg,Gr%jeg
                do i=iMin,iMax
                    xyz = Gr%xyzcc (i,j,k,:)
                    dV  = Gr%volume(i,j,k)
                    r = norm2(xyz-xyz0)
                    bs1 = Jxyz(i,j,k,XDIR)*(xyz(YDIR)-xyz0(YDIR))
                    bs2 = Jxyz(i,j,k,YDIR)*(xyz(XDIR)-xyz0(XDIR))
                    bScl = B0*dV/(4*PI)

                    dBz = -(bs1 - bs2)/(r**3.0)
                    Dst = Dst + bScl*dBz
                enddo ! i loop
            enddo
        enddo !k loop

    end subroutine EstDST

end module voltio

