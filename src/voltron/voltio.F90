!Routines to write restart/outputs from Voltron
module voltio
    use gamapp
    use volttypes
    use mixio
    use clocks
    use innermagsphere
    use wind

    implicit none

    integer, parameter, private :: MAXVOLTIOVAR = 20
    logical, private :: isConInit = .false.
    real(rp), private ::  oMJD = 0.0
    character(len=strLen), private :: vh5File

    contains

    subroutine consoleOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(in) :: vApp

        !Using console output from Gamera
        call consoleOutput(gApp%Model,gApp%Grid,gApp%State)

        !Using console output from Voltron
        call consoleOutputVOnly(vApp,gApp,gApp%Model%MJD0)

    end subroutine consoleOutputV

    subroutine consoleOutputVOnly(vApp,gApp,MJD0)
        class(voltApp_T), intent(in) :: vApp
        class(gamApp_T) , intent(in) :: gApp
        real(rp), intent(in) :: MJD0

        real(rp) :: cpcp(2) = 0.0

        real(rp) :: dpT,dtWall,cMJD,dMJD,simRate

        integer :: iYr,iDoY,iMon,iDay,iHr,iMin
        real(rp) :: rSec
        character(len=strLen) :: utStr
        real(rp) :: dD,dP,Dst

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

        !Get Dst estimate
        call EstDST(gApp%Model,gApp%Grid,gApp%State,Dst)

        if (vApp%isLoud) then
            write(*,*) ANSIBLUE
            write(*,*) 'VOLTRON'
            write (*,'(a,a)')                    '      UT   = ', trim(utStr)
            write (*, '(a,1f8.3,a)')             '      tilt = ' , dpT, ' [deg]'
            write (*, '(a,2f8.3,a)')             '      CPCP = ' , cpcp(NORTH), cpcp(SOUTH), ' [kV, N/S]'
            write (*, '(a, f8.3,a)')             '    BSDst  ~ ' , Dst, ' [nT]'
            write (*, '(a,1f7.3,a)')             '      Running @ ', simRate*100.0, '% of real-time'
            
            write (*, *) ANSIRESET, ''
        endif

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
        call fOutputVOnly(vApp,gApp)

    end subroutine fOutputV

    subroutine fOutputVOnly(vApp,gApp)
        class(voltApp_T), intent(inout) :: vApp
        class(gamApp_T) , intent(inout) :: gApp

        !Write ReMIX data
        call writeMix(vApp%remixApp%ion,vApp%IO%nOut,mjd=vApp%MJD,time=vApp%time)

        !Write inner mag IO if needed
        if (vApp%doDeep) then
            call InnerMagIO(vApp,vApp%IO%nOut)
        endif

        call WriteVolt(vApp,gApp,vApp%IO%nOut)

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

    !Output voltron data
    subroutine WriteVolt(vApp,gApp,nOut)
        class(voltApp_T), intent(inout) :: vApp
        type(gamApp_T)  , intent(inout) :: gApp
        integer, intent(in) :: nOut

        integer :: Ni,Nj,Nk,Njp,Nkp
        integer :: i,j,k
        character(len=strLen) :: gStr
        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars
        real(rp) :: cpcp(2)

        real(rp), dimension(:,:,:,:), allocatable :: inEijk,inExyz,gJ,Veb
        real(rp), dimension(:,:,:), allocatable :: psi
        real(rp), dimension(NDIM) :: Exyz,Bdip,xcc

        !Get data
        call getCPCP(vApp%mix2mhd%mixOutput,cpcp)

        associate(Gr => gApp%Grid)
        !Cell-centers w/ ghosts
        Nj = gApp%Grid%Nj
        Nk = gApp%Grid%Nk
        Ni = vApp%mix2mhd%PsiShells

        !Cell centers
        Njp = gApp%Grid%Njp
        Nkp = gApp%Grid%Nkp

        allocate(inEijk(Ni+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,1:NDIM))
        allocate(inExyz(Ni  ,Gr%jsg:Gr%jeg  ,Gr%ksg:Gr%keg  ,1:NDIM))

        !Active cell centers
        allocate(gJ (Ni,Gr%js:Gr%je,Gr%ks:Gr%ke,1:NDIM))
        allocate(Veb(Ni,Gr%js:Gr%je,Gr%ks:Gr%ke,1:NDIM))
        allocate(psi(Ni,Gr%js:Gr%je,Gr%ks:Gr%ke)) !Cell-centered potential

        call Ion2MHD(gApp%Model,gApp%Grid,vApp%mix2mhd%gPsi,inEijk,inExyz,vApp%mix2mhd%rm2g)

        !Subtract dipole before calculating current
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,Bdip,xcc,Exyz)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=1,Ni
                    psi(i,j,k) = 0.125*( vApp%mix2mhd%gPsi(i+1,j  ,k) + vApp%mix2mhd%gPsi(i+1,j  ,k+1) &
                                       + vApp%mix2mhd%gPsi(i  ,j+1,k) + vApp%mix2mhd%gPsi(i  ,j+1,k+1) &
                                       + vApp%mix2mhd%gPsi(i+1,j+1,k) + vApp%mix2mhd%gPsi(i+1,j+1,k+1) &
                                       + vApp%mix2mhd%gPsi(i  ,j  ,k) + vApp%mix2mhd%gPsi(i  ,j  ,k+1) )
                    xcc = Gr%xyzcc(i,j,k,:)
                    Bdip = MagsphereDipole(xcc,gApp%Model%MagM0)
                    Exyz = inExyz(i,j,k,:)
                    Veb(i,j,k,:) = cross(Exyz,Bdip)/dot_product(Bdip,Bdip)
                enddo
            enddo
        enddo

        write(gStr,'(A,I0)') "Step#", nOut
        !Reset IO chain
        call ClearIO(IOVars)

        call AddOutVar(IOVars,"Ex",inExyz(:,Gr%js:Gr%je,Gr%ks:Gr%ke,XDIR))
        call AddOutVar(IOVars,"Ey",inExyz(:,Gr%js:Gr%je,Gr%ks:Gr%ke,YDIR))
        call AddOutVar(IOVars,"Ez",inExyz(:,Gr%js:Gr%je,Gr%ks:Gr%ke,ZDIR))

        !Add inner currents
        gJ = 0.0
        gJ(Ni,:,:,XDIR:ZDIR) =  vApp%mhd2mix%gJ(1,:,:,XDIR:ZDIR) !Just assuming 1 shell
        call AddOutVar(IOVars,"Jx",gJ(:,:,:,XDIR))
        call AddOutVar(IOVars,"Jy",gJ(:,:,:,YDIR))
        call AddOutVar(IOVars,"Jz",gJ(:,:,:,ZDIR))

        call AddOutVar(IOVars,"Vx",Veb(:,:,:,XDIR))
        call AddOutVar(IOVars,"Vy",Veb(:,:,:,YDIR))
        call AddOutVar(IOVars,"Vz",Veb(:,:,:,ZDIR))

        call AddOutVar(IOVars,"psi",psi)
        call AddOutVar(IOVars,"cpcpN",cpcp(1))
        call AddOutVar(IOVars,"cpcpS",cpcp(2))

        call WriteVars(IOVars,.true.,vh5File,gStr)

        end associate
    end subroutine WriteVolt

    !Initialize Voltron-unique IO
    subroutine InitVoltIO(vApp,gApp)
        class(voltApp_T), intent(inout) :: vApp
        type(gamApp_T)  , intent(inout) :: gApp

        character(len=strLen) :: RunID
        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars
        logical :: fExist, isRestart

        integer :: Ni,Nj,Nk,Ng

        isRestart = gApp%Model%isRestart
        RunID = trim(gApp%Model%RunID)

        !Create filename
        vh5File = trim(RunID) // ".volt.h5" !Voltron output
        fExist = CheckFile(vh5File)
        write(*,*) 'Voltron outputting to ',trim(vh5File)

        if ( (.not. isRestart) .or. (isRestart .and. (.not. fExist)) ) then
            !Not a restart or it is a restart and no file
            call CheckAndKill(vh5File) !For non-restart but file exists

            !Reset IO chain
            call ClearIO(IOVars)

            !Identify part of grid that we want
            Nj = gApp%Grid%Njp+1
            Nk = gApp%Grid%Nkp+1
            Ni = vApp%mix2mhd%PsiShells+1
            Ng = 4 !Number of ghosts

            call AddOutVar(IOVars,"X",gApp%Grid%xyz(1-Ng:Ni-Ng,1:Nj,1:Nk,XDIR))
            call AddOutVar(IOVars,"Y",gApp%Grid%xyz(1-Ng:Ni-Ng,1:Nj,1:Nk,YDIR))
            call AddOutVar(IOVars,"Z",gApp%Grid%xyz(1-Ng:Ni-Ng,1:Nj,1:Nk,ZDIR))

            call AddOutVar(IOVars,"UnitsID","VOLTRON")
            call WriteVars(IOVars,.true.,vh5File)
        endif

    end subroutine InitVoltIO

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
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
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

    !Calculate relative difference between source/MHD
    subroutine IMagDelta(Model,Gr,State,dD,dP)
        type(Model_T), intent(in)  :: Model
        type(Grid_T) , intent(in)  :: Gr
        type(State_T), intent(in)  :: State
        real(rp)     , intent(out) :: dD,dP

        real(rp) :: Dsrc,Dmhd,Psrc,Pmhd
        real(rp) :: dV
        real(rp), dimension(NVAR) :: pCon,pW
        logical :: doInD,doInP,doIngest
        integer :: i,j,k

        dD = 0.0
        dP = 0.0
        if (.not. Model%doSource) return

        !Zero out accumulators
        Dsrc = 0.0
        Dmhd = 0.0
        Psrc = 0.0
        Pmhd = 0.0
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie
                    dV  = Gr%volume(i,j,k)
                    doInD = (Gr%Gas0(i,j,k,IMDEN,BLK)>TINY)
                    doInP = (Gr%Gas0(i,j,k,IMPR ,BLK)>TINY)
                    doIngest = doInD .or. doInP
                    
                    if (.not. doIngest) cycle
                    pCon = State%Gas(i,j,k,:,BLK)
                    call CellC2P(Model,pCon,pW)

                    if (doInD) then
                        Dsrc = Dsrc + dV*Gr%Gas0(i,j,k,IMDEN,BLK)
                        Dmhd = Dmhd + dV*pW(DEN)
                    endif

                    if (doInP) then
                        Psrc = Psrc + dV*Gr%Gas0(i,j,k,IMPR,BLK)
                        Pmhd = Pmhd + dV*pW(PRESSURE)
                    endif
                    
                enddo
            enddo
        enddo
        if (Dsrc>TINY) dD = Dmhd/Dsrc
        if (Psrc>TINY) dP = Pmhd/Psrc
        
    end subroutine IMagDelta

end module voltio

