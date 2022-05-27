!Routines to write restart/outputs from Voltron
module voltio
    use gamapp
    use volttypes
    use cmiutils
    use mixio
    use clocks
    use innermagsphere
    use wind
    use dyncoupling
    use dstutils
    
    implicit none

    integer , parameter, private :: MAXVOLTIOVAR = 50
    real(rp), parameter, private :: dtWallMax = 1.0 !How long between timer resets[hr]
    logical , private :: isConInit = .false.
    real(rp), private ::  oMJD = 0.0
    integer , private :: oTime = 0.0
    real(rp), private :: gamWait = 0.0
    character(len=strLen), private :: vh5File

    contains

    subroutine consoleOutputV(vApp,gApp)
        class(gamApp_T) , intent(inout) :: gApp
        class(voltApp_T), intent(inout) :: vApp

        !Using console output from Gamera
        call consoleOutput(gApp%Model,gApp%Grid,gApp%State)

        !Using console output from Voltron
        call consoleOutputVOnly(vApp,gApp,gApp%Model%MJD0)

    end subroutine consoleOutputV

    subroutine consoleOutputVOnly(vApp,gApp,MJD0)
        class(voltApp_T), intent(inout) :: vApp
        class(gamApp_T) , intent(in) :: gApp
        real(rp), intent(in) :: MJD0

        real(rp) :: cpcp(2) = 0.0

        real(rp) :: dpT,dtWall,cMJD,dMJD,simRate

        integer :: nTh,curCount,countMax
        real(rp) :: clockRate
        character(len=strLen) :: utStr
        real(rp) :: DelD,DelP,dtIM,BSDst0,AvgBSDst,DPSDst,symh,BSSMRs(4)

        !Augment Gamera console output w/ Voltron stuff
        call getCPCP(vApp%mix2mhd%mixOutput,cpcp)

        dpT = vApp%tilt%evalAt(vApp%time)*180.0/PI

        !Figure out some perfromance info
        cMJD = T2MJD(vApp%time,MJD0) !Current MJD

        if (isConInit) then
            !Console output has been initialized
            dMJD = cMJD - oMJD !Elapsed MJD since first recorded value
            call system_clock(curCount,clockRate,countMax)
            dtWall = (curCount - oTime)/clockRate
            if(dtWall < 0) dtWall = dtWall + countMax / clockRate
            simRate = dMJD*24.0*60.0*60.0/dtWall !Model seconds per wall second
            gamWait = 0.8*gamWait + 0.2*readClock('GameraSync')/(readClock(1)+TINY) ! Weighted average to self-correct
        else
            simRate = 0.0
            oMJD = cMJD
            call system_clock(count=oTime)
            isConInit = .true.
            dtWall = 0.0
            gamWait = 0.0
        endif

        if ( (simRate<0) .or. (abs(dtWall/3600.0) >= dtWallMax) ) then
            ! Partially reset counters so that the values don't become so large they don't change
            oMJD = cMJD - 0.1*dMJD
            oTime = curCount - 0.1*dtWall*clockRate
            if(oTime < 0) oTime = oTime + countMax
        endif
        
        !Get MJD info
        call mjd2utstr(cMJD,utStr)

        !Get Dst estimate: DPS, center of earth, MLT avg of equatorial stations
        call EstDST(gApp%Model,gApp%Grid,gApp%State,BSDst0,AvgBSDst,BSSMRs,DPSDst)

        vApp%BSDst = AvgBSDst
        
        !Get symh from input time series
        symh = vApp%symh%evalAt(vApp%time)

        !Get imag ingestion info
        if (vApp%doDeep) then
            call IMagDelta(gApp%Model,gApp%Grid,gApp%State,DelD,DelP,dtIM)
        endif

        if (vApp%isLoud) then
            write(*,*) ANSIBLUE
            write(*,*) 'VOLTRON'
            write (*,'(a,a)')                    '      UT   = ', trim(utStr)
            write (*, '(a,1f8.3,a)')             '      tilt = ' , dpT, ' [deg]'
            write (*, '(a,2f8.3,a)')             '      CPCP = ' , cpcp(NORTH), cpcp(SOUTH), ' [kV, N/S]'
            write (*, '(a, f8.3,a)')             '    Sym-H  = ' , symh  , ' [nT]'
            write (*, '(a, f8.3,a)')             '    BSDst  ~ ' , AvgBSDst , ' [nT]'
            write (*, '(a,4f8.2,a)')             '           dSMRs  ~ ' , BSSMRs-AvgBSDst, ' [nT, 12/18/00/06]'
            !write (*, '(a,4f8.2,a)')             '   BSSMRs  ~ ' , BSSMRs, ' [nT, 12/18/00/06]'

            if (vApp%doDeep .and. (vApp%time>0.0)) then
                write (*, '(a, f8.3,a)')             '   DPSDst  ~ ' , DPSDst, ' [nT]'
                write (*, '(a)'                 )    '   IMag Ingestion'
                write (*, '(a,1f7.2,a,1f7.2,a)' )    '       D/P = ', 100.0*DelD,'% /',100.0*DelP,'%'
                write (*, '(a,1f7.2,a)'         )    '        dt = ', dtIM, ' [s]'
                
                !write (*,'(a,1f8.3,I6,a)')           '      xTrc = ', vApp%rTrc,vApp%nTrc, ' [r/n]'
            endif
            write (*, '(a,1f7.1,a)' ) '   Spent ', gamWait*100.0, '% of time waiting for Gamera'
            if (simRate>TINY) then
                if (vApp%isSeparate) then
                    nTh = NumOMP()
                    write (*, '(a,1f8.3,a,I0,a)')             '    Running @ ', simRate*100.0, '% of real-time (',nTh,' threads)'  
                else
                    write (*, '(a,1f8.3,a)'     )             '    Running @ ', simRate*100.0, '% of real-time'
                endif
            endif
            
            write(*,'(a)',advance="no") ANSIRESET!, ''
        endif

        !Write inner mag console IO if needed
        if (vApp%doDeep) then
            call vApp%imagApp%doConIO(vApp%MJD,vApp%time)
        endif

        !Setup for next output
        vApp%IO%tsNext = vApp%ts + vApp%IO%tsOut
        
        if (vApp%doDynCplDT) then
            call UpdateCouplingCadence(vApp)
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
        call resOutput(gApp%Model,gApp%Grid,gApp%oState,gApp%State)

        !Write Voltron restart data
        call resOutputVOnly(vApp,gApp)

    end subroutine resOutputV

    subroutine resOutputVOnly(vApp, gApp)
        class(voltApp_T), intent(inout) :: vApp
        class(gamApp_T) , intent(inout) :: gApp

        if (vApp%writeFiles) then
            call writeMIXRestart(vApp%remixApp%ion,vApp%IO%nRes,mjd=vApp%MJD,time=vApp%time)
            !Write inner mag restart
            if (vApp%doDeep) then
                call vApp%imagApp%doRestart(vApp%IO%nRes,vApp%MJD,vApp%time)
            endif
            call writeVoltRestart(vApp,gApp)
        endif

        vApp%IO%tRes = vApp%IO%tRes + vApp%IO%dtRes
        vApp%IO%nRes = vApp%IO%nRes + 1

    end subroutine resOutputVOnly

    subroutine writeVoltRestart(vApp,gApp)
        class(voltApp_T), intent(in) :: vApp
        class(gamApp_T) , intent(in) :: gApp

        character(len=strLen) :: ResF,lnResF
        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars

        write (ResF, '(A,A,I0.5,A)') trim(gApp%Model%RunID), ".volt.Res.", vApp%IO%nRes, ".h5"
        call CheckAndKill(ResF)

        call ClearIO(IOVars)

        !Main attributes
        call AddOutVar(IOVars,"nOut",vApp%IO%nOut)
        call AddOutVar(IOVars,"nRes",vApp%IO%nRes)
        call AddOutVar(IOVars,"ts"  ,vApp%ts)
        call AddOutVar(IOVars,"MJD" ,vApp%MJD)
        call AddOutVar(IOVars,"time",vApp%time)

        !Coupling info
        call AddOutVar(IOVars,"ShallowT",vApp%ShallowT)
        call AddOutVar(IOVars,"DeepT"   ,vApp%DeepT)
        call AddOutVar(IOVars,"gBAvg", vApp%mhd2Mix%gBAvg)
        
        call WriteVars(IOVars,.false.,ResF)
        !Create link to latest restart
        write (lnResF, '(A,A,A,A)') trim(gApp%Model%RunID), ".volt.Res.", "XXXXX", ".h5"
        call MapSymLink(ResF,lnResF)

    end subroutine writeVoltRestart

    subroutine readVoltronRestart(vApp,xmlInp)
        class(voltApp_T), intent(inout) :: vApp
        type(XML_Input_T), intent(inout) :: xmlInp

        character(len=strLen) :: ResF,resID,nStr
        type(IOVAR_T), dimension(MAXVOLTIOVAR) :: IOVars
        logical :: fExist
        integer :: nRes,n0

        call xmlInp%Set_Val(resID,"/Kaiju/gamera/restart/resID","msphere")
        call xmlInp%Set_Val(nRes,"/Kaiju/gamera/restart/nRes" ,-1)

        !Get number string
        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') nRes
        endif

        write (ResF, '(A,A,A,A)') trim(resID), ".volt.Res.", trim(nStr), ".h5"
        write(*,*) 'Reading Voltron restart from ', trim(ResF)
        inquire(file=ResF,exist=fExist)
        if (.not. fExist) then
            !Error out and leave
            write(*,*) 'Unable to open input voltron restart file, exiting'
            stop
        endif

        call ClearIO(IOVars)
        call AddInVar(IOVars,"gBAvg")

        call AddInVar(IOVars,"nOut"    ,vTypeO=IOINT)
        call AddInVar(IOVars,"nRes"    ,vTypeO=IOINT)
        call AddInVar(IOVars,"ts"      ,vTypeO=IOINT)
        call AddInVar(IOVars,"MJD"     ,vTypeO=IOREAL)
        call AddInVar(IOVars,"time"    ,vTypeO=IOREAL)
        call AddInVar(IOVars,"ShallowT",vTypeO=IOREAL)
        call AddInVar(IOVars,"DeepT"   ,vTypeO=IOREAL)


        !Get data
        call ReadVars(IOVars,.false.,ResF)

        vApp%IO%nOut  = GetIOInt(IOVars,"nOut")
        vApp%IO%nRes  = GetIOInt(IOVars,"nRes") + 1
        vApp%ts       = GetIOInt(IOVars,"ts")
        vApp%MJD      = GetIOReal(IOVars,"MJD")
        vApp%time     = GetIOReal(IOVars,"time")
        vApp%ShallowT = GetIOReal(IOVars,"ShallowT")
        vApp%DeepT    = GetIOReal(IOVars,"DeepT")

        !Check to see if gB0 is present
        n0 = FindIO(IOVars,"gBAvg")
        if (IOVars(n0)%isDone) then
            if (.not. allocated(vApp%mhd2Mix%gBAvg)) then
                !Allocate this if necessary
                allocate( vApp%mhd2Mix%gBAvg(IOVars(n0)%dims(1),IOVars(n0)%dims(2),IOVars(n0)%dims(3),IOVars(n0)%dims(4)) )
            endif
            call IOArray4DFill(IOVars,"gBAvg",vApp%mhd2Mix%gBAvg)
        else
            write(*,*) "gBAvg not found in Voltron restart, assuming dipole ..."
        endif

    end subroutine readVoltronRestart


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
        
        if(vApp%writeFiles) then
            !Write ReMIX data
            call writeMix(vApp%remixApp%ion,vApp%IO%nOut,mjd=vApp%MJD,time=vApp%time)

            !Write inner mag IO if needed
            if (vApp%doDeep) then
                call vApp%imagApp%doIO(vApp%IO%nOut,vApp%MJD,vApp%time)
            endif

            call WriteVolt(vApp,gApp,vApp%IO%nOut)
        endif

        vApp%IO%tOut = vApp%IO%tOut + vApp%IO%dtOut
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
        real(rp) :: cpcp(2),symh

        real(rp), dimension(:,:,:,:), allocatable :: inEijk,inExyz,gJ,gB,Veb
        real(rp), dimension(:,:,:), allocatable :: psi,D,Cs
        real(rp), dimension(NDIM) :: Exyz,Bdip,xcc
        real(rp) :: Csijk,Con(NVAR)
        real(rp) :: BSDst0,AvgBSDst,DPSDst,BSSMRs(4)

        !Get data
        call getCPCP(vApp%mix2mhd%mixOutput,cpcp)

        !Get symh from input time series
        symh = vApp%symh%evalAt(vApp%time)

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
        allocate(gB (Ni,Gr%js:Gr%je,Gr%ks:Gr%ke,1:NDIM))
        allocate(Veb(Ni,Gr%js:Gr%je,Gr%ks:Gr%ke,1:NDIM))
        allocate(psi(Ni,Gr%js:Gr%je,Gr%ks:Gr%ke)) !Cell-centered potential
        allocate(D  (Ni,Gr%js:Gr%je,Gr%ks:Gr%ke))
        allocate(Cs (Ni,Gr%js:Gr%je,Gr%ks:Gr%ke))

        call Ion2MHD(gApp%Model,gApp%Grid,vApp%mix2mhd%gPsi,inEijk,inExyz,vApp%mix2mhd%rm2g)

        D  = 0.0
        Cs = 0.0

        !Subtract dipole before calculating current
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,Bdip,xcc,Exyz,Con,Csijk)
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
                    if (i == Ni) then
                        Con = gApp%State%Gas(JpSt,j,k,:,BLK)
                        D (i,j,k) = Con(DEN)
                        call CellPress2Cs(gApp%Model,Con,Csijk)
                        Cs(i,j,k) = Csijk
                    endif
                enddo
            enddo
        enddo

        !Get Dst estimate: DPS, center of earth, MLT avg of equatorial stations
        call EstDST(gApp%Model,gApp%Grid,gApp%State,BSDst0,AvgBSDst,BSSMRs,DPSDst)
        vApp%BSDst = AvgBSDst

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

        gB = 0.0
        gB(Ni,:,:,XDIR:ZDIR) = gApp%State%Bxyz(JpSt,Gr%js:Gr%je,Gr%ks:Gr%ke,XDIR:ZDIR)
        call AddOutVar(IOVars,"dBx",gB(:,:,:,XDIR))
        call AddOutVar(IOVars,"dBy",gB(:,:,:,YDIR))
        call AddOutVar(IOVars,"dBz",gB(:,:,:,ZDIR))

        call AddOutVar(IOVars,"D" ,D (:,:,:))
        call AddOutVar(IOVars,"Cs",Cs(:,:,:))


        call AddOutVar(IOVars,"Vx",Veb(:,:,:,XDIR))
        call AddOutVar(IOVars,"Vy",Veb(:,:,:,YDIR))
        call AddOutVar(IOVars,"Vz",Veb(:,:,:,ZDIR))

        call AddOutVar(IOVars,"psi",psi)
        !---------------------
        !Do attributes

        call AddOutVar(IOVars,"cpcpN",cpcp(1))
        call AddOutVar(IOVars,"cpcpS",cpcp(2))

        call AddOutVar(IOVars,"BSDst" ,AvgBSDst)
        call AddOutVar(IOVars,"BSDst0",BSDst0)
        call AddOutVar(IOVars,"DPSDst",DPSDst)
        call AddOutVar(IOVars,"BSSMR12",BSSMRs(1))
        call AddOutVar(IOVars,"BSSMR18",BSSMRs(2))
        call AddOutVar(IOVars,"BSSMR00",BSSMRs(3))
        call AddOutVar(IOVars,"BSSMR06",BSSMRs(4))

        call AddOutVar(IOVars,"SymH",symh)

        call AddOutVar(IOVars,"time" ,vApp%time)
        call AddOutVar(IOVars,"MJD"  ,vApp%MJD)
        call AddOutVar(IOVars,"timestep",vApp%ts)

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

    !Calculate relative difference between source/MHD
    subroutine IMagDelta(Model,Gr,State,dD,dP,dtIM)
        type(Model_T), intent(in)  :: Model
        type(Grid_T) , intent(in)  :: Gr
        type(State_T), intent(in)  :: State
        real(rp)     , intent(out) :: dD,dP,dtIM

        real(rp) :: Dsrc,Dmhd,Psrc,Pmhd,dtP
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
        dtP  = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,dV,doInD,doInP,doIngest) &
        !$OMP private(pCon,pW) &
        !$OMP reduction(+:Dsrc,Dmhd,Psrc,Pmhd,dtP)
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
                        dtP  = dtP  + dV*pW(PRESSURE)*Gr%Gas0(i,j,k,IMTSCL,BLK)
                    endif
                    
                enddo
            enddo
        enddo
        if (Dsrc>TINY) dD = Dmhd/Dsrc
        if (Psrc>TINY) dP = Pmhd/Psrc
        if (Pmhd>TINY) then
            dtIM = Model%Units%gT0*dtP/Pmhd
        else
            dtIM = 0.0
        endif

    end subroutine IMagDelta

end module voltio

