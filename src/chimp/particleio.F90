module particleio

    use ioH5
    use chmpdefs
    use chmpunits
    use tptypes
    use ebtypes
    use tputils
    use files
    implicit none

    character(len=strLen) :: tpOutF

    integer, parameter :: MAXTPVS = 30
    logical :: doMuOut = .false.

    contains

    subroutine initPio(Model,tpState)
        type(chmpModel_T), intent(in) :: Model
        type(tpState_T), intent(in) :: tpState

        type(IOVAR_T), dimension(MAXTPVS) :: IOVars

        !write(tpOutF,'(2a)') trim(adjustl(Model%RunID)),'.h5part'
        write(tpOutF,'(a,a,I0.6,a)') trim(adjustl(Model%RunID)),'.',Model%Nblk,'.h5part'

        ! root info is copied from restart file
        if (Model%isRestart) return

        !Add root info for h5p output file
        call CheckAndKill(tpOutF)
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"Np",tpState%Np)
        call AddOutVar(IOVars,"NpT",tpState%NpT)
        call AddOutVar(IOVars,"m0",Model%m0)
        call AddOutVar(IOVars,"q0" ,Model%q0)
        call WriteVars(IOVars,.true.,tpOutF)
        call ClearIO(IOVars)


    end subroutine initPio

    !Do output slice for TP data
    subroutine writeTP(Model,ebState,tpState,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(tpState_T), intent(in)   :: tpState
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXTPVS) :: IOVars
        integer :: i,Np
        integer, dimension(:), allocatable :: OCb
        real(rp), dimension(:), allocatable :: Kev,Mu
        real(rp), dimension(:,:), allocatable :: Qeq,tpLL
        real(rp) :: K2kev,eb2nT
        real(rp) :: tpEq(NVAREQ)
        logical :: isG,doTPTrc

        associate( TPs=>tpState%TPs )

        !Set scaling factors
        K2kev = mec2*1.0e+3 !K in mu -> keV
        eb2nT = G2nT/ebScl !B in mu -> nT
        !Get various metrics
        Np = tpState%Np

        allocate(Kev(Np),Mu(Np))
        allocate(OCb(Np))
        OCb = 0
        allocate(Qeq(Np,NVAREQ))

        if (Model%doLL) then
            allocate(tpLL(Np,2))
            tpLL = 0.0
        endif

        !$OMP PARALLEL DO default(shared) &
        !$OMP& private(tpEq,isG,doTPTrc) &
        !$OMP& schedule(dynamic)
        do i=1,Np
            Kev(i) = prt2kev(Model,TPs(i))
            if (doMuOut .or. TPs(i)%isGC) then
                Mu(i) = (K2kev/eb2nT)*prtMu(TPs(i),Model%t,Model,ebState)
            else
                Mu(i) = 0.0
            endif

            !Save most recent recorded EQX data
            Qeq(i,:) = TPs(i)%Qeq
            doTPTrc = TPs(i)%isIn .and. (Model%doEQProj .or. Model%doTrc)

            if (doTPTrc) then
                !Trace this particles local field line
                tpEq = tpProject(TPs(i),Model%t,Model,ebState,isG,OCb(i))
                if (Model%doEQProj .and. isG) then
                    !If we're doing EQ-projections and this TP projection was good
                    Qeq(i,:) = tpEq !Save projection
                endif
                !Lat-lon projection (northern hemi)
                if (Model%doLL) then
                    call Map2NH(Model,ebState,TPs(i)%Q(XPOS:ZPOS),Model%t,tpLL(i,1),tpLL(i,2))
                endif
            endif
        enddo

        call ClearIO(IOVars)
        call AddOutVar(IOVars,"id",1.0_dp*tpState%TPs(:)%id)
        call AddOutVar(IOVars,"x",TPs(:)%Q(XPOS))
        call AddOutVar(IOVars,"y",TPs(:)%Q(YPOS))
        call AddOutVar(IOVars,"z",TPs(:)%Q(ZPOS))
        call AddOutVar(IOVars,"K",Kev)
        call AddOutVar(IOVars,"Mu",Mu,uStr="keV/nT")
        call AddOutVar(IOVars,"A" ,rad2deg*TPs(:)%alpha)
        
        !Equatorial values
        call AddOutVar(IOVars,"xeq"  ,Qeq(:,EQX))
        call AddOutVar(IOVars,"yeq"  ,Qeq(:,EQY))
        call AddOutVar(IOVars,"Aeq"  ,rad2deg*Qeq(:,EQALP))
        call AddOutVar(IOVars,"Teq"  ,oTScl  *Qeq(:,EQTIME))
        call AddOutVar(IOVars,"Keq"  ,        Qeq(:,EQKEV))
        call AddOutVar(IOVars,"ebKeq",        Qeq(:,EQKEB))

        if (Model%doLL) then
            call AddOutVar(IOVars,"lat"  ,rad2deg*tpLL(:,1))
            call AddOutVar(IOVars,"lon"  ,rad2deg*tpLL(:,2))
        endif

        !Topology (always outputting)
        call AddOutVar(IOVars,"OCb",1.0_dp*OCb)

        !Wave-particle values
        if (Model%doWPI) then
            call AddOutVar(IOVars,"dAwpi",rad2deg*TPs(:)%dAwpi)
            call AddOutVar(IOVars,"dKwpi",TPs(:)%dKwpi)
            call AddOutVar(IOVars,"Nwpi",TPs(:)%Nwpi)
        endif

        !Logical values
        call AddOutVar(IOVars,"isIn",merge(1.0_rp,0.0_rp,TPs(:)%isIn))

        !Birthdays (if streaming)
        if (Model%doStream) then
            call AddOutVar(IOVars,"T0p",oTScl*TPs(:)%T0p)
        endif
        
        !Scalar attributes
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"MJD",ioTabMJD(ebState%ebTab,Model%t))
        call AddOutVar(IOVars,"TimeValue",oTScl*Model%t)
        call WriteVars(IOVars,.true.,tpOutF,gStr)
        call ClearIO(IOVars)

        end associate
    end subroutine writeTP

    !Get info from an h5p file
    subroutine TPinfo(fIn,Np,Nt,dtStp,T0,m0,q0)
        character(len=strLen), intent(in) :: fIn
        integer, intent(out) :: Np,Nt
        real(rp), intent(out) :: dtStp,T0,m0,q0

        real(rp) :: T1,tpScl
        type(IOVAR_T), dimension(MAXTPVS) :: IOVars
        character(len=strLen) :: gStr0,gStr1
        call ClearIO(IOVars)
        Nt = NumSteps(fIn)
        call AddInVar(IOVars,"Np",vTypeO=IOINT)
        call AddInVar(IOVars,"m0",vTypeO=IOREAL)
        call AddInVar(IOVars,"q0",vTypeO=IOREAL)

        !Get first pass at data
        call ReadVars(IOVars,.false.,fIn)
        Np = int(IOVars(1)%data(1))
        m0 =     IOVars(2)%data(1)
        q0 =     IOVars(3)%data(1)
        call ClearIO(IOVars)

        !TP timing data is always in seconds
        !Use specific scaling (may disagree) with MHD data scaling
        tpScl = 1.0*vc_cgs/L0 !in2s = 1
        !Now get timestep
        if (Nt>1) then
            gStr0 = "Step#0"
            gStr1 = "Step#1"
            call AddInVar(IOVars,"time",vTypeO=IOREAL)
            call ReadVars(IOVars,.false.,fIn,gStr0)
            T0 = tpScl*IOVars(1)%data(1)
            call ReadVars(IOVars,.false.,fIn,gStr1)
            T1 = tpScl*IOVars(1)%data(1)
            dtStp = (T1-T0)
        else
            T0 = 0.0
            dtStp = 0.0
            write(*,*) 'Not enough steps found ...'
            stop
        endif
    end subroutine TPinfo

    !Write restart dump of TP data to "ResF"
    subroutine writeTPres(Model,ebState,tpState,ResF)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(tpState_T), intent(in)   :: tpState
        character(len=*), intent(in)  :: ResF

        type(IOVAR_T), dimension(MAXTPVS) :: IOVars
        real(rp), dimension(:,:), allocatable :: Q,Qeq,ijk0
        integer :: i,Np,minID,maxID

        associate( TPs=>tpState%TPs )

        Np = tpState%Np
        allocate(Q(Np,NVARTP))
        allocate(Qeq(Np,NVAREQ))
        allocate(ijk0(Np,NDIM))

        do i=1,Np
            Q(i,:) = TPs(i)%Q
            Qeq(i,:) = TPs(i)%Qeq
            ijk0(i,:) = TPs(i)%ijk0
        end do

        minID = minval(TPs(:)%id)
        maxID = maxval(TPs(:)%id)

        call StampIO(ResF)

        !Reset IO chain
        call ClearIO(IOVars)

        !Main attributes
        call AddOutVar(IOVars,"nOut",Model%IO%nOut)
        call AddOutVar(IOVars,"nRes",Model%IO%nRes)
        call AddOutVar(IOVars,"ts"  ,Model%ts)
        call AddOutVar(IOVars,"t"   ,Model%t)

        call AddOutVar(IOVars,"Np",tpState%Np)
        call AddOutVar(IOVars,"NpT",tpState%NpT)
        call AddOutVar(IOVars,"m0",Model%m0)
        call AddOutVar(IOVars,"q0",Model%q0)
        !Time steps
        call AddOutVar(IOVars,"ddt",TPs(:)%ddt)

        call AddOutVar(IOVars,"id",1.0_dp*TPs(:)%id)
        call AddOutVar(IOVars,"minID",minID)
        call AddOutVar(IOVars,"maxID",maxID)
        call AddOutVar(IOVars,"Ngc",1.0_dp*TPs(:)%Ngc)
        call AddOutVar(IOVars,"Nfo",1.0_dp*TPs(:)%Nfo)
        !Dynamic varaiables
        call AddOutVar(IOVars,"Q",Q)
        call AddOutVar(IOVars,"alpha",TPs(:)%alpha)
        !Equatorial values
        call AddOutVar(IOVars,"Qeq",Qeq)
        !Last known location in eb grid
        call AddOutVar(IOVars,"ijk0",ijk0)
        !Topology
        call AddOutVar(IOVars,"OCb",1.0_dp*TPs(:)%OCb)
        !Wave-particle interaction diagnostics
        call AddOutVar(IOVars,"dAwpi",TPs(:)%dAwpi)
        call AddOutVar(IOVars,"dKwpi",TPs(:)%dKwpi)
        call AddOutVar(IOVars,"Nwpi",TPs(:)%Nwpi)
        call AddOutVar(IOVars,"xj",TPs(:)%xj)
        call AddOutVar(IOVars,"yj",TPs(:)%yj)
        !Logical values
        call AddOutVar(IOVars,"isIn",merge(1.0_rp,0.0_rp,TPs(:)%isIn))
        call AddOutVar(IOVars,"isGC",merge(1.0_rp,0.0_rp,TPs(:)%isGC))
        call AddOutVar(IOVars,"isInit",merge(1.0_rp,0.0_rp,TPs(:)%isInit))
        !Birthdays/Deaths
        call AddOutVar(IOVars,"T0p",TPs(:)%T0p)
        call AddOutVar(IOVars,"Tfp",TPs(:)%Tfp)

        !Write out, force real precision not IO precision
        call WriteVars(IOVars,.false.,ResF)

        end associate
    end subroutine writeTPres

    subroutine readTPrestart(Model,ebState,tpState,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(inout)   :: ebState
        type(tpState_T), intent(inout)   :: tpState
        type(XML_Input_T), intent(in) :: inpXML

        type(IOVAR_T), dimension(MAXTPVS) :: IOVars
        real(rp), dimension(:), allocatable   :: ids,Ngc,Nfo,OCb,isIn,isInit,isGC
        real(rp), dimension(:,:), allocatable :: Q,Qeq,ijk0
        logical :: doSP !whether input is single precision
        integer :: nRes
        integer :: i,Np,minID,maxID
        character(len=strLen) :: resID,nStr,jStr,inH5
        integer :: nvar,dims(2)

        doSP = .false. !Restarts are always double precision

        ! Generate file name
        call inpXML%Set_Val(resID,"restart/resID",trim(Model%RunID))
        call inpXML%Set_Val(nRes,"restart/nRes" ,-1)

        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write (nStr,'(I0.5)') nRes
        endif
        write (jStr,'(I0)') Model%Nblk
        inH5 = trim(resID)// "." // trim(jStr) // ".Res." // trim(nStr) // ".h5"
        call CheckFileOrDie(inH5,"Restart file not found ...")
        write(*,*) 'Reading restart from : ', trim(inH5)

        !Reset IO chain
        call ClearIO(IOVars)

        call AddInVar(IOVars,"nOut",vTypeO=IOINT)
        call AddInVar(IOVars,"nRes",vTypeO=IOINT)
        call AddInVar(IOVars,"ts"  ,vTypeO=IOINT)
        call AddInVar(IOVars,"t"   ,vTypeO=IOREAL)
        call AddInVar(IOVars,"Np"  ,vTypeO=IOINT)
        call AddInVar(IOVars,"NpT" ,vTypeO=IOINT)
        call AddInVar(IOVars,"m0"  ,vTypeO=IOREAL)
        call AddInVar(IOVars,"q0"  ,vTypeO=IOREAL)
        call AddInVar(IOVars,"minID"  ,vTypeO=IOINT)
        call AddInVar(IOVars,"maxID"  ,vTypeO=IOINT)

        call AddInVar(IOVars,"ddt")
        call AddInVar(IOVars,"id")
        call AddInVar(IOVars,"Ngc")
        call AddInVar(IOVars,"Nfo")
        call AddInVar(IOVars,"Q")
        call AddInVar(IOVars,"alpha")
        call AddInVar(IOVars,"Qeq")
        call AddInVar(IOVars,"ijk0")
        call AddInVar(IOVars,"OCb")
        call AddInVar(IOVars,"dAwpi")
        call AddInVar(IOVars,"dKwpi")
        call AddInVar(IOVars,"Nwpi")
        call AddInVar(IOVars,"xj")
        call AddInVar(IOVars,"yj")
        call AddInVar(IOVars,"T0p")
        call AddInVar(IOVars,"Tfp")
        call AddInVar(IOVars,"isIn")
        call AddInVar(IOVars,"isGC")
        call AddInVar(IOVars,"isInit")

        call ReadVars(IOVars,doSP,inH5)

        ! set variables
        Model%IO%nOut = GetIOInt(IOVars,"nOut")
        Model%IO%nRes = GetIOInt(IOVars,"nRes")
        Model%nOut    = Model%IO%nOut ! to allow backward compatability
        Model%ts      = GetIOInt(IOVars,"ts")
        Model%t       = GetIOReal(IOVars,"t")
        ! updating output times
        Model%tOut    = Model%t
        Model%IO%tOut = Model%t
        Model%IO%tRes = Model%t

        !species info
        Model%m0      = GetIOReal(IOVars,"m0")
        Model%q0      = GetIOReal(IOVars,"q0")

        !min and max IDs if want to reset TP ids in the future
        minID         = GetIOReal(IOVars,"minID")
        maxID         = GetIOReal(IOVars,"maxID")

        tpState%NpT   = GetIOInt(IOVars,"NpT")
        Np            = GetIOInt(IOVars,"Np")
        tpState%Np    = Np

        allocate(tpState%TPs(Np))
        allocate(ids(Np),Ngc(Np),Nfo(Np),OCb(Np))
        allocate(isIn(Np),isGC(Np),isInit(Np))
        allocate(Q(Np,NVARTP))
        allocate(Qeq(Np,NVAREQ))
        allocate(ijk0(Np,NDIM))

        associate( TPs=>tpState%TPs )

        ! integer arrays saved as reals
        call IOArray1DFill(IOVars,"id", ids)
        call IOArray1DFill(IOVars,"Ngc",Ngc)
        call IOArray1DFill(IOVars,"Nfo",Nfo)
        call IOArray1DFill(IOVars,"OCb",OCb)

        TPs(:)%id  = nint(ids)
        TPs(:)%Ngc = nint(Ngc)
        TPs(:)%Nfo = nint(Nfo)
        TPs(:)%OCb = nint(OCb)

        ! logical arrays
        call IOArray1DFill(IOVars,"isIn",isIn)
        call IOArray1DFill(IOVars,"isGC",isGC)
        call IOArray1DFill(IOVars,"isInit",isInit)
        TPs(:)%isIn   = .false.
        TPs(:)%isGC   = .false.
        TPs(:)%isInit = .false.
        TPs(:)%isIn   = isIn > 0.5
        TPs(:)%isGC   = isGC > 0.5
        TPs(:)%isInit = isInit > 0.5

        ! 1D arrays
        call IOArray1DFill(IOVars,"ddt",  TPs(:)%ddt)
        call IOArray1DFill(IOVars,"alpha",TPs(:)%alpha)
        call IOArray1DFill(IOVars,"T0p",  TPs(:)%T0p)
        call IOArray1DFill(IOVars,"Tfp",  TPs(:)%Tfp)
        call IOArray1DFill(IOVars,"dAwpi",TPs(:)%dAwpi)
        call IOArray1DFill(IOVars,"dKwpi",TPs(:)%dKwpi)
        call IOArray1DFill(IOVars,"Nwpi", TPs(:)%Nwpi)
        call IOArray1DFill(IOVars,"xj",   TPs(:)%xj)
        call IOArray1DFill(IOVars,"yj",   TPs(:)%yj)

        ! 2D arrays
        call IOArray2DFill(IOVars,"Q",   Q)
        call IOArray2DFill(IOVars,"Qeq", Qeq)
        call IOArray2DFill(IOVars,"ijk0",ijk0)

        do i=1,Np
            TPs(i)%Q    = Q(i,:)
            TPs(i)%Qeq  = Qeq(i,:)
            TPs(i)%ijk0 = ijk0(i,:)
        end do

        end associate
    end subroutine readTPrestart

end module particleio
