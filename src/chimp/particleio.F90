module particleio

    use ioH5
    use chmpdefs
    use chmpunits
    use tptypes
    use ebtypes
    use tputils
    use files
    use ebtabutils
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
            call AddOutVar(IOVars,"dKwpi",K2kev*TPs(:)%dKwpi)
        endif

        !Logical values
        call AddOutVar(IOVars,"isIn",merge(1.0_rp,0.0_rp,TPs(:)%isIn))

        !Birthdays (if streaming)
        if (Model%doStream) then
            call AddOutVar(IOVars,"T0p",oTScl*TPs(:)%T0p)
        endif
        
        !Scalar attributes
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"MJD",MJDAt(ebState%ebTab,Model%t))
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

end module particleio
