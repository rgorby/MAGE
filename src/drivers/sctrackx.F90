!Driver for pulling space craft trajectory out of data
program sctrackx
    use clocks
    use chmpdefs
    use starter
    use chmpio
    use ebtypes
    use chmpfields
  
    implicit none

    integer, parameter :: MAXIOVAR = 35

    !Spacecraft trajectory type
    type SCTrack_T
        !Everything stored in code units
        integer :: NumP !Number of points
        real(rp), dimension(:), allocatable :: X,Y,Z,T,MJDs
        real(rp), dimension(:,:), allocatable :: Q,E,B !MHD vars
        real(rp), dimension(:)  , allocatable :: Beq,OCb
        real(rp), dimension(:,:), allocatable :: xyz2EQ,xyz2NH !Projection variables
        real(rp), dimension(:), allocatable :: inDom !Lazy real-valued boolean
        logical :: doSmooth
        integer :: Ns=0
    end type SCTrack_T

    !Main data structures
    type(chmpModel_T) :: Model
    type(ebState_T)   :: ebState
    type(XML_Input_T) :: inpXML
    type(SCTrack_T)   :: SCTrack

    integer :: n,topo,nS,nE
    integer , dimension(NDIM) :: ijkG !ijk location guess
    real(rp), dimension(NDIM) :: xyz,Et,Bt,xyzEQ,xyzMin
    real(rp), dimension(NVARMHD) :: Qt
    real(rp) :: R0,R,mlat,mlon,Beq
    logical  :: isIn,doParInT

    !Setup timers
    call initClocks()

    !----------------------------
    !Initialize model and fields
    call goApe(Model,ebState,iXML=inpXML)

    !----------------------------
    !Read SC trajectory data
    call GetTrack(inpXML,SCTrack)

    call GetBounds(inpXML,SCTrack,nS,nE)

    !Find inner radius of grid
    R0 = norm2(ebState%ebGr%xyz(1,1,1,XDIR:ZDIR))
    write(*,*) 'Chopping out values inside of R0 = ', R0

    !Loop over trajectory positions/times
    do n=1,SCTrack%NumP
        !Updates
        Model%t = SCTrack%T(n)
        Model%ts = Model%ts+1

        call Tic("Omega")

    !Update fields to current time
        call Tic("Step")
        
        call updateFields(Model,ebState,Model%t)
        call Toc("Step")
    !Evaluate at specific point on trajectory
        call Tic("Eval")
        xyz = [SCTrack%X(n),SCTrack%Y(n),SCTrack%Z(n)]

        if (n == 1) then
            !Do locate on first try
            call locate(xyz,ijkG,Model,ebState%ebGr,isIn)
        endif
        
        R = norm2(xyz)

        if (R>R0) then
            !Inside domain
            SCTrack%inDom(n) = 1.0
            call ebFields (xyz,Model%t,Model,ebState,Et,Bt,ijkO=ijkG)
            Qt = mhdInterp(xyz,Model%t,Model,ebState,ijkO=ijkG)

            call getEquatorProjection(Model,ebState,xyz,Model%t,xyzEQ)
            call getMagEQ(Model,ebState,xyz,Model%t,xyzMin,Beq,topo)
            call Map2NH(Model,ebState,xyz,Model%t,mlat,mlon)

        else
            !Outside domain
            SCTrack%inDom(n) = 0.0
            Et = 0.0
            Bt = 0.0
            Qt = 0.0
            xyzEQ = 0.0
            mlat = 0.0; mlon = 0.0
            xyzMin = 0.0
            Beq = 0.0
            topo = 0
        endif

        SCTrack%MJDs(n) = MJDAt(ebState%ebTab,Model%t)

        !Store values
        SCTrack%Q(n,:) = Qt
        SCTrack%B(n,:) = Bt
        SCTrack%E(n,:) = Et
        SCTrack%xyz2NH(n,:) = (180.0/PI)*[mlat,mlon]
        SCTrack%xyz2EQ(n,:) = [xyzEQ(XDIR),xyzEQ(YDIR)]
        SCTrack%Beq(n) = Beq
        SCTrack%OCb(n) = 1.0*topo !Lazily storing integer as real
        call Toc("Eval")

        call Toc("Omega")
    enddo
    
    write(*,*) 'Done trajectory cut, outputting ...'
    !Done scraping trajectory, now output data
    call OutputTrack(Model,SCTrack,nS,nE)

    contains

        subroutine OutputTrack(Model,SCTrack,nS,nE)
            type(chmpModel_T), intent(in) :: Model
            type(SCTrack_T)  , intent(inout)   :: SCTrack
            character(len=strLen) :: H5Out
            type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
            integer, intent(in) :: nS,nE

            integer :: n
            if (doParInT) then
                write(H5Out,'(a,a,I0.4,a)') trim(adjustl(Model%RunID)),'.',Model%Nblk,'.sc.h5'
            else
                write(H5Out,'(2a)') trim(adjustl(Model%RunID)),'.sc.h5'
            endif

            call CheckAndKill(H5Out)

            if (SCTrack%doSmooth .and. (.not. doParInT) ) then
                write(*,*) 'Smoothing data ...'
                do n=1,NVARMHD
                    call SmoothTS(SCTrack%Q(:,n),SCTrack%inDom,SCTrack%NumP,SCTrack%Ns)
                enddo
                do n=1,NDIM
                    call SmoothTS(SCTrack%B(:,n),SCTrack%inDom,SCTrack%NumP,SCTrack%Ns)
                    call SmoothTS(SCTrack%E(:,n),SCTrack%inDom,SCTrack%NumP,SCTrack%Ns)
                enddo
            endif
            
            !Setup output chain
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"nS",nS)
            call AddOutVar(IOVars,"nE",nE)
            call AddOutVar(IOVars,"NumP",SCTrack%NumP)

            call AddOutVar(IOVars,"X"    ,SCTrack%X       (nS:nE),uStr="SM-Re")
            call AddOutVar(IOVars,"Y"    ,SCTrack%Y       (nS:nE),uStr="SM-Re")
            call AddOutVar(IOVars,"Z"    ,SCTrack%Z       (nS:nE),uStr="SM-Re")
            call AddOutVar(IOVars,"T"    ,oTScl*SCTrack%T (nS:nE),uStr="s")
            call AddOutVar(IOVars,"MJDs" ,SCTrack%MJDs    (nS:nE),uStr="days")
            call AddOutVar(IOVars,"inDom",SCTrack%inDom   (nS:nE),uStr="BOOLEAN")

            !Output field variables
            call AddOutVar(IOVars,"Bx",oBScl*SCTrack%B(nS:nE,XDIR),uStr="nT")
            call AddOutVar(IOVars,"By",oBScl*SCTrack%B(nS:nE,YDIR),uStr="nT")
            call AddOutVar(IOVars,"Bz",oBScl*SCTrack%B(nS:nE,ZDIR),uStr="nT")

            call AddOutVar(IOVars,"Beq",oBScl*SCTrack%Beq(nS:nE),uStr="nT")
            call AddOutVar(IOVars,"OCb",oBScl*SCTrack%OCb(nS:nE),uStr="Topo")

            call AddOutVar(IOVars,"Ex",oEScl*SCTrack%E(nS:nE,XDIR),uStr="mV/m")
            call AddOutVar(IOVars,"Ey",oEScl*SCTrack%E(nS:nE,YDIR),uStr="mV/m")
            call AddOutVar(IOVars,"Ez",oEScl*SCTrack%E(nS:nE,ZDIR),uStr="mV/m")

            !Output fluid variables
            call AddOutVar(IOVars,"D"       ,SCTrack%Q(nS:nE,DEN),uStr="#/cc")
            call AddOutVar(IOVars,"P"       ,SCTrack%Q(nS:nE,PRESSURE),uStr="nPa")
            call AddOutVar(IOVars,"Vx",oVScl*SCTrack%Q(nS:nE,VELX),uStr="km/s")
            call AddOutVar(IOVars,"Vy",oVScl*SCTrack%Q(nS:nE,VELY),uStr="km/s")
            call AddOutVar(IOVars,"Vz",oVScl*SCTrack%Q(nS:nE,VELZ),uStr="km/s")

            !Output projection variables
            call AddOutVar(IOVars,"xeq" ,SCTrack%xyz2EQ(nS:nE,1),uStr="SM-Re")
            call AddOutVar(IOVars,"yeq" ,SCTrack%xyz2EQ(nS:nE,2),uStr="SM-Re")
            call AddOutVar(IOVars,"MLAT",SCTrack%xyz2NH(nS:nE,1),uStr="deg")
            call AddOutVar(IOVars,"MLON",SCTrack%xyz2NH(nS:nE,2),uStr="deg")

            !Let loose (do double precision)
            call WriteVars(IOVars,.false.,H5Out)

        end subroutine OutputTrack

        !Get bounds if doing parallel in time
        subroutine GetBounds(inpXML,SCTrack,nS,nE)
            type(XML_Input_T), intent(in ) :: inpXML
            type(SCTrack_T)  , intent(in ) :: SCTrack
            integer          , intent(out) :: nS,nE

            integer :: NumB,dN

            call inpXML%Set_Val(NumB,'parintime/NumB',NumB)
            if (NumB > 1) then
                doParInT = .true.
                if ( (Model%Nblk>NumB) .or. (Model%Nblk<1) ) then
                    write(*,*) "This block outside of acceptable bounds"
                    write(*,*) "Block = ",Model%Nblk
                    write(*,*) "Bounds = ",1,NumB
                    write(*,*) "Bailing ..."
                    stop
                endif

                !Now set nS,nE for this part of work
                dN = nint(1.0*SCTrack%NumP/NumB)
                nS = (Model%Nblk-1)*dN + 1 !Always start at 1
                nE = nS+dN-1

                if (NumB == 1) nS = 1
                if (NumB == Model%Nblk) nE = SCTrack%NumP

                write(*,*) '// In Time'
                write(*,*) '   Doing ', nS,nE, ' out of ', SCTrack%NumP
            else
                doParInT = .false.
                !Not // in time
                nS = 1
                nE = SCTrack%NumP
            endif

            if (SCTrack%doSmooth .and. doParInT) then
                write(*,*) 'Warning, smoothing and // in time not implemented (why not do that?)'
                stop
            endif

        end subroutine GetBounds

        subroutine GetTrack(inpXML,SCTrack)
            type(XML_Input_T), intent(in)      :: inpXML
            type(SCTrack_T)  , intent(inout)   :: SCTrack

            type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
            character(len=strLen) :: H5In
            integer :: Nt

            !Start by getting file and checking it exists
            call inpXML%Set_Val(H5In,"trajectory/H5Traj","sctrack.h5")
            call CheckFileOrDie(H5In,"Trajectory file not found ...")

            call inpXML%Set_Val(SCTrack%doSmooth,"trajectory/doSmooth",.false.)
            if (SCTrack%doSmooth) then
                call inpXML%Set_Val(SCTrack%Ns,"trajectory/Ns",2)
            endif

            !Setup input chain
            call ClearIO(IOVars)
            call AddInVar(IOVars,"X")
            call AddInVar(IOVars,"Y")
            call AddInVar(IOVars,"Z")
            call AddInVar(IOVars,"T")

            call ReadVars(IOVars,.false.,H5In) !Don't use io precision

            Nt = IOVars(1)%N
            SCTrack%NumP = Nt
            
            allocate(SCTrack%X(Nt))
            allocate(SCTrack%Y(Nt))
            allocate(SCTrack%Z(Nt))
            allocate(SCTrack%T(Nt))
            allocate(SCTrack%MJDs(Nt))
            allocate(SCTrack%inDom(Nt))

            allocate(SCTrack%Beq(Nt))
            allocate(SCTrack%OCb(Nt))

            allocate(SCTrack%Q(Nt,NVARMHD))
            allocate(SCTrack%E(Nt,NDIM))
            allocate(SCTrack%B(Nt,NDIM))
            allocate(SCTrack%xyz2EQ(Nt,2))
            allocate(SCTrack%xyz2NH(Nt,2))
            
            call IOArray1DFill(IOVars,"X",SCTrack%X)
            call IOArray1DFill(IOVars,"Y",SCTrack%Y)
            call IOArray1DFill(IOVars,"Z",SCTrack%Z)
            call IOArray1DFill(IOVars,"T",SCTrack%T)

            !Scale input time
            SCTrack%T = inTScl*SCTrack%T

        end subroutine GetTrack

        !Do 3-pt smoothing window in place on time series
        subroutine SmoothTS(Q,inDom,Nt,Ns)
            real(rp), intent(inout) :: Q(Nt)
            real(rp), intent(in)    :: inDom(Nt)
            integer, intent(in) :: Nt,Ns

            real(rp), dimension(:), allocatable :: Qs
            logical , dimension(:), allocatable :: isIn

            integer :: n
            real(rp) :: w

            allocate(Qs(Nt))
            allocate(isIn(Nt))

            isIn = (inDom>0.5)

            Qs = Q

            do n=1+Ns,Nt-Ns
                if (any(isIn(n-Ns:n+Ns))) then
                    Qs(n) = sum(Q(n-Ns:n+Ns),mask=isIn(n-Ns:n+Ns))/count(isIn(n-Ns:n+Ns))
                else
                    Qs(n) = Q(n)
                endif
            enddo

            Q = Qs

        end subroutine SmoothTS
end program sctrackx