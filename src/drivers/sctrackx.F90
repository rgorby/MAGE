!Driver for pulling space craft trajectory out of data
program sctrackx
    use clocks
    use chmpdefs
    use starter
    use chmpio
    use ebtypes
    use chmpfields
  
    implicit none

    integer, parameter :: MAXIOVAR = 15

    !Spacecraft trajectory type
    type SCTrack_T
        !Everything stored in code units
        integer :: NumP !Number of points
        real(rp), dimension(:), allocatable :: X,Y,Z,T,MJDs
        real(rp), dimension(:,:), allocatable :: Q,E,B !MHD vars
    end type SCTrack_T

    !Main data structures
    type(chmpModel_T) :: Model
    type(ebState_T)   :: ebState
    type(XML_Input_T) :: inpXML
    type(SCTrack_T)   :: SCTrack

    integer :: n
    real(rp), dimension(NDIM) :: xyz,Et,Bt
    real(rp), dimension(NVARMHD) :: Qt

    !Setup timers
    call initClocks()

    !----------------------------
    !Initialize model and fields
    call goApe(Model,ebState,iXML=inpXML)

    !----------------------------
    !Read SC trajectory data
    call GetTrack(inpXML,SCTrack)

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
        call ebFields(xyz,Model%t,Model,ebState,Et,Bt)
        Qt = mhdInterp(xyz,Model%t,Model,ebState)
        SCTrack%MJDs(n) = MJDAt(ebState%ebTab,Model%t)

        !Store values
        SCTrack%Q(n,:) = Qt
        SCTrack%B(n,:) = Bt
        SCTrack%E(n,:) = Et
        call Toc("Eval")

        call Toc("Omega")
    enddo
    
    write(*,*) 'Done trajectory cut, outputting ...'
    !Done scraping trajectory, now output data
    call OutputTrack(Model,SCTrack)

    contains

        subroutine OutputTrack(Model,SCTrack)
            type(chmpModel_T), intent(in) :: Model
            type(SCTrack_T)  , intent(inout)   :: SCTrack
            character(len=strLen) :: H5Out
            type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

            write(H5Out,'(2a)') trim(adjustl(Model%RunID)),'.sc.h5'
            call CheckAndKill(H5Out)

            !Setup output chain
            call ClearIO(IOVars)
            
            call AddOutVar(IOVars,"X",SCTrack%X,uStr="SM-Re")
            call AddOutVar(IOVars,"Y",SCTrack%Y,uStr="SM-Re")
            call AddOutVar(IOVars,"Z",SCTrack%Z,uStr="SM-Re")
            call AddOutVar(IOVars,"T",oTScl*SCTrack%T,uStr="s")
            call AddOutVar(IOVars,"MJDs",SCTrack%MJDs)

            !Output field variables
            call AddOutVar(IOVars,"Bx",oBScl*SCTrack%B(:,XDIR),uStr="nT")
            call AddOutVar(IOVars,"By",oBScl*SCTrack%B(:,YDIR),uStr="nT")
            call AddOutVar(IOVars,"Bz",oBScl*SCTrack%B(:,ZDIR),uStr="nT")

            call AddOutVar(IOVars,"Ex",oEScl*SCTrack%E(:,XDIR),uStr="mV/m")
            call AddOutVar(IOVars,"Ey",oEScl*SCTrack%E(:,YDIR),uStr="mV/m")
            call AddOutVar(IOVars,"Ez",oEScl*SCTrack%E(:,ZDIR),uStr="mV/m")

            !Output fluid variables
            call AddOutVar(IOVars,"D" ,SCTrack%Q(:,DEN),uStr="#/cc")
            call AddOutVar(IOVars,"P" ,SCTrack%Q(:,PRESSURE),uStr="nPa")
            call AddOutVar(IOVars,"Vx",oVScl*SCTrack%Q(:,VELX),uStr="km/s")
            call AddOutVar(IOVars,"Vy",oVScl*SCTrack%Q(:,VELY),uStr="km/s")
            call AddOutVar(IOVars,"Vz",oVScl*SCTrack%Q(:,VELZ),uStr="km/s")

            !Let loose
            call WriteVars(IOVars,.true.,H5Out)

        end subroutine OutputTrack

        subroutine GetTrack(inpXML,SCTrack)
            type(XML_Input_T), intent(in)      :: inpXML
            type(SCTrack_T)  , intent(inout)   :: SCTrack

            type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
            character(len=strLen) :: H5In
            integer :: Nt

            !Start by getting file and checking it exists
            call inpXML%Set_Val(H5In,"trajectory/H5Traj","sctrack.h5")
            call CheckFileOrDie(H5In,"Trajectory file not found ...")

            !Setup input chain
            call ClearIO(IOVars)
            call AddInVar(IOVars,"X")
            call AddInVar(IOVars,"Y")
            call AddInVar(IOVars,"Z")
            call AddInVar(IOVars,"T")

            call ReadVars(IOVars,.false.,H5In) !Don't use io precision

            Nt = IOVars(1)%N
            SCTrack%NumP = Nt
            !write(*,*) 'N = ', Nt
            
            allocate(SCTrack%X(Nt))
            allocate(SCTrack%Y(Nt))
            allocate(SCTrack%Z(Nt))
            allocate(SCTrack%T(Nt))
            allocate(SCTrack%MJDs(Nt))

            allocate(SCTrack%Q(Nt,NVARMHD))
            allocate(SCTrack%E(Nt,NDIM))
            allocate(SCTrack%B(Nt,NDIM))

            call IOArray1DFill(IOVars,"X",SCTrack%X)
            call IOArray1DFill(IOVars,"Y",SCTrack%Y)
            call IOArray1DFill(IOVars,"Z",SCTrack%Z)
            call IOArray1DFill(IOVars,"T",SCTrack%T)

            !Scale input time
            SCTrack%T = inTScl*SCTrack%T

        end subroutine GetTrack
end program sctrackx