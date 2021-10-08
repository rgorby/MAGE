!Routines to handle SST inner magnetosphere model -- LonLat model projecting to ionosphere, like RCM
module sstLLimag
    use volttypes
    use ebtypes
    use ebinit
    use ioh5
    use files
    use earthhelper
    use mixdefs
    use mixgeom

    implicit none

    type, extends(innerMagBase_T) :: empData_T ! replace "emp" (specific to equator) with a more descriptive "emp" for "empirical"

        type(ebTab_T)   :: ebTab
        logical :: doStatic = .true.
        integer :: Nt,Np
        real(rp), dimension(:,:), allocatable :: X,Y
        real(rp) :: empT1,empT2 !Times of two data slices
        integer  :: empN1,empN2 !Indices of two data slices
        real(rp), dimension(:,:,:), allocatable :: empW1,empW2
        real(rp) :: rDeep   ! where we're ingesting
        type(mixGrid_T) :: sstG   ! remix-style grid
        real(rp), dimension(:,:), allocatable :: sstP  ! pressure interpolated to sstG

        contains

        ! over-ride base functions
        procedure :: doInit => initSST
        procedure :: doAdvance=> advanceSST
        !procedure :: doEval => evalSST
        !procedure :: doIO => 
        !procedure :: doRestart => 

    end type empData_T

    contains


    !Initialize Empirical Map data
    subroutine initSST(imag,iXML,isRestart,vApp)
        class(empData_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart !Do you even care? VGM: NO
        type(voltApp_T), intent(inout) :: vApp

        character(len=strLen) :: empFile
        integer :: i,j,n1,n2
        integer :: dims(2)
        integer, parameter :: NIOVAR = 3
        type(IOVAR_T), dimension(NIOVAR) :: IOVars

        !Get name of file holding emp data
        call iXML%Set_Val(empFile,"empmap/empFile","ts07.h5")
        call CheckFileOrDie(empFile,"Error opening Empirical Map data")

        imag%ebTab%bStr = empFile        
        !Scrape info from file (don't use CHIMP time scaling)
        call rdTab(imag%ebTab,iXML,empFile,doTSclO=.false.)
        if (imag%ebTab%N>1) then
            imag%doStatic = .false.
        endif

        imag%Np = imag%ebTab%dNi
        imag%Nt = imag%ebTab%dNj
        
        imag%rDeep = vApp%rTrc  ! probably, will remain unused

        allocate(imag%X(1:imag%Np+1,1:imag%Nt+1))
        allocate(imag%Y(1:imag%Np+1,1:imag%Nt+1))

        allocate(imag%empW1(1:imag%Np,1:imag%Nt,NVARIMAG))
        allocate(imag%empW2(1:imag%Np,1:imag%Nt,NVARIMAG))

        allocate(imag%sstP(1:imag%Np,1:imag%Nt))        

    !Read grid
        call ClearIO(IOVars)
        call AddInVar(IOVars,"X")
        call AddInVar(IOVars,"Y")

    ! Note, X == Phi, Y == Theta from the sst2h5ion.py script

        call ReadVars(IOVars,.false.,empFile) !Don't use io precision
        dims = [imag%Np+1,imag%Nt+1]

        imag%X = reshape(IOVars(1)%data,dims)
        imag%Y = reshape(IOVars(2)%data,dims)

    !Initialize data
        ! this returns n1=1, n2=2 for t0<=0
        ! if imag%doStatic, it returns n1=n2=1

        ! TODO: can we just do AdvanceSST(imag,vApp,vApp%time) here? 
        call findSlc(imag%ebTab,vApp%time,n1,n2)

        call rdEmpMap(imag,n1,imag%empW1)
        imag%empN1 = n1
        imag%empT1 = imag%ebTab%times(n1)

        call rdEmpMap(imag,n2,imag%empW2)
        imag%empN2 = n2
        imag%empT2 = imag%ebTab%times(n2)

    ! init pressure, set to 1 nPa arbitrarily
        imag%sstP = 1.  ! doesn't matter

        ! start with the first time slice
        call EvalSST(imag,imag%empT1)

        ! initialize the remix-style grid to interpolate to RCM
        ! note, -1 in both dimensions 
        ! -1 in Phi cuts off the periodic point
        ! -1 in Theta removes the last line so we can interpolate from the cell-centered data
        call init_grid_fromTP(imag%sstG,imag%Y(1:imag%Np,1:imag%Nt),imag%X(1:imag%Np,1:imag%Nt),isSolverGrid=.false.) 

    end subroutine initSST

    !Update emp map state for current time
    subroutine AdvanceSST(imag,vApp,tAdv)
        class(empData_T), intent(inout) :: imag
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv

        integer :: n1,n2

        !This results in EvalSST not being called, so no time interpolation happens while T1 < tAdv < T2
        !if ( (tAdv >= imag%empT1) .and. (tAdv <= imag%empT2) ) then
            !Nothing to do here
        !    return
        !endif

        !If tAdv is out of range of loaded slices, need to update
        if (tAdv < imag%empT1 .or. tAdv > imag%empT2) then
            call findSlc(imag%ebTab,tAdv,n1,n2)
            if (imag%empN1 /= n1) then
                !Read slice
                call rdEmpMap(imag,n1,imag%empW1)
                imag%empN1 = n1
                imag%empT1 = imag%ebTab%times(n1)
            endif

            if (imag%empN2 /= n2) then
                !Read slice
                call rdEmpMap(imag,n2,imag%empW2)
                imag%empN2 = n2
                imag%empT2 = imag%ebTab%times(n2)
            endif
        endif

        !Always want to run this
        call EvalSST(imag,tAdv)

    end subroutine AdvanceSST

    ! interpolate sst values to corners for representation on the remix-style grid and interpolation to rcm

    subroutine EvalSST(empData,t)  
        type(empData_T), intent(inout) :: empData
        real(rp) :: t,w1,w2
        real(rp), dimension(:,:), allocatable :: P

        ! first interpolate in time
        call tWeights(empData,t,w1,w2)
        P = w1*empData%empW1(:,:,1) + w2*empData%empW2(:,:,1)

        ! now interpolate in space (cell-centers to corners)
        ! P is dynamically allocated above (Fortran 2003 standard) to (Np,Nt) size

        associate(Np=>empData%Np, Nt=>empData%Nt)
        empData%sstP(2:Np,2:Nt) = 0.25*( P(2:Np,1:Nt-1)+P(2:Np,2:Nt)+P(1:Np-1,1:Nt-1)+P(1:Np-1,2:Nt) )
            
        ! fix periodic
        empData%sstP(1,2:Nt) = 0.25*( P(Np,1:Nt-1)+P(Np,2:Nt)+P(1,1:Nt-1)+P(1,2:Nt) )

        ! fix pole
        empData%sstP(:,1) = empData%sstP(:,2)

        ! now sstP has the size (1:Np,1:Nt)
        end associate

    end subroutine EvalSST    

    subroutine rdEmpMap(empData,n,W)
        type(empData_T), intent(inout) :: empData
        integer, intent(in) :: n
        real(rp), dimension(:,:,:), intent(inout) :: W

        integer, parameter :: NIOVAR = 2
        type(IOVAR_T), dimension(NIOVAR) :: IOVars
        integer :: dims(2)

        write(*,*) 'Reading file/group = ', &
             trim(empData%ebtab%bStr),'/',trim(empData%ebTab%gStrs(n))

        call ClearIO(IOVars)
        call AddInVar(IOVars,"P")

        call ReadVars(IOVars,.false.,empData%ebTab%bStr,empData%ebTab%gStrs(n))

        dims = [empData%Np,empData%Nt]
        W(:,:,1) = reshape(IOVars(1)%data,dims)
    end subroutine rdEmpMap


    ! Below is the version of the doEval function intended for standalone (without rcm)  use of sst pressures
    ! I wrote this function copy/pasting stuff from rcm doEval but never tested it because it's never been used
    !
    ! !Evaluate emp map at a given point
    ! !Returns density (#/cc) and pressure (nPa)
    ! subroutine EvalSST(imag,x1,x2,t,imW,isEdible)
    !     class(empData_T), intent(inout) :: imag
    !     real(rp), intent(in) :: x1,x2,t
    !     real(rp), intent(out) :: imW(NVARIMAG)
    !     logical, intent(out) :: isEdible

    !     integer, dimension(2) :: ij0

    !     real(rp) :: D,P,x0,y0
    !     integer :: ij0(2),i0,j0
    !     real(rp) :: w1,w2

    !     associate(empData => imag, lat=>x1, lon=>x2)

    !     !Set defaults
    !     imW = 0.0
    !     isEdible = .false.

    !     colat = PI/2 - lat

    !     ! Do 1st short cut tests
    !     ! see if we're inside SST grid, which is always the case since it covers the whole hemisphere
    !     ! just keeping it here for consistency with RCM
    !     isEdible =  (colat >= empData%Y(1,1)) .and. (colat <= empData%Y(1,empData%Nt)) &
    !                 .and. (lat > TINY) .and. (t>0)

    !     if (.not. isEdible) return

    !     !If still here, find mapping (i,j) on RCM grid of point
    !     call GetSSTLoc(lat,lon,ij0)

    !     if (isEdible) then
    !         call tWeights(empData,t,w1,w2)
    !         !D = psphD(r) !Gallagher plasmasphere
    !         D = GallagherRP(r,phi)   ! FIXME REPLACE Gallagher with the lon-lat version

    !         P = w1*empData%empW1(i0,j0,1) + w2*empData%empW2(i0,j0,1)
    !     else
    !         D = 0.0
    !         P = 0.0
    !     endif

    !     imW(IMDEN)  = D
    !     imW(IMPR)   = P
    !     imW(IMTSCL) = 0.0 !Rely on coupling timescale
    !     imW(IMX1)   = (180.0/PI)*lat
    !     imW(IMX2)   = (180.0/PI)*lon
        
    !     end associate

    !     contains

    !     subroutine GetSSTLoc(lat,lon,ij0)
    !         real(rp), intent(in) :: lat,lon
    !         integer, intent(out) :: ij0(2)

    !         integer :: iX,jX,iC
    !         real(rp) :: colat,dp,dcol,dI,dJ

    !         associate(gcolat=>empData%Y(1,:),glong=>empData%Y(:,1), &
    !                   nLat=>empData%Nt,nLon=>empData%Np)

    !         !Assuming constant lon spacing
    !         dp = glong(2) - glong(1)

    !         !Get colat point
    !         colat = PI/2 - lat
    !         iC = findloc(gcolat >= colat,.true.,dim=1) - 1
    !         dcol = gcolat(iC+1)-gcolat(iC)
    !         dI = (colat-gcolat(iC))/dcol
    !         if (dI <= 0.5) then
    !             iX = iC
    !         else
    !             iX = iC+1
    !         endif

    !         !Get lon point
    !         dJ = lon/dp
    !         if ( (dJ-floor(dJ)) <= 0.5 ) then
    !             jX = floor(dJ)+1
    !         else
    !             jX = floor(dJ)+2
    !         endif

    !         !Impose bounds just in case
    !         iX = max(iX,1)
    !         iX = min(iX,nLat)
    !         jX = max(jX,1)
    !         jX = min(jX,nLon)

    !         ij0 = [jX,iX]  ! note, flipping the order compared to RCM, since SST data have (lon,lat) order

    !         end associate
    !     end subroutine GetSSTLoc

    ! end subroutine EvalSST

    !Get weights for given time
    subroutine tWeights(empdata,t,w1,w2)
        type(empData_T), intent(inout) :: empData
        real(rp), intent(in)  :: t
        real(rp), intent(out) :: w1,w2
        real(rp) :: dt
        if (empData%doStatic) then
            w1 = 1.0
            w2 = 0.0
            return
        endif

        if (t > empData%empT2) then
            w1 = 0.0
            w2 = 1.0
        else if (t < empData%empT1) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = empData%empT2-empData%empT1
            if (dt>TINY) then
                w1 = (empData%empT2-t)/dt
                w2 = (t - empData%empT1)/dt
            else
                w1 = 0.5
                w2 = 0.5
            endif
        endif !Weights
    end subroutine

end module sstLLimag
