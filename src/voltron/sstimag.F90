!Routines to handle SST inner magnetosphere model
module sstimag
    use volttypes
    use ebtypes
    use ebinit
    use ioh5
    use files
    use earthhelper

    implicit none

    type eqData_T
        type(ebTab_T)   :: ebTab
        logical :: doStatic = .true.
        integer :: Nr,Np
        real(rp), dimension(:,:), allocatable :: X,Y,xxc,yyc
        real(rp) :: eqT1,eqT2 !Times of two data slices
        integer  :: eqN1,eqN2 !Indices of two data slices
        real(rp), dimension(:,:,:), allocatable :: eqW1,eqW2
        

    end type eqData_T
    
    type(eqData_T), private :: eqData

    contains

    !Initialize EQ Map data
    subroutine initSST(iXML,isRestart)
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart !Do you even care?

        character(len=strLen) :: eqFile
        integer :: i,j,n1,n2
        integer :: dims(2)
        integer, parameter :: NIOVAR = 3
        type(IOVAR_T), dimension(NIOVAR) :: IOVars

        !Get name of file holding eq data
        call iXML%Set_Val(eqFile,"eqmap/eqFile","ts07.h5")
        call CheckFileOrDie(eqFile,"Error opening EQ Map data")

        eqData%ebTab%bStr = eqFile        
        !Scrape info from file (don't use CHIMP time scaling)
        call rdTab(eqData%ebTab,iXML,eqFile,doTSclO=.false.)
        if (eqData%ebTab%N>1) then
            eqData%doStatic = .false.
        endif

        eqData%Nr = eqData%ebTab%dNi
        eqData%Np = eqData%ebTab%dNj

        allocate(eqData%X(1:eqData%Nr+1,1:eqData%Np+1))
        allocate(eqData%Y(1:eqData%Nr+1,1:eqData%Np+1))

        allocate(eqData%eqW1(1:eqData%Nr,1:eqData%Np,NVARIMAG))
        allocate(eqData%eqW2(1:eqData%Nr,1:eqData%Np,NVARIMAG))

    !Read grid
        call ClearIO(IOVars)
        call AddInVar(IOVars,"X")
        call AddInVar(IOVars,"Y")

        call ReadVars(IOVars,.false.,eqFile) !Don't use io precision
        dims = [eqData%Nr+1,eqData%Np+1]

        eqData%X = reshape(IOVars(1)%data,dims)
        eqData%Y = reshape(IOVars(2)%data,dims)
    !Create helper grids
        allocate(eqData%xxc(1:eqData%Nr,1:eqData%Np))
        allocate(eqData%yyc(1:eqData%Nr,1:eqData%Np))

        !$OMP PARALLEL DO default(shared)
        do i=1,eqData%Nr
            do j=1,eqData%Np
                eqData%xxc(i,j) = 0.25*( eqData%X(i,j) + eqData%X(i,j+1) + eqData%X(i+1,j) + eqData%X(i+1,j+1) )
                eqData%yyc(i,j) = 0.25*( eqData%Y(i,j) + eqData%Y(i,j+1) + eqData%Y(i+1,j) + eqData%Y(i+1,j+1) )
            enddo
        enddo

    !Initialize data
        n1 = 1
        n2 = 2
        if (eqData%doStatic) n2 = 1

        call rdEQMap(eqData%ebTab,n1,eqData%eqW1)
        eqData%eqN1 = n1
        eqData%eqT1 = eqData%ebTab%times(n1)

        call rdEQMap(eqData%ebTab,n2,eqData%eqW2)
        eqData%eqN2 = n2
        eqData%eqT2 = eqData%ebTab%times(n2)

    end subroutine initSST

    !Update eq map state for current time
    subroutine AdvanceSST(t)
        real(rp), intent(in) :: t

        integer :: n1,n2
        if ( (t >= eqData%eqT1) .and. (t <= eqData%eqT2) ) then
            !Nothing to do here
            return
        endif

        !Otherwise we need to update
        call findSlc(eqData%ebTab,t,n1,n2)
        if (eqData%eqN1 /= n1) then
            !Read slice
            call rdEQMap(eqData%ebTab,n1,eqData%eqW1)
            eqData%eqN1 = n1
            eqData%eqT1 = eqData%ebTab%times(n1)
        endif

        if (eqData%eqN2 /= n2) then
            !Read slice
            call rdEQMap(eqData%ebTab,n2,eqData%eqW2)
            eqData%eqN2 = n2
            eqData%eqT2 = eqData%ebTab%times(n2)
        endif

    end subroutine AdvanceSST

    subroutine rdEQMap(ebTab,n,W)
        type(ebTab_T), intent(in) :: ebTab
        integer, intent(in) :: n
        real(rp), dimension(:,:,:), intent(inout) :: W

        integer, parameter :: NIOVAR = 2
        type(IOVAR_T), dimension(NIOVAR) :: IOVars
        integer :: dims(2)

        write(*,*) 'Reading file/group = ', trim(ebtab%bStr),'/',trim(ebTab%gStrs(n))
        call ClearIO(IOVars)
        call AddInVar(IOVars,"P")

        call ReadVars(IOVars,.false.,ebTab%bStr,ebTab%gStrs(n))

        dims = [eqData%Nr,eqData%Np]
        W(:,:,1) = reshape(IOVars(1)%data,dims)
    end subroutine rdEQMap

    !Evaluate eq map at a given point
    !Returns density (#/cc) and pressure (nPa)
    subroutine EvalSST(r,phi,t,imW)
        real(rp), intent(in) :: r,phi,t
        real(rp), intent(out) :: imW(NVARIMAG)

        real(rp) :: D,P,x0,y0
        integer :: ij0(2),i0,j0
        real(rp) :: w1,w2

        x0 = r*cos(phi)
        y0 = r*sin(phi)

        ij0 = minloc( (eqData%xxc-x0)**2.0 + (eqData%yyc-y0)**2.0 )
        
        D = psphD(r) !Gallagher plasmasphere

        i0 = ij0(IDIR)
        j0 = ij0(JDIR)

        call tWeights(t,w1,w2)
        P = w1*eqData%eqW1(i0,j0,1) + w2*eqData%eqW2(i0,j0,1)
        if (t<=0) then
            D = 0.0
            P = 0.0
        endif
        imW(IMDEN) = D
        imW(IMPR)  = P
        
    end subroutine EvalSST

    !Get weights for given time
    subroutine tWeights(t,w1,w2)
        real(rp), intent(in)  :: t
        real(rp), intent(out) :: w1,w2
        real(rp) :: dt
        if (eqData%doStatic) then
            w1 = 1.0
            w2 = 0.0
            return
        endif

        if (t > eqData%eqT2) then
            w1 = 1.0
            w2 = 0.0
        else if (t < eqData%eqT1) then
            w1 = 0.0
            w2 = 1.0
        else
            dt = eqData%eqT2-eqData%eqT1
            if (dt>TINY) then
                w1 = (eqData%eqT2-t)/dt
                w2 = (t - eqData%eqT1)/dt
            else
                w1 = 0.5
                w2 = 0.5
            endif
        endif !Weights
    end subroutine

end module sstimag
