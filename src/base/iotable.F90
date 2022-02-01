
!Table for holding input/output data slice information/decomposition/etc

    !Main routines

module iotable
    use kdefs
    use ioH5
    use xml_input

    implicit none

    integer, parameter, private :: IOTABVARS = 10

    !Data necessary to update fields, ie time->field data file mapping
    !File: bStr + gStr (base+group)
    type ioTab_T
        integer :: N !Number of time slices
        character(len=strLen) :: bStr
        character(len=strLen), dimension(:), allocatable :: gStrs
        real(rp), dimension(:), allocatable :: times
        real(rp), dimension(:), allocatable :: MJDs

        !Information for decomposed data
        logical :: isMPI = .false.,hasMJD=.false.
        integer :: Ri,Rj,Rk
        integer :: dNi,dNj,dNk
    end type ioTab_T

    contains

    !Read times for input data slices, convert times to code units
    !Figure out grid sizes
    subroutine InitIOTab(ioTab,inpXML,ioFile)
        type(ioTab_T), intent(inout)      :: ioTab
        type(XML_Input_T), intent(in)     :: inpXML
        character(len=strLen), intent(in) :: ioFile

        integer :: s0,sE,Nstp,i,Nd,dims(NDIM)
        character(len=strLen) :: gStr
        real(rp), allocatable, dimension(:) :: Ts
        type(IOVAR_T), dimension(IOTABVARS) :: inIOs


        call StepInfo(ioFile,s0,sE,Nstp)

        write(*,'(a,a,a,I0,a,I0,a,I0,a)') '<', trim(ioFile), ': Found ',Nstp,' timeslices, ', s0, ' to ', sE,'>'

        ioTab%N = Nstp
        allocate(ioTab%times(Nstp))
        allocate(ioTab%MJDs (Nstp))
        allocate(ioTab%gStrs(Nstp))
        allocate(Ts(Nstp))

        call StepTimes(ioFile,s0,sE,Ts)
        call StepMJDs (ioFile,s0,sE,ioTab%MJDs)
        if (maxval(ioTab%MJDs)>TINY) then
            ioTab%hasMJD = .true.
            write(*,*) 'Found MJD data ...'
        else
            ioTab%hasMJD = .false.
        endif

        do i=1,Nstp
            write(gStr,'(A,I0)') "Step#", s0+i-1

            ioTab%gStrs(i) = gStr
            ioTab%times(i) = Ts(i)
        enddo !stp loop
        
        !Get grid size info
        call ClearIO (inIOs)
        call AddInVar(inIOs,"X")
        call ReadVars(inIOs,.false.,ioFile) !Use IO precision
        Nd = inIOs(1)%Nr
        
        dims = inIOs(1)%dims(1:Nd)
        ioTab%dNi = dims(IDIR)-1
        ioTab%dNj = dims(JDIR)-1
        ioTab%dNk = dims(KDIR)-1

    end subroutine InitIOTab


    !Finds bounding slices from ebTab file
    !NOTE: findloc isn't supported by most gfortran versions, so this is a lazy workaround
    !When gfortran gets its act together, just comment out last few lines using findloc
    subroutine GetTabSlc(ioTab,t,i1,i2)
        type(ioTab_T), intent(in) :: ioTab
        real(rp), intent(in) :: t
        integer, intent(out) :: i1,i2

        integer :: n

        !Old code
        ! i1 = findloc(ebTab%times .le. t,.true.,dim=1,back=.true.)
        ! i2 = findloc(ebTab%times .gt. t,.true.,dim=1)
        ! i1 = max(1,i1)
        ! i2 = min(ebTab%N,i2)
    
        !Work-around code        
        do n=1,ioTab%N
            if (ioTab%times(n)>t) exit
        enddo
        i1 = n-1
        i2 = i1+1
        i1 = max(1,i1)
        
        if (i2 == i1) i2=i1+1 !Possible if none of the tab slices are in range
        i2 = min(ioTab%N,i2)
        
    end subroutine GetTabSlc
end module iotable
