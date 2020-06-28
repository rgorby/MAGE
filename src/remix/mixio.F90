  ! meant to replace mixio.F90 by using Kareems pipeline in base/ioH5.F90

module mixio
  use ioH5
  use mixdefs
  use mixtypes
  use files
  implicit none
  
  integer, parameter :: MAXIOVAR = 50
  type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
  character(len=strLen) :: h5File
  character(len=strLen), dimension(nVars) :: mixVarNames
  character(len=strLen), dimension(nVars) :: mixUnitNames   

contains

  ! since Fortran doesn't allow arrays of strings of variable length
  ! use this function to define these arrays
  ! note, this needs to be called before use (e.g., in writeMIX and readMIX)
  subroutine initMIXNames()
    ! NOTE: these have to be in the same order as the
    ! variable enumerator in mixdefs
    mixVarNames(POT)           = "Potential"
    mixUnitNames(POT)          = "kV"
    mixVarNames(FAC)           = "Field-aligned current"
    mixUnitNames(FAC)          = "muA/m**2"
    mixVarNames(SIGMAP)        = "Pedersen conductance"
    mixUnitNames(SIGMAP)       = "S"
    mixVarNames(SIGMAH)        = "Hall conductance"
    mixUnitNames(SIGMAH)       = "S"
    mixVarNames(SOUND_SPEED)   = "Sound speed"
    mixUnitNames(SOUND_SPEED)  = "cm/s"
    mixVarNames(DENSITY)       = "Density"    
    mixUnitNames(DENSITY)      = "g/cm^3"
    mixVarNames(AVG_ENG)       = "Average energy"
    mixUnitNames(AVG_ENG)      = "keV"
    mixVarNames(NUM_FLUX)      = "Number flux"
    mixUnitNames(NUM_FLUX)     = "1/cm^2 s"
    mixVarNames(NEUTRAL_WIND)  = "Neutral wind"
    mixUnitNames(NEUTRAL_WIND) = "cm/s"
    mixVarNames(EFIELD)        = "Electric field"
    mixUnitNames(EFIELD)       = "mV/m" 
  end subroutine initMIXNames

  subroutine initMIXIO(I,RunID,isRestart)
    type(mixIon_T),dimension(:),intent(in) :: I
    character(len=*),optional, intent(in) :: RunID   ! these two things are needed when we're coupled with Gamera
    logical,optional, intent(in) :: isRestart
    
    logical :: fExist
    real(rp), dimension(:,:),allocatable :: xc,yc
    integer :: n0

    !Setup remix output file
    h5File = trim(RunID) // ".mix.h5"

    inquire(file=h5File,exist=fExist)

    ! in case it's a restart, we either already have the file, in which case we don't do anything
    ! or we need to create it because the previous one got moved
    if ( (.not. isRestart).or.(isRestart.and.(.not.fExist)) ) then
       call CheckAndKill(h5File)

       !Reset IO chain
       call ClearIO(IOVars)

       call genOutGrid(I(NORTH)%G%x,I(NORTH)%G%y,xc,yc)
       
       ! save grid only for north
       call AddOutVar(IOVars,"X",xc,uStr="Ri")
       call AddOutVar(IOVars,"Y",yc,uStr="Ri")

       call AddOutVar(IOVars,"UnitsID","ReMIX")
                
       !Write out the chain (to root)
       call WriteVars(IOVars,.true.,h5File)
    endif

    ! Finally, set var and unit names
    ! for use by all subsequent functions
    ! NOTE: this assumes that initMIXIO ALWAYS gets called in the beginning
    ! e.g., as part of MIX initialization
    call initMIXNames()

  end subroutine initMIXIO

  subroutine writeMIX(I,Step,mjd,time)
    type(mixIon_T),dimension(:),intent(in) :: I    
    integer, intent(in) :: Step
    real(rp), optional, intent(in) :: time, mjd
    character(len=strLen) :: vStr

    integer :: v,h,n0
    character(len=strLen) :: gStr,uStr,hStr
    logical :: doDump = .true.
    real(rp) :: cpcp = 0.0

    !Reset IO chain
    call ClearIO(IOVars)

    ! I don't know how to handle the IO when the I object has more than one hemisphere
    ! because they're labeled by appending the hemisphere stirng to the data set name
    ! so kill it for now and worry about it when a case like this actually appears
    if (size(I)>2) then
       write(*,*) "writeMIX: Wrong hemisphere identifier. Stopping..."
       stop
    end if

    do h=1,size(I)
       ! hemisphere should be set up properly by now but still check just in case
       if (I(h)%St%hemisphere.eq.NORTH) then
          hStr = "NORTH"          
       else if (I(h)%St%hemisphere.eq.SOUTH) then
          hStr = "SOUTH"          
       else
          write(*,*) "writeMIX: Wrong hemisphere identifier. Stopping..."
          stop
       end if
       
       do v=1,nVars
          select case (v)
          case (POT)
             doDump = .true.
          case (FAC)
             doDump = .true.             
          case (SIGMAP)
             doDump = .true.             
          case (SIGMAH)
             doDump = .true.             
          case (SOUND_SPEED)
             doDump = .true.             
          case (DENSITY)
             doDump = .true.             
          case (AVG_ENG)
             doDump = .true.             
          case (NUM_FLUX)
             doDump = .true.             
          case (NEUTRAL_WIND) 
             doDump = .false.
          case (EFIELD)
             ! we never compute it
             doDump = .false.
          case DEFAULT
             doDump = .false. ! only dump the variables explicitely set to be dumped above
          end select

          ! NOTE: assuming initMIXNames got called before
          vStr = trim(mixVarNames(v)) // " "//trim(hStr)
          uStr = trim(mixUnitNames(v))
          
          if (doDump) then
             call AddOutVar(IOVars,vStr,I(h)%St%Vars(:,2:,v))
             ! inelegantly specifying the units       
             n0 = FindIO(IOVars,vStr)
             IOVars(n0)%unitStr = uStr
          endif
       enddo
    enddo

    ! now add time
    if (present(time)) call AddOutVar(IOVars,"time",time)
    if (present(mjd))  call AddOutVar(IOVars,"MJD",mjd)
    ! also add tilt
    call AddOutVar(IOVars,"tilt",I(NORTH)%St%tilt)

    ! add cpcp
    call AddOutVar(IOVars,"nCPCP",maxval(I(NORTH)%St%Vars(:,:,POT))-minval(I(NORTH)%St%Vars(:,:,POT)))
    call AddOutVar(IOVars,"sCPCP",maxval(I(SOUTH)%St%Vars(:,:,POT))-minval(I(SOUTH)%St%Vars(:,:,POT)))    
    
    !Write out the chain (to root)
    write(gStr,'(A,I0)') "Step#", Step
    call WriteVars(IOVars,.true.,h5File,gStr)
  end subroutine writeMIX

  subroutine readMIX(inH5,Step,mixIOobj)
    character(len=*), intent(in) :: inH5
    integer, intent(in) :: Step
    type(mixIO_T), intent(inout) :: mixIOobj
    
    integer,parameter :: hmsphrs(2) = [NORTH,SOUTH]
    real(rp),dimension(:,:),allocatable :: xc,yc
    integer :: h, v, n0
    integer :: dims(2)    
    
    character(len=strLen) :: gStr,hStr,uStr,vStr

    ! filling in var and unit names
    call initMIXNames()    

    call CheckFileOrDie(inH5,"Input H5 file does not exist.")

    !Reset IO chain
    call ClearIO(IOVars)
    
    ! read grid corners from root
    call AddInVar(IOVars,"X")
    call AddInVar(IOVars,"Y")
    call ReadVars(IOVars,.true.,inH5)

    dims = IOVars(1)%dims(1:2)

    if (.not.allocated(xc)) allocate(xc(dims(1),dims(2)))
    if (.not.allocated(yc)) allocate(yc(dims(1),dims(2)))
    xc = reshape(IOVars(1)%data,dims)
    yc = reshape(IOVars(2)%data,dims)

    ! convert to original mix grid
    ! and fill in mixIOobj%x,y
    call genInGrid(xc,yc,mixIOobj%x,mixIOobj%y)

    ! now read from step

    !Reset IO chain
    call ClearIO(IOVars)
    
    call AddInVar(IOVars,"time")
    call AddInVar(IOVars,"MJD")
    call AddInVar(IOVars,"tilt")    

    if ( (.not.(size(hmsphrs).eq.2)) ) then
       write(*,*) 'Code is only implemented to do two hemispheres (north,south)'
       stop
    end if 

    do h=1,size(hmsphrs)
       if (h.eq.NORTH) then
          hStr = "NORTH"          
       else if (h.eq.SOUTH) then
          hStr = "SOUTH"          
       else
          write(*,*) "readMIX: Wrong hemisphere identifier. Stopping..."
          stop
       end if

       do v=1,nVars
          vStr = trim(mixVarNames(v)) // " "// trim(hStr)
          call AddInVar(IOVars,vStr)
       end do
    end do 
      
    write(gStr,'(A,I0)') "Step#", Step
    call ReadVars(IOVars,.true.,inH5,trim(gSTr)) ! note, this checks if step exists

    ! finally fill in the mixIO object for passing to calling program
    ! (mixIOobj%x,y already filled in above by genInGrid

    ! time & mjd
    mixIOobj%time = IOVars(1)%data(1)
    mixIOobj%mjd  = IOVars(2)%data(1)

    ! allow for no tilt in the restart file
    ! for backward compatibility
    if (IOVars(3)%isDone) then
       mixIOobj%tilt  = IOVars(3)%data(1)
    else
       mixIOobj%tilt  = 0
    end if

    ! allocate as necessary
    ! NOTE: we're assuming readMIX is not called in a loop
    ! and is just used for 1-step calculation
    ! thus allocating here
    ! also, remember dims below was taken from the x,y arrays
    ! those have one extra point in phi but THE SAME size in theta as our target array here
    ! (although the actual stored arrays had -1 point in the theta direction as well)
    ! since we cut out the pole but extrapolated the low lat boundary
    if (.not.allocated(mixIOobj%Vars)) allocate(mixIOobj%Vars(dims(1)-1,dims(2),nVars,size(hmsphrs)))
    mixIOobj%Vars = 0 ! and initialize to zero

    do h=1,size(hmsphrs)
       if (h.eq.NORTH) then
          hStr = "NORTH"          
       else if (h.eq.SOUTH) then
          hStr = "SOUTH"          
       end if

       do v=1,nVars
          vStr = trim(mixVarNames(v)) // " "// trim(hStr)
          n0 = FindIO(IOVars,vStr)

          ! check whether the variable exists in the file
          ! since we didn't necessarily dump all of them in writeMIX
          if (IOVars(n0)%isDone) then
             mixIOobj%Vars(:,2:dims(2),v,h) = reshape(IOVars(n0)%data,dims-1)
             ! fix pole
             mixIOobj%Vars(:,1,v,h) = sum(mixIOobj%Vars(:,2,v,h))/size(mixIOobj%Vars(:,2,v,h))
          end if
       end do
    end do

  end subroutine readMIX

  subroutine genOutGrid(x,y,xc,yc)  
    real(rp), dimension(:,:),intent(in) :: x,y
    real(rp), dimension(:,:),allocatable,intent(out) :: xc,yc ! with corners 1/2-cell shifted from original
    real(rp), dimension(:,:),allocatable :: xtmp,ytmp ! temporary arrays
    integer, dimension(2) :: dims
    integer :: Np, Nt

    dims = shape(x)
    Np = dims(1); Nt = dims(2)

    if (.not.allocated(xc)) allocate(xc(Np+1,Nt))
    if (.not.allocated(yc)) allocate(yc(Np+1,Nt))
    if (.not.allocated(xtmp)) allocate(xtmp(Np+2,Nt))
    if (.not.allocated(ytmp)) allocate(ytmp(Np+2,Nt))

    ! first fix up periodic
    ! overlap two lines and roll back one
    xtmp(2:Np+1,:) = x
    xtmp(1,:)    = x(Np,:)
    xtmp(Np+2,:) = x(1,:)
    ytmp(2:Np+1,:) = y
    ytmp(1,:)    = y(Np,:)
    ytmp(Np+2,:) = y(1,:)

    ! move to cell centers
    xc(:,1:Nt-1) = 0.25*(xtmp(2:Np+2,2:Nt)+xtmp(1:Np+1,2:Nt)+xtmp(2:Np+2,1:Nt-1)+xtmp(1:Np+1,1:Nt-1))
    yc(:,1:Nt-1) = 0.25*(ytmp(2:Np+2,2:Nt)+ytmp(1:Np+1,2:Nt)+ytmp(2:Np+2,1:Nt-1)+ytmp(1:Np+1,1:Nt-1))

    ! finally extrapolate below equatorward boundary
    xc(:,Nt) = 2*xc(:,Nt-1)-xc(:,Nt-2)
    yc(:,Nt) = 2*yc(:,Nt-1)-yc(:,Nt-2)      
  end subroutine genOutGrid

  ! mirror of genOutGrid
  ! take staggered grid from what's stored in the H5 file
  ! and turn it into the original remix grid
  subroutine genInGrid(xc,yc,x,y)
    real(rp), dimension(:,:),intent(in) :: xc,yc ! with corners 1/2-cell shifted from original    
    real(rp), dimension(:,:),allocatable,intent(out) :: x,y

    real(rp), dimension(:,:),allocatable :: xtmp,ytmp ! temporary arrays
    integer, dimension(2) :: dims
    integer :: Np, Nt

    dims = shape(xc)
    Np = dims(1)-1; Nt = dims(2)

    if (.not.allocated(x)) allocate(x(Np,Nt))
    if (.not.allocated(y)) allocate(y(Np,Nt))

    x(:,2:Nt) = 0.25*(xc(2:Np+1,2:Nt)+xc(1:Np,2:Nt)+xc(2:Np+1,1:Nt-1)+xc(1:Np,1:Nt-1))
    y(:,2:Nt) = 0.25*(yc(2:Np+1,2:Nt)+yc(1:Np,2:Nt)+yc(2:Np+1,1:Nt-1)+yc(1:Np,1:Nt-1))

    ! fix pole
    x(:,1) = 0
    y(:,1) = 0
  end subroutine genInGrid


  
end module mixio
