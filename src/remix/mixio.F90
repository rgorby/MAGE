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

contains

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
       call AddOutVar(IOVars,"X",xc)
       call AddOutVar(IOVars,"Y",yc)

       ! inelegantly specifying the units
       n0 = FindIO(IOVars,"X")
       IOVars(n0)%unitStr = "Ri"
       n0 = FindIO(IOVars,"Y")
       IOVars(n0)%unitStr = "Ri"
          
       !Write out the chain (to root)
       call WriteVars(IOVars,.true.,h5File)
    endif

  contains
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
             vStr = "Potential"//" "//hStr
             uStr = "kV"
             doDump = .true.
          case (FAC)
             vStr = "Field-aligned current"//" "//hStr
             uStr = "muA/m**2"
             doDump = .true.             
          case (SIGMAP)
             vStr = "Pedersen conductance"//" "//hStr
             uStr = "S"
             doDump = .true.             
          case (SIGMAH)
             vStr = "Hall conductance"//" "//hStr
             uStr = "S"
             doDump = .true.             
          case (SOUND_SPEED)
             vStr = "Sound speed"//" "//hStr
             uStr = "cm/s"
             doDump = .true.             
          case (DENSITY)
             vStr = "Density"//" "//hStr
             uStr = "g/cm^3"
             doDump = .true.             
          case (AVG_ENG)
             vStr = "Average energy"//" "//hStr
             uStr = "keV"
             doDump = .true.             
          case (NUM_FLUX)
             vStr = "Number flux"//" "//hStr
             uStr = "1/cm^2 s"
             doDump = .true.             
          case (NEUTRAL_WIND) 
             vStr = "Neutral wind"//" "//hStr
             uStr = "cm/s"
             doDump = .false.
          case (EFIELD)
             ! we never compute it
             vStr = "Electric field"//" "//hStr
             uStr = "mV/m" !???
             doDump = .false.
          case DEFAULT
             doDump = .false. ! only dump the variables explicitely set to be dumped above
          end select

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

    ! add cpcp
    call AddOutVar(IOVars,"nCPCP",maxval(I(NORTH)%St%Vars(:,:,POT))-minval(I(NORTH)%St%Vars(:,:,POT)))
    call AddOutVar(IOVars,"sCPCP",maxval(I(SOUTH)%St%Vars(:,:,POT))-minval(I(SOUTH)%St%Vars(:,:,POT)))    
    
    !Write out the chain (to root)
    write(gStr,'(A,I0)') "Step#", Step
    call WriteVars(IOVars,.true.,h5File,gStr)
  end subroutine writeMIX

end module mixio
