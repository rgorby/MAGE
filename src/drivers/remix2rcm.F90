  ! standalone remix copied and simplified from omega
  ! hardcoding the current source to the test in Merkin&Lyon 2010

program MIX
  use mixtypes
  use mixdefs
  use mixmain
  use mixio
  
  implicit none

  ! Input deck
  character(len=strLen) :: inpXML

  integer,parameter :: hmsphrs(2) = [NORTH,SOUTH]
  real(rp),parameter :: tilt=0.23456
  type(mixApp_T) :: remixApp

  ! rcm stuff
  integer :: i,j
  integer :: isize = 200, jsize=100
  real(rp) :: rcmLowLat = 80., rcmHighLat = 20.  ! in degrees colatitude
  real(rp), dimension(:,:), allocatable :: rcmt, rcmp
  real(rp), dimension(:,:), allocatable :: rcmPsi
  type(mixGrid_T) :: rcmG
  type(Map_T) :: rcmMap

  call readArgs(inpXML)
  call init_mix(remixApp%ion,hmsphrs,inpXML,'mixtest',.false.)
  call fill_fac(remixApp%ion)
  call run_mix(remixApp%ion,tilt)
  call writeMIX(remixApp%ion,0,1344._rp,1555._rp)

  call init_uniform(rcmG,jsize,isize,rcmLowLat*pi/180._rp,rcmHighLat*pi/180._rp,isSolverGrid=.false.)
  call mix_set_map(remixApp%ion(NORTH)%G,rcmG,rcmMap)
  call mix_map_grids(rcmMap,remixApp%ion(NORTH)%St%Vars(:,:,POT),rcmPsi)

  write(*,*) minval(remixApp%ion(NORTH)%St%Vars(:,:,POT)),maxval(remixApp%ion(NORTH)%St%Vars(:,:,POT))  
  write(*,*) minval(rcmPsi),maxval(rcmPsi)
  open(unit=10, file="rcm.dat")
  do i=1,isize
     do j=1,jsize
        write(10,'(3f10.2)') rcmG%p(j,i),rcmG%t(j,i),rcmPsi(j,i)
     enddo
  enddo
  close(10)
  
  contains

  subroutine readArgs(inpXML)
    character(len=*),intent(inout) :: inpXML
    integer :: Narg
    logical :: fExist
    
    ! Input deck
    Narg = command_argument_count()
    if (Narg .eq. 0) then
       write(*,*) 'No input deck specified, defaulting to Input.xml'
       inpXML = "Input.xml"
    else
       call get_command_argument(1,inpXML)
    endif
    
    write(*,*) 'Reading input deck from ', trim(inpXML)
    inquire(file=inpXML,exist=fExist)
    if (.not. fExist) then
       write(*,*) 'Error opening input deck, exiting ...'
       write(*,*) ''
       stop
    endif

  end subroutine readArgs
  
  subroutine fill_fac(I)
    type(mixIon_T),dimension(:),intent(inout) :: I
    
    ! current setup
    real(rp) :: thetaMin, thetaDelta
    integer :: h,ii,jj

    thetaMin = 22.0_rp*PI/180.
    thetaDelta = 12.0_rp*PI/180.
    
    do h=1,size(I)
       do ii=1,I(h)%G%Np
          do jj=1,I(h)%G%Nt
             if(I(h)%G%t(ii,jj) .ge. thetaMin .and. I(h)%G%t(ii,jj) .le. (thetaMin+thetaDelta)) then
                I(h)%St%Vars(ii,jj,FAC) = sin(I(h)%G%t(ii,jj)) * sin(I(h)%G%p(ii,jj))
!             I(h)%St%Vars(ii,jj,FAC) = Mollify(I(h)%G%t(ii,jj)-thetaMin-0.5*thetaDelta,thetaDelta) * sin(I(h)%G%p(ii,jj))
             else
                I(h)%St%Vars(ii,jj,FAC) = 0.0_rp
             endif
          end do
       end do
    end do
  end subroutine fill_fAC
  
end program MIX
