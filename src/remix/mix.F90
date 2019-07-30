  ! standalone remix copied and simplified from omega
  ! hardcoding the current source to the test in Merkin&Lyon 2010

program MIX
  use mixtypes
  use mixdefs
  use mixmain
  use mixconductance
  use mixio
  
  
  implicit none

  ! Input deck
  character(len=strLen) :: inpXML
  integer :: Narg
  logical :: fExist

  integer,parameter :: hmsphrs(1) = [NORTH]
  integer, parameter :: MAXMIXIOVAR = 10
  
  type(mixIon_T),dimension(1) :: ion  ! just north
  type(mixConductance_T) :: conductance

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

  call init_mix(ion,hmsphrs,conductance,inpXML)
  call fill_fac(ion)
  call run_mix(ion,0._rp,conductance)
  call writeMIX('mixtest.h5',ion,hmsphrs)

  write(*,*) 'Min/Max potential',minval(ion(1)%St%Vars(:,:,POT)),maxval(ion(1)%St%Vars(:,:,POT))

contains

  subroutine fill_fac(I)
    type(mixIon_T),dimension(:),intent(inout) :: I
    
    ! current setup
    real(rp) :: thetaMin, thetaDelta
    integer :: h,ii,jj

    thetaMin = 22.0_rp*PI/180.
    thetaDelta = 12.0_rp*PI/180.
    
    do h=1,size(ion)
       do ii=1,ion(h)%G%Np
          do jj=1,ion(h)%G%Nt
             if(ion(h)%G%t(ii,jj) .ge. thetaMin .and. ion(h)%G%t(ii,jj) .le. (thetaMin+thetaDelta)) then
                ion(h)%St%Vars(ii,jj,FAC) = sin(ion(h)%G%t(ii,jj)) * sin(ion(h)%G%p(ii,jj))
!             ion(h)%St%Vars(ii,jj,FAC) = Mollify(ion(h)%G%t(ii,jj)-thetaMin-0.5*thetaDelta,thetaDelta) * sin(ion(h)%G%p(ii,jj))
             else
                ion(h)%St%Vars(ii,jj,FAC) = 0.0_rp
             endif
          end do
       end do
    end do
  end subroutine fill_fAC
  
end program MIX
