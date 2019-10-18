  ! read new style remix files and produce a solution
  ! useful for testing different conductance models

program remix2remix
  use mixtypes
  use mixdefs
  use mixmain
  use mixio

  use xml_input  

  implicit none

  ! Input deck
  character(len=strLen) :: inpXML,inH5,outH5
  type(mixIO_T) :: mixIOobj
  type(mixIon_T),dimension(:),allocatable :: I
  ! assume the file has both hemispheres
  ! if not change it here and recompile
  ! could put in logic to figure it out from data set names in inH5
  ! but feels like too much trouble for a fairly specific application
  integer,parameter :: hmsphrs(2) = [NORTH,SOUTH]
  integer :: Step = 0  ! read from xml eventually
  

  call readArgs(inpXML,inH5,outH5)
  call readMIX(trim(inH5),Step,mixIOobj)
  call init_mix(I,hmsphrs,inpXML,outH5,.false.,mixIOobj)
  call run_mix(I,mixIOobj%tilt)
  call writeMIX(I,0,mixIOobj%mjd,mixIOobj%time)

contains
  
  subroutine readArgs(inpXML,inH5,outH5)
    character(len=*),intent(inout) :: inpXML  ! input deck
    character(len=*),intent(out) :: outH5   ! where to write (this
                                              ! is not the file from
                                              ! which we're restarting

    character(len=*),intent(out) :: inH5      ! construct this based on runid from xml file

    integer :: Narg
    logical :: fExist
    
    ! Input deck
    Narg = command_argument_count()

    if (.not.(Narg .eq. 3)) then
       write(*,*) 'Usage: remix2remix <input deck file> <input hdf file> <output hdf file>'
       stop
    else
       call get_command_argument(1,inpXML)
       call get_command_argument(2,inH5)
       call get_command_argument(3,outH5)       
    end if
    
    write(*,*) 'Reading input deck from ', trim(inpXML)
    inquire(file=inpXML,exist=fExist)
    if (.not. fExist) then
       write(*,*) 'Error opening input deck, exiting ...'
       write(*,*) ''
       stop
    endif

  end subroutine readArgs
  
end program remix2remix
