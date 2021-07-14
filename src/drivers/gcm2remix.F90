  ! read new style remix files and produce a solution
  ! useful for testing different conductance models

program remix2remix
  use mixtypes
  use mixdefs
  use mixmain
  use mixio
  use gcmtypes
  use gcminterp

  use xml_input  
  use files

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
  integer :: Step 
  type(XML_Input_T) :: xmlInp  
  type(gcm_T) :: gcmT

  call readArgs(inpXML)

  ! First set input deck reader
  xmlInp = New_XML_Input(trim(inpXML),'Kaiju/REMIX',.true.)

  ! get input H5 file name
  call xmlInp%Set_Val(inH5,"remix2remix/inH5","inH5.h5")  
  ! get output H5 file name
  call xmlInp%Set_Val(outH5,"remix2remix/outH5","outH5.h5")  
  ! get time step
  call xmlInp%Set_Val(Step,"remix2remix/Step",5)
  !inH5 = "inH5.h5"
  write(*,*) trim(inH5)
  call CheckFileOrDie(inH5,"Couldn't find input h5 file. Exiting...")
  call CheckAndKill(outH5)

  call readMIX(trim(inH5),Step,mixIOobj)
  write(*,*) 'Init Mix: ',trim(inpXML)
  call init_mix(I,hmsphrs,inpXML,outH5,.false.,mixIOobj)
  call init_gcm_mix(gcmT,I)
  call mapGCM2MIX(gcmT,I)
  write(*,*) 'Run Mix'
  call run_mix(I,mixIOobj%tilt,gcm=gcmT)
  write(*,*) 'Write Mix'
  call writeMIX(I,0,mixIOobj%mjd,mixIOobj%time)

contains
  
  subroutine readArgs(inpXML)
    character(len=*),intent(inout) :: inpXML  ! input deck
    integer :: Narg
    logical :: fExist
    
    ! Input deck
    Narg = command_argument_count()

    if (.not.(Narg .eq. 1)) then
       write(*,*) 'Usage: remix2remix <input deck file>'
       stop
    else
       call get_command_argument(1,inpXML)
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
