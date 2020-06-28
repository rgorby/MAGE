module mixmain
  use xml_input
  use clocks
  use mixtypes
  use mixparams
  use mixgeom
  use mixconductance
  use mixstate
  use mixsolver
  use mixio

  implicit none

  contains

    subroutine init_mix(I,hmsphrs,optFilename,RunID,isRestart,mixIOobj)
      type(mixIon_T),dimension(:),allocatable,intent(inout) :: I ! I for ionosphere (is an array of 1 or 2 elements for north and south) or it can be artibrarily many, e.g., for different solves done in loop
      integer, dimension(:), intent(in) :: hmsphrs       
      character(len=*), optional, intent(in) :: optFilename
      character(len=*),optional, intent(in) :: RunID   ! these two things are needed when we're coupled with Gamera
      logical,optional, intent(in) :: isRestart
      type(mixIO_T),optional, intent(in) :: mixIOobj      

      integer :: h

      if (.not.allocated(I)) allocate(I(size(hmsphrs)))
      
      I%rad_iono_m = RIonE*1.e+6 ! Default to RIonE. To change, overwrite directly after the init_mix call

      do h=1,size(I)
         if(present(optFilename)) then
            call initMIXParams(I(h)%P, optFilename)
         else
            call initMIXParams(I(h)%P)
         endif

         if (present(mixIOobj)) then
            call init_grid(I(h),mixIOobj)
            call init_state(I(h)%G,I(h)%St,mixIOobj%vars(:,:,:,h)) ! passing only the hemisphere variables            
         else
            call init_grid(I(h))
            call init_state(I(h)%G,I(h)%St)            
         end if
         
         call conductance_init(I(h)%conductance,I(h)%P,I(h)%G)

         I(h)%St%hemisphere = hmsphrs(h)

         ! check that hemisphere makes sense.
         if ((I(h)%St%hemisphere.ne.NORTH).and.(I(h)%St%hemisphere.ne.SOUTH)) then
            write(*,*) 'Hemisphere is set to an unallowable value: ',I(h)%St%hemisphere
            write(*,*) 'Stopping...'
            stop
         end if

         ! initialize solver for each ionosphere instance
         call init_solver(I(h)%P,I(h)%G,I(h)%S)

         ! flip the sign of cosd for south in this init_mix function
         ! -- which is the only place where it makes sense.  Could've
         ! done it in the mixgeom.F90 functions but they don't share
         ! the info of the entire mixIon_T object.
         if (I(h)%St%hemisphere.eq.SOUTH) I(h)%G%cosd = -I(h)%G%cosd
      end do

      ! initialize the mix dump file
      call initMIXIO(I,RunID,isRestart)
    end subroutine init_mix

    subroutine get_potential(I)
      type(mixIon_T),intent(inout) :: I

      !I%St%Vars(:,:,POT) = reshape(I%S%solution,[I%G%Np,I%G%Nt])*RionE**2*1.D3 ! in kV
      I%St%Vars(:,:,POT) = reshape(I%S%solution,[I%G%Np,I%G%Nt])*(I%rad_iono_m*1.e-6)**2*1.D3 ! in kV
    end subroutine get_potential

    subroutine run_mix(I,tilt,doModelOpt)
      type(mixIon_T),dimension(:),intent(inout) :: I 
      real(rp),intent(in) :: tilt
      logical, optional, intent(in) :: doModelOpt  ! allow to change on the fly whether we use conductance model

      logical :: doModel=.true.   ! always default to xml input deck unless doModelOpt is present and on
      integer :: h,NumH

      if (present(doModelOpt)) doModel = doModelOpt

      NumH = size(I)

      !Removing OMP loop due to weird runtime errors
      do h=1,NumH
        if (I(h)%St%hemisphere.eq.NORTH) then
          I(h)%St%tilt = tilt
        else
          I(h)%St%tilt = -tilt
        end if

        if (doModel) then
          I(h)%conductance%const_sigma = I(h)%P%const_sigma         
        else
          I(h)%conductance%const_sigma = .true.            
        end if

        call Tic("MIX-COND")
        call conductance_total(I(h)%conductance,I(h)%G,I(h)%St)
        call Toc("MIX-COND")
        call Tic("MIX-SOLVE")
        call run_solver(I(h)%P,I(h)%G,I(h)%St,I(h)%S)
        call Toc("MIX-SOLVE")
        call Tic("MIX-POT")
        call get_potential(I(h))
        call Toc("MIX-POT")
      end do
    end subroutine run_mix



end module mixmain
