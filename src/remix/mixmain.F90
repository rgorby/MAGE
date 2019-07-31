module mixmain
  use xml_input

  use mixtypes
  use mixparams
  use mixgeom
  use mixconductance
  use mixstate
  use mixsolver



  implicit none

  contains

    subroutine init_mix(I,hmsphrs,optFilename)
      type(mixIon_T),dimension(:),intent(inout) :: I ! I for ionosphere (is an array of 1 or 2 elements for north and south) or it can be artibrarily many, e.g., for different solves done in loop
      integer, dimension(:), intent(in) :: hmsphrs ! array of integers marking hemispheres for the I object array.
      character(len=*), optional, intent(in) :: optFilename

      integer :: h ! h for hemisphere

      ! check that hmsphs and I have the same size
      if (size(I).ne.size(hmsphrs)) then
         write(*,*) "The sizes of the mixIon_T array and hemispheres array are different. Stopping..."
         stop
      end if

      do h=1,size(I)
         if(present(optFilename)) then
            call initMIXParams(I(h)%P, optFilename)
         else
            call initMIXParams(I(h)%P)
         endif
         ! FIXME: replace with a function pointer allowing an arbitrary grid specification, e.g., init_grid=>init_uniform
         call init_uniform(I(h)%G,I(h)%P%Np,I(h)%P%Nt,I(h)%P%LowLatBoundary*pi/180._rp,.true.)
         call init_state(I(h)%G,I(h)%St) 
         call conductance_init(I(h)%conductance,I(h)%P,I(h)%G)

         ! copy hemisphere and tilt parameters from the xml file only
         ! if hmsphrs array sets them to -1. This allows flexibility
         ! in allowing the possibility to define the hemisphere in the
         ! xml file for some applications, in which case they should
         ! initialize the hmsphrs array to -1.
         I(h)%St%hemisphere = hmsphrs(h)
         if (I(h)%St%hemisphere.eq.-1) I(h)%St%hemisphere = I(h)%P%hemisphere

         ! check that hemisphere makes sense by now. note, mixparams
         ! has already checked that what's read from the xml file is
         ! either NORTH or SOUTH. So the check here is only against
         ! hmsphrs array not set correctly (to -1 or NORTH or SOUTH)
         ! in the calling driver program.
         if ((I(h)%St%hemisphere.ne.NORTH).and.(I(h)%St%hemisphere.ne.SOUTH)) then
            write(*,*) 'Hemisphere is set to an unallowable value: ',I(h)%St%hemisphere
            write(*,*) 'Stopping...'
            stop
         end if

         ! tilt will be used from the param file or should be set
         ! in calling the run_mix function.
         if (I(h)%St%hemisphere.eq.NORTH) then
            I(h)%St%tilt = I(h)%P%tilt
         else 
            I(h)%St%tilt = -I(h)%P%tilt
         end if

         ! initialize solver for each ionosphere instance
         call init_solver(I(h)%P,I(h)%G,I(h)%S)

         ! flip the sign of cosd for south in this init_mix function
         ! -- which is the only place where it makes sense.  Could've
         ! done it in the mixgeom.F90 functions but they don't share
         ! the info of the entire mixIon_T object.
         if (I(h)%St%hemisphere.eq.SOUTH) I(h)%G%cosd = -I(h)%G%cosd
      end do
    end subroutine init_mix

    subroutine get_potential(I)
      type(mixIon_T),intent(inout) :: I

      I%St%Vars(:,:,POT) = reshape(I%S%solution,[I%G%Np,I%G%Nt])*RionE**2*1.D3 ! in kV
    end subroutine get_potential

    subroutine run_mix(I,tilt)
      type(mixIon_T),dimension(:),intent(inout) :: I 
      real(rp),intent(in) :: tilt

      integer :: h

      do h=1,size(I)
         ! note, if tilt is set to -9999. in the calling function the
         ! value from the param.xml file will be used.  if it is not
         ! set in the xml file, the xml reader will have set it to default
         ! value by now.
         if (tilt.ne.-9999._rp) then
            if (I(h)%St%hemisphere.eq.NORTH) then
               I(h)%St%tilt = tilt
            else
               I(h)%St%tilt = -tilt
            end if
         end if

         call conductance_total(I(h)%conductance,I(h)%G,I(h)%St)
         call run_solver(I(h)%P,I(h)%G,I(h)%St,I(h)%S)
         call get_potential(I(h))
      end do
    end subroutine run_mix

end module mixmain
