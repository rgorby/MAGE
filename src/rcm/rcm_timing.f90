 module rcm_timing_module

!  Use  Rcm_mod_subs, only : iprec,rprec,rcmdir 
  use rcm_precision

  integer(iprec), allocatable,dimension(:),save :: rcm_timing

  contains

      subroutine AddToList(element,list)
          IMPLICIT NONE

          integer(iprec) :: i, idim
          integer(iprec), intent(in) :: element
          integer(iprec), dimension(:), allocatable, intent(inout) :: list
          integer(iprec), dimension(:), allocatable :: clist


          if(allocated(list)) then
              idim = size(list)
              allocate(clist(idim+1))
              do i=1,idim          
                clist(i) = list(i)
              end do
              clist(idim+1) = element

              deallocate(list)
              call move_alloc(clist, list)

          else
              allocate(list(1))
              list(1) = element
          end if

          !K: Commenting out output
          !write(*,*)' addtolist adding t=',element
          !write(*,*)' rcm_time=',list
      return

      end subroutine AddToList

      subroutine read_rcm_timing(list)
        IMPLICIT NONE
        INCLUDE 'rcmdir.h'
         integer(iprec) :: i, idim
         integer(iprec), dimension(:), allocatable, intent(inout) :: list
         integer(iprec) :: iin,tin

         open(unit=20,file='RCMfiles/rcm_timing.dat',status='old',err=21)
         idim = 0
         do
         read(20,*,end=20)tin,iin
         idim = idim + 1
         end do
     20  close(20)

         allocate(list(idim))

         open(unit=20,file='RCMfiles/rcm_timing.dat',status='old')
         do i=1,idim
           read(20,*,end=20)tin,iin
           list(i) = tin
         end do
         close(20)

         return

     21  STOP 'read error in read_rcm_timing'

        end subroutine read_rcm_timing

        subroutine write_rcm_timing(list)
        IMPLICIT NONE
        INCLUDE 'rcmdir.h'
        integer(iprec) :: i, idim
        integer(iprec), dimension(:),intent(in) :: list

        write(*,*)' writing rcm timing '

         open(unit=20,file=rcmdir//'rcm_timing.dat',status='replace',err=21)
         idim = size(list)

         do i=1,idim
          write(20,*)list(i),i
         end do
         close(20)

         return

     21  STOP 'write error in write_rcm_timing'

         end subroutine write_rcm_timing

      subroutine find_record(itime,list,record)

       IMPLICIT NONE
       integer(iprec), dimension(:), intent(in) :: list
       integer(iprec), intent(in) :: itime
       integer(iprec), intent(out) :: record
       integer(iprec) :: idim,i

       idim = size(list)

       do i=1,idim
           if(itime == list(i))then
                record = i
                write(*,'(2(a,i5))')' At t=',itime,' record =',record
                return
        end if
        end do

        write(*,*) 'find_record: error in find record, for itime=',itime
        record = -1
        return

       end subroutine find_record


  end module rcm_timing_module


