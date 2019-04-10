! routines for pressure ingestion

module ingestpress
  use types
  use gamutils
  use ioH5

  implicit none

  character(len=strLen) :: fname="ts07d.h5" ! get from xml eveuntually
  real(rp),parameter :: zeq=0.05  ! distance within which we consider the point to be at equator
  real(rp),parameter :: rscale = 3/1.6  ! scale TS07d radius to move away from MHD boundary. This is a temporary fudge to make things work for now (05/16/2018)

  contains

  ! get equatorial pressure from whatever source (right now assuming empirical)
  subroutine getEqPress(Model,Grid,State)
    type(Model_T), intent(in) :: Model
    type(Grid_T), intent(in) :: Grid
    type(State_T), intent(inout) :: State
    
    integer :: i,j,k

    real(rp), dimension(:,:,:),allocatable :: eqp
    real(rp) :: pout

    call readTS07d(fname,eqp)
    
    !$OMP PARALLEL DO default (shared) collapse(2) &
    !$OMP private(i,j,k)
     do k=Grid%ks, Grid%ke
       do j=Grid%js, Grid%je
          do i=Grid%is, Grid%ie
             ! check that we're at equator
             ! we might not be, i.e.:
             ! 1. FLs traced to outer boundary of the tracing domain
             ! 2. FLs traced from seed points outside of tracing domain (assigned position [-999.]*3 by chimp)
             if (abs(State%eqMap(i,j,k,3)).le.zeq) then 

                ! this converts the empirical pressure [nPa] to code units
                State%eqPres(i,j,k) = eqPLoc(State%eqMap(i,j,k,1:2),eqp)/1.67e-2  ! TODO: generalize unit conversion

                ! add plasmaspheric density (always on if heating)
                if (Model%doPsphere) then
                   State%eqDen(i,j,k) = gallagherD(norm2(State%eqMap(i,j,k,1:2)))    ! gallagher returns /cc which is our standard units
                end if
             end if
          end do
       end do
    end do

  contains 
    subroutine readTS07d(h5fname,eqp)
      ! reading in empirical pressure from h5 file
      ! similarly to wsa.F90 interface
      character(len=*), intent(in) :: h5fname
      real(rp), dimension(:,:,:),allocatable,intent(inout) :: eqp
      logical :: fExist
      integer :: i,nvar,dims(2)
      integer, parameter :: NIOVAR = 3
      type(IOVAR_T), dimension(NIOVAR) :: IOVars

      !Reset IO chain
      call ClearIO(IOVars)

      inquire(file=h5fname,exist=fExist)
      if (.not. fExist) then
         !Error out and leave
         write(*,*) 'Unable to open TS07d pressure file, exiting'
         stop
      endif

      !Setup input chain
      call AddInVar(IOVars,"x")
      call AddInVar(IOVars,"y")
      call AddInVar(IOVars,"p")

      call ReadVars(IOVars,.false.,h5fname) !Don't use io precision
      dims=IOVars(1)%dims(1:2)   ! NOTE order: (nt,nr)
      if (.not.allocated(eqp)) allocate(eqp(dims(1),dims(2),NIOVAR))

      do i=1,NIOVAR 
         eqp(:,:,i) = reshape(IOvars(i)%data,dims)
      end do

      !scale (temporary!!! -- set to 1. at top eventually)
      eqp(:,:,1:2) = rscale*eqp(:,:,1:2)

      ! make sure the pressure is not negative
      eqp(:,:,3) = max(TINY,eqp(:,:,3))
    end subroutine readTS07d

    ! find equatorial location in eqp closest to xin (which is a 2D (x,y) vector)
    ! and return the corresponding pressure
    function eqPLoc(xin,eqp) result(pout)
      real(rp), intent(in) :: xin(2)
      real(rp), dimension(:,:,:), intent(in) :: eqp

      real(rp) :: pout
      integer :: ind(2)

      ind = minloc(sqrt((eqp(:,:,1)-xin(1))**2+(eqp(:,:,2)-xin(2))**2))
      pout = eqp(ind(1),ind(2),3)
    end function eqPLoc

    ! VGM stole from LTR-para/RCM/src/tomhd.F90
    ! originally written by Frank Toffoletto for LFM-RCM coupling
    function gallagherD(L) result(density)
      ! approx based on gallagher et al, figure 1
      ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 105, NO. A8, PAGES 18,819-18,833, AUGUST 1, 2000
      ! returns values in ples/cc
      
      implicit none
      real(rp),intent(in) :: L
      real(rp)  :: density
      real(rp), parameter :: L0 = 4.5
      real(rp), parameter :: alpha = 10.
      real(rp), parameter :: a1 = -0.25
      real(rp), parameter :: b1 = 2.4
      real(rp), parameter :: a2 = -0.5
      real(rp), parameter :: b2 = 4.5
      real ::f,q
      
      f = 0.5*(1.0+tanh(alpha*(L-L0)))
      q = f*(a1*L + b1) + (1.0-f)*(a2*L + b2)
      density = 10.**q
    end function gallagherD


    ! subroutine readTS07d(fname)
    !   character(len=*),intent(in) :: fname

    !   logical :: fExist
    !   character(len=strLen) :: header
    !   integer :: nr,nt,i,k
    !   real(rp),dimension(:,:),allocatable :: xeq,yeq,peq

    !   inquire(file=trim(fname),exist=fExist)
    !   if (.not. fExist) then
    !      write(*,*) 'Error opening TS07d pressure file, exiting ...'
    !      write(*,*) ''
    !      stop
    !   endif

    !   open(unit=111,file=trim(fname),status='old')
    !   ! red header
    !   read(111,"(A)") header
    !   call getNrNt(header,nr,nt)

    !   if (.not.allocated(xeq)) allocate(xeq(nr,nt))
    !   if (.not.allocated(yeq)) allocate(yeq(nr,nt))
    !   if (.not.allocated(peq)) allocate(peq(nr,nt))

    !   do k=1,nt
    !      do i=1,nr
    !         read(111,"(E9.3,E9.3,E9.3)") xeq(i,k),yeq(i,k),peq(i,k)
    !         print *,xeq(i,k),yeq(i,k),peq(i,k)
    !      end do
    !   end do

    !   close(111)

    ! end subroutine readTS07d

    ! subroutine getNrNt(header,nr,nt)
    !   character(len=*),intent(inout) :: header
    !   integer,intent(out) :: nr,nt
    !   integer :: is,ie   ! first and last index of '=' in the header
    !   character(len=strLen) :: nrstr,ntstr

    !   header = trim(header)
    !   ! find Nr and Nt
    !   is = scan(header,"=")
    !   nrstr= header(is+1:)
    !   ie = scan(nrstr,",")
    !   nrstr= nrstr(2:ie-1)

    !   ie = scan(header,"=",.true.)
    !   ntstr= header(ie+1:)

    !   read(nrstr,"(I)") nr
    !   read(ntstr,"(I)") nt

    ! end subroutine getNrNt
  end subroutine getEqPress

end module ingestpress
