! program to convert rcmu2.dat file to tecplot (as written by torcm.f90)
! 5/05 frt
      program rcmu22tecplot
      implicit none

      integer :: idim,jdim,kdim
      integer :: i,j,k,kk
      integer :: itime,itime0,itimef
      integer, allocatable, dimension(:) :: ikflav
      integer :: ilast

      real :: time,start_time_min,time_min,eps,r,pi,dpi
      real, allocatable, dimension(:) :: alam,etac
      real, allocatable, dimension(:,:) :: xi,yi,zi,rmin,pmin
      real, allocatable, dimension(:,:) :: xe,ye,ze
      real, allocatable, dimension(:,:) :: vm,v,birk,bmin,pressure,pvg
      real, allocatable, dimension(:,:) :: dens,ti,te
      real, allocatable, dimension(:,:,:) :: eeta,veff
      integer, allocatable, dimension(:,:) :: open
      real,parameter :: boltz = 1.38e-23 ! boltzmann constant
      real,parameter :: ev = 1.6e-19 ! electron volt

      logical :: threed
      logical :: onefile,opened

      character (len=9) :: chartime
      character (len=45) :: tecoutfile
      character (len=1) :: ans

      pi = acos(-1.0)
      dpi = pi/180.

      write(*,'(a,$)')' enter the time to start, end reading: '
      read(5,*)itime0,itimef
      write(*,'(a,$)')' enter start time on min: '
      read(5,*)start_time_min

      write(*,'(a,$)')' do you want 3d data information (y/n)?: '
      threed=.false.
      read(5,'(a1)')ans
      if(ans.eq.'y'.or.ans.eq.'Y')threed=.true.
      write(*,'(a,$)')' do you want it all in one file (y/n)?: '
      onefile=.false.
      read(5,'(a1)')ans
      if(ans.eq.'y'.or.ans.eq.'Y')onefile=.true.
      opened = .false.

      open(unit=1,file='rcmu2.dat',status='old',action='read',&
      form='unformatted')

      itime = 0

      do
      if(itime <= itimef)then

      read(1,end=10)idim,jdim,kdim

      write(*,*)' idim =',idim,'jdim =',jdim,' kdim =',kdim

! only allocate once
      if( .not.allocated(alam)) then
      allocate (alam(kdim),etac(kdim),ikflav(kdim))
      allocate (xi(idim,jdim),yi(idim,jdim),zi(idim,jdim))
      allocate (rmin(idim,jdim),pmin(idim,jdim))
      allocate (xe(idim,jdim),ye(idim,jdim),ze(idim,jdim))
      allocate (vm(idim,jdim),v(idim,jdim),birk(idim,jdim))
      allocate (dens(idim,jdim),ti(idim,jdim),te(idim,jdim))
      allocate (eeta(idim,jdim,kdim),bmin(idim,jdim))
      allocate (open(idim,jdim))
      allocate (veff(idim,jdim,kdim))
      allocate (pressure(idim,jdim))
      allocate (pvg(idim,jdim))
      end if

      read(1) itime, alam, xi, yi, zi, rmin, pmin, &
              open, vm, pressure, dens, bmin, ti, te, eeta
      write(*,*)' reading in time =',itime
      if(itime > itimef)stop

!      write(*,*)'alam: ',alam
! output file
      if(itime >= itime0)then
              time_min = itime/60.+start_time_min
      call min2hr(time_min,chartime)

      if(.not.opened.and.onefile)then
      if(threed)then

       tecoutfile = adjustr('tec3d') // chartime // '.dat'
              write(*,*)' writing file =',tecoutfile

      else

       tecoutfile = adjustr('tec') // chartime // '.dat'
              write(*,*)' writing file =',tecoutfile

       endif
       endif
! compute xe,ye
       do j=1,jdim
        do i=1,idim
         xe(i,j) = rmin(i,j)*cos(pmin(i,j))
         ye(i,j) = rmin(i,j)*sin(pmin(i,j))
        end do
       end do

       if(.not.opened) open(unit=2,file=tecoutfile,status='unknown',&
                            recl=200)
! add corotation
       do j=1,jdim
        do i=1,idim
         if(vm(i,j) > 0)then
          r = sqrt(xe(i,j)**2+ye(i,j)**2)  
         end if
        end do
       end do
! now set the pressure
         do i=ilast,idim
         do kk=1,kdim
          pressure(i,j) = pressure(i,j) + 1.67e-35*abs(alam(kk))*&
                          eeta(i,j,kk)*vm(i,j)**2.5
          veff(i,j,k) = v(i,j) + vm(i,j)*alam(k)
         end do
         end do
         where(vm > 0) 
          pvg = pressure*vm**(-2.5)*1.0e9 ! converts to nPa/(re/nt)^5/3
         elsewhere
          pressure = 0.0
          pvg = 0.0
         end where
! process the data
         do j=1,jdim
          ilast = 1
          do i=idim,2,-1
           if(vm(i-1,j) < 0)then
           ilast = i
           exit
           end if
          end do
!         write(*,*)' for j= ',j,' ilast =',ilast
!         write(*,*)' vm ilast =',vm(ilast,j)

          eps = 1.0e-3
          do i=ilast-1,1,-1
           xe(i,j) = xe(i+1,j)*(1 + eps)
           ye(i,j) = ye(i+1,j)*(1 + eps)
          end do
         end do
! now output the data

       if(threed)then
       if(.not.opened)then
       write(2,*)' VARIABLES =, "xe(Re)" "ye(Re)" "xi(Re)" '&
               ,'"yi(Re)" "zi(Re)"  "pressure(Pa)"'&
               ,'"pV<sup><greek>g</greek></sup>(nPa(re/nt)^5/3)" '&
               ,'"eeta" "alam"'&
               ,'"dens(cc)" "bmin(nT)" "ti(eV)" "te(eV)" "vm"'
       end if
       write(2,*)' ZONE, T ="',itime,'" ,I=',idim,' ,J=',jdim ,' ,K=',kdim &
                 ,', DATAPACKING=POINT'
               else
       if(.not.opened)then
       write(2,*)' VARIABLES =, "xe(Re)" "ye(Re)" "xi(Re)" "yi(Re)" '&
               ,'"zi(Re)" "pressure(Pa)"'&
               ,'"pV<sup><greek>g</greek></sup>(nPa(re/nt)^5/3)" '&
               ,'"eeta" "alam" "dens(cc)" "bmin(nT)" "ti(eV)" "te(eV)"'&
               ,'"vm"'
       end if
       write(2,*)' ZONE, T ="',itime,'" ,I=',idim,' ,J=',jdim ,& 
                 ', DATAPACKING=POINT'
       end if

       write(2,*)' AUXDATA TIME ="',chartime,'"'

       opened = .true.

        if(threed)then
        do k=1,kdim
        do j=1,jdim

         do i=1,idim
          write(2,21)xe(i,j),ye(i,j),xi(i,j),yi(i,j),zi(i,j),&
                     pressure(i,j),&
                     pvg(i,j),eeta(i,j,k),alam(k),&
                     dens(i,j)/1.0e6,bmin(i,j),&
                     ti(i,j)*boltz/ev,te(i,j)*boltz/ev,vm(i,j)
        21 format(14(g12.6,1x))             
          end do
         end do
        end do
        else

        do j=1,jdim
         k = 1
         do i=1,idim
          write(2,11)xe(i,j),ye(i,j),xi(i,j),yi(i,j),zi(i,j),&
                     pressure(i,j),&
                     pvg(i,j),eeta(i,j,k),alam(k),&
                     dens(i,j)/1.0e6,bmin(i,j),&
                     ti(i,j)*boltz/ev,te(i,j)*boltz/ev,vm(i,j)
        11 format(13(g12.6,1x))             
          end do
         end do
        end if

        if(.not.onefile)then
                opened=.false.
                close(2)
        else
                opened=.true.
        endif

! system dependent, comment out         
!       if(.not.opened)then
!       CALL System ('/Applications/Tec100/bin/preplot '//tecoutfile)
!       CALL System ('rm '//tecoutfile)
!       end if
      end if

      else
              exit
      end if

      end do
        if(onefile)then
         close(2)
! system dependent, comment out         
!        CALL System ('/Applications/Tec100/bin/preplot '//tecoutfile)
!        CALL System ('rm '//tecoutfile)
        end if

    10 stop

      end program rcmu22tecplot
!-----------------------------------
      subroutine min2hr(time,hours)
! returns a string of hr:min:sec from an input in minutes
! 2/04 frt
      implicit none
      real :: time
      real :: time_hours
      real :: time_minutes
      real :: time_seconds
      character (len=*) :: hours
      character (len=3) :: char_time_hours
      character (len=2) :: char_time_minutes,char_time_seconds
! time is in minutes

      time_hours = floor(time/60.)
      time_minutes = floor(time - 60*time_hours)
      time_seconds = 60*(time -floor(time))

      write(char_time_hours,'(i3.3)')int(time_hours)
      write(*,*)' char_time_hours =',char_time_hours
!     if(int(time_hours) < 10)then
!             char_time_hours = '00'//adjustr(char_time_hours)
!     endif
!     write(*,*)' char_time_hours =',char_time_hours
      write(char_time_minutes,'(i2.2)')int(time_minutes)
      write(char_time_seconds,'(i2.2)')int(time_seconds)
      write(*,*)' char_time_seconds =',char_time_seconds
       
      hours = adjustr(char_time_hours) //':'//char_time_minutes//':' &
      //adjustl(char_time_seconds)

      return
      end subroutine min2hr

