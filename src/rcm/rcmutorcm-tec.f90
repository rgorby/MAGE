      program rcmutorcm
! attempt to add eeta channels
!      use kdefs
      USE Rcm_mod_subs, ONLY : iprec,rprec

      IMPLICIT NONE

      integer(iprec) :: idim,jdim,kdim
      integer(iprec) :: i,j,k,kk
      integer(iprec) :: itime,itime0,itimef
      integer(iprec), allocatable, dimension(:) :: ikflav
      integer(iprec) :: ks,kf

      real(rprec) :: time,start_time_min,time_min,eps,r
      integer(iprec), allocatable, dimension(:) :: imin_j
      real(rprec), allocatable, dimension(:) :: alam,etac
      real(rprec), allocatable, dimension(:,:) :: xi,yi,zi
      real(rprec), allocatable, dimension(:,:) :: xmin,ymin,zmin,rmin,pmin,beta_avg
      real(rprec), allocatable, dimension(:,:) :: vm,v,birk,bmin,pressure,pvg,dsob3,pressure_mhd,pressure_rcm,t_mhd
      real(rprec), allocatable, dimension(:,:) :: dens_mhd,densi, dense, ti, te, pressurei,pressuree,dens_corr
      real(rprec), allocatable, dimension(:,:,:) :: eeta,veff,partial_pressures,partial_densities
      integer(iprec), allocatable, dimension(:,:) :: iopen
      real(rprec), parameter :: mass_proton = 1.67e-27
      real(rprec), parameter :: mass_electron = 9.1e-31
      real(rprec), parameter :: coulomb = 1.6022e-19
      real(rprec), parameter :: radius_earth = 6380.e3
      real(rprec), parameter :: nt = 1.0e-9
      real(rprec) :: pressure_factor, density_factor

      logical :: threed,tecplot360
      logical :: onefile,opened
      logical :: firstpass =.true.

      character (len=9) :: chartime
      character (len=45) :: tecoutfile
      character (len=1) :: ans
      CHARACTER (LEN=05) :: char_kvalue='xxxxx'

      pressure_factor = 2./3.*coulomb/(radius_earth)*nt
  
      density_factor = nt/radius_earth

      write(6,*)' NOTE: this version writes out ze'

      write(*,'(a,$)')' enter the time to start, end reading: '
      read(5,*)itime0,itimef
      write(*,'(a,$)')' enter start time on min: '
      read(5,*)start_time_min

      write(*,'(a,$)')' do you want 3d data information (y/n)?: '
      read(5,'(a1)')ans
      if(ans.eq.'y'.or.ans.eq.'Y')threed=.true.

      write(*,'(a,$)')' do you want it all in one file (y/n)?: '
      onefile=.false.
      read(5,'(a1)')ans
      if(ans.eq.'y'.or.ans.eq.'Y')onefile=.true.
      opened = .false.

!      write(*,'(a,$)')' do you want it for tecplot360(y/n)?: '
!      tecplot360=.false.
!      read(5,'(a1)')ans
!      if(ans.eq.'y'.or.ans.eq.'Y')i
       tecplot360=.true.

      open(unit=1,file='rcmu_torcm.dat',status='old',action='read',&
      form='unformatted')

      itime = 0

      do
      if(itime <= itimef)then

      read(1,end=10)idim,jdim,kdim

      write(*,*)'idim =',idim,' jdim=',jdim,' kdim =',kdim

      write(*,'(a,i6)')' reading in time =',itime
! only allocate once
      if( .not.allocated(alam)) then
      allocate (alam(kdim),etac(kdim),ikflav(kdim))
      allocate (xi(idim,jdim),yi(idim,jdim),zi(idim,jdim))
      allocate (xmin(idim,jdim),ymin(idim,jdim),zmin(idim,jdim))
      allocate (vm(idim,jdim),v(idim,jdim))
      allocate (rmin(idim,jdim),pmin(idim,jdim))
      allocate (birk(idim,jdim),dsob3(idim,jdim))
      allocate (eeta(idim,jdim,kdim),bmin(idim,jdim))
      allocate (veff(idim,jdim,kdim))
      allocate (partial_pressures(idim,jdim,kdim))
      allocate (partial_densities(idim,jdim,kdim))
      allocate (pressure_mhd(idim,jdim))
      allocate (pressure_rcm(idim,jdim))
      allocate (pressure(idim,jdim))
      allocate (pressuree(idim,jdim))
      allocate (pressurei(idim,jdim))
      allocate (dens_mhd(idim,jdim))
      allocate (densi(idim,jdim))
      allocate (dense(idim,jdim))
      allocate (t_mhd(idim,jdim))
      allocate (ti(idim,jdim))
      allocate (te(idim,jdim))
      allocate (beta_avg(idim,jdim))
      allocate (pvg(idim,jdim))
      allocate (iopen(idim,jdim))
      allocate (imin_j(jdim))
      end if

      read(1) itime, alam, xi, yi, zi, rmin, pmin, &
              iopen, vm, pressure_mhd, dens_mhd, bmin, ti, te, beta_avg,v,eeta

      if(firstpass)then
      if(threed)then
      do k=1,kdim
       write(6,*)' alam(',k,')', alam(k)
      end do
              write(6,*)' enter the start and end k vals threed output'
              write(6,*)'enter 0 0 for all the data'
              read(5,*)ks,kf
              if(ks == 0 .and. kf ==0)then
                      ks = 1
                      kf = kdim
              endif
      endif
      END IF
      firstpass=.false.

! output file
      if(itime >= itime0)then
              time_min = itime/60.+start_time_min
      call min2hr(time_min,chartime)

      if(.not.opened.and.onefile)then

       if(threed)then
        tecoutfile = adjustr('tec3dtorcm') // chartime // '.dat'
              write(*,*)' writing file =',tecoutfile
       else
        tecoutfile = adjustr('tectorcm') // chartime // '.dat'
              write(*,*)' writing file =',tecoutfile
       endif

       endif

       ! compute xmin,ymin
       do j=1,jdim
        do i=1,idim
         xmin(i,j) = rmin(i,j)*cos(pmin(i,j))
         ymin(i,j) = rmin(i,j)*sin(pmin(i,j))
        end do
       end do

       if(.not.opened) open(unit=2,file=tecoutfile,status='unknown',recl=500)
! add corotation
!       do j=1,jdim
!        do i=1,idim
!         if(vm(i,j) > 0)then
!          r = sqrt(xmin(i,j)**2+ymin(i,j)**2)  
!          !v(i,j) = v(i,j) - 92400./r
!
!         do k=1,kdim
!          veff(i,j,k) = v(i,j) + alam(k)*vm(i,j) -92400./r
!         end do
!         end if
!        end do
!       end do

! now output the files

       if(threed)then
       if(.not.opened)then
       write(2,*)' VARIABLES =, "xmin(Re)" "ymin(Re)" "zmin(Re)" "xi(Re)"'&
               ,' "yi(Re)" "zi(Re)" ' &
               ,' "v(V)" "vm" "bmin(nT)" "mhd pressure(Pa)" "mhd density(cc)"'&
               ,' "pV<sup><greek>g</greek></sup>"'&
               ,' "ion density(cc)" "electron density(cc)" '&
               ,' "ion pressure(Pa)" "electron pressure(Pa)" "rcm pressure(Pa)" '&
               ,' "mhd temperature(eV)" "ion temperature(eV)" "electron temperature(eV)"'&
               ,' "eeta"  "veff" "partial_pressures(Pa)" "partial_densities(cc)"'

       end if
! 2d
       else
       if(.not.opened)then
       write(2,*)' VARIABLES =, "xmin(Re)" "ymin(Re)" "zmin(Re)" "xi(Re)" '&
               ,' "yi(Re)" "zi(Re)" ' &
               ,' "v(V)" "vm" "bmin(nT)" '&
               ,' "mhd pressure(Pa)" "mhd density(cc)"'&
               ,' "pV<sup><greek>g</greek></sup>"'&
               ,' "ion density(cc)" "electron density(cc)" "rcm density(cc)"'&
               ,' "ion pressure(Pa)" "electron pressure(Pa)" "mhd pressure(Pa)" '&
               ,' "mhd temperature(eV)" "ion temperature(eV)" "electron temperature(eV)"'&
               ,' "eeta0" "veff0" "partial_pressures(Pa)" "partial_densities(cc)"'
       end if

       end if

       pressure = 0.0
       pressuree = 0.0
       pressurei = 0.0
       partial_pressures = 0.0
       partial_densities = 0.0
       densi = 0.0
       dense = 0.0
       ti = 0.0
       te = 0.0
       opened = .true.
        do j=1,jdim
! process the data
! imin_j is the last closed fieldline        
         imin_j(j) = 2
         do i=idim,2,-1
          if(vm(i,j) < 0)then
           imin_j(j) = i + 1
          exit
          end if
         end do

         eps = 1.0e-3
         do i=imin_j(j)-1,1,-1
          xmin(i,j) = xmin(i+1,j)*(1 + eps)
          ymin(i,j) = ymin(i+1,j)*(1 + eps)
         end do
! now set the pressure
         do i=imin_j(j),idim
         do kk=1,kdim
          pressure(i,j) = pressure(i,j) + pressure_factor*abs(alam(kk))*&
                          eeta(i,j,kk)*vm(i,j)**2.5
          veff(i,j,kk) = v(i,j) + vm(i,j)*alam(kk)-92400./rmin(i,j)
! partial pressures
          partial_pressures(i,j,kk) =  pressure_factor*abs(alam(kk))*&
                          eeta(i,j,kk)*vm(i,j)**2.5
          partial_densities(i,j,kk) =  density_factor*eeta(i,j,kk)*vm(i,j)**1.5
          if(alam(kk) > 0.0)then
           densi(i,j) = densi(i,j) + density_factor*eeta(i,j,kk)*vm(i,j)**1.5
           pressurei(i,j) = pressurei(i,j) + pressure_factor*abs(alam(kk))*&
                          eeta(i,j,kk)*vm(i,j)**2.5
          else
           dense(i,j) = dense(i,j) + density_factor*eeta(i,j,kk)*vm(i,j)**1.5
           pressuree(i,j) = pressuree(i,j) + pressure_factor*abs(alam(kk))*&
                          eeta(i,j,kk)*vm(i,j)**2.5
          end if
         end do

         pressure_rcm = pressurei+pressuree

         end do
         end do
! pv^gamma
         where(vm > 0) 
          pvg = pressure*vm**(-2.5)*1.0e9 ! convert to wolf units
         end where
! ion temperature
         where(vm > 0)         
                 where(densi > 0)
                 ti = pressurei/densi/coulomb
                 end where
                 where(dense > 0)
                 te = pressuree/dense/coulomb
                 end where
         end where
! mhd temperature
        t_mhd = pressure_mhd/dens_mhd/coulomb
! reset to avoid Inf and NaN
         where (vm < 0.0)
                 pressure = 0.0
                 pressurei = 0.0
                 pressuree = 0.0
                 densi = 0.0
                 dense = 0.0
                 ti = 0.0
                 te = 0.0
                 pvg = 0.0
         end where

! output the data
        if(threed)then
        do k=ks,kf
        WRITE (char_kvalue,'(A2,I3.3)') 'K=', k
        if(tecplot360)then
        if(k==ks)then
        WRITE (2,*) 'ZONE T="RCM-'//char_kvalue//&
                    '" I=', idim, ', J=',jdim, ',DATAPACKING=BLOCK, SOLUTIONTIME=',itime
        else
        WRITE (2,*) 'ZONE T="RCM-'//char_kvalue//&
                    '" I=', idim, ', J=',jdim, ',DATAPACKING=BLOCK, VARSHARELIST=([1-20]), SOLUTIONTIME=',itime
        end if

        else
        if(k==ks)then
        WRITE (2,*) 'ZONE T="RCM-'//char_kvalue//&
                    '" I=', idim, ', J=',jdim, ',  DATAPACKING=BLOCK'
        else
        WRITE (2,*) 'ZONE T="RCM-'//char_kvalue//&
    '" I=', idim, ', J=',jdim, ',  DATAPACKING=BLOCK,  VARSHARELIST=([1-20])'
        end if

        end if
         write(2,'(A,F9.3,A)') 'AUXDATA Alamc='//'"', alam(k),'"'
         write(2,'(a,a,a)')' AUXDATA TIME ="',chartime,'"'

         if(k == ks)then
         DO j=1,jdim; write(2,'(15es14.5)')(xmin(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(ymin(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(zmin(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(xi(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(yi(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(zi(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')( v(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(vm(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(bmin(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pressure_mhd(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(dens_mhd(i,j)/1.0e6,i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pvg(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(densi(i,j)/1.0e6,i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(dense(i,j)/1.0e6,i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pressurei(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pressuree(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pressure_rcm(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(t_mhd(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(ti(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(te(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(eeta(i,j,k),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(veff(i,j,k),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')&
                        (partial_pressures(i,j,k),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')&
                        (partial_densities(i,j,k)/1.0e6,i=1,idim);END DO
        else
         DO j=1,jdim; write(2,'(15es14.5)')(eeta(i,j,k),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(veff(i,j,k),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')&
                        (partial_pressures(i,j,k),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')&
                        (partial_densities(i,j,k)/1.0e6,i=1,idim);END DO
        end if

        end do

        else
! 2d        
        WRITE (2,*) 'ZONE T="RCM-'//&
                    '" I=', idim, ', J=',jdim, ', &
                            DATAPACKING=BLOCK, SOLUTIONTIME=',itime
         write(2,'(a,a,a)')' AUXDATA TIME ="',chartime,'"'
         DO j=1,jdim; write(2,'(15es14.5)')(xmin(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(ymin(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(zmin(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(xi(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(yi(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(zi(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')( v(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(vm(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(bmin(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pressure_mhd(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(dens_mhd(i,j)/1.0e6,i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pvg(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(densi(i,j)/1.0e6,i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(dense(i,j)/1.0e6,i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pressurei(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pressuree(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(pressure_rcm(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(t_mhd(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(ti(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(te(i,j),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(eeta(i,j,1),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')(veff(i,j,1),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')&
                              (partial_pressures(i,j,1),i=1,idim);END DO
         DO j=1,jdim; write(2,'(15es14.5)')&
                        (partial_densities(i,j,1)/1.0e6,i=1,idim);END DO
        end if

        if(.not.onefile)then
                opened=.false.
                close(2)
        else
                opened=.true.
        endif

        if(.not.opened)then

        end if
      end if

      else
              exit

      end if

      end do
        if(onefile)then
         close(2)
        end if

    10 stop

      end program rcmutorcm

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
!     write(*,*)' char_time_hours =',char_time_hours
!     if(int(time_hours) < 10)then
!             char_time_hours = '00'//adjustr(char_time_hours)
!     endif
!     write(*,*)' char_time_hours =',char_time_hours
      write(char_time_minutes,'(i2.2)')int(time_minutes)
      write(char_time_seconds,'(i2.2)')int(time_seconds)
!     write(*,*)' char_time_seconds =',char_time_seconds
       
      hours = adjustr(char_time_hours) //'-'//char_time_minutes//'-' &
      //adjustl(char_time_seconds)
      write(*,'(a,g14.6,a,a)')' time (min) =',time,' h-m-s:',hours

      return
      end subroutine min2hr

