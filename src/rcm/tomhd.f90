      SUBROUTINE tomhd (RM,rec, ierr)
      USE rcm_precision
      USE Rcm_mod_subs, ONLY : isize, jsize, kcsize,jwrap, nptmax, &
                              colat, aloct, v, birk, &
                              bmin, xmin, ymin, zmin, vm, pmin, rmin, &
                              birk_avg, v_avg, eeta, eeta_avg, alamc, ikflavc,&
                              pi, read_array, label, LUN, &
                              boundary, bndloc, pressrcm,&
                              Read_grid, Read_plasma,Get_boundary, &
                              xmass, densrcm,imin_j,rcmdir
      USE constants, ONLY : mass_proton,radius_earth_m,nt,ev,pressure_factor,density_factor
      USE rice_housekeeping_module
      Use rcm_mhd_interfaces
! 
!==============================================================
! purpose:
!      To read rcm data at a specified record (rec) and place 
!      into a common block '/press2_eq/' 2d pressure information
!                           now replaced - frt 10.03
!  It also writes to files:
!    'rcmu.dat' = unformatted time-averaged rcm data for plotting
!
! inputs:
!   rec = record to read rcm data
!
!
!  program to read rcm output files and write out results in
!  this version also spits out pressure contours in the equatorial
!  plane for input to setupdg.f which setups the pressure for the
! multi channel version 02/00 frt
!
!   4/18/95             rws
!   9/18/95       frt
!   9/1/98 this version takes into account the the rcm record
!           can differ to the record
!     7/09 -restructured to use modules and allow to transfer
!           the LFM grid - frt      
!     2/19 -modified version to connect to gamera - frt
!
!
!==============================================================
!
      IMPLICIT NONE
      type(rcm_mhd_T),intent(inout) :: RM
      INTEGER(iprec), INTENT (IN) :: rec
      INTEGER(iprec), INTENT (OUT) :: ierr
!
      INTEGER(iprec) :: iunit, iunit1, iunit2, iout, ig, jg
      INTEGER(iprec) :: kkk, i, j, L, k, itime, n
      INTEGER(iprec), PARAMETER :: idim=isize, jdim=jsize, jdimr=jsize
      REAL(rprec), PARAMETER :: vmmax = 2000.0
      REAL(rprec) :: bbi, bbj, di, dj,  &
                     emp, vmp, ptemp, rhotemp,ctemp
      REAL(rprec) :: xval,yval,pval,rval
      LOGICAL :: L_flag
!
      real(rprec), parameter:: pressmin0 =1.0e-12
!
      real(rprec) :: xg,yg
      real(rprec) :: xn,yn,zn
      real(rprec) :: xtrans,ytrans,ztrans,rtrans
      real(rprec) :: rinter
      real(rprec) :: peq,denseq
      real(rprec) :: max_xmin,min_xmin,max_ymin,min_ymin
      real(rprec), PARAMETER :: rtrans_min = 2.0 ! min radius to nudge back
      real(rprec), PARAMETER :: empmin = 1.0E+10
      
      integer(iprec) ::  in_tot,ierr2,mask
      integer(iprec) :: istart,iend,iit,jstart,jend,jit,kstart,kend,kit
      integer(iprec) :: ier,ifound
      integer(iprec) :: nxp,nyp
      integer(iprec) :: bnd_max
      integer(iprec), parameter :: i_offset = 1, imin_grid = 2

      real(rprec) :: dens_plasmasphere
      logical, parameter :: use_plasmasphere = .true.
      logical, parameter :: use_ionosphere_interpolation = .true.

      save ig,jg,bbi,bbj,di,dj,iout,ierr2
      save pval,rval,ptemp,rhotemp
      save xn,yn,zn,xtrans,ytrans,ztrans
      save ifound,mask

      INTEGER (iprec) :: itime0
      REAL (rprec), ALLOCATABLE :: v0(:,:), birk0(:,:), vm0(:,:), bmin0(:,:),&
                    xmin0(:,:), ymin0(:,:), rmin0(:,:), pmin0(:,:), eeta0(:,:,:),&
                    bi0(:),bj0(:),etab0(:),v_avg0(:,:),birk_avg0(:,:),eeta_avg0(:,:,:)

      LOGICAL,PARAMETER :: avoid_boundaries = .false.
      INTEGER(iprec) :: im,ipl,jm,jpl,km,kpl

!      INCLUDE 'rcmdir.h'

      ierr = 0

      ! Logic of this thing is such that rec must be >= 2 (RCM had to run at least once):

      write(6,*)'tomhd, rec=',rec
      IF (rec < 2) STOP ' tomhd called before RCM was called, aborting...'

      ! At this point, we assume that RCM arrays are all populated
      ! (plasma, b-field, bndloc, etc.):

   IF (L_write_rcmu) then

      IF (.NOT.allocated(v0) ) then  
         ALLOCATE (v0(isize,jsize), birk0(isize,jsize), vm0(isize,jsize),&
                   bmin0(isize,jsize),xmin0(isize,jsize), ymin0(isize,jsize),&
                   rmin0(isize,jsize),pmin0(isize,jsize), &
                   birk_avg0(isize,jsize), v_avg0(isize,jsize), &
                   eeta0(isize,jsize,kcsize), eeta_avg0(isize,jsize,kcsize))
      END IF
 
      ! save current RCM arrays so that another set can be retrieved then will put them back.
       v0 = v
       birk0 = birk
       vm0   = vm
       bmin0 = bmin
       xmin0 = xmin
       ymin0 = ymin
       rmin0 = rmin
       pmin0 = pmin
       eeta0 = eeta
  
       v_avg0  = v_avg
       birk_avg0 = birk_avg
       eeta_avg0 = eeta_avg
       itime0   = label%intg(6)

       L = rec - 1

         CALL Read_array (rcmdir//'rcmv'   , L, label, ARRAY_2D = v)
         CALL Read_array (rcmdir//'rcmbirk', L, label, ARRAY_2D = birk)
         CALL Read_array (rcmdir//'rcmvm',   L, label, ARRAY_2D = vm)
         CALL Read_array (rcmdir//'rcmbmin', L, label, ARRAY_2D = bmin)
         CALL Read_array (rcmdir//'rcmxmin', L, label, ARRAY_2D = xmin)
         CALL Read_array (rcmdir//'rcmymin', L, label, ARRAY_2D = ymin)
         CALL Read_array (rcmdir//'rcmzmin', L, label, ARRAY_2D = zmin)
         rmin = SQRT (xmin**2+ymin**2)
         WHERE (xmin == 0.0 .AND. ymin == 0.0)
             pmin = 0.0
         ELSEWHERE
             pmin = ATAN2 (ymin, xmin)
         END WHERE
         WHERE (pmin < 0.0) pmin = pmin + 2.0*pi
         CALL Read_array (rcmdir//'rcmeeta', L, label, ARRAY_3D = eeta)

  ! change on 6/19/2013 (Stan Sazykin): turn off reading plasma edges info as it is not used
  !      CALL Read_array (rcmdir//'rcmbi', L, label, ARRAY_1D = bi)
  !      CALL Read_array (rcmdir//'rcmbj', L, label, ARRAY_1D = bj)
  !      CALL Read_array (rcmdir//'rcmetab', L, label, ARRAY_1D = etab)
  !      CALL Read_array (rcmdir//'rcmitrack', L, label, ARRAY_1D = itrack)
  !      CALL Read_array (rcmdir//'rcmmpoint', L, label, ARRAY_1D = mpoint)
  !      CALL Read_array (rcmdir//'rcmnpoint', L, label, ARRAY_1D = npoint)
  !  end of change 6/19/2013 

         CALL Read_array (rcmdir//'rcmvavg'   , L, label, ARRAY_2D = v_avg)
         CALL Read_array (rcmdir//'rcmbirkavg', L, label, ARRAY_2D = birk_avg)
         CALL Read_array (rcmdir//'rcmeetaavg', L, label, ARRAY_3D = eeta_avg)
     
         itime = label%intg(6)
         WRITE (*,'(A,I9.9,A,I5.5)') 'TOMHD: Read RCM, T=',itime,', REC=',L

         call write_rcmu (L,0_iprec)

         ! now read again at the next record
         L = rec
         
         CALL Read_array (rcmdir//'rcmv'   , L, label, ARRAY_2D = v)
         CALL Read_array (rcmdir//'rcmbirk', L, label, ARRAY_2D = birk)
         CALL Read_array (rcmdir//'rcmvm',   L, label, ARRAY_2D = vm)
         CALL Read_array (rcmdir//'rcmbmin', L, label, ARRAY_2D = bmin)
         CALL Read_array (rcmdir//'rcmxmin', L, label, ARRAY_2D = xmin)
         CALL Read_array (rcmdir//'rcmymin', L, label, ARRAY_2D = ymin)
         CALL Read_array (rcmdir//'rcmzmin', L, label, ARRAY_2D = zmin)
         rmin = SQRT (xmin**2+ymin**2)
         WHERE (xmin == 0.0 .AND. ymin == 0.0)
             pmin = 0.0
         ELSEWHERE
             pmin = ATAN2 (ymin, xmin)
         END WHERE
         WHERE (pmin < 0.0) pmin = pmin + 2.0*pi
         CALL Read_array (rcmdir//'rcmeeta', L, label, ARRAY_3D = eeta)

         CALL Read_array (rcmdir//'rcmvavg'   , L, label, ARRAY_2D = v_avg)
         CALL Read_array (rcmdir//'rcmbirkavg', L, label, ARRAY_2D = birk_avg)
         CALL Read_array (rcmdir//'rcmeetaavg', L, label, ARRAY_3D = eeta_avg)
     
         itime = label%intg(6)
         WRITE (*,'(A,I9.9,A,I5.5)') 'TOMHD: Read RCM, T=',itime,', REC=',L

     call write_rcmu (L,0_iprec)

   END IF

!     Compute rcm pressure and density (on the ionospheric RCM grid):
!
      DO j = 1, jdim
      DO i = 1, idim
         pressrcm (i,j) = 0.0
         densrcm (i,j)  = 0.0
         IF (vm(i,j) < 0.0) CYCLE
         DO k = 1, kcsize
            pressrcm(i,j) = pressrcm(i,j) + &
                 pressure_factor*ABS(alamc(k))*eeta_avg(i,j,k)*vm(i,j)**2.5 ! in pascals

!           normalize everything to the mass_proton, otherwise answer is below
!           floating point minimum answer and gets zero in ples/m^3
!           FIXME: This version is mass weighted, not sure why.
           if(alamc(k) >0.0)then ! only add the ion contribution
            densrcm(i,j) = densrcm(i,j) + &
               density_factor/mass_proton*xmass(ikflavc(k))*eeta_avg(i,j,k)*vm(i,j)**1.5
           end if
         END DO
         if(use_plasmasphere)then
! add a simple plasmasphere model based on carpenter 1992 or gallagher 2002 in ples/cc
!           call carpenter(rmin(i,j),dens_plasmasphere)
            call gallagher(rmin(i,j),dens_plasmasphere)
            densrcm(i,j) = densrcm(i,j) + dens_plasmasphere*1.0e6
         end if
      END DO
      END DO
 
      max_xmin = maxval(xmin)
      max_ymin = maxval(ymin)
      min_xmin = minval(xmin)
      min_ymin = minval(ymin)
 
!     now update the pressure and density in the mhd code
! use the lfm grid to transfer the rcm information       

      RM%Prcm = pressrcm(:,jwrap:jdim)
      RM%Nrcm   = densrcm(:,jwrap:jdim)

! if the locations are within 1 grid point of the boundary, then set the mask to zero


      RETURN
      END SUBROUTINE tomhd
!
!-----------------------------------------------------------
      subroutine writeu_2d(idim,jdim,imax,jmax,file,array,x,y,rec)
! routine to writeout array unformatted at some interval
! 8/08? frt
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      implicit none
      integer(iprec) :: idim,jdim
      integer(iprec) :: imax,jmax
      real(rprec) :: array(idim,jdim)
      real(rprec) :: x(idim),y(jdim)

      integer(iprec) :: rec,i
      character (LEN=*) :: file

!     if(rec.eq.1)then
      if(rec.le.2)then
       open(80,file=file,status='unknown',form='unformatted')
       write(80)idim,jdim,imax,jmax
       write(80)rec,array,x,y
       close(80)
      else
       open(80,file=file,status='unknown',form='unformatted' &
            ,position='append')
       write(80)rec,array,x,y
       close(80)
      end if

      return
      end subroutine  writeu_2d
 !---------------------------------------------------

      subroutine write_rcmu (record,offset)
! routine to write out rcm data in one unformated data file
! based on rcmu.dat  10/07 frt
      USE RCM_MOD_SUBS
      IMPLICIT NONE
! offset adds a time to itime to offset it
     INTEGER(iprec), INTENT(IN) :: record,offset
     INTEGER(iprec) :: itime

     IF (record == 1) THEN ! first time we write, record = 1
         OPEN (LUN, FILE=rcmdir//'rcmu.dat', status = 'replace', &
               FORM = 'UNFORMATTED')
         write(*,*)' Creating new rcmu.dat file'
     ELSE
         OPEN (LUN, FILE=rcmdir//'rcmu.dat', status = 'old', &
               FORM = 'UNFORMATTED', POSITION='APPEND')
     END IF

!     write(*,*)'---rcm time=',label%intg(6),' offset =',offset

     itime = label%intg(6) + offset

     write(6,*)'  label%intg(6) =', label%intg(6) ,' offset =',offset
     write(6,'(a,i3,a,i3)')'--->writing rcmu.dat at record =',record,' time =',itime

     WRITE (LUN) itime, isize, jsize, kcsize                         
     WRITE (LUN) alamc, etac, ikflavc                                
     WRITE (LUN) xmin,ymin,zmin, &
            sin(colat)*cos(aloct), &
            sin(colat)*sin(aloct),cos(colat), &
            vm,v_avg,eeta,birk_avg,bmin !,dsob3 
                                                             
     close(LUN)                                          
                                                           
     return                                                          
     end subroutine write_rcmu
 !------------------------------------------
         subroutine gallagher(L,density)
! approx based on gallagher et al, figure 1
! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 105, NO. A8, PAGES 18,819-18,833, AUGUST 1, 2000
! returns values in ples/cc
!         USE Rcm_mod_subs, ONLY : rprec
         USE rcm_precision, ONLY : rprec
         implicit none
         real(rprec),intent(in) :: L
         real(rprec),intent(out) :: density
         real(rprec), parameter :: L0 = 4.5
         real(rprec), parameter :: alpha = 10.
         real(rprec), parameter :: a1 = -0.25
         real(rprec), parameter :: b1 = 2.4
         real(rprec), parameter :: a2 = -0.5
         real(rprec), parameter :: b2 = 4.5
         real ::f,q

         f = 0.5*(1.0+tanh(alpha*(L-L0)))
         q = f*(a1*L + b1) + (1.0-f)*(a2*L + b2)
         density = 10.**q

         return
         end subroutine gallagher


