      SUBROUTINE tomhd (RM, ierr)
      USE earthhelper, ONLY : GallagherXY
      USE rcm_precision
      USE Rcm_mod_subs, ONLY : isize, jsize, kcsize,jwrap, nptmax, &
                              colat, aloct, v, birk, &
                              bmin, xmin, ymin, zmin, vm, pmin, rmin, &
                              birk_avg, v_avg, eeta, eeta_avg, &
                              alamc, ikflavc, &
                              pi, &
                              boundary, bndloc, pressrcm, &
                              Read_grid, &
                              xmass, densrcm,denspsph,imin_j,rcmdir, &
                              eflux,eavg,ie_el
      USE constants, ONLY : mass_proton,nt,ev!,radius_earth_m,pressure_factor,density_factor
      USE rice_housekeeping_module
      Use rcm_mhd_interfaces
! 
!==============================================================
! purpose:
!  To convert RCM information (eta), to MHD information (p,n) 
!  It also writes to files (optional):
!    'rcmu.dat' = unformatted time-averaged rcm data for analysis
!
! inputs:
!
!   4/18/95             rws
!   9/18/95       frt
!   9/1/98 this version takes into account the the rcm record
!           can differ to the record
!     7/09 -restructured to use modules and allow to transfer
!           the LFM grid - frt      
!     2/19 -modified version to connect to gamera - frt
!     5/20 -removed use of record numbers in rcm bookkeeping
!          - also adds output from plasmasphere model by Shanshan Bao
!     5/20/20 - removed idim,jdim -frt
!
!
!==============================================================
!
      IMPLICIT NONE
      type(rcm_mhd_T),intent(inout) :: RM
      INTEGER(iprec), INTENT (OUT) :: ierr
!
      INTEGER(iprec) :: iunit, iunit1, iunit2, iout, ig, jg
      INTEGER(iprec) :: kkk, i, j, L, k, itime, n
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

      INTEGER (iprec) :: itime0,jp,iC
      real(rp) :: dRad,RadC,rIJ,dRxy

      LOGICAL,PARAMETER :: avoid_boundaries = .false.
      INTEGER(iprec) :: im,ipl,jm,jpl,km,kpl

      !AMS 04-22-2020
      real(rprec) :: pressure_factor,density_factor

      pressure_factor = 2./3.*ev/RM%planet_radius*nt
      density_factor = nt/RM%planet_radius

      ierr = 0

      ! At this point, we assume that RCM arrays are all populated
      ! (plasma, b-field, bndloc, etc.):

   IF (L_write_rcmu) then

         rmin = SQRT (xmin**2+ymin**2+zmin**2)
         WHERE (xmin == 0.0 .AND. ymin == 0.0)
             pmin = 0.0
         ELSEWHERE
             pmin = ATAN2 (ymin, xmin)
         END WHERE
         WHERE (pmin < 0.0) pmin = pmin + 2.0*pi

         rmin = SQRT (xmin**2+ymin**2)
         WHERE (xmin == 0.0 .AND. ymin == 0.0)
             pmin = 0.0
         ELSEWHERE
             pmin = ATAN2 (ymin, xmin)
         END WHERE
         WHERE (pmin < 0.0) pmin = pmin + 2.0*pi
         WRITE (*,'(A,I9.9,A,I5.5)') 'TOMHD: Read RCM, T=',itime,', REC=',L

     call write_rcmu (L,0_iprec)

   END IF

!     Compute rcm pressure and density (on the ionospheric RCM grid):
!
      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j,k,dens_plasmasphere)
      DO j = 1, jsize
       DO i = 1, isize
        pressrcm (i,j) = 0.0
        densrcm  (i,j) = 0.0
        denspsph (i,j) = 0.0

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
        if (use_plasmasphere) then
          if (dp_on) then 
            ! use plasmasphere channel eeta_avg(:,:,1) sbao 03/2020
            denspsph(i,j) = density_factor/mass_proton*xmass(2)*eeta_avg(i,j,1)*vm(i,j)**1.5
          else
            ! add a simple plasmasphere model based on carpenter 1992 or gallagher 2002 in ples/cc
            dens_plasmasphere = GallagherXY(xmin(i,j),ymin(i,j))
            denspsph(i,j) = dens_plasmasphere*1.0e6
          endif
        endif

       END DO
      END DO

  !    write(*,*) "rcm/tomhd.f90: eeta_avg(50,50,50)=",eeta_avg(50,50,50)
 
      max_xmin = maxval(xmin)
      max_ymin = maxval(ymin)
      min_xmin = minval(xmin)
      min_ymin = minval(ymin)
 
!     now update the pressure and density in the mhd code  

      RM%Prcm    = pressrcm(:,jwrap:jsize)
      RM%Nrcm    = densrcm (:,jwrap:jsize)
      RM%Npsph   = denspsph(:,jwrap:jsize)
      RM%flux    = eflux   (:,jwrap:jsize,:)
      RM%eng_avg = eavg    (:,jwrap:jsize,:)
      RM%fac     = birk    (:,jwrap:jsize)
      
      RM%toMHD = .false.
      !dRad = ellBdry%dRadMHD*radius_earth_m
      dRad = ellBdry%dRadMHD*RM%planet_radius

      do j=jwrap,jsize
        jp = j-jwrap+1
        iC = imin_j(j)
        RadC = norm2(RM%X_bmin(iC,jp,1:2))-dRad
        do i=iC+1,isize
          rIJ = norm2(RM%X_bmin(i,jp,1:2))
          dRxy = norm2(RM%X_bmin(i,jp,1:2)) - norm2(RM%X_bmin(i+1,jp,1:2))
          !write(*,*) 'RadC/rIJ = ',RadC/radius_earth_m,rIj/radius_earth_m
          if ( (rIJ<=RadC) .and. (dRxy>0) ) exit
        enddo
        !write(*,*) 'i/iC = ',i,iC

        RM%toMHD(i:,jp) = .true.
      enddo
      

      RETURN
      END SUBROUTINE tomhd
!
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

