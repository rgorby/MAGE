!==================================================================      
      SUBROUTINE Torcm (RM, itimei, ierr, icontrol) 
      USE rcm_precision
      USE Rcm_mod_subs, ONLY : isize,jsize, jwrap, kcsize, iesize, &
                               vm, bmin, xmin, ymin, pmin, rmin,v, & 
                               alamc, etac, ikflavc, fudgec, eeta, eeta_pls0, &
                               imin_j, bndloc, vbnd,               &
                               colat, aloct, bir, sini,            &
                               ibnd_type,rcmdir
     USE conversion_module
     USE rice_housekeeping_module
     USE constants, only: big_vm,tiote
     Use rcm_mhd_interfaces


! NOTE: This version fixes the rcm boundary condition at rec=1
!===================================================================
!
! purpose: to setup rcm for a run, converts mhd information to 
!          rcm information
!
! inputs:
!     rec_in   = record to be updated with new values
!               rec == 1: resets the rcm and initializes the files
!                rec > 1: only boundary conditions are updated
!    itimei = time to write to rcm records
!
! outputs:
!       ierr = error flag if there is a problem
!
!===================================================================
!
!  last update:  05.18.90
!                04.13.95
!                      96
!                20.01.00 frt mc version
!                    8.03 frt computes integral average pressure
!                06.01.04 frt cleaned up version
!                09.10.12 frt added code to handle tilt
!                19.05.20 frt removed use of records in rcm 
!                         frt remove use of idim,jdim,kdim
!
!===================================================================

      IMPLICIT NONE
      type(rcm_mhd_T),intent(inout) :: RM
      INTEGER(iprec), INTENT (IN) :: itimei,icontrol
      INTEGER(iprec), INTENT (IN OUT) :: ierr

      REAL(rprec) :: xp,yp,red
      REAL(rprec) :: eetabnd_max
      REAL(rprec), PARAMETER :: offseti = 0.0  
      REAL(rprec), PARAMETER :: ti0 = 5.0E+8
      REAL(rprec), PARAMETER :: small_eta = 1.0E+4
      REAL(rprec), PARAMETER :: xline = -25.0
      REAL(rprec), PARAMETER :: lmin = 2.10
!      REAL(rprec), PARAMETER :: big_vm = -1.0e5      
      
      INTEGER(iprec) :: kin,jmid,ibnd, inew, iold
      INTEGER(iprec) :: jm,jp,ii,itmax
      INTEGER(iprec) :: min0,i,j,k,n,ns
      INTEGER(iprec), PARAMETER :: n_smooth = 5
      INTEGER(iprec), PARAMETER :: LUN =100
      LOGICAL,PARAMETER :: use_ellipse = .true.
      LOGICAL,PARAMETER :: set_boundary_with_mach = .false.
      LOGICAL,PARAMETER :: set_boundary_with_beta = .false.
      REAL(rprec), PARAMETER :: max_beta = 1.0 ! max averaged beta to set the boundary

      ierr = 0

      IF (ierr < 0) RETURN
       x0_sm = x0; y0_sm = y0; z0_sm = z0;

!      if(rcm_tilted)then
! temporarily assign  geo grid as x0 and convert to sm
! this preserves the default behaviour when tilt is off
!       call geo2sm(isize,jsize,x0,y0,z0,x0_sm,y0_sm,z0_sm,1)
!      end if

!      write(6,*)' x0 y0 z0',x0(1,1),y0(1,1),z0(1,1)
!      write(6,*)' x0_sm y0_sm z0_sm',x0_sm(1,1),y0_sm(1,1),z0_sm(1,1)

      ! Save current RCM arrays before modifying them.
      ! T=0 or restart is a special case:

      ! If T>0, we already have old array values
      ! T=0 case means we have to defer this to after
      ! magnetic field arrays are processed.

      IF (icontrol==RCMADVANCE.or.icontrol==RCMRESTART) then  
        bndloc_old = bndloc
        imin_j_old = imin_j
        if(minval(bndloc) <= 0)then
          write(6,*) ' TORCM: boundary problem'
          write(6,*)' bndloc:',bndloc
          stop
        end if
      END IF  

      ! get B-field lines starting from RCM ionospheric grid
      ! points, compute flux-tube volume of field lines that are
      ! closed, and mark open ones with a mask array OPEN
      ! (1=open, 0-closed). For open field lines, FTV is set to big_vm:

      CALL Calc_ftv (RM,big_vm,ierr)
      IF (ierr < 0) RETURN

      ! Now, set RCM high-latitude grid boundary. Initially,
      ! for each MLT, we find the grid point with highest
      ! latitude and still on closed field lines; the point
      ! next to it (down in latitude) is set to be boundary.
      ! Result of this subsection is to populate arrays BNDLOC
      ! and IMIN_J with values:

      do j=1,jsize
         bndloc(j) = 2  ! if everything else fails, can use this...
         do i=isize,2,-1
            if(iopen(i,j) >= 0)then
               bndloc(j) = i + 1
               exit
            endif
         end do
       ! reset imin_j
       imin_j = ceiling(bndloc)
         IF (L_write_vars_debug) then
            write(*,*)' bndy ',bndloc(j),j,vm(imin_j(j),j)
         END IF
      end do

      ! fits an ellipse to the boundary
      if (use_ellipse)then
        CALL Set_ellipse(isize,jsize,rmin,pmin,vm,big_vm,bndloc,iopen)
        imin_j = ceiling(bndloc)
      end if

!      IF(set_boundary_with_mach)then
!       write(6,*)' boundary before set_bndy_w_mach'
!       write(6,*) bndloc
!       CALL Set_bndy_w_mach(isize,jsize,rmin,pmin,vm,big_vm,bndloc,iopen)
!       write(6,*)' boundary after set_bndy_w_mach'
!       write(6,*)bndloc
!      END IF

      IF(set_boundary_with_beta)then
       do j=1,jsize
        do i=ceiling(bndloc(j)),isize-1
          IF(beta_average(i,j) > max_beta)then
            bndloc(j) = i+1
            vm(1:i,j) = big_vm
            iopen(1:i,j) = 0
          END IF
        end do
       end do
       imin_j = ceiling(bndloc)
      END IF

      ! smooth boundary location 
      do ns=1,n_smooth
         if (doRCMVerbose) write(6,*)' smoothing rcm boundary, ns =', ns
         call smooth_boundary_location(isize,jsize,jwrap,bndloc)
         call reset_rcm_vm(isize,jsize,bndloc,big_vm,imin_j,vm,iopen) ! adjust Imin_j
      end do

      if (n_smooth < 1) call reset_rcm_vm(isize,jsize,bndloc,big_vm,imin_j,vm,iopen)

! reset mapping points on open field lines
       do j=1,jsize
        do i=1,imin_j(j)-1
          rmin(i,j) = 0.0
          pmin(i,j) = 0.0
        end do
       end do
! find the inner boundary of the lfm in the rcm grid, key on pressure
       do j=1,jsize
         inner_bndy(j) = isize ! default
        do i=imin_j(j),isize
         if(press(i,j) < 1.0e-15)then
           inner_bndy(j) = i+1
           exit
         endif
        end do
       end do
       
       IF (L_write_vars_debug) then
          write(*,*)' location of lfm inner boundary in the rcm is:'
          write(*,*)inner_bndy
       END IF

!---->Set new EETA from lfm code pressure. 
!     On open field lines, values of ETA will be zero:

      CALL Gettemp (ierr)

      IF (ierr < 0) RETURN
      
      CALL Read_alam (kcsize, alamc, ikflavc, fudgec, almdel, almmax, almmin, iesize, ierr)
      IF (ierr < 0) RETURN

      !CALL Press2eta       ! this populates EETA_NEW array
      CALL Press2eta(RM%planet_radius)       ! this populates EETA_NEW array

      if(maxval(eeta_new) <=0)then
       write(6,*)' something is wrong in the new eeta arrays'
       Stop
      end if

      ! It must come after calls to calc_ftv and press2eta.
      ! this is where we must set initial conditions:

      IF (icontrol==RCMCOLDSTART) THEN
         write(6,*)' TORCM: initializing the RCM arrays at t=',itimei
         bndloc_old = bndloc
         imin_j_old = imin_j
      ! initialize the dynamic plasmasphere   sbao 03282020
         call set_plasmasphere(isize,jsize,kcsize,xmin,ymin,vm,eeta_new,imin_j)
       ! eeta_new(:,:,1) = 0.0  ! to see the refilling model alone
         eeta       = eeta_new  ! this is initial conditions on plasma
      END IF

       ! reset the static part of the plasmasphere sbao 07292020
         if (staticR > 2.0) then
          do j=1,jsize
           do i=imin_j(j),isize
             if(rmin(i,j) < 2.0 .and. vm(i,j) > 0.0)then
                  eeta (i,j,1) = eeta_pls0 (i,j)
            end if
           end do
           end do
         end if
 
      ! just in case:
      imin_j     = CEILING(bndloc)
      imin_j_old = CEILING(bndloc_old)

      if (doRCMVerbose) then
        write(*,*)'imin_j',imin_j
        write(*,*)'imin_j_old',imin_j_old
      endif

   ! Set boundary conditions on plasma (EETA) for all MLT's and energy levels:
   !$OMP PARALLEL DO default(shared) &
   !$OMP schedule(dynamic) &
   !$OMP private(j,k,inew,iold)
    DO k=1,kcsize
      DO j=1,jsize

        inew = imin_j(j)
        iold = imin_j_old(j)

        ! values on the (new) boundary and outside are from MHD:

        IF (inew < iold) then

          ! There are newly-acquired points inside the modeling region,
          ! assign MHD-produced values to them:

          eeta (inew+1:iold,j,k) = eeta_new(inew+1:iold,j,k)

        END IF

      END DO
    END DO


! get the boundary eeta local time
      do j=1,jsize
        do k=1,kcsize
          eetabnd(j,k) = eeta(imin_j(j),j,k)
          if (eetabnd(j,k) <= 0.0 .and. doRCMVerbose) then
            write(6,'(a,i3,1x,i3)')' warning: eetabnd <= 0 at j,k =',j,k
            write(6,'(a,i3,a,i3,a,g14.5)')' eetabnd(',j,',',k,')=',eetabnd(j,k)
            write(6,'(a,i3,a,i3,a,g14.5)')' vm(',imin_j(j),',',j,')=',vm(imin_j(j),j)
          end if
        end do

      end do
      
! this is off for now
! set eeta to start with and empty magnetosphere
!     if(itimei==0)then
!     do j=1,jsize
!      do i=imin_j(j)+1,isize
!       do k = 1,ksize
!        eeta(i,j,k) = 1.0
!       end do
!      end do
!     end do
!     end if
      

! smooth eeta at the boundary
      CALL Smooth_eta_at_boundary(isize,jsize,kcsize,jwrap,eeta,imin_j)

! now reset eeta outside the rcm to be eeta at the boundary
      do j=1,jsize
        do i=1,imin_j(j)-1
          eeta(i,j,:) = eeta(imin_j(j),j,:)
        end do
      end do
! import ionosphere
   call Ionosphere_toRCM(RM)

!-----------------------Write to rcmu_torcm.dat-----------------
! this file is for plotting and/or debugging. RCM itselt does not need it.

      IF (L_write_rcmu_torcm ) then
          IF (itimei==0) THEN
              OPEN (LUN, FILE=rcmdir//'rcmu_torcm.dat',form='unformatted', STATUS = 'replace')
              WRITE (LUN) isize,jsize,kcsize
          ELSE
              OPEN (LUN,file=rcmdir//'rcmu_torcm.dat',form='unformatted', status='old',position='append')
          END IF
          WRITE (LUN) itimei, alamc, x0, y0, z0, rmin, pmin,&
                      iopen, vm, press, den, bmin, ti, te, beta_average,v,eeta_new
          CLOSE (LUN)
      END IF
!----------------- end write to rcmu_torcm.dat-------------------


      ! Complete preparing B-field arrays for RCM. We already have
      ! vm, bndloc, rmin, and pmin. Still need xmin and Ymin:
      ! xmin,ymin comes from the MHD code

!      xmin = rmin * COS (pmin)
!      ymin = rmin * SIN (pmin)
      DO j = 1, jsize
        xmin(1:imin_j(j)-1,j) = xmin(imin_j(j),j) 
        ymin(1:imin_j(j)-1,j) = ymin(imin_j(j),j) 
      END DO

      DO j = 1, jsize
        DO i = imin_j(j),isize
          IF (vm(i,j) <= 0.0) STOP 'vm problem in TORCM'
        END DO
      END DO

      RETURN
      END SUBROUTINE Torcm
!----------------------------------------------------------
      SUBROUTINE Calc_ftv (RM,big_vm,ierr) 
      USE conversion_module
      USE rcm_precision
      USE constants, ONLY : mass_proton,gamma,one_over_gamma,mu0,radius_earth_m,nt
      USE RCM_mod_subs,ONLY : isize,jsize,kcsize,bmin, vm, bir,sini,rmin,pmin,&
                              xmin,ymin,zmin,vbnd,pi,jwrap
      USE rice_housekeeping_module
      use rcm_mhd_interfaces

      IMPLICIT NONE
      type(rcm_mhd_T), intent(in) :: RM
      REAL(rprec), INTENT (IN) :: big_vm
      INTEGER(iprec), INTENT (OUT) :: ierr
!
!===================================================================
!
! calc_ftv: gets flux tube volumes, mapping parameters (r,p)
!            be, sini, open/closed flag, and pressure on each 2D
!
! inputs: 
!      isize = dimension in latitude
!      jsize = dimension in local time
!     kcsize = channel dimensions (not used here)
!  x0_sm,y0_sm,z0_sm = location of grid points in ionosphere in sm coordinates
!   big_vm  = value to set open field lines
!
!  outputs:
!        be = magnitude of b field in the equatorial plane in nT
!        vm = flux tube volume**(-2/3)
!             note: vm for open fieldlines is flagged with a
!             large value (big_vm) (input)
!      sini = dip angle in the ionosphere from ground
!      rmin = equatorial r location of mapping in Re
!      pmin = equatorial lt location of mapping in radians
!    iopen  = 2d array flagging open/undefined (0/1) or closed (-1) pnts
!             on rcm grid
!    press  = pressure on field line, data from the mhd code in Pa
!    den    = density on field line, data from the mhd code in ple/cc 
!
!===================================================================
!
!  this routine maps either to the magnetopause        
!  or to the conjugate hemisphere and generates        
!  a potential distribution on a grid in the polar cap 
!  it does this by tracing field lines                 
!  using the trac tracer.                             
!  uses gehavo for the field line tracing              
!  created 03/90 last mod 08/95                       
!  01/19 - frt -this version gets all the values from the MHD code
!
      INCLUDE 'rcmdir.h'
!
      character(LEN=20) grid
      character(LEN=80) header
      integer(iprec),parameter :: imax = 2000
      real(rprec):: xx(imax), yy(imax), zz(imax) ,bbb(imax)
      real(rprec):: ftv(imax)

      real(rprec):: xs,ys,zs,xf,yf,zf,xe,ye,ze
      real(rprec):: dir,gla,glo,pot,poten,ftv1
      real(rprec):: ds,pp1,pp2,pp3
      real(rprec):: press1,dens1,x1,y1,z1,beta1
      real(rprec):: bxe,bye,bze,bbe

      real, parameter :: press_min =0.0 !1.e-11 ! Pa

      integer(iprec) :: jj,ival,i,j,iiit,ii,numm,k
      integer(iprec) :: i0,j0,i2
      integer(iprec), parameter:: use_avg = 1
      REAL(rprec), PARAMETER :: den_min = 0.0e6 ! in ples/m^3
      real(rprec) :: rade
      REAL(rprec) :: rdist1,rdist2
      REAL(rprec) :: dx,dy,dz,bf,bx_ion_sm,by_ion_sm,bz_ion_sm
      REAL(rprec) :: bx_ion_geo,by_ion_geo,bz_ion_geo,radius_ion,bradial_ion

      press (:,jwrap:jsize) = RM%Pave (:,      :)
      den (:,jwrap:jsize) = RM%Nave (:,      :)

      xmin (:,jwrap:jsize) = RM%x_bmin (:,      :,1)/RM%planet_radius
      ymin (:,jwrap:jsize) = RM%x_bmin (:,      :,2)/RM%planet_radius
      zmin (:,jwrap:jsize) = RM%x_bmin (:,      :,3)/RM%planet_radius

      bmin (:,jwrap:jsize) = RM%bmin (:,      :)/nt ! in nT
      beta_average (:,jwrap:jsize) = RM%beta_average (:,      :)
      iopen (:,jwrap:jsize) = RM%iopen (:,      :)
! wrap
      do j=1,jwrap-1
       press (:,         j)  = press (:,jsize-jwrap+j)
       den (:,           j)  = den (:,jsize-jwrap+j)
       xmin (:,          j)  = xmin (:,jsize-jwrap+j)
       ymin (:,          j)  = ymin (:,jsize-jwrap+j)
       zmin (:,          j)  = zmin (:,jsize-jwrap+j)
       bmin (:,          j)  = bmin (:,jsize-jwrap+j)
       beta_average (:,  j)  = beta_average (:,jsize-jwrap+j)
       iopen (:,         j)  = iopen (:,jsize-jwrap+j)
      end do

      ! now compute vm and find the boundary
      do j=jwrap,jsize
       vbnd(j) =1
       vm(:,j) = big_vm
       do i=isize,1,-1
        if( iopen(i,j).ge.0)then
                vbnd(j) = i
                exit
        else
                vm(i,j) = 1.0/(RM%vol(i,j-jwrap+1)*nt)**0.667 ! (nt/re)^0.667
        endif

       end do
      end do

      do j=1,jwrap-1
       vm (:, j) = vm (:,jsize-jwrap+j)
       vbnd(j)   = vbnd(jsize-jwrap+j)
      end do

      ! compute rmin,pmin
      rmin = sqrt(xmin**2 + ymin**2 + zmin**2)
      !rmin = sqrt(xmin**2 + ymin**2)

      pmin = atan2(ymin,xmin)

      ierr = 0

      RETURN 
      END SUBROUTINE Calc_ftv
!--------------------------------------------------
!
      SUBROUTINE Gettemp ( ierr)
!
! routine to compute an estimate for temperature given
! the pressure assume that alam*vm = energy = kt
! 2/00 frt
! bug fix to fac 2/19 frt
!
! inputs:
! isize,jsize - rcm grid dimensions
!       r,p - equatorial mapping location of a grid point in Re and rad
!     press - press in Pa
!        vm - flux tube volume^(-2/3) in (nt/Re/^(-2/3)   
!      open - open/closed flag (-1 = closed)
!       den - ple density in ple/m^3
! outputs:
!      ti,te - temperaure if ions and electrons in kelvin
!      ierr  - error flag
!
!-------------------------------------------------
      USE Rcm_mod_subs, ONLY : isize,jsize
      USE rcm_precision
      USE conversion_module
      USE CONSTANTS, ONLY: boltz,tiote
      IMPLICIT NONE
!
      REAL(rprec), PARAMETER :: den0=0.5E+6
      REAL(rprec), PARAMETER :: fac = tiote/(1.+tiote)
!
      INTEGER(iprec) :: i,j,ierr
!
!    set the temperature:

      DO j=1,jsize
       DO i=1,isize
        IF (iopen(i,j) == -1 .and. den(i,j)>0.) THEN
         ti(i,j) = fac*press(i,j)/den(i,j)/boltz
         te(i,j) = ti(i,j)/tiote
        ELSE
         ti(i,j) = 0.0
         te(i,j) = 0.0
        END IF
       END DO
      END DO
!
      ierr = 0

      if(maxval(press) <=0.)then
       write(6,*)' maxval pressure < 0 in gettemp'
       ierr = -1
      end if
!
      RETURN
      END SUBROUTINE Gettemp
!

!====================================================================
!
      SUBROUTINE Press2eta(planet_radius) 

! ====================================================================
!
!  purpose: 
!     convert mhd/lfm pressure quantities to rcm quantities
! 
!  input:  
!     isize     number of grid points in latitudinal direction
!     jsize      number of grid points in longitudinal direction
!     vm(isize,jsize)       (flux tube vloume)^(-2/3) at each grid point
!     den(i,j)             density in particles/m^3
!     tempi(isize,jsize)     ion number temperature on that field line in K
!     tempe(isize,jsize)     electron temperature on that field line in K
!     iopen(isize,jsize)     label to define if fieldline is open (0/1)
!                          or closed (-1)
!     press(isize,jsize)     total pressure on that field line in Pa
!     kcsize            number of energy channels
!     alam(kcsize)      the energy invarant of each energy channel
!     almdel(kcsize)    the width  of each energy channel
!     almmax(kcsize)    the max alam of each energy channel
!     almmin(kcsize)    the min alam of each energy channel
!     iflav(kcsize)     the species of each energy channel
!     numspe(3)        number of channels of each species, e-, H+, O+
!
!  output:
!     eeta_new(isize,jsize,kcsize)       flux tube content at this time
! 
!  other:    
!     mass of proton      mass of ion (kg)
!     mass of electron    mass of electron (kg)
!     boltz    Boltzman constant
!     eV       electron volts
!     trans    conversion factor (re/nt )
!
! ====================================================================
! notes:
!    This subroutine makes use of the Erf (actually the double 
!    precision version derf, erf returns erronous results), which is
!    not part of the standard fortran functions, but is very common.
!
!    1/19/2000 - frt
!    May 19, 2020 - removed idim,jdim,kdim - frt
!
! --------------------------------------------------------------------
!
      USE conversion_module
      USE rcm_precision
      USE RCM_mod_subs, ONLY : ikflavc,vm,alamc,pi,isize,jsize,kcsize
      USE CONSTANTS, ONLY : boltz,ev,nt!,pressure_factor,radius_earth_m
      IMPLICIT NONE
      real(rprec):: xmin,xmax
      real(rprec):: ptemp,eta_correction
!
!      real(rprec):: Erf
!      EXTERNAL Erf
!
      real(rprec),parameter :: eps=1.0e-30
      real(rprec):: sqrtpi,factor,fac0,fac1,t
      integer(iprec) :: i,j,k
      
      real(rprec), intent(in) :: planet_radius
      real(rprec) pressure_factor
      pressure_factor = 2./3.*ev/planet_radius*nt
      
      sqrtpi = SQRT(pi)
    
! factor is a modification to the temperature to get a reasonable
! population in the higher energy channels.  Usually 100 suffices,
! set it to 1 for production runs.

      factor = 1.0
      if(factor.gt.1.)then
       write(*,*)' increasing temp in press2eta by a factor', factor
       ti = factor * ti
       te = factor * te
      end if
!
      DO j = 1, jsize
      DO i = 1, isize
!        IF (iopen(i,j) /= -1 .AND. vm(i,j) > 0.0) THEN
!        IF (iopen(i,j) /= -1 .OR. press(i,j) <= 0.0 ) THEN
         IF (iopen(i,j) /= -1 ) THEN
            eeta_new (i,j,:) = 0.0
         ELSE
            DO k = 1, kcsize
               IF (ikflavc(k) == 1) THEN  ! electrons
                  t = te (i,j)
               ELSE  IF (ikflavc(k) == 2) THEN ! ions (protons)
                  t = ti (i,j)
               ELSE
                  STOP 'ILLEGAL IKFLAVC(K) IN PRESS2ETA'
               END IF
               if (t > 0.)then
                fac0 = planet_radius/nt*den(i,j)/((vm(i,j))**1.5)
                xmax = SQRT(ev*ABS(almmax(k))*vm(i,j)/boltz/t)
                xmin = SQRT(ev*ABS(almmin(k))*vm(i,j)/boltz/t)
                fac1 = (Erf(xmax)-Erf(xmin)) -2.0/sqrtpi* &
                       (xmax*EXP(-xmax**2)-xmin*EXP(-xmin**2))
                fac1 = MAX (fac1,eps)
                eeta_new(i,j,k)= fac0 * fac1
               else
                eeta_new(i,j,k) = 0.0
               endif
            END DO
         END IF
      END DO
      END DO
! now check to see if we get the original pressure back

      DO j=1, jsize
       DO i=1, isize
        IF (iopen(i,j) < 0) then
         ptemp = 0.0
         DO k=1,kcsize
          ptemp = ptemp + &
           pressure_factor*ABS(alamc(k))*eeta_new(i,j,k)*vm(i,j)**2.5 
         END DO
         eta_correction = 0.0
         IF (ptemp/= 0)eta_correction = press(i,j)/ptemp
!         write(100,*)i,j,eta_correction
         DO k=1,kcsize
          eeta_new(i,j,k) = eeta_new(i,j,k)*eta_correction
         END DO
         ptemp = 0.0
         DO k=1,kcsize
          ptemp = ptemp + &
           pressure_factor*ABS(alamc(k))*eeta_new(i,j,k)*vm(i,j)**2.5 
         END DO
         IF (ABS(ptemp-press(i,j)) > 0.01*press(i,j))then
          Write(6,*)' Warning, pressures do not match at i,j =',i,j
          Write(6,*)' input pressure = ',press(i,j)
          Write(6,*)' RCM pressure = ',ptemp
          write(6,*)' vm =', vm(i,j)
!         PAUSE
         END IF
        END IF
       END DO
      END DO

!
      RETURN
      END SUBROUTINE Press2eta
 
!--------------atran10
 
      subroutine atran10(xs,ys,zs,gla,glo)
      USE rcm_precision
      USE Rcm_mod_subs, ONLY : pi
      implicit none
      real(rprec) :: xs,ys,zs,gla,glo,rrr,rho
      real(rprec) :: atan5 
                         
!      pi = acos(-1.0)

      gla = 0.0
      glo = 0.0    
      rrr = sqrt(xs**2+ys**2+zs**2)
      if(rrr.eq.0.0)return 

      xs = xs/rrr
      ys = ys/rrr
      zs = zs/rrr

      rho = sqrt(xs**2+ys**2)   

      glo = 0.0
      if(rho.ne.0.)then
      glo = atan5(-ys,xs)    
      glo = glo/pi*180. 
      if(glo.lt.0.0)glo=glo+360.
      end if   

      gla = atan5(rho,zs)
      gla = 90.-gla*180./pi

      return
      end        
! -----------------atan5
      function atan5(y,x)
      USE rcm_precision
      USE Rcm_mod_subs, ONLY : pi
      implicit none

      real(rprec) :: atan5,y,x
!      pi = acos(-1.0)
      atan5 = 0.00
      if((x.ne.0.0).or.(y.ne.0.0))   atan5 = atan2(y,x)
 
      return
      end

!===================================================================
    SUBROUTINE Set_ellipse(idim,jdim,rmin,pmin,vm,big_vm,bndloc,iopen)

! routine that fits an ellipse and resets the modeling
! boundary to be inside the ellipse
! 6/03 frt
! inputs:
!	idim,jdim - size of the 2d rcm arrays (lat, long)
!	xe,ye - equatorial mapping point of field line from rcm grid
!	vm - computed flux tube volume ^(-2/3) set to big_vm if open
!            or outside the ellipse boundary (output)
! 	bndloc/imin_j - boundary location (output)
!	a1 - dayside end of ellipse
!	a2 - nightside location of the ellipse
!	b - semi minor axis (y) of ellipse

!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      USE rice_housekeeping_module, ONLY : ellBdry
      implicit none
      integer(iprec) :: idim,jdim
      real(rprec) :: rmin(idim,jdim), pmin(idim,jdim)
      real(rprec) :: xe(idim,jdim), ye(idim,jdim)
      real(rprec) :: vm(idim,jdim)
      integer(iprec) :: iopen(idim,jdim)
      real(rprec) :: bndloc(jdim)
      real(rprec) :: big_vm,a1,a2,a,b,x0,ell
      real(rprec) :: xP,xM,yMax
      integer(iprec) :: i,j
!  x0 = (a1 + a2)/2
!   a = (a1 - a2)/2
!
!                   b
!           |       |
!           |       |
! x<a1------0-------x0-------------a2
!           |       |
!           |       |
!                   b
!

      xe = rmin * cos(pmin)
      ye = rmin * sin(pmin)

!K: Replacing these hard-coded values with ellipse type set by XML file
      a1 = ellBdry%xSun
      a2 = ellBdry%xTail
      b  = ellBdry%yDD

      if (ellBdry%isDynamic) then
        !Tune to current equatorial bounds
        xP = maxval(xe)
        xM = minval(xe)
        yMax = maxval(abs(ye))
        !Enforce max's from XML ellipse
        a1 = min(a1,xP)
        a2 = max(a2,xM)
        b  = min(b ,yMax)        
      endif
      
      x0 = (a1 + a2)/2.
      a  = (a1 - a2)/2.
      do j=1,jdim
        do i=ceiling(bndloc(j)),idim-1
! now check to see if the point is outside the ellipse, if so
! reset open and bndloc        
          ell = ((xe(i,j)-x0)/a)**2+(ye(i,j)/b)**2
          if(ell > 1.) then
            bndloc(j) = i+1
            iopen(i,j) = 0
            vm(i,j) = big_vm
          end if
        end do
      end do 

    end subroutine Set_ellipse


!------------------------------------

      SUBROUTINE Smooth_eta_at_boundary(idim,jdim,kdim,jwrap,eeta,imin_j)
! this routine attempts to smooth out high frequency noise at the boundary
! of the rcm 
! written 2/06 frt
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      IMPLICIT NONE
      INTEGER(iprec) :: idim,jdim,kdim,jwrap
      INTEGER(iprec) :: imin_j(jdim)
      INTEGER(iprec) :: i,j,k,jm,jmm,jp,jpp
      REAL(rprec) :: eeta(idim,jdim,kdim)
      REAL(rprec) :: eetas2d(jdim,kdim)
! these are the smoothing weights
      REAL(rprec), PARAMETER :: a1 = 1.0  
      REAL(rprec), PARAMETER :: a2 = 1.0  
      REAL(rprec), PARAMETER :: a3 = 2.0  
      REAL(rprec), PARAMETER :: a4 = 1.0  
      REAL(rprec), PARAMETER :: a5 = 1.0  
! now do the smoothing

      do k=1,kdim
       do j=1,jdim

! 1 <=> jdim -jwrap +1
! jdim <=> jwrap
       jmm = j - 2
       if(jmm < 1)jmm = jdim - jwrap - 1       
       jm  = j - 1
       if(jm < 1) jm = jdim - jwrap       
       jpp = j + 2
       if(jpp > jdim)jpp = jwrap + 2
       jp  = j + 1
       if(jp > jdim) jp = jwrap + 1
! now smooth
       eetas2d(j,k) = &
          ( a1*eeta(imin_j(jmm),jmm,k) + &
            a2*eeta(imin_j(jm ),jm ,k) + &
            a3*eeta(imin_j(j  ),j  ,k) + &
            a4*eeta(imin_j(jp ),jp ,k) + &
            a5*eeta(imin_j(jpp),jpp,k) )/(a1+a2+a3+a4+a5) 

           end do
          end do
! now reset the boundary values
        do k=1,kdim
         do j=1,jdim
          eeta(imin_j(j),j,k) = eetas2d(j,k)
         end do      
        end do      
       return
END SUBROUTINE Smooth_eta_at_boundary

      SUBROUTINE Smooth_boundary_location(idim,jdim,jwrap,bndloc)
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rice_housekeeping_module
      IMPLICIT NONE
      INTEGER(iprec), INTENT(IN) :: idim,jdim,jwrap
      REAL(rprec), INTENT(IN OUT) :: bndloc(jdim)

      INTEGER(iprec) :: i,j,jp,jm
      REAL(rprec) :: bndloc_new(jdim)
      REAL(rprec), PARAMETER :: am = 1, a0 =2, ap = 1

     
      if (L_write_vars_debug) then
       write(6,*)' smooth_boundary_location, old boundary'
       write(6,*)bndloc
      end if
! 1 <=> jdim -jwrap +1
! jdim <=> jwrap
      do j=1,jdim
       jm  = j - 1
       if(jm < 1) jm = jdim - jwrap 
       jp  = j + 1
       if(jp > jdim) jp = jwrap + 1

       bndloc_new(j) = (am*bndloc(jm) + a0*bndloc(j) + ap*bndloc(jp))/(am+a0+ap)

      end do

      bndloc(:) = max(bndloc_new(:),bndloc(:))

      if (L_write_vars_debug) then
       write(6,*)' smooth_boundary_location, new boundary'
       write(6,*)bndloc
      end if

      return
    
      END SUBROUTINE Smooth_boundary_location

!
      subroutine set_plasmasphere(idim,jdim,kdim,xmin,ymin,vm,eeta,imin_j)
! subroutine set_plasmasphere(idim,jdim,kdim,rmin,pmin,vm,eeta,imin_j)
! crude routine to set a plasmasphere model in the rcm
! alam(1) should be set to a small value (0.01)
! 2/07 frt
! Use the gallagher model for initial condition, density in ple/m^3 and store it in eeta_pls0
! alam(1) is set to be 0
! sbao 03/25

      USE Rcm_mod_subs, ONLY : eeta_pls0 
      USE rcm_precision
      USE earthhelper, ONLY : GallagherXY
      USE constants, ONLY: density_factor
      USE rice_housekeeping_module, ONLY: InitKp
      IMPLICIT NONE

      integer(iprec) :: idim,jdim,kdim
      real(rprec) :: dens_gal = 0.0
      integer(iprec) :: imin_j(jdim)
      real(rprec) :: vm(idim,jdim),xmin(idim,jdim),ymin(idim,jdim)
      real(rprec) :: eeta(idim,jdim,kdim)

      integer(iprec) :: i,j,k
!      real(rprec), parameter :: radius_ps = 5.0
!      real(rprec), parameter :: dens_ps = 10.0e6 ! ple/m^3

!      do j=1,jdim
!       do i=imin_j(j),idim
!        if(rmin(i,j) < radius_ps .and. vm(i,j) > 0.0)then
!        eeta(i,j,1) = dens_ps/(1.5695e-16*vm(i,j)**1.5)
!        end if
!       end do
!      end do

        do j=1,jdim
        do i=imin_j(j),idim
                if(vm(i,j) > 0.0)then
                  dens_gal = GallagherXY(xmin(i,j),ymin(i,j),InitKp)*1.0e6
                ! add to the existing eeta - frt
                  eeta(i,j,1) = eeta(i,j,1) + dens_gal/(density_factor*vm(i,j)**1.5)
                  eeta_pls0(i,j) = eeta(i,j,1)
                end if
        end do
        end do


      return

      end subroutine set_plasmasphere

!-------------------------------------
      subroutine print_max(nx,ny,nz,label,array,x,y,z)
! routine to find, and print the max value of a quantity on a 3D array
! and print its location and value
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      IMPLICIT NONE
      INTEGER(iprec), INTENT(IN) :: nx,ny,nz
      REAL(rprec), INTENT(IN) :: x(nx),y(ny),z(nz)
      REAL(rprec), INTENT(IN) :: array(nx,ny,nz)
      CHARACTER(*), INTENT(IN) :: label

      INTEGER(iprec) :: location(3)
      REAL(rprec) :: max_array
!     INTEGER :: inlen

      max_array = maxval(array)
      location = maxloc(array)
!     inlen = 1 + int(log10(real(max(nx,ny,nz)))) ! length of the arrays

      write(6,'(a,a,a,g14.6,a,i4,a,f10.3,a,i4,a,f10.3,a,i4,a,f10.3)') &
                   '  max value of ',label,' is ',max_array,&
                   ' at: x(',location(1),')=',x(location(1)),&
                   ', y(',location(2),')=',y(location(2)),&
                   ', z(',location(3),')=',z(location(3))

      return
      end subroutine print_max
!-------------------------------------
      subroutine print_max_2d(nx,ny,label,array,x,y)
! routine to find, and print the max value of a quantity on a 3D array
! and print its location and value
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      IMPLICIT NONE
      INTEGER(iprec), INTENT(IN) :: nx,ny
      REAL(rprec), INTENT(IN) :: x(nx),y(ny)
      REAL(rprec), INTENT(IN) :: array(nx,ny)
      CHARACTER(*), INTENT(IN) :: label

      INTEGER(iprec) :: location(2)
      REAL(rprec) :: max_array
!     INTEGER :: inlen

      max_array = maxval(array)
      location = maxloc(array)
!     inlen = 1 + int(log10(real(max(nx,ny,nz)))) ! length of the arrays

      write(6,'(a,a,a,g14.6,a,i4,a,f10.3,a,i4,a,f10.3,a,i4,a,f10.3)') &
                   '  max value of ',label,' is ',max_array,&
                   ' at: x(',location(1),')=',x(location(1)),&
                   ', y(',location(2),')=',y(location(2))

      return
      end subroutine print_max_2d
!-------------------------------------
      subroutine print_max_integer(nx,ny,nz,label,array,x,y,z)
! routine to find, and print the max value of a quantity on a integer 3D array
! and print its location and value
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      IMPLICIT NONE
      INTEGER(iprec), INTENT(IN) :: nx,ny,nz
      REAL(rprec), INTENT(IN) :: x(nx),y(ny),z(nz)
      INTEGER(iprec), INTENT(IN) :: array(nx,ny,nz)
      CHARACTER(*), INTENT(IN) :: label

      INTEGER(iprec) :: location(3)
      INTEGER(iprec) :: max_array
!     INTEGER :: inlen

      max_array = maxval(array)
      location = maxloc(array)
!     inlen = 1 + int(log10(real(max(nx,ny,nz)))) ! length of the arrays

      write(6,'(a,a,a,i5,a,i4,a,f10.3,a,i4,a,f10.3,a,i4,a,f10.3)') &
                   '  max value of ',label,' is ',max_array,&
                   ' at: x(',location(1),')=',x(location(1)),&
                   ', y(',location(2),')=',y(location(2)),&
                   ', z(',location(3),')=',z(location(3))

      return
      end subroutine print_max_integer
!-------------------------------------
      subroutine print_min(nx,ny,nz,label,array,x,y,z)
! routine to find, and print the min value of a quantity on a 3D array
! and print its location and value
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      IMPLICIT NONE
      INTEGER(iprec), INTENT(IN) :: nx,ny,nz
      REAL(rprec), INTENT(IN) :: x(nx),y(ny),z(nz)
      REAL(rprec), INTENT(IN) :: array(nx,ny,nz)
      CHARACTER(*), INTENT(IN) :: label

      INTEGER(iprec) :: location(3)
      REAL(rprec) :: min_array
!     INTEGER :: inlen

      min_array = minval(array)
      location = minloc(array)
!     inlen = 1 + int(log10(real(max(nx,ny,nz)))) ! length of the arrays

      write(6,'(a,a,a,g14.6,a,i4,a,f10.3,a,i4,a,f10.3,a,i4,a,f10.3)') &
                   '  min value of ',label,' is ',min_array,&
                   ' at: x(',location(1),')=',x(location(1)),&
                   ', y(',location(2),')=',y(location(2)),&
                   ', z(',location(3),')=',z(location(3))

      return     
      end subroutine print_min
!-------------------------------------
      subroutine print_min_2d(nx,ny,label,array,x,y)
! routine to find, and print the min value of a quantity on a 3D array
! and print its location and value
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      IMPLICIT NONE
      INTEGER(iprec), INTENT(IN) :: nx,ny
      REAL(rprec), INTENT(IN) :: x(nx),y(ny)
      REAL(rprec), INTENT(IN) :: array(nx,ny)
      CHARACTER(*), INTENT(IN) :: label

      INTEGER(iprec) :: location(2)
      REAL(rprec) :: min_array
!     INTEGER :: inlen

      min_array = minval(array)
      location = minloc(array)
!     inlen = 1 + int(log10(real(max(nx,ny,nz)))) ! length of the arrays

      write(6,'(a,a,a,g14.6,a,i4,a,f10.3,a,i4,a,f10.3,a,i4,a,f10.3)') &
                   '  min value of ',label,' is ',min_array,&
                   ' at: x(',location(1),')=',x(location(1)),&
                   ', y(',location(2),')=',y(location(2))

      return     
      end subroutine print_min_2d
!-------------------------------------
      subroutine print_min_integer(nx,ny,nz,label,array,x,y,z)
! routine to find, and print the max value of a quantity on a 3D array
! and print its location and value
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      IMPLICIT NONE
      INTEGER(iprec), INTENT(IN) :: nx,ny,nz
      REAL(rprec), INTENT(IN) :: x(nx),y(ny),z(nz)
      INTEGER(iprec), INTENT(IN) :: array(nx,ny,nz)
      CHARACTER(*), INTENT(IN) :: label

      INTEGER(iprec) :: location(3)
      INTEGER(iprec) :: min_array
!     INTEGER :: inlen

      min_array = minval(array)
      location = minloc(array)
!     inlen = 1 + int(log10(real(max(nx,ny,nz)))) ! length of the arrays

      write(6,'(a,a,a,i5,a,i4,a,f10.3,a,i4,a,f10.3,a,i4,a,f10.3)') &
                   '  min value of ',label,' is ',min_array,&
                   ' at: x(',location(1),')=',x(location(1)),&
                   ', y(',location(2),')=',y(location(2)),&
                   ', z(',location(3),')=',z(location(3))

      return
      end subroutine print_min_integer
!------------------------------------------      
      subroutine reset_rcm_vm(idim,jdim,bndloc,big_vm,imin_j,vm,iopen)
! this routine resets imin_j, vm, and open based on a newly set bndloc      
!      USE Rcm_mod_subs, ONLY: rprec,iprec
      USE rcm_precision
      implicit none
      integer(iprec), intent(in) :: idim,jdim
      integer(iprec) :: i,j
      integer(iprec),intent(inout) :: imin_j(jdim),iopen(idim,jdim)
      real(rprec), intent(in) :: bndloc(jdim)
      real(rprec), intent(in) :: big_vm
      real(rprec), intent(inout) :: vm(idim,jdim)

       imin_j = CEILING(bndloc)

       do j=1,jdim
        do i=1,imin_j(j)-1
          vm(i,j) = big_vm
          iopen(i,j) = 0
          end do
       end do
      return

      end subroutine reset_rcm_vm
! converts geo to sm for iflag = 1 or sm to geo for iflag =-1
! ! frt 7/12 
!       subroutine geo2sm(idim,jdim,xgeo,ygeo,zgeo,xsm,ysm,zsm,iflag)
!       USE Rcm_mod_subs, ONLY : iprec,rprec 
!       USE mhd_scalars, only: iaScalars, YEAR,MONTH,DAY,HOUR,MINUTE,SECOND
!       implicit none
!       INTEGER, INTENT(IN) :: iflag,idim,jdim
!       CHARACTER(LEN=3) :: months(12)=(/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT', 'NOV','DEC'/)
!       INTEGER(iprec) :: iyear,imonth,iday,ihour,iminute,isecond,idoy,i,j
!       REAL(rprec) :: xgeo(idim,jdim),ygeo(idim,jdim),zgeo(idim,jdim)
!       REAL(rprec) ::xsm(idim,jdim),ysm(idim,jdim),zsm(idim,jdim)
!       REAL*8 :: dxgeo,dygeo,dzgeo,dxsm,dysm,dzsm,dxgsw,dygsw,dzgsw

!      iyear   = iaScalars(YEAR)
!      imonth  = iaScalars(MONTH)
!      iday    = iaScalars(DAY)
!      ihour   = iaScalars(HOUR)
!      iminute = iaScalars(MINUTE)
!      isecond = iaScalars(SECOND)

! !    write(*,*)'---time information---'
! !    write(*,*) '   year =',iyear
! !    write(*,*) '  month =',months(imonth)
! !    write(*,*) '    day =',iday
! !    write(*,*) '   hour =',ihour
! !    write(*,*) 'minutes =',iminute
! !    write(*,*) 'seconds =',isecond

!       call date_doy(0,iyear,months(imonth),imonth,iday,idoy)

! !    write(*,*) 'day of year =',idoy
!      call recalc(iyear,idoy,ihour,iminute,isecond)

!      if(iflag == 1)then
!      ! geopack requires double precision
!       do j=1,jdim
!        do i=1,idim
!           ! geo -> sm
!           dxgeo = dble(xgeo(i,j)); dygeo = dble(ygeo(i,j)); dzgeo = dble(zgeo(i,j))
!           call geogsw_08(dxgeo,dygeo,dzgeo,dxgsw,dygsw,dzgsw,1)
!           call smgsw_08(dxsm,dysm,dzsm,dxgsw,dygsw,dzgsw,-1)
!           xsm(i,j) = sngl(dxsm); ysm(i,j) = sngl(dysm); zsm(i,j) = sngl(dzsm)
!         end do
!        end do
!       elseif(iflag==-1)then
!       do j=1,jdim
!        do i=1,idim
!          ! sm -> geo
!          dxsm = dble(xsm(i,j)); dysm = dble(ysm(i,j)); dzsm = dble(zsm(i,j))
!          call smgsw_08(dxsm,dysm,dzsm,dxgsw,dygsw,dzgsw,1)
!          call geogsw_08(dxgeo,dygeo,dzgeo,dxgsw,dygsw,dzgsw,-1)
!           xgeo(i,j) = sngl(dxgeo); ygeo(i,j) = sngl(dygeo); zgeo(i,j) = sngl(dzgeo)
!         end do
!        end do
!       else
!          write(*,*)'Error in geo2sm, wrong iflag'
!          return
!       endif

!      end subroutine geo2sm
!
!======================================
      subroutine allocate_conversion_arrays(isize,jsize,kcsize)
! used to allocate memory for the exchange arrays      
! 7/09 frt
      use conversion_module
!      use rcm_mod_subs, only : iprec
      USE rcm_precision, only : iprec
      implicit none
      integer(iprec),intent(in) :: isize,jsize,kcsize
      integer(iprec) :: idim,jdim,kdim
! if the arrays are allocated, then return
      if(allocated(x0))return

      idim = isize
      jdim = jsize
      kdim = kcsize

      write(*,*)' Allocating conversion arrays'

      ! 1d arrays
      allocate(imin_j_old(jdim))
      allocate(inner_bndy(jdim))
      allocate(bndloc_old(jdim))
      allocate(almmin(kdim))
      allocate(almmax(kdim))
      allocate(almdel(kdim))
      ! 2d arrays
      allocate(x0_sm(idim,jdim))
      allocate(y0_sm(idim,jdim))
      allocate(z0_sm(idim,jdim))
      allocate(x0(idim,jdim))
      allocate(y0(idim,jdim))
      allocate(z0(idim,jdim))
      allocate(den(idim,jdim))
      allocate(press(idim,jdim))
      allocate(deno(idim,jdim))
      allocate(presso(idim,jdim))
      allocate(te(idim,jdim))
      allocate(ti(idim,jdim))
      allocate(to(idim,jdim))
      allocate(beta_average(idim,jdim))
      allocate(eetabnd(jdim,kdim))
      allocate(iopen(idim,jdim))
      ! 3d arrays
      allocate(eeta_new(idim,jdim,kdim))
    
     return

     end subroutine allocate_conversion_arrays 

     
