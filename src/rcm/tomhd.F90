      SUBROUTINE tomhd (RM, ierr)
      USE earthhelper, ONLY : GallagherXY,DP2kT
      USE rcm_precision
      USE Rcm_mod_subs, ONLY : isize, jsize, kcsize,jwrap, nptmax, &
                              colat, aloct, v, birk, &
                              bmin, xmin, ymin, zmin, vm, pmin, rmin, &
                              birk_avg, v_avg, eeta, eeta_avg, &
                              alamc, ikflavc, &
                              boundary, bndloc, pressrcm, &
                              xmass, densrcm,denspsph,imin_j,rcmdir, &
                              eflux,eavg,ie_el
      USE constants, ONLY : mass_proton,mass_electron,nt,ev
      USE rice_housekeeping_module
      Use rcm_mhd_interfaces
      USE kdefs, ONLY : TINY
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

      INTEGER(iprec) :: i, j, L, k, itime,n,klow
      real(rprec) :: max_xmin,min_xmin,max_ymin,min_ymin
      real(rprec) :: dens_plasmasphere
      real(rp) :: sclmass(RCMNUMFLAV) !xmass prescaled to proton
      real(rp) :: kevRCM,kevMHD
      real(rp), parameter :: kRatMax = 0.5
      !AMS 04-22-2020
      real(rprec) :: pressure_factor,density_factor

      pressure_factor = 2./3.*ev/RM%planet_radius*nt
      density_factor = nt/RM%planet_radius

      RM%MaxAlam = maxval(alamc)
      
      ierr = 0

      !Set lowest RC channel
      if (use_plasmasphere) then
        klow = 2
      else
        klow = 1
      endif


      !Set scaled mass by hand here to avoid precision issues
      sclmass(RCMELECTRON) = mass_electron/mass_proton
      sclmass(RCMPROTON) = 1.0

!     Compute rcm pressure and density (on the ionospheric RCM grid):
      !Tweaking scaling for better precision and testing eeta instead of avg (K: 8/20)
      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j,k,dens_plasmasphere,kevRCM,kevMHD)
      DO j = 1, jsize
        DO i = 1, isize
          pressrcm (i,j) = 0.0
          densrcm  (i,j) = 0.0
          denspsph (i,j) = 0.0

          IF (vm(i,j) < 0.0) CYCLE
          DO k = klow, kcsize
            !Pressure calc in pascals
            pressrcm(i,j) = pressrcm(i,j) + pressure_factor*ABS(alamc(k))*eeta(i,j,k)*vm(i,j)**2.5

            !Density calc
            if (alamc(k) > 0.0) then ! only add the ion contribution
              densrcm(i,j) = densrcm(i,j) + density_factor*sclmass(ikflavc(k))*eeta(i,j,k)*vm(i,j)**1.5
            endif
          ENDDO !k

          if (use_plasmasphere) then
            if (dp_on) then 
              ! use plasmasphere channel eeta_avg(:,:,1) sbao 03/2020
              denspsph(i,j) = density_factor*sclmass(RCMPROTON)*eeta(i,j,1)*vm(i,j)**1.5
            else
              ! add a simple plasmasphere model based on carpenter 1992 or gallagher 2002 in ples/cc
              dens_plasmasphere = GallagherXY(xmin(i,j),ymin(i,j))
              denspsph(i,j) = dens_plasmasphere*1.0e6
            endif !dp_on
          endif !use_plasmasphere

          !Do some checking on how well resolved energy channel stuff is
          kevRCM = abs(alamc(kcsize))*vm(i,j)*1.0e-3 !Max keV of RCM channels here
          kevMHD = DP2kT(densrcm(i,j)*rcmNScl,pressrcm(i,j)*rcmPScl) !Get keV from RCM moments
          
          if (kevMHD > kRatMax*kevRCM) then
            !Effective "MHD" temperature, P=nkT_{MHD} is above kRatMax the max RCM channel energy
            !This is probably bad for resolving the distribution so we do some shady cooling here

            !Rescale eeta's to clamp P_{RCM}
            eeta(i,j,klow:) = (kRatMax*kevRCM/kevMHD)*eeta(i,j,klow:)
            pressrcm(i,j)   = (kRatMax*kevRCM/kevMHD)*pressrcm(i,j)

            !write(*,*) 'kevRCM, kevMHD = ', kevRCM,kevMHD,densrcm(i,j)*rcmNScl,pressrcm(i,j)*rcmPScl
          endif

        ENDDO !i
      ENDDO !j

      if (any(isnan(densrcm)) .or. any(isnan(pressrcm))) then
        write(*,*) 'NaN in RCM-DEN/P at tomhd'
        stop
      endif
 
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



      ! At this point, we assume that RCM arrays are all populated
      ! (plasma, b-field, bndloc, etc.):

      IF (L_write_rcmu) then

        rmin = SQRT (xmin**2+ymin**2+zmin**2)
        WHERE (xmin == 0.0 .AND. ymin == 0.0)
          pmin = 0.0
        ELSE WHERE
          pmin = ATAN2 (ymin, xmin)
        END WHERE
        WHERE (pmin < 0.0) pmin = pmin + 2.0*pi

        rmin = SQRT (xmin**2+ymin**2)
        WHERE (xmin == 0.0 .AND. ymin == 0.0)
          pmin = 0.0
        ELSE WHERE
          pmin = ATAN2 (ymin, xmin)
        END WHERE
        WHERE (pmin < 0.0) pmin = pmin + 2.0*pi
        WRITE (*,'(A,I9.9,A,I5.5)') 'TOMHD: Read RCM, T=',itime,', REC=',L

        !call write_rcmu (L,0_iprec)

      END IF

      
      RETURN
      END SUBROUTINE tomhd


