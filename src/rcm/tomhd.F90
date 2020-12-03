MODULE tomhd_mod
  USE rcm_precision
  USE rice_housekeeping_module
  USE rcm_mhd_interfaces
  USE kdefs, ONLY : TINY
  USE constants, ONLY : mass_proton,mass_electron,nt,ev,tiote,boltz
  USE Rcm_mod_subs, ONLY : isize, jsize, kcsize,jwrap,alamc,ikflavc
  USE earthhelper, ONLY : GallagherXY,DP2kT

  implicit none

  real(rp), private :: density_factor = 0.0 !module private density_factor using planet radius
  real(rp), private :: pressure_factor = 0.0
  integer , private :: klow = 1

  logical, parameter, private :: doClamp=.true. !Whether to clamp poorly resolved pressure
  !real(rprec), parameter :: kRatMax = 0.9 !Ratio of kT_{Avg} =P_{RCM}/n_{RCM} over kT_{RCM,Max}

  contains

    SUBROUTINE tomhd (RM, ierr)
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

      USE Rcm_mod_subs, ONLY : bmin,birk,xmin,ymin,zmin,vm,eeta,eeta_avg, &
                               bndloc,pressrcm,densrcm,denspsph,imin_j,eflux,eavg


      IMPLICIT NONE
      type(rcm_mhd_T),intent(inout) :: RM
      INTEGER(iprec), INTENT (OUT) :: ierr

      !Always set p/d_factors
      pressure_factor = 2./3.*ev/RM%planet_radius*nt
      density_factor = nt/RM%planet_radius

      !Set lowest RC channel
      if (use_plasmasphere) then
        klow = 2
      else
        klow = 1
      endif

      !Do some checks on the RCM eta distribution
      if (doClamp) call ClampEta(eeta,eeta_avg,vm)

      !Now pick which eta (instant vs. avg) and calculate moments
      if (doAvg2MHD) then
        call eeta2DP(eeta_avg,vm,densrcm,denspsph,pressrcm)
      else
        call eeta2DP(eeta    ,vm,densrcm,denspsph,pressrcm)
      endif

      RM%MaxAlam = maxval(alamc)
      
!     now update the pressure and density in the mhd code  
      RM%Prcm    = pressrcm(:,jwrap:jsize)
      RM%Nrcm    = densrcm (:,jwrap:jsize)
      RM%Npsph   = denspsph(:,jwrap:jsize)
      RM%flux    = eflux   (:,jwrap:jsize,:)
      RM%eng_avg = eavg    (:,jwrap:jsize,:)
      RM%fac     = birk    (:,jwrap:jsize)

      ierr = 0
    END SUBROUTINE tomhd


    !Do some safety stuff to eta w/ temperature is high
    SUBROUTINE ClampEta(eta,eta_avg,vm)
      REAL(rprec), intent(inout), dimension(isize,jsize,kcsize)  :: eta,eta_avg
      REAL(rprec), intent(in)  :: vm(isize,jsize)

      REAL(rprec), dimension(isize,jsize) :: Drc,Dpp,Prc
      integer :: i,j
      REAL(rprec) :: kevRCM,kevMHD,wMax
      REAL(rprec), dimension(kcsize) :: etaMax

      !Get moments from eta
      call eeta2DP(eta,vm,Drc,Dpp,Prc)

      !Now loop over i,j and check "MHD" vs. top RCM energy
      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j,kevRCM,kevMHD,wMax,etaMax)
      DO j = 1, jsize
        DO i = 1, isize
          IF (vm(i,j) < 0.0) CYCLE
          kevRCM = abs(alamc(kcsize))*vm(i,j)*1.0e-3 !Max keV of RCM channels here
          kevMHD = DP2kT(Drc(i,j)*rcmNScl,Prc(i,j)*rcmPScl) !Get keV from RCM moments
          wMax = kevMHD/kevRCM !Want this to be small for well-resolved pressure
          wMax = min(wMax,1.0) !Ensure <= 1

          !Blend w/ Maxwellian
          call DP2eeta(Drc(i,j),Prc(i,j),vm(i,j),etaMax)

          eta(i,j,klow:) = (1.0-wMax)*eta(i,j,klow:) + wMax*etaMax(klow:)

          ! if (kevMHD >= kRatMax*kevRCM) then
          !   !Effective "MHD" temperature, P=nkT_{MHD} is above kRatMax the max RCM channel energy
          !   !This is probably bad for resolving the distribution so we do some shady cooling here

          !   !Rescale eeta's to clamp P_{RCM}
          !   eta    (i,j,klow:) = (kRatMax*kevRCM/kevMHD)*eta    (i,j,klow:)
          !   eta_avg(i,j,klow:) = (kRatMax*kevRCM/kevMHD)*eta_avg(i,j,klow:)
          ! endif

        ENDDO
      ENDDO

    END SUBROUTINE ClampEta

    !Convert given single density/pressure to eeta
    SUBROUTINE DP2eeta(Drc,Prc,vm,eta)
      USE conversion_module, ONLY : almmax,almmin,erfexpdiff
      REAL(rprec), intent(in)  :: Drc,Prc,vm
      REAL(rprec), intent(out) :: eta(kcsize)

      REAL(rprec), PARAMETER :: fac = tiote/(1.+tiote)

      REAL(rprec) :: Tk,ti,te,A0,prcmI,prcmE,pmhdI,pmhdE
      REAL(rprec) :: xp,xm,pcon,psclI,psclE
      INTEGER(iprec) :: k
      logical :: isIon

      eta = 0.0
      if ( (vm<0) .or. (Drc<TINY) ) return

      !Get ion/electron temperature
      ti = fac*Prc/Drc/boltz
      te = ti/tiote

      A0 = (Drc/density_factor)/(vm**1.5)
      prcmI = 0.0 !Cumulative ion pressure
      prcmE = 0.0 !Cumulative electron pressure

      !Loop over non-zero channels
      do k=klow,kcsize
        !Get right temperature
        IF (ikflavc(k) == RCMELECTRON) THEN  ! electrons
          Tk = te
          isIon = .false.
        ELSE IF (ikflavc(k) == RCMPROTON) THEN ! ions (protons)
          Tk = ti
          isIon = .true.
        ENDIF

        xp = SQRT(ev*ABS(almmax(k))*vm/boltz/Tk)
        xm = SQRT(ev*ABS(almmin(k))*vm/boltz/Tk)
        !Use quad prec calc of erf/exp differences
        eta(k) = erfexpdiff(A0,xp,xm)
        

        !Pressure contribution from this channel
        pcon = pressure_factor*ABS(alamc(k))*eta(k)*vm**2.5

        if (isIon) then
          prcmI = prcmI + pcon
        else
          prcmE = prcmE + pcon
        endif

      enddo !k loop

      !Now rescale eeta channels to conserve pressure integral between MHD/RCM
      !In particular, we separately conserve ion/electron contribution to total pressure
      pmhdI = Prc*tiote/(1.0+tiote) !Desired ion pressure
      pmhdE = Prc*  1.0/(1.0+tiote) !Desired elec pressure

      psclI = pmhdI/prcmI
      psclE = pmhdE/prcmE

      !Loop over channels and rescale      
      do k=klow,kcsize
        IF (ikflavc(k) == RCMELECTRON) THEN  ! electrons
          eta(k) = psclE*eta(k)
        ELSE IF (ikflavc(k) == RCMPROTON) THEN ! ions (protons)
          eta(k) = psclI*eta(k)
        ENDIF
      enddo

    END SUBROUTINE DP2eeta

    !Convert given eeta to density (RC/plasmasphere) and pressure
    SUBROUTINE eeta2DP(eta,vm,Drc,Dpp,Prc)
      USE Rcm_mod_subs, ONLY : xmin,ymin,zmin
      IMPLICIT NONE
      REAL(rprec), intent(in)  :: eta(isize,jsize,kcsize)
      REAL(rprec), intent(in)  :: vm(isize,jsize)
      REAL(rprec), intent(out), dimension(isize,jsize) :: Drc,Dpp,Prc

      
      integer :: i,j,k
      real(rprec) :: sclmass(RCMNUMFLAV) !xmass prescaled to proton

      !Set lowest RC channel
      if (use_plasmasphere) then
        klow = 2
      else
        klow = 1
      endif

      !Set scaled mass by hand here to avoid precision issues
      sclmass(RCMELECTRON) = mass_electron/mass_proton
      sclmass(RCMPROTON) = 1.0

      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j,k)
      DO j = 1, jsize
        DO i = 1, isize
          Prc(i,j) = 0.0
          Drc(i,j) = 0.0
          Dpp(i,j) = 0.0

          IF (vm(i,j) < 0.0) CYCLE
          DO k = klow, kcsize
            !Pressure calc in pascals
            Prc(i,j) = Prc(i,j) + pressure_factor*ABS(alamc(k))*eta(i,j,k)*vm(i,j)**2.5

            !Density calc (ring current)
            if (alamc(k) > 0.0) then ! only add the ion contribution
              Drc(i,j) = Drc(i,j) + density_factor*sclmass(ikflavc(k))*eta(i,j,k)*vm(i,j)**1.5
            endif
          ENDDO !k loop

          if (use_plasmasphere) then
            if (dp_on) then 
              ! use plasmasphere channel eeta_avg(:,:,1) sbao 03/2020
              Dpp(i,j) = density_factor*sclmass(RCMPROTON)*eta(i,j,1)*vm(i,j)**1.5
            else
              ! add a simple plasmasphere model based on carpenter 1992 or gallagher 2002 in ples/cc
              Dpp(i,j) = (1.0e6)*GallagherXY(xmin(i,j),ymin(i,j))
              
              
            endif !dp_on
          endif !use_plasmasphere
        ENDDO !i
      ENDDO !j

    END SUBROUTINE eeta2DP

!       SUBROUTINE tomhd (RM, ierr)
!       USE earthhelper, ONLY : GallagherXY,DP2kT
      
!       USE Rcm_mod_subs, ONLY : isize, jsize, kcsize,jwrap, nptmax, &
!                               colat, aloct, v, birk, &
!                               bmin, xmin, ymin, zmin, vm, pmin, rmin, &
!                               birk_avg, v_avg, eeta, eeta_avg, &
!                               alamc, ikflavc, &
!                               boundary, bndloc, pressrcm, &
!                               xmass, densrcm,denspsph,imin_j,rcmdir, &
!                               eflux,eavg,ie_el
      
      
      
      
! ! 
! !==============================================================
! ! purpose:
! !  To convert RCM information (eta), to MHD information (p,n) 
! !  It also writes to files (optional):
! !    'rcmu.dat' = unformatted time-averaged rcm data for analysis
! !
! ! inputs:
! !
! !   4/18/95             rws
! !   9/18/95       frt
! !   9/1/98 this version takes into account the the rcm record
! !           can differ to the record
! !     7/09 -restructured to use modules and allow to transfer
! !           the LFM grid - frt      
! !     2/19 -modified version to connect to gamera - frt
! !     5/20 -removed use of record numbers in rcm bookkeeping
! !          - also adds output from plasmasphere model by Shanshan Bao
! !     5/20/20 - removed idim,jdim -frt
! !
! !
! !==============================================================
! !
!       IMPLICIT NONE
!       type(rcm_mhd_T),intent(inout) :: RM
!       INTEGER(iprec), INTENT (OUT) :: ierr

!       INTEGER(iprec) :: i, j, L, k, itime,n,klow
!       real(rprec) :: max_xmin,min_xmin,max_ymin,min_ymin
!       real(rprec) :: dens_plasmasphere
!       real(rprec) :: sclmass(RCMNUMFLAV) !xmass prescaled to proton
!       real(rprec) :: kevRCM,kevMHD
!       real(rprec), parameter :: kRatMax = 0.9
!       !AMS 04-22-2020
!       real(rprec) :: pressure_factor,density_factor

!       !Array to hold eeta to ingest into MHD, should do this as pointer but requires adding "target" various places
!       real(rprec), dimension(isize,jsize,kcsize) :: eeta2mhd

!       !Pick which eeta to use, avg or final after advance
!       if (doAvg2MHD) then
!         eeta2mhd = max(eeta_avg,0.0)
!       else
!         eeta2mhd = max(eeta    ,0.0)
!       endif

!       pressure_factor = 2./3.*ev/RM%planet_radius*nt
!       density_factor = nt/RM%planet_radius

!       RM%MaxAlam = maxval(alamc)
      
!       ierr = 0

!       !Set lowest RC channel
!       if (use_plasmasphere) then
!         klow = 2
!       else
!         klow = 1
!       endif

!       !Set scaled mass by hand here to avoid precision issues
!       sclmass(RCMELECTRON) = mass_electron/mass_proton
!       sclmass(RCMPROTON) = 1.0

! !     Compute rcm pressure and density (on the ionospheric RCM grid):
!       !Tweaking scaling for better precision and testing eeta instead of avg (K: 8/20)
!       !$OMP PARALLEL DO default(shared) &
!       !$OMP schedule(dynamic) &
!       !$OMP private(i,j,k,dens_plasmasphere,kevRCM,kevMHD)
!       DO j = 1, jsize
!         DO i = 1, isize
!           pressrcm (i,j) = 0.0
!           densrcm  (i,j) = 0.0
!           denspsph (i,j) = 0.0

!           IF (vm(i,j) < 0.0) CYCLE
!           DO k = klow, kcsize
!             !Pressure calc in pascals
!             pressrcm(i,j) = pressrcm(i,j) + pressure_factor*ABS(alamc(k))*eeta2mhd(i,j,k)*vm(i,j)**2.5

!             !Density calc
!             if (alamc(k) > 0.0) then ! only add the ion contribution
!               densrcm(i,j) = densrcm(i,j) + density_factor*sclmass(ikflavc(k))*eeta2mhd(i,j,k)*vm(i,j)**1.5
!             endif
!           ENDDO !k

!           if (use_plasmasphere) then
!             if (dp_on) then 
!               ! use plasmasphere channel eeta_avg(:,:,1) sbao 03/2020
!               denspsph(i,j) = density_factor*sclmass(RCMPROTON)*eeta2mhd(i,j,1)*vm(i,j)**1.5
!             else
!               ! add a simple plasmasphere model based on carpenter 1992 or gallagher 2002 in ples/cc
!               dens_plasmasphere = GallagherXY(xmin(i,j),ymin(i,j))
!               denspsph(i,j) = dens_plasmasphere*1.0e6
!             endif !dp_on
!           endif !use_plasmasphere

!           !Do some checking on how well resolved energy channel stuff is
!           kevRCM = abs(alamc(kcsize))*vm(i,j)*1.0e-3 !Max keV of RCM channels here
!           kevMHD = DP2kT(densrcm(i,j)*rcmNScl,pressrcm(i,j)*rcmPScl) !Get keV from RCM moments
          
!           if (kevMHD >= kRatMax*kevRCM) then
!             !Effective "MHD" temperature, P=nkT_{MHD} is above kRatMax the max RCM channel energy
!             !This is probably bad for resolving the distribution so we do some shady cooling here

!             !Rescale eeta's to clamp P_{RCM} (doing directly to eeta, not copyed eeta2mhd)
!             eeta(i,j,klow:) = (kRatMax*kevRCM/kevMHD)*eeta(i,j,klow:)
!             pressrcm(i,j)   = (kRatMax*kevRCM/kevMHD)*pressrcm(i,j)
!             !write(*,*) 'kevRCM, kevMHD = ', kevRCM,kevMHD,densrcm(i,j)*rcmNScl,pressrcm(i,j)*rcmPScl
!           endif

!         ENDDO !i
!       ENDDO !j

!       if (any(isnan(densrcm)) .or. any(isnan(pressrcm))) then
!         write(*,*) 'NaN in RCM-DEN/P at tomhd'
!         stop
!       endif
 
!       max_xmin = maxval(xmin)
!       max_ymin = maxval(ymin)
!       min_xmin = minval(xmin)
!       min_ymin = minval(ymin)
 
! !     now update the pressure and density in the mhd code  
!       RM%Prcm    = pressrcm(:,jwrap:jsize)
!       RM%Nrcm    = densrcm (:,jwrap:jsize)
!       RM%Npsph   = denspsph(:,jwrap:jsize)
!       RM%flux    = eflux   (:,jwrap:jsize,:)
!       RM%eng_avg = eavg    (:,jwrap:jsize,:)
!       RM%fac     = birk    (:,jwrap:jsize)



!       ! At this point, we assume that RCM arrays are all populated
!       ! (plasma, b-field, bndloc, etc.):

!       IF (L_write_rcmu) then

!         rmin = SQRT (xmin**2+ymin**2+zmin**2)
!         WHERE (xmin == 0.0 .AND. ymin == 0.0)
!           pmin = 0.0
!         ELSE WHERE
!           pmin = ATAN2 (ymin, xmin)
!         END WHERE
!         WHERE (pmin < 0.0) pmin = pmin + 2.0*pi

!         rmin = SQRT (xmin**2+ymin**2)
!         WHERE (xmin == 0.0 .AND. ymin == 0.0)
!           pmin = 0.0
!         ELSE WHERE
!           pmin = ATAN2 (ymin, xmin)
!         END WHERE
!         WHERE (pmin < 0.0) pmin = pmin + 2.0*pi
!         WRITE (*,'(A,I9.9,A,I5.5)') 'TOMHD: Read RCM, T=',itime,', REC=',L

!         !call write_rcmu (L,0_iprec)

!       END IF

      
!       RETURN
!       END SUBROUTINE tomhd


END MODULE tomhd_mod