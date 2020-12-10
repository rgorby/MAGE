MODULE tomhd_mod
  USE rcm_precision
  USE rice_housekeeping_module
  USE rcm_mhd_interfaces
  USE kdefs, ONLY : TINY
  USE constants, ONLY : mass_proton,mass_electron,nt,ev,tiote,boltz
  USE Rcm_mod_subs, ONLY : isize, jsize, kcsize,jwrap,alamc,ikflavc
  USE earthhelper, ONLY : GallagherXY,DP2kT
  USE etautils
  
  implicit none

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

      REAL(rprec), dimension(isize,jsize) :: Pircm,Percm

      !Always set p/d_factors
      call SetFactors(RM%planet_radius)

      !Do some checks on the RCM eta distribution
      if (doRelax) call RelaxEta(eeta,eeta_avg,vm,RM)

      !Now pick which eta (instant vs. avg) and calculate moments
      if (doAvg2MHD) then
        call rcm2moments(eeta_avg,vm,densrcm,denspsph,pressrcm,Pircm,Percm)
      else
        call rcm2moments(eeta    ,vm,densrcm,denspsph,pressrcm,Pircm,Percm)
      endif

      RM%MaxAlam = maxval(alamc)
      
!     now update the pressure and density in the mhd code  
      RM%Prcm    = pressrcm(:,jwrap:jsize)
      RM%Nrcm    = densrcm (:,jwrap:jsize)
      RM%Npsph   = denspsph(:,jwrap:jsize)
      RM%flux    = eflux   (:,jwrap:jsize,:)
      RM%eng_avg = eavg    (:,jwrap:jsize,:)
      RM%fac     = birk    (:,jwrap:jsize)
      RM%Percm   = Percm   (:,jwrap:jsize)
      ierr = 0

    END SUBROUTINE tomhd

    !Do some safety stuff to eta w/ temperature is high
    SUBROUTINE RelaxEta(eta,eta_avg,vm,RCMApp)
      USE Rcm_mod_subs, ONLY : rmin
      IMPLICIT NONE

      REAL(rprec), intent(inout), dimension(isize,jsize,kcsize)  :: eta,eta_avg
      REAL(rprec), intent(in)  :: vm(isize,jsize)
      type(rcm_mhd_T),intent(in) :: RCMApp

      REAL(rprec), dimension(isize,jsize) :: Drc,Dpp,Prc,Lb,Tb,Pion,Pele
      integer :: i,j,jp,klow,k
      REAL(rprec), dimension(kcsize) :: etaMax,etaNew,etaOld
      REAL(rprec) :: TauDP,wDP,wgt
      
      !Set lowest RC channel
      if (use_plasmasphere) then
        klow = 2
      else
        klow = 1
      endif

      !Get moments from eta
      call rcm2moments(eta,vm,Drc,Dpp,Prc,Pion,Pele)

      !Map Lb from RCM-MHD grid to RCM grid, Lb = Tube length [m]
      call EmbiggenWrap(RCMApp%Lb,Lb)
      !Convert to km
      Lb = RCMApp%planet_radius*Lb*1.0e-3

      call EmbiggenWrap(RCMApp%Tb,Tb)

    !Now loop over i,j and relax in energy space by blending w/ Maxwellian
      !Get weights for Maxwellian part of blend
      !eta_{R} = w * eta_{Maxwellian} + (1-w) * eta
      !Thermal bounce: w = dtCpl/tau, tau = sound wave bounce period
      !Grid bounce: w = kT_{RCM}/kT_{Top}, ratio of effective temperature to largest grid energy

      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j,jp,TauDP,wDP,wgt,k) &
      !$OMP private(etaMax,etaNew,etaOld)
      DO j = 1, jsize
        !i,j is index in RCM grid
        !i,jp is index in RCM-MHD grid
        if (j>=jwrap) then
          jp = j-jwrap+1
        else
          !j (RCM) => jsize-jwrap+j (RCM) => jsize-jwrap+j-jwrap+1 (RCM-MHD)
          jp = jsize-jwrap+j-jwrap+1
        endif

        DO i = 1, isize
          IF (vm (i,j) < 0.0) CYCLE
          IF (Drc(i,j) < TINY) CYCLE

        !Get Maxwellian to blend with
          call DPP2eta(Drc(i,j),Pion(i,j),Pele(i,j),vm(i,j),etaMax,doRescaleO=.true.)
        !Get timescale to blend over
          TauDP = DriftPeriod(Drc(i,j),Prc(i,j),rmin(i,j),RCMApp%Bmin(i,jp),RCMApp%radcurv(i,jp),RCMApp%planet_radius)
          wDP = RCMApp%dtCpl/TauDP !Drift period
          call ClampWeight(wDP)
          !Choose which weight to use
          wgt = wDP 

        !Now do blending
          etaOld = eta(i,j,:)
          do k=1,kcsize
            if (alamc(k) > TINY) then
              etaNew(k) = (1-wgt)*etaOld(k) + wgt*etaMax(k)
            else
              etaNew(k) = etaOld(k)
            endif
          enddo
          eta(i,j,:) = etaNew

        ENDDO
      ENDDO !j loop

      contains

        !L [Rp], Rc [Rp], Rp[m]
        function DriftPeriod(n,P,L,Bmin,Rc,Rp) result(TauD)
          REAL(rprec), intent(in) :: n,P,L,Bmin,Rc,Rp
          REAL(rprec) :: TauD
          REAL(rprec) :: keV,Bnt,Vd

          keV = DP2kT(n*rcmNScl,P*rcmPScl) !Temp in keV
          Bnt = Bmin*1.0e+9 !B in nT

          !Using, Vd = 156*K/Rc/B, K[keV],Rc[Re],B[nT]

          Vd = (1.0e+3)*156.0*keV/Rc/Bnt !m/s
          TauD = 2.0*PI*L*Rp/Vd
          !write(*,*) 'Tau = ', TauD,DipoleDriftPeriod(n,P,L)
        end function DriftPeriod

        function DipoleDriftPeriod(n,P,L) result(TauD)
          REAL(rprec), intent(in) :: n,P,L
          REAL(rprec) :: TauD
          REAL(rprec) :: keV

          keV = DP2kT(n*rcmNScl,P*rcmPScl) !Temp in keV
          !Using, Td = 700/K/L hrs
          TauD = (60.0*60.0)*700/(keV*L)
        end function DipoleDriftPeriod

        !Return sound wave bounce period [s], take n/P in RCM units and L [km]
        function CsBounce(n,P,L) result(TauCS)
          REAL(rprec), intent(in) :: n,P,L
          REAL(rprec) :: TauCS
          REAL(rprec) :: TiEV,CsMKS

          integer, parameter :: nBounce = 8 !Number of bounces to equilibrate

          TiEV = (1.0e+3)*DP2kT(n*rcmNScl,P*rcmPScl) !Temp in eV
          !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
          CsMKS = 9.79*sqrt((5.0/3)*TiEV)
          TauCS = 2*nBounce*L/CsMKS
        end function CsBounce

        function GridWeight(n,P,vm,alamax) result(wgt)
          REAL(rprec), intent(in) :: n,P,vm,alamax
          REAL(rprec) :: wgt
          REAL(rprec) :: kevMHD,kevRCM

          kevMHD = DP2kT(n*rcmNScl,P*rcmPScl) !Get keV from RCM moments
          kevRCM = abs(alamax)*vm*1.0e-3 !keV of max RCM energy channel

          wgt = kevMHD/kevRCM
          call ClampWeight(wgt)
        end function GridWeight

        !Clamps weight in [0,1]
        subroutine ClampWeight(wgt)
          REAL(rprec), intent(inout) :: wgt
          if (wgt<0.0) wgt = 0.0
          if (wgt>1.0) wgt = 1.0
        end subroutine ClampWeight

    END SUBROUTINE RelaxEta


    !Convert given eeta to density (RC/plasmasphere) and pressure
    SUBROUTINE rcm2moments(eta,vm,Drc,Dpp,Prc,Pion,Pele)
      USE Rcm_mod_subs, ONLY : xmin,ymin,zmin
      IMPLICIT NONE
      REAL(rprec), intent(in)  :: eta(isize,jsize,kcsize)
      REAL(rprec), intent(in)  :: vm(isize,jsize)
      REAL(rprec), intent(out), dimension(isize,jsize) :: Drc,Dpp,Prc
      REAL(rprec), intent(out), dimension(isize,jsize), optional :: Pion,Pele
      INTEGER (iprec) :: i,j
      LOGICAL :: doIE

      Drc = 0.0
      Dpp = 0.0
      Prc = 0.0
      if (present(Pion) .and. present(Pele)) then
        Pion = 0.0
        Pele = 0.0
        doIE = .true.
      else
        doIE = .false.
      endif

      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j)
      DO j = 1, jsize
        DO i = 1, isize

          call eta2DP(eta(i,j,:),vm(i,j),Drc(i,j),Dpp(i,j),Prc(i,j))
          if (doIE) then
            call IntegratePressureIE(eta(i,j,:),vm(i,j),Pion(i,j),Pele(i,j)) !Get separated pressures
          endif
        ENDDO
      ENDDO !J loop

    END SUBROUTINE rcm2moments

END MODULE tomhd_mod