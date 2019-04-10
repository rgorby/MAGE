module mixconductance
  use mixdefs
  use mixtypes

  implicit none

  integer :: euv_model_type, et_model_type
  real(rp) :: alpha, beta, R, F107,pedmin,hallmin,sigma_ratio,ped0
  logical :: const_sigma, do_ramp, apply_cap

  ! auxilary variables
  real(rp) :: PI2, ang65, ang100, pref, href, shall
  real(rp) :: speder, pedslope, pedslope2, hallslope,sigmap65, sigmah65, sigmap100
  real(rp), dimension(:,:), allocatable :: zenith, coszen
  real(rp), dimension(:,:), allocatable :: euvSigmaP, euvSigmaH
  real(rp), dimension(:,:), allocatable :: deltaSigmaP, deltaSigmaH
  real(rp), dimension(:,:), allocatable :: E0, phi0, deltaE, aRes
  real(rp), dimension(:,:), allocatable :: rampFactor
  real(rp), dimension(:,:), allocatable :: engFlux

  contains
    subroutine conductance_init(Params,G)
      type(mixParams_T), intent(in) :: Params
      type(mixGrid_T), intent(in) :: G

      ! define these module-wide variables so we don't have to pass the Params object to all the conductance functions
      ! this is similar to how it was done in MIX
      euv_model_type = Params%euv_model_type
      et_model_type = Params%et_model_type
      alpha = Params%alpha
      beta = Params%beta
      R = Params%R
      F107 = Params%F107
      pedmin = Params%pedmin
      hallmin = Params%hallmin
      sigma_ratio = Params%sigma_ratio
      ped0 = Params%ped0
      const_sigma = Params%const_sigma
      do_ramp = Params%do_ramp
      apply_cap = Params%apply_cap

      PI2       = pi/2.0D0
      ang65     = pi/180.0D0*65.0D0
      ang100    = pi*5.0D0/9.0D0
      pref      = 2.0D0*250.0D0**(-0.666666)
      href      = 1.0D0/(1.8D0*sqrt(250.0D0))
      shall     = 1.8D0*sqrt(f107)
      speder    = 0.5D0*f107**(0.666666)
      pedslope  = 0.24D0*pref*speder*rad2deg
      pedslope2 = 0.13D0*pref*speder*rad2deg
      hallslope = 0.27D0*href*shall*rad2deg;
      sigmap65  = speder*cos(ang65)**0.666666
      sigmah65  = shall*cos(ang65)
      sigmap100 = sigmap65-(ang100-ang65)*pedslope

      if (.not. allocated(zenith)) allocate(zenith(G%Np,G%Nt))
      if (.not. allocated(coszen)) allocate(coszen(G%Np,G%Nt))
      if (.not. allocated(euvSigmaP)) allocate(euvSigmaP(G%Np,G%Nt))
      if (.not. allocated(euvSigmaH)) allocate(euvSigmaH(G%Np,G%Nt))
      if (.not. allocated(deltaSigmaP)) allocate(deltaSigmaP(G%Np,G%Nt))
      if (.not. allocated(deltaSigmaH)) allocate(deltaSigmaH(G%Np,G%Nt))

      if (.not. allocated(rampFactor)) allocate(rampFactor(G%Np,G%Nt))
      if (.not. allocated(ares)) allocate(ares(G%Np,G%Nt))
      if (.not. allocated(deltaE)) allocate(deltaE(G%Np,G%Nt))
      if (.not. allocated(E0)) allocate(E0(G%Np,G%Nt))
      if (.not. allocated(phi0)) allocate(phi0(G%Np,G%Nt))
      if (.not. allocated(engFlux)) allocate(engFlux(G%Np,G%Nt))
    end subroutine conductance_init

    subroutine conductance_euv(G,St)
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      zenith = PI2 - ( asin(G%x) + St%tilt )
      ! An alternative (correct) definition of the zenith angle (for Moen-Brekke)
      coszen = G%x*cos(St%tilt)+sqrt(1.-G%x**2-G%y**2)*sin(St%tilt) ! as it should be
      zenith = acos(coszen)

      select case ( euv_model_type )
         case (AMIE)
            where (zenith <= ang65) 
               euvSigmaP = speder*cos(zenith)**0.666666
               euvSigmaH = shall*cos(zenith)
            elsewhere (zenith <= ang100)
               euvSigmaP = sigmap65 - pedslope*(zenith - ang65)
               euvSigmaH = sigmah65 - hallslope*(zenith - ang65)
            elsewhere (zenith > ang100)
               euvSigmaP = sigmap100 - pedslope2*(zenith-ang100)
               euvSigmaH = sigmah65 - hallslope*(zenith-ang65)
            end where
         case (MOEN_BREKKE) !!! Needs testing
            ! This works only for the dayside (zenith <= pi/2)
            ! Set it to pedMin (hallMin) on the nightside (may rethink it later)
            where (coszen >=0.0) 
!               euvSigmaP = f107**0.49*( 0.34*cos(zenith)+0.93*sqrt(cos(zenith)) )
!               euvSigmaH = f107,0.53*( 0.81*cos(zenith)+0.54*sqrt(cos(zenith)) )
               euvSigmaP = f107**0.49*( 0.34*coszen+0.93*sqrt(coszen) )
               euvSigmaH = f107**0.53*( 0.81*coszen+0.54*sqrt(coszen) )
            elsewhere
               euvSigmaP = pedmin
               euvSigmaH = hallmin
            end where
         case default
            stop "The EUV model type entered is not supported."
      end select

!      where (euvSigmaP<pedmin) euvSigmaP=pedmin
!      where (euvSigmaH<hallmin) euvSigmaH=hallmin
      euvSigmaP = max(euvSigmaP,pedmin)
      euvSigmaH = max(euvSigmaH,hallmin)
    end subroutine conductance_euv


    subroutine conductance_fedder95(G,St)
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      
      real(rp) :: Redge, Rmin, Rmin2, Rmax, rfac
      real(rp) :: signOfY, signOfJ
      real(rp) :: Rout = 6.D0, Rin = 1.2D0
      real(rp) :: rhoFactor = 3.3D-24*0.5D0

      if (St%hemisphere==NORTH) then
         signOfY = -1  ! note, factor2 (dawn-dusk asymmetry is not
                           ! implemented since factor2 in the old
                           ! fedder95 code was removed, i.e., set to
                           ! 1, anyway). I think Mike did this when he
                           ! implemented his ramp function.
         signOfJ = -1  
      elseif (St%hemisphere==SOUTH) then
         signOfY = 1
         signOfJ = 1
      else
         stop 'Wrong hemisphere label. Stopping...'
      endif

      ! fills in rampFactor
      if (do_ramp) then
         call conductance_ramp(G,20.0D0*pi/180.D0,30.0D0*pi/180.0D0,0.02D0)
      else 
         rampFactor = 1.0D0
      end if

      E0 = alpha*mp*heFrac*erg2kev*St%Vars(:,:,SOUND_SPEED)**2*RampFactor
      phi0 = sqrt(kev2erg)/(heFrac*mp)**1.5D0*beta*St%Vars(:,:,DENSITY)*sqrt(E0)*RampFactor

      ! resistence out of the ionosphere is 2*rout resistence into the
      ! ionosphere is 2*rin outward current is positive
      where ( signOfJ*St%Vars(:,:,FAC) >=0. ) 
         aRes = 2.D0*Rout
      elsewhere
         aRes = 2.D0*Rin
      end where

      ! Density floor to limit characteristic energy.  See Wiltberger et al. 2009 for details.
      where (St%Vars(:,:,DENSITY) < rhoFactor*euvSigmaP) 
         St%Vars(:,:,DENSITY) = rhoFactor*euvSigmaP
      end where
      deltaE = (heFrac*mp)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*R*aRes*signOfJ*(St%Vars(:,:,FAC)*1.e-6)*sqrt(E0)/St%Vars(:,:,DENSITY)

      ! limit the max potential energy drop to 20 [keV]
      deltaE = min(20.D0,deltaE)

      ! floor on total energy
      St%Vars(:,:,AVG_ENG) = max(E0 + deltaE,1.D-8)

      where  ( deltaE > 0. )
         St%Vars(:,:,NUM_FLUX) = phi0*(8.D0-7.D0*exp(-deltaE/7.D0/E0))
      elsewhere 
         St%Vars(:,:,NUM_FLUX) = phi0*exp(deltaE/E0)
      end where
    end subroutine conductance_fedder95

    subroutine conductance_aurora(G,St)
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      ! first call fedder to fill in AVG_ENERGY and NUM_FLUX
      call conductance_fedder95(G,St)

      engFlux = kev2erg*St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)  ! Energy flux in ergs/cm^2/s
      deltaSigmaP = 40.D0*St%Vars(:,:,AVG_ENG)*sqrt(engFlux)/(16.D0+St%Vars(:,:,AVG_ENG)**2);
      deltaSigmaH = 0.45D0*deltaSigmaP*St%Vars(:,:,AVG_ENG)**0.85D0/(1.D0+0.0025D0*St%Vars(:,:,AVG_ENG)**2)
    end subroutine conductance_aurora

    subroutine conductance_total(G,St)
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      if (const_sigma) then
         St%Vars(:,:,SIGMAP) = ped0
         St%Vars(:,:,SIGMAH) = 0.D0
      else
         call conductance_euv(G,St)
         call conductance_aurora(G,St)
      
         St%Vars(:,:,SIGMAP) = sqrt( euvSigmaP**2 + deltaSigmaP**2) 
         St%Vars(:,:,SIGMAH) = sqrt( euvSigmaH**2 + deltaSigmaH**2)
      endif

      ! Apply cap
      if ((apply_cap).and.(.not.const_sigma)) then
         St%Vars(:,:,SIGMAP) = max(pedmin,St%Vars(:,:,SIGMAP))
         St%Vars(:,:,SIGMAH) = min(max(hallmin,St%Vars(:,:,SIGMAH)),&
              St%Vars(:,:,SIGMAP)*sigma_ratio)
       endif
    end subroutine conductance_total

    subroutine conductance_ramp(G,rPolarBound,rEquatBound,rLowLimit)
      type(mixGrid_T), intent(in) :: G      
      real(rp), intent(in) :: rPolarBound,rEquatBound,rLowLimit
      
      where (G%r < rPolarBound)
         rampFactor = 1.0D0
      elsewhere ( (G%r > rPolarBound).and.(G%r <= rEquatBound) )
         rampFactor = 1.0D0+(rLowLimit - 1.0D0)*(G%r - rPolarBound)/(rEquatBound-rPolarBound)
      elsewhere ! G%r > rEquatBound
         rampFactor = rLowLimit
      end where
    end subroutine conductance_ramp
    
end module mixconductance
