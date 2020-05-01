module mixconductance
  use mixdefs
  use mixtypes
  use earthhelper
  use gcmtypes
  use gcminterp
  
  implicit none

  real(rp), dimension(:,:), allocatable :: tmpD,tmpC ! used for chilling in Fedder95. Declare it here so we can allocate in init.

  contains
    subroutine conductance_init(conductance,Params,G)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixParams_T), intent(in) :: Params
      type(mixGrid_T), intent(in) :: G

      ! define these module-wide variables so we don't have to pass the Params object to all the conductance functions
      ! this is similar to how it was done in MIX
      conductance%euv_model_type = Params%euv_model_type
      conductance%et_model_type = Params%et_model_type
      conductance%alpha = Params%alpha
      conductance%beta = Params%beta
      conductance%R = Params%R
      conductance%F107 = Params%F107
      conductance%pedmin = Params%pedmin
      conductance%hallmin = Params%hallmin
      conductance%sigma_ratio = Params%sigma_ratio
      conductance%ped0 = Params%ped0
      conductance%const_sigma = Params%const_sigma
      conductance%doRamp = Params%doRamp
      conductance%doChill = Params%doChill
      conductance%doStarlight = Params%doStarlight      
      conductance%doMR = Params%doMR      
      conductance%apply_cap = Params%apply_cap
      conductance%doGCM = Params%doGCM
      conductance%doGCM2way = Params%doGCM2way

      conductance%PI2       = pi/2.0D0
      conductance%ang65     = pi/180.0D0*65.0D0
      conductance%ang100    = pi*5.0D0/9.0D0
      conductance%pref      = 2.0D0*250.0D0**(-0.666666)
      conductance%href      = 1.0D0/(1.8D0*sqrt(250.0D0))
      conductance%shall     = 1.8D0*sqrt(conductance%f107)
      conductance%speder    = 0.5D0*conductance%f107**(0.666666)
      conductance%pedslope  = 0.24D0*conductance%pref*conductance%speder*rad2deg
      conductance%pedslope2 = 0.13D0*conductance%pref*conductance%speder*rad2deg
      conductance%hallslope = 0.27D0*conductance%href*conductance%shall*rad2deg;
      conductance%sigmap65  = conductance%speder*cos(conductance%ang65)**0.666666
      conductance%sigmah65  = conductance%shall*cos(conductance%ang65)
      conductance%sigmap100 = conductance%sigmap65-(conductance%ang100-conductance%ang65)*conductance%pedslope

      if (.not. allocated(conductance%zenith)) allocate(conductance%zenith(G%Np,G%Nt))
      if (.not. allocated(conductance%coszen)) allocate(conductance%coszen(G%Np,G%Nt))
      if (.not. allocated(conductance%euvSigmaP)) allocate(conductance%euvSigmaP(G%Np,G%Nt))
      if (.not. allocated(conductance%euvSigmaH)) allocate(conductance%euvSigmaH(G%Np,G%Nt))
      if (.not. allocated(conductance%deltaSigmaP)) allocate(conductance%deltaSigmaP(G%Np,G%Nt))
      if (.not. allocated(conductance%deltaSigmaH)) allocate(conductance%deltaSigmaH(G%Np,G%Nt))

      if (.not. allocated(conductance%rampFactor)) allocate(conductance%rampFactor(G%Np,G%Nt))
      if (.not. allocated(conductance%ares)) allocate(conductance%ares(G%Np,G%Nt))
      if (.not. allocated(conductance%deltaE)) allocate(conductance%deltaE(G%Np,G%Nt))
      if (.not. allocated(conductance%E0)) allocate(conductance%E0(G%Np,G%Nt))
      if (.not. allocated(conductance%phi0)) allocate(conductance%phi0(G%Np,G%Nt))
      if (.not. allocated(conductance%engFlux)) allocate(conductance%engFlux(G%Np,G%Nt))

      if (.not. allocated(tmpD)) allocate(tmpD(G%Np,G%Nt))
      if (.not. allocated(tmpC)) allocate(tmpC(G%Np,G%Nt))  

      write(*,*) "doGCM?",conductance%doGCM
      write(*,*) "doGCM2way?",conductance%doGCM2way    

    end subroutine conductance_init

    subroutine conductance_euv(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      conductance%zenith = conductance%PI2 - ( asin(G%x) + St%tilt )
      ! An alternative (correct) definition of the zenith angle (for Moen-Brekke)
      conductance%coszen = G%x*cos(St%tilt)+sqrt(1.-G%x**2-G%y**2)*sin(St%tilt) ! as it should be
      conductance%zenith = acos(conductance%coszen)

      select case ( conductance%euv_model_type )
         case (AMIE)
            where (conductance%zenith <= conductance%ang65) 
               conductance%euvSigmaP = conductance%speder*cos(conductance%zenith)**0.666666
               conductance%euvSigmaH = conductance%shall*cos(conductance%zenith)
            elsewhere (conductance%zenith <= conductance%ang100)
               conductance%euvSigmaP = conductance%sigmap65 - conductance%pedslope*(conductance%zenith - conductance%ang65)
               conductance%euvSigmaH = conductance%sigmah65 - conductance%hallslope*(conductance%zenith - conductance%ang65)
            elsewhere (conductance%zenith > conductance%ang100)
               conductance%euvSigmaP = conductance%sigmap100 - conductance%pedslope2*(conductance%zenith-conductance%ang100)
               conductance%euvSigmaH = conductance%sigmah65 - conductance%hallslope*(conductance%zenith-conductance%ang65)
            end where
         case (MOEN_BREKKE) !!! Needs testing
            ! This works only for the dayside (zenith <= pi/2)
            ! Set it to pedMin (hallMin) on the nightside (may rethink it later)
            where (conductance%coszen >=0.0) 
!               euvSigmaP = f107**0.49*( 0.34*cos(zenith)+0.93*sqrt(cos(zenith)) )
!               euvSigmaH = f107,0.53*( 0.81*cos(zenith)+0.54*sqrt(cos(zenith)) )
               conductance%euvSigmaP = conductance%f107**0.49*( 0.34*conductance%coszen+0.93*sqrt(conductance%coszen) )
               conductance%euvSigmaH = conductance%f107**0.53*( 0.81*conductance%coszen+0.54*sqrt(conductance%coszen) )
            elsewhere
               conductance%euvSigmaP = conductance%pedmin
               conductance%euvSigmaH = conductance%hallmin
            end where
         case default
            stop "The EUV model type entered is not supported."
      end select

      if (conductance%doStarlight) then ! Makes sense to turn off apply_cap if this is on
         ! add background conductance quadratically instead of a sharp cutoff
         ! first need to remove negative conductances resulting form AMIE for strong tilt
         ! otherwise, the quadratic addition results in conductances growing toward nightside.
         conductance%euvSigmaP = max(conductance%euvSigmaP,0.)
         conductance%euvSigmaH = max(conductance%euvSigmaH,0.)
         
         conductance%euvSigmaP = sqrt(conductance%euvSigmaP**2 + conductance%pedmin**2)
         conductance%euvSigmaH = sqrt(conductance%euvSigmaH**2 + conductance%hallmin**2)
      else
         ! otherwise, default to the standard way of applying the floor
         ! I'll leave it as default for now (doStarlight=.false.) but we can reconsider later
         conductance%euvSigmaP = max(conductance%euvSigmaP,conductance%pedmin)
         conductance%euvSigmaH = max(conductance%euvSigmaH,conductance%hallmin)
      end if
      
    end subroutine conductance_euv


    subroutine conductance_fedder95(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
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
      if (conductance%doRamp) then
         call conductance_ramp(conductance,G,20.0D0*pi/180.D0,30.0D0*pi/180.0D0,0.02D0)
      else 
         conductance%rampFactor = 1.0D0
      end if

      if (conductance%doChill) then
         ! MHD density replaced with gallagher where it's lower      
         ! and temperature changed correspondingly
          tmpD = max(G%D0*Mp_cgs,St%Vars(:,:,DENSITY))
          tmpC = St%Vars(:,:,SOUND_SPEED)*sqrt(St%Vars(:,:,DENSITY)/tmpD)
       else
          tmpD = St%Vars(:,:,DENSITY)
          tmpC = St%Vars(:,:,SOUND_SPEED)
      end if
      
      conductance%E0 = conductance%alpha*Mp_cgs*heFrac*erg2kev*tmpC**2*conductance%RampFactor
      conductance%phi0 = sqrt(kev2erg)/(heFrac*Mp_cgs)**1.5D0*conductance%beta*tmpD*sqrt(conductance%E0)*conductance%RampFactor
      ! resistence out of the ionosphere is 2*rout resistence into the
      ! ionosphere is 2*rin outward current is positive
      where ( signOfJ*St%Vars(:,:,FAC) >=0. ) 
         conductance%aRes = 2.D0*Rout
      elsewhere
         conductance%aRes = 2.D0*Rin
      end where
   
      ! Density floor to limit characteristic energy.  See Wiltberger et al. 2009 for details.
      where (tmpD < rhoFactor*conductance%euvSigmaP) 
         tmpD = rhoFactor*conductance%euvSigmaP
      end where
      conductance%deltaE = (heFrac*Mp_cgs)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*conductance%R*conductance%aRes*signOfJ*(St%Vars(:,:,FAC)*1.e-6)*sqrt(conductance%E0)/tmpD

      ! limit the max potential energy drop to 20 [keV]
      conductance%deltaE = min(20.D0,conductance%deltaE)

      ! floor on total energy
      St%Vars(:,:,AVG_ENG) = max(conductance%E0 + conductance%deltaE,1.D-8)

      where  ( conductance%deltaE > 0. )
         St%Vars(:,:,NUM_FLUX) = conductance%phi0*(8.D0-7.D0*exp(-conductance%deltaE/7.D0/conductance%E0))
      elsewhere 
         St%Vars(:,:,NUM_FLUX) = conductance%phi0*exp(conductance%deltaE/conductance%E0)
      end where

    end subroutine conductance_fedder95

    subroutine conductance_aurora(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      ! note, this assumes that fedder has been called prior
      conductance%engFlux = kev2erg*St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)  ! Energy flux in ergs/cm^2/s
      conductance%deltaSigmaP = 40.D0*St%Vars(:,:,AVG_ENG)*sqrt(conductance%engFlux)/(16.D0+St%Vars(:,:,AVG_ENG)**2);
      conductance%deltaSigmaH = 0.45D0*conductance%deltaSigmaP*St%Vars(:,:,AVG_ENG)**0.85D0/(1.D0+0.0025D0*St%Vars(:,:,AVG_ENG)**2)

      ! correct for multiple reflections if you're so inclined
      if (conductance%doMR) call conductance_mr(conductance,St)
    end subroutine conductance_aurora

    ! George Khazanov's multiple reflection(MR) corrections
    subroutine conductance_mr(conductance,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixState_T), intent(in) :: St

      conductance%deltaSigmaP = (2.16-0.87*exp(-0.16*St%Vars(:,:,AVG_ENG)))*conductance%deltaSigmaP
      conductance%deltaSigmaH = (1.87-0.54*exp(-0.16*St%Vars(:,:,AVG_ENG)))*conductance%deltaSigmaH
    end subroutine conductance_mr

    subroutine conductance_total(conductance,G,St,gcm,h)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      type(gcm_T),optional,intent(in) :: gcm
      integer,optional,intent(in) :: h

      ! always call fedder to fill in AVG_ENERGY and NUM_FLUX
      ! even if const_sigma, we still have the precip info that way

      ! compute EUV though because it's used in fedder
      call conductance_euv(conductance,G,St)
      call conductance_fedder95(conductance,G,St)
      
      if (present(gcm)) then
         write(*,*) 'going to apply!'
         call apply_gcm2mix(gcm,St,h)
         St%Vars(:,:,SIGMAP) = max(conductance%pedmin,St%Vars(:,:,SIGMAP))
         !St%Vars(:,:,SIGMAH) = min(max(conductance%hallmin,St%Vars(:,:,SIGMAH)),&
         !     St%Vars(:,:,SIGMAP)*conductance%sigma_ratio)
      else if (conductance%const_sigma) then
         St%Vars(:,:,SIGMAP) = conductance%ped0
         St%Vars(:,:,SIGMAH) = 0.D0
      else
         call conductance_aurora(conductance,G,St)
      
         St%Vars(:,:,SIGMAP) = sqrt( conductance%euvSigmaP**2 + conductance%deltaSigmaP**2) 
         St%Vars(:,:,SIGMAH) = sqrt( conductance%euvSigmaH**2 + conductance%deltaSigmaH**2)
      endif

      ! Apply cap
      if ((conductance%apply_cap).and.(.not. conductance%const_sigma).and.(.not. present(gcm))) then
         St%Vars(:,:,SIGMAP) = max(conductance%pedmin,St%Vars(:,:,SIGMAP))
         St%Vars(:,:,SIGMAH) = min(max(conductance%hallmin,St%Vars(:,:,SIGMAH)),&
              St%Vars(:,:,SIGMAP)*conductance%sigma_ratio)
      endif

    end subroutine conductance_total

    subroutine conductance_ramp(conductance,G,rPolarBound,rEquatBound,rLowLimit)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G      
      real(rp), intent(in) :: rPolarBound,rEquatBound,rLowLimit
      
      where (G%r < rPolarBound)
         conductance%rampFactor = 1.0D0
      elsewhere ( (G%r > rPolarBound).and.(G%r <= rEquatBound) )
         conductance%rampFactor = 1.0D0+(rLowLimit - 1.0D0)*(G%r - rPolarBound)/(rEquatBound-rPolarBound)
      elsewhere ! G%r > rEquatBound
         conductance%rampFactor = rLowLimit
      end where
    end subroutine conductance_ramp

  end module mixconductance
