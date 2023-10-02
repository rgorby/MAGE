module mixconductance
  use mixdefs
  use mixtypes
  use planethelper
  use gcmtypes
  use gcminterp
  use math
  use euvhelper
  use auroralhelper
  use kai2geo
  use rcmdefs, ONLY : tiote_RCM
  
  implicit none

  real(rp), dimension(:,:), allocatable, private :: tmpD,tmpC ! used for chilling in Fedder95. Declare it here so we can allocate in init.
  real(rp), dimension(:,:), allocatable, private :: JF0,RM,RRdi ! used for zhang15
  real(rp), dimension(:,:), allocatable, private :: tmpE,tmpF ! used for smoothing precipitation avg_eng and num_flux
  real(rp), dimension(:,:), allocatable, private :: beta_RCM,alpha_RCM,gtype_RCM ! two-dimensional beta based on RCM fluxes.

  !Replacing some hard-coded inline values (bad) w/ module private values (slightly less bad)
  real(rp), parameter, private :: maxDrop = 20.0 !Hard-coded max potential drop [kV]
  real(rp), parameter, private :: eTINY = 1.D-8 ! Floor of average energy [keV]
  real(rp), private :: RinMHD = 0.0 !Rin of MHD grid (0 if not running w/ MHD)
  real(rp), private :: MIXgamma
  logical , private :: doDrift = .false. !Whether to add drift term from Zhang

  contains

    subroutine conductance_init(conductance,Params,G)
      ! Initialize precipitation and conductance variables.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixParams_T)     , intent(in)    :: Params
      type(mixGrid_T)       , intent(in)    :: G

      ! define these module-wide variables so we don't have to pass the Params object to all the conductance functions
      ! this is similar to how it was done in MIX
      conductance%euv_model_type    = Params%euv_model_type
      conductance%et_model_type     = Params%et_model_type
      conductance%alpha             = Params%alpha
      conductance%beta              = Params%beta
      conductance%R                 = Params%R
      conductance%F107              = Params%F107
      conductance%pedmin            = Params%pedmin
      conductance%hallmin           = Params%hallmin
      conductance%sigma_ratio       = Params%sigma_ratio
      conductance%ped0              = Params%ped0
      conductance%const_sigma       = Params%const_sigma
      conductance%doRamp            = Params%doRamp
      conductance%doChill           = Params%doChill
      conductance%doStarlight       = Params%doStarlight      
      conductance%doMR              = Params%doMR      
      conductance%doAuroralSmooth   = Params%doAuroralSmooth      
      conductance%apply_cap         = Params%apply_cap
      conductance%aurora_model_type = Params%aurora_model_type

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

      if (.not. allocated(conductance%avgEng)) allocate(conductance%avgEng(G%Np,G%Nt))
      if (.not. allocated(conductance%drift)) allocate(conductance%drift(G%Np,G%Nt))      
      if (.not. allocated(conductance%AuroraMask)) allocate(conductance%AuroraMask(G%Np,G%Nt))      
      if (.not. allocated(conductance%PrecipMask)) allocate(conductance%PrecipMask(G%Np,G%Nt))    

      ! these arrays are global and should not be! reallocate them
      if(allocated(tmpD)) deallocate(tmpD)
      if(allocated(tmpC)) deallocate(tmpC)
      if(allocated(JF0)) deallocate(JF0)
      if(allocated(RM)) deallocate(RM)
      if(allocated(RRdi)) deallocate(RRdi)
      if(allocated(tmpE)) deallocate(tmpE)
      if(allocated(tmpF)) deallocate(tmpF)
      if(allocated(beta_RCM)) deallocate(beta_RCM)
      if(allocated(alpha_RCM)) deallocate(alpha_RCM)
      if(allocated(gtype_RCM)) deallocate(gtype_RCM)

      allocate(tmpD(G%Np,G%Nt))
      allocate(tmpC(G%Np,G%Nt))  
      allocate(JF0(G%Np,G%Nt))      
      allocate(RM(G%Np,G%Nt))      
      allocate(RRdi(G%Np,G%Nt))      
      allocate(tmpE(G%Np+4,G%Nt+4)) ! for boundary processing.
      allocate(tmpF(G%Np+4,G%Nt+4))
      allocate(beta_RCM(G%Np,G%Nt))
      allocate(alpha_RCM(G%Np,G%Nt))
      allocate(gtype_RCM(G%Np,G%Nt))

      call SetMIXgamma(Params%gamma)
      ! MHD inner boundary, used to calc mirror ratio.
      RinMHD = Params%RinMHD
      ! Te/Tmhd
      alpha_RCM = 1.0/(tiote_RCM+1.0)
      ! Loss cone rate
      beta_RCM  = conductance%beta
      ! RCM grid weight: 1. Totally on closed RCM; 0. Totally outside RCM.
      ! if conductance_IM_GTYPE is not called, gtype_RCM has all zero, MHD values have a weight of 1.
      gtype_RCM = 0.0

    end subroutine conductance_init

    subroutine conductance_total(conductance,G,St,gcm,h)
      ! The main subroutine in module mixconductance.
      ! It derives auroral precipitation and conductance caused by both solar EUV and aurora.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      type(gcm_T),optional,intent(in) :: gcm
      integer,optional,intent(in) :: h

      call GenMirrorRatio(G,St)

      ! Compute EUV though because it's used in fedder
      call conductance_euv(conductance,G,St)
      select case ( conductance%aurora_model_type )
         case (FEDDER)
            call conductance_fedder95(conductance,G,St)
         case (ZHANG)
            doDrift = .true.
            call conductance_zhang15(conductance,G,St)
         case (LINMRG)
            call conductance_linmrg(conductance,G,St)
         case default
            stop "The auroral precipitation model type entered is not supported."
      end select

      ! Correct for multiple reflections if you're so inclined
      if (conductance%doMR) call conductance_mr(conductance,G,St)

      ! Smooth precipitation energy and flux before calculating conductance.
      if (conductance%doAuroralSmooth) call conductance_auroralsmooth(St,G,conductance)
      
      if (present(gcm)) then
         ! Use GCM conductance.
         call apply_gcm2mix(gcm,St,h)
         St%Vars(:,:,SIGMAP) = max(conductance%pedmin,St%Vars(:,:,SIGMAP))
         St%Vars(:,:,SIGMAH) = max(conductance%hallmin,St%Vars(:,:,SIGMAH))
      else if (conductance%const_sigma) then
         ! Use constant conductance.
         St%Vars(:,:,SIGMAP) = conductance%ped0
         St%Vars(:,:,SIGMAH) = 0.D0
      else
         ! Use vector addition of auroral and EUV conductance.
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

    subroutine GenMirrorRatio(G,St,doIGRFO)
      ! Calculate mirror ratio array RM(G%Np,G%Nt).
      ! NOTE: Leaving this to be done every time at every lat/lon to accomodate improved model later
      type(mixGrid_T) , intent(in) :: G
      type(mixState_T), intent(in) :: St
      logical,optional,intent(in) :: doIGRFO
      real(rp) :: mlat,mlon
      integer  :: i,j
      logical  :: doIGRF
      
      if(present(doIGRFO)) then
         doIGRF = doIGRFO
      else ! default is using IGRF, i.e, default is unequal split.
         doIGRF = .true. 
      endif

      if (RinMHD > 0) then
         !Calculate actual mirror ratio
         !NOTE: Should replace this w/ actual inner boundary field strength

         !$OMP PARALLEL DO default(shared) &
         !$OMP private(i,j,mlat,mlon)
         do j=1,G%Nt
            do i=1,G%Np
               mlat = PI/2 - G%t(i,j)
               mlon = G%p(i,j)

               if(.not.doIGRF) then
                  RM(i,j) = MirrorRatio(mlat,RinMHD)
               elseif (St%hemisphere==NORTH) then
                  RM(i,j) = IGRFMirrorRatio(+mlat,+mlon,RinMHD)
               else
                  !Southern, always a right-hand system based on the local pole.
                  !SM phi (mlon) goes in clock-wise as opposed to counter-clockwise if looking down on the southern pole from above.
                  RM(i,j) = IGRFMirrorRatio(-mlat,-mlon,RinMHD)
               endif
            enddo
         enddo !j loop
      else
         !Set mirror ratio everywhere no matter the inner boundary to 10
         !Note: This may have been okay for simulating magnetospheres 40 years ago, but it's 2021 now
         RM = 10.0
      endif
    end subroutine GenMirrorRatio

    subroutine conductance_euv(conductance,G,St)
      ! Derive solar EUV conductance: euvSigmaP, euvSigmaH.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      real(rp) :: ang65, ang100, pref, href, shall
      real(rp) :: speder, pedslope, pedslope2, hallslope,sigmap65, sigmah65, sigmap100

      conductance%zenith = PI/2 - ( asin(G%x) + St%tilt )
      ! An alternative (correct) definition of the zenith angle (for Moen-Brekke)
      conductance%coszen = G%x*cos(St%tilt)+sqrt(1.-G%x**2-G%y**2)*sin(St%tilt) ! as it should be
      conductance%zenith = acos(conductance%coszen)

      select case ( conductance%euv_model_type )
         case (AMIE)
            ang65     = pi/180.0*65.0
            ang100    = pi*5.0/9.0
            pref      = 2.0*250.0**(-2.0/3.0)
            href      = 1.0/(1.8*sqrt(250.0))
            shall     = 1.8*sqrt(conductance%f107)
            speder    = 0.5*conductance%f107**(2.0/3.0)
            pedslope  = 0.24*pref*speder*rad2deg
            pedslope2 = 0.13*pref*speder*rad2deg   
            hallslope = 0.27*href*shall*rad2deg;
            sigmap65  = speder*cos(ang65)**(2.0/3.0)
            sigmah65  = shall*cos(ang65)
            sigmap100 = sigmap65-(ang100-ang65)*pedslope

            where (conductance%zenith <= ang65) 
               conductance%euvSigmaP = speder*cos(conductance%zenith)**(2.0/3.0)
               conductance%euvSigmaH = shall *cos(conductance%zenith)
            elsewhere (conductance%zenith <= ang100)
               conductance%euvSigmaP = sigmap65 - pedslope *(conductance%zenith - ang65)
               conductance%euvSigmaH = sigmah65 - hallslope*(conductance%zenith - ang65)
            elsewhere (conductance%zenith > ang100)
               conductance%euvSigmaP = sigmap100 - pedslope2*(conductance%zenith-ang100)
               conductance%euvSigmaH = sigmah65  - hallslope*(conductance%zenith-ang65)
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
         case (LOMPE) 
            conductance%euvSigmaP = SigP_EUV_LOMPE(conductance%zenith,conductance%f107)
            conductance%euvSigmaH = SigH_EUV_LOMPE(conductance%zenith,conductance%f107)

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
      ! Derive electron precipitation energy flux and avg energy using Fedder95 formula.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      
      real(rp) :: signOfY, signOfJ
      real(rp) :: Rout = 6.D0, Rin = 1.2D0
      real(rp) :: rhoFactor = 3.3D-24*0.5D0
      real(rp) :: rPolarBound = 20.0D0*pi/180.D0, rEquatBound = 30.0D0*pi/180.0D0, rLowLimit = 0.02D0
      real(rp) :: D,Cs,Pe,Ne,dE,phi0,E0,aRes,dEc
      integer :: i,j

      if (St%hemisphere==NORTH) then
         signOfY = -1
         signOfJ = -1  
      elseif (St%hemisphere==SOUTH) then
         signOfY = 1
         signOfJ = 1
      else
         stop 'Wrong hemisphere label. Stopping...'
      endif

      ! fills in rampFactor
      if (conductance%doRamp) then
         call conductance_ramp(conductance,G,rPolarBound,rEquatBound,rLowLimit)
      else 
         conductance%rampFactor = 1.0D0
      endif
      
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,D,Cs,dE,phi0,E0,aRes,dEc)
      do j=1,G%Nt
         do i=1,G%Np
            if (conductance%doChill) then
               ! MHD density replaced with gallagher where it's lower      
               ! and temperature changed correspondingly
               D  = max(G%D0(i,j)*Mp_cgs, St%Vars(i,j,DENSITY)) ! [g/cm^3]
               Cs = St%Vars(i,j,SOUND_SPEED)*sqrt(St%Vars(i,j,DENSITY)/D) ! [cm/s]
            else
               D  = St%Vars(i,j,DENSITY)
               Cs = St%Vars(i,j,SOUND_SPEED)
            endif
            ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
            E0 = conductance%alpha*Mp_cgs*heFrac*erg2kev*Cs**2*conductance%rampFactor(i,j)
            conductance%E0  (i,j) = E0
            ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
            phi0 = sqrt(kev2erg)/(heFrac*Mp_cgs)**1.5D0*conductance%beta*D*sqrt(E0)*conductance%rampFactor(i,j)
            conductance%phi0(i,j) = phi0
            
            ! resistence out of the ionosphere is 2*rout resistence into the
            ! ionosphere is 2*rin outward current is positive
            if( signOfJ*St%Vars(i,j,FAC) >=0. ) then
              aRes = 2.D0*Rout
            else
              aRes = 2.D0*Rin
            endif

            ! Density floor to limit characteristic energy.  See Wiltberger et al. 2009 for details.
            D = min(D,rhoFactor*conductance%euvSigmaP(i,j))
            ! Original form for ref: 
            ! (heFrac*Mp_cgs)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*conductance%R*aRes*signOfJ*(St%Vars(:,:,FAC)*1.D-6)*sqrt(E0)/D
            dE = (heFrac*Mp_cgs)**1.5D0/eCharge*1.D-4*sqrt(erg2kev)*conductance%R*aRes*signOfJ*(St%Vars(i,j,FAC)*1.D-6)*sqrt(E0)/D
            conductance%deltaE(i,j) = dE
            dEc = min(dE,maxDrop)

            St%Vars(i,j,AVG_ENG) = max( E0+dEc, eTINY)
            if(dEc>0.) then
              St%Vars(i,j,NUM_FLUX) = phi0*(8.D0-7.D0*exp(-dEc/7.D0/E0))
            else
              St%Vars(i,j,NUM_FLUX) = phi0*exp(dEc/E0)
            endif
            St%Vars(i,j,DELTAE) = conductance%deltaE(i,j) ! [kV]
         enddo ! i
      enddo ! j

    end subroutine conductance_fedder95

    subroutine conductance_zhang15(conductance,G,St)
      ! Derive electron precipitation energy flux and avg energy using the nonlinear Fridman-Lemaire relation [Zhang et al., 2014JA020615].
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      
      real(rp) :: signOfY, signOfJ
      real(rp) :: D,Cs,Pe,Ne,J2eF0,eV2kT,dE,phi0,kT
      integer :: i,j
      
      if (St%hemisphere==NORTH) then
         signOfY = -1
         signOfJ = -1  
      elseif (St%hemisphere==SOUTH) then
         signOfY = 1
         signOfJ = 1
      else
         stop 'Wrong hemisphere label. Stopping...'
      endif
      
      if (doDrift) then
        ! Use artificial drift to get dawn-preferred diffuse electron precipitation
        ! when only using MHD information to derive it.
        call conductance_auroralmask(conductance,G,signOfY)
      else
        ! conductance%drift should be turned off when using RCM for diffuse
        conductance%drift = 1.0
      endif
      
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,D,Cs,Pe,Ne,J2eF0,eV2kT,dE,phi0,kT)
      do j=1,G%Nt
         do i=1,G%Np
            if (conductance%doChill) then
               ! MHD density replaced with gallagher where it's lower      
               ! and temperature changed correspondingly
               D  = max(G%D0(i,j)*Mp_cgs, St%Vars(i,j,DENSITY)) ! [g/cm^3]
               Cs = St%Vars(i,j,SOUND_SPEED)*sqrt(St%Vars(i,j,DENSITY)/D) ! [cm/s]
            else
               D  = St%Vars(i,j,DENSITY)
               Cs = St%Vars(i,j,SOUND_SPEED)
            endif
            ! electron pressure from MHD side in [Pa]
            ! 0.1 is to convert [g/cm^3]*[cm/s]^2=[g/cm/s^2] to [Pa].
            Pe = 0.1/MIXgamma*alpha_RCM(i,j)*D*(Cs**2)
            ! electron number density from MHD side in [/m^3].
            Ne = D/(Mp_cgs*heFrac)*1.0D6
            ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
            kT = Pe/Ne/kev2J
            conductance%E0  (i,j) = 2.0*kT
            ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
            phi0 = beta_RCM(i,j)* sqrt(Pe*Ne/(2.0D-3*PI*Me_cgs))*1.0D-4
            conductance%phi0(i,j) = phi0

            ! Note JF0 may go inf where phi0/Pe_MHD/Ne_MHD is zero.
            J2eF0 = min( 1.D-4*signOfJ*(St%Vars(i,j,FAC)*1.D-6)/(eCharge*phi0), RM(i,j)*0.99 )
            ! NonLinear Fridman-Lemaire relation: 
            ! eV = 2*kB*Te + eV*(1-exp(-eV/(RM-1)/kB/Te))/(1-(1-1/RM)*exp(-eV/(RM-1)/kB/Te)), when
            ! 1<=J/e/F0<=RM.
            if (J2eF0>1.0) then
               dE = kT*(RM(i,j)-1.D0)*dlog((RM(i,j)-1.D0)/(RM(i,j)-J2eF0))
               conductance%deltaE(i,j) = dE
               St%Vars(i,j,Z_NFLUX) = J2eF0*phi0
            else
               conductance%deltaE(i,j) = 0.0
               St%Vars(i,j,Z_NFLUX) = phi0*conductance%drift(i,j)
            endif

            ! limit the max potential energy drop to maxDrop (20 [keV])
            ! Use capped dE to calculate the mean energy.
            ! eV2kT = exp(-eV/(RM-1)/kB/Te)
            eV2kT = exp( -min(dE,maxDrop)/(kT*(RM(i,j)-1.D0)) )
            ! Floor the mean energy.
            ! Note in the original code of Zhang, EAVG is simply 2*kT+dE, a good enough approximation.
            St%Vars(i,j,Z_EAVG) = max( 2.0*kT+min(dE,maxDrop)*(1.D0-eV2kT)/(1.D0-(1.D0-1.D0/RM(i,j))*eV2kT), eTINY)

            ! Apply Zhang15 precipitation to main arrays
            St%Vars(i,j,DELTAE)   = conductance%deltaE(i,j) ! [kV]
            St%Vars(i,j,AVG_ENG)  = St%Vars(i,j,Z_EAVG)     ! [keV]
            St%Vars(i,j,NUM_FLUX) = St%Vars(i,j,Z_NFLUX)    ! [#/cm^2/s]
         enddo ! i
      enddo ! j

    end subroutine conductance_zhang15
 
    subroutine conductance_linmrg(conductance,G,St)
      ! Derive mono-diffuse electron precipitation where mono is based on linearized FL relation,
      ! diffuse is derived from RCM. The merging code was duplicated from kmerge. 
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      integer :: i,j
      real(rp) :: wC1,wC2,wRCM,gtype
      real(rp) :: mhd_eavg,mhd_nflx,mhd_eflx,mhd_delE
      real(rp) :: rcm_eavg,rcm_nflx,rcm_eflx
      real(rp) :: mix_eavg,mix_nflx,mix_eflx
      logical :: isRCM,isMHD,isMIX
      St%Vars(:,:,AUR_TYPE) = 0

      ! Two thresholds of rcm grid type between which both MHD and RCM precipitation will be merged.
      wC1 = 0.15
      wC2 = 1.0-wC1

      !Get RCM grid weighting: 1=RCM and 0=MHD
      call conductance_IM_GTYPE(G,St)

      ! derive spatially varying beta using RCM precipitation and thermal fluxes. Need IM_GTYPE.
      call conductance_beta(G,St)

      ! Derive mono using the linearized FL relation.
      call conductance_linmono(conductance,G,St)

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,isRCM,isMHD,isMIX,wRCM,gtype) &
      !$OMP private(mhd_eavg,mhd_nflx,mhd_eflx,mhd_delE) &
      !$OMP private(rcm_eavg,rcm_nflx,rcm_eflx         ) &
      !$OMP private(mix_eavg,mix_nflx,mix_eflx         )
      do j=1,G%Nt
         do i=1,G%Np
            !Start by figuring out where we are
            !gtype = St%Vars(i,j,IM_GTYPE)
            gtype = gtype_RCM(i,j)         
            isRCM = gtype >= wC2
            isMHD = gtype <= wC1
            if ( (.not. isRCM) .and. (.not. isMHD) ) then
               isMIX = .true.
            endif

            !Grab values
            mhd_eavg = St%Vars(i,j,AVG_ENG)
            mhd_nflx = St%Vars(i,j,NUM_FLUX)
            mhd_delE = conductance%deltaE(i,j)
            mhd_eflx = St%Vars(i,j,AVG_ENG)*St%Vars(i,j,NUM_FLUX)*kev2erg
            rcm_nflx = St%Vars(i,j,IM_ENFLX)
            rcm_eflx = St%Vars(i,j,IM_EFLUX)
            rcm_eavg = rcm_eflx/(rcm_nflx*kev2erg) !Back to keV

            if (isMHD) then
               St%Vars(i,j,AVG_ENG ) = mhd_eavg
               St%Vars(i,j,NUM_FLUX) = mhd_nflx
            else if (isRCM) then
               if (rcm_nflx <= TINY) then
                  rcm_eavg = eTINY
               endif
               St%Vars(i,j,AVG_ENG ) = rcm_eavg
               St%Vars(i,j,NUM_FLUX) = rcm_nflx
               conductance%deltaE(i,j) = 0.0
            else
               !Mixture
               wRCM = RampUp(gtype,wC1,wC2-wC1)
               if ( (rcm_nflx > TINY) ) then
                  !Mix both
                  mix_nflx = wRCM*rcm_nflx + (1-wRCM)*mhd_nflx
                  mix_eflx = wRCM*rcm_eflx + (1-wRCM)*mhd_eflx
               else
                  !RCM data wasn't good so just use MHD
                  mix_nflx = mhd_nflx
                  mix_eflx = mhd_eflx
               endif
               mix_eavg = mix_eflx/(mix_nflx*kev2erg)

               St%Vars(i,j,AVG_ENG ) = mix_eavg
               St%Vars(i,j,NUM_FLUX) = mix_nflx
            endif
         enddo
      enddo

    end subroutine conductance_linmrg

    subroutine conductance_linmono(conductance,G,St)
      ! Linearized Fridman-Lemaire relation for mono.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      
      real(rp) :: signOfY, signOfJ
      integer :: i,j
      real(rp) :: D,Cs,Pe,Ne,J2eF0,eV2kT,dE,phi0,kT
      real(rp) :: Ne_floor

      if (St%hemisphere==NORTH) then
         signOfY = -1
         signOfJ = -1  
      elseif (St%hemisphere==SOUTH) then
         signOfY = 1
         signOfJ = 1
      else
         stop 'Wrong hemisphere label. Stopping...'
      endif

      Ne_floor = 0.03e6 ! minimum Ne in [/m^3] when evaluating the linearized FL relation.

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,D,Cs,Pe,Ne,J2eF0,eV2kT,dE,phi0,kT)
      do j=1,G%Nt
         do i=1,G%Np
            if (conductance%doChill) then
               ! MHD density replaced with gallagher where it's lower      
               ! and temperature changed correspondingly
               D  = max(G%D0(i,j)*Mp_cgs, St%Vars(i,j,DENSITY)) ! [g/cm^3]
               Cs = St%Vars(i,j,SOUND_SPEED)*sqrt(St%Vars(i,j,DENSITY)/D) ! [cm/s]
            else
               D  = St%Vars(i,j,DENSITY)
               Cs = St%Vars(i,j,SOUND_SPEED)
            endif
            ! electron pressure from MHD side in [Pa]
            ! 0.1 is to convert [g/cm^3]*[cm/s]^2=[g/cm/s^2] to [Pa].
            Pe = 0.1/MIXgamma*alpha_RCM(i,j)*D*(Cs**2)
            ! electron number density from MHD side in [/m^3].
            Ne = D/(Mp_cgs*heFrac)*1.0D6
            ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
            kT = Pe/Ne/kev2J
            conductance%E0  (i,j) = 2.0*kT
            ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
            phi0 = beta_RCM(i,j)* sqrt(Pe*Ne/(2.0D-3*PI*Me_cgs))*1.0D-4
            conductance%phi0(i,j) = phi0

            ! Note JF0 may go inf where phi0/Pe_MHD/Ne_MHD is zero.
            J2eF0 = signOfJ*(St%Vars(i,j,FAC)*1.0D-6)/(eCharge*phi0*1.0D4)
            ! Linearize the Fridman-Lemaire relation: 
            ! eV = kB*Te*(RM-1)*ln((RM-1)/(RM-J/e/F0)) when 1<=J/e/F0<=RM.
            ! eV ~ kB*Te*(J/e/F0-1) when RM>>J/e/F0
            if (J2eF0>1.0 .and. Ne>Ne_floor) then
               dE = kT*(J2eF0-1.0)
               conductance%deltaE(i,j) = dE
               eV2kT = min(dE,maxDrop)/kT + 1.0
               St%Vars(i,j,Z_NFLUX) = phi0*eV2kT
               St%Vars(i,j,Z_EAVG) = kT*(eV2kT+1.0/eV2kT)
            else
               ! When there is no acceleration, MHD flux represents the diffuse precipitation.
               conductance%deltaE(i,j) = 0.0
               St%Vars(i,j,Z_NFLUX) = phi0
               St%Vars(i,j,Z_EAVG) = kT 
            endif
            St%Vars(i,j,DELTAE)   = conductance%deltaE(i,j)
            St%Vars(i,j,AVG_ENG)  = max(St%Vars(i,j,Z_EAVG),eTINY)   ! [keV]
            St%Vars(i,j,NUM_FLUX) = St%Vars(i,j,Z_NFLUX)  ! [#/cm^2/s]
         enddo ! i
      enddo ! j
      
    end subroutine conductance_linmono

    subroutine conductance_mr(conductance,G,St)
      ! George Khazanov's multiple reflection(MR) corrections
      ! Modify the diffuse precipitation energy and mean energy based on equation (5) in Khazanov et al. [2019JA026589]
      ! Kc = 3.36-exp(0.597-0.37*Eavg+0.00794*Eavg^2)
      ! Eavg_c = 0.073+0.933*Eavg-0.0092*Eavg^2
      ! Nflx_c = Eflx_c/Eavg_c = Eflx*Kc/Eavg_c = Nflx*Eavg*Kc/Eavg_c
      ! where Eavg is the mean energy in [keV] before MR, Eavg_c is the modified mean energy. 
      ! Kc is the ratio between modified enflux and unmodified enflux.
      ! No longer use the conductance modification factor in Khazanov et al. [2018SW001837] (eq 7).
      type(mixConductance_T), intent(in) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      real(rp) :: eavg,nflx,eavgMR,nflxMR
      integer :: i,j

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,eavg,nflx,eavgMR,nflxMR)
      do j=1,G%Nt
        do i=1,G%Np
          if(conductance%deltaE(i,j)<=0.0) then
            ! Only modify diffuse precipitation. May modify mono precipitation when a formula is available.
            eavg = St%Vars(i,j,AVG_ENG)
            nflx = St%Vars(i,j,NUM_FLUX)
            call AugmentMR(eavg,nflx,eavgMR,nflxMR)
            St%Vars(i,j,AVG_ENG ) = eavgMR
            St%Vars(i,j,NUM_FLUX) = nflxMR
          endif
        enddo
      enddo
    end subroutine conductance_mr

    subroutine AugmentMR(eavg,nflx,eavgMR,nflxMR)
      ! Cleaned up routine to calculate MR-augmented eavg and nflx
      real(rp), intent(in)  :: eavg  ,nflx
      real(rp), intent(out) :: eavgMR,nflxMR
      real(rp) :: eflux,Kc

      if ( (eavg<=30.0) .and. (eavg>=0.5) ) then
         ! Kc becomes negative when E_avg > 47-48 keV.
         Kc = 3.36 - exp(0.597-0.37*eavg+0.00794*eavg**2)
         eflux = eavg*nflx
         eavgMR = 0.073+0.933*eavg-0.0092*eavg**2 ! [keV]
         nflxMR = Kc*eflux/eavgMR
      else
         !Nothing to do
         eavgMR = eavg
         nflxMR = nflx
      endif
    end subroutine AugmentMR

    subroutine conductance_auroralsmooth(St,G,conductance)
      type(mixState_T), intent(inout) :: St
      type(mixGrid_T), intent(in) :: G      
      type(mixConductance_T), intent(inout) :: conductance
      integer :: i, j
      logical :: smthDEPonly = .true.
      integer :: smthE
      smthE = 1 ! 1. smooth EFLUX; 2. smooth EAVG

      tmpC = 0.D0
      tmpD = 0.D0
      tmpE = 0.D0
      tmpF = 0.D0

      ! Test only smoothing diffuse precipitation.
      ! Diffuse mask is St%Vars(:,:,Z_NFLUX)<0.
      ! Get AVG_ENG*mask and NUM_FLUX*mask for smoothing.
      ! Fill AVG_ENG and NUM_FLUX with smoothed where mask is true.
      ! this domain of conductance is never used. Here it's used for diffuse precipitation mask. Be careful if not in RCMONO mode.
      conductance%PrecipMask = 1.D0 
      if(smthDEPonly) then
        where(St%Vars(:,:,Z_NFLUX)>=0)
          conductance%PrecipMask = 0.D0
        end where
        ! back up mono
        tmpC = St%Vars(:,:,AVG_ENG)
        tmpD = St%Vars(:,:,NUM_FLUX)
      endif

      ! get expanded diffuse with margin grid (ghost cells)
      if(smthE==1) then
        call conductance_margin(G,St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)*conductance%PrecipMask, tmpE)
      else
        call conductance_margin(G,St%Vars(:,:,AVG_ENG)*conductance%PrecipMask, tmpE)
      end if
      call conductance_margin(G,St%Vars(:,:,NUM_FLUX)*conductance%PrecipMask, tmpF)

      !$OMP PARALLEL DO default(shared) collapse(2) &
      !$OMP private(i,j)
      do j=1,G%Nt
        do i=1,G%Np
          ! SmoothOperator55(Q,isGO)
          St%Vars(i,j,NUM_FLUX)=SmoothOperator55(tmpF(i:i+4,j:j+4))
          if(smthE==1) then
            St%Vars(i,j,AVG_ENG) = max(SmoothOperator55(tmpE(i:i+4,j:j+4))/St%Vars(i,j,NUM_FLUX),1.D-8)
          else
            St%Vars(i,j,AVG_ENG) = max(SmoothOperator55(tmpE(i:i+4,j:j+4)),1.D-8)
          end if
        end do
      end do

      if(smthDEPonly) then
        ! combine smoothed diffuse and unsmoothed mono
        where(St%Vars(:,:,Z_NFLUX)>=0)
          St%Vars(:,:,AVG_ENG)  = tmpC
          St%Vars(:,:,NUM_FLUX) = tmpD
        end where
      endif
    end subroutine conductance_auroralsmooth

    subroutine conductance_aurora(conductance,G,St)
      ! Use Robinson formula to get Pedersen and Hall conductance from electron precipitation.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      conductance%engFlux = kev2erg*St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)  ! Energy flux in ergs/cm^2/s
      conductance%deltaSigmaP = SigmaP_Robinson(St%Vars(:,:,AVG_ENG),conductance%engFlux)
      conductance%deltaSigmaH = SigmaH_Robinson(St%Vars(:,:,AVG_ENG),conductance%engFlux)
    end subroutine conductance_aurora

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

    subroutine conductance_auroralmask(conductance,G,signOfY)
      ! Artificial auroral mask used in Zhang15 to represent dawnward shift of diffuse electron precipitation.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G      
      real(rp), intent(in) :: signOfY
      real(rp) :: Rio, al0, alp, Radi, order, Rady
      
      Rio = 1.02
      al0 = -5.0*deg2rad
      alp = 28*deg2rad
      alp = min(28*deg2rad,alp)
      Rady = Rio*sin(alp-al0)
      Radi = Rady**2 ! Rio**2*(1-cos(alp-al0)*cos(alp-al0))
      order = 2.0
      
      RRdi = 0.D0
      RRdi = (G%y-0.03*signOfY)**2 + ( G%x/cos(al0) - Rio*cos(alp-al0)*tan(al0) )**2
      where(RRdi < Radi)
         conductance%AuroraMask = cos((RRdi/Radi)**order*PI/2)+0.D0
      elsewhere
         conductance%AuroraMask = 0.D0
      end where
      where(abs(G%y)<Rady)
         conductance%drift = 1.D0 + 0.5*G%y*signOfY/Rady
      elsewhere
         conductance%drift = 1.D0
      end where
    end subroutine conductance_auroralmask

    subroutine conductance_IM_GTYPE(G,St)
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(in) :: St
      logical :: isAnc(G%Np,G%Nt)
      integer :: i,j

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j)
      do j=1,G%Nt ! use open BC for lat.
        do i=1,G%Np ! use periodic BC for lon.
          if(St%Vars(i,j,IM_GTYPE)>0.01 .and. St%Vars(i,j,IM_GTYPE)<0.99) then
            isAnc(i,j) = .false.
          else
            isAnc(i,j) = .true.
          endif
        enddo ! i
      enddo ! j
      gtype_RCM = St%Vars(:,:,IM_GTYPE) ! supposed to be between 0 and 1.
      call conductance_smooth(G,gtype_RCM,isAnc)

      gtype_RCM = min(gtype_RCM,1.0)
      gtype_RCM = max(gtype_RCM,0.0)
    end subroutine conductance_IM_GTYPE

    subroutine conductance_beta(G,St)
      ! Use RCM precipitation and source population to derive the loss cone rate beta.
      ! Assume beta is one in the polar cap and smooth across RCM boundary.
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      logical :: isAnc(G%Np,G%Nt)
      real(rp) :: temp(G%Np,G%Nt)
      real(rp) :: phi0_rcm
      integer :: i,j

      St%Vars(:,:,IM_BETA) = 1.0
      isAnc = .false.
      ! set the first circle around pole as anchore.
      isAnc(:,1) = .true.

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,phi0_rcm)
      do j=1,G%Nt
         do i=1,G%Np
            ! Total RCM thermal flux includes the trapped and precipitated in both NH and SH.
            ! 1.0e-4 is to convert to [#/cm^2/s]
            ! sqrt([Pa]*[#/m^3]/[kg]) = sqrt([#/m^4/s^2]) = 1e-4*[#/cm^2/s]
            phi0_rcm = sqrt(St%Vars(i,j,IM_EPRE)*St%Vars(i,j,IM_EDEN)/(Me_cgs*1e-3*2*pi))*1.0e-4 + St%Vars(i,j,IM_ENFLX)*2.0
            if(phi0_rcm>TINY) then
               St%Vars(i,j,IM_BETA) = St%Vars(i,j,IM_ENFLX)/phi0_rcm
               isAnc(i,j) = .true. ! set points with valid rcm beta as anchore.
            elseif(St%Vars(i,j,IM_GTYPE) > 0.5) then
               ! In the low lat, if there is no meaningful RCM precipitation,
               ! set beta=0 to freeze other precipitation mechanism.
               ! also make it anchor points to avoid smoothing.
               St%Vars(i,j,IM_BETA) = 0.0
               isAnc(i,j) = .true.
            endif
         enddo
      enddo
      temp = St%Vars(:,:,IM_BETA) ! supposed to be between 0 and 1.
      call conductance_smooth(G,temp,isAnc)
      St%Vars(:,:,IM_BETA) = temp

      ! IM_BETA is for output. beta_RCM is used in calculation.
      beta_RCM = min(St%Vars(:,:,IM_BETA), 1.0)
    end subroutine conductance_beta

    subroutine conductance_smooth(Gr,Q,isAnchor)
      ! Do smoothing window on ReMIX grid quantity
      ! Skip certain points
      type(mixGrid_T), intent(in) :: Gr
      real(rp), intent(inout) :: Q(Gr%Np,Gr%Nt)
      logical, intent(in) :: isAnchor(Gr%Np,Gr%Nt)
      real(rp) :: temp(Gr%Np,Gr%Nt)
      real(rp) :: thres,mad,Ttmp
      integer :: i,j,it,im1,ip1,jm1,jp1,MaxIter

      thres = 0.025
      MaxIter = 15
      call FixPole(Gr,Q)

      do it=1,MaxIter
        mad = 0.D0 ! max abs difference from last iteration.
        temp = Q

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,jm1,jp1,im1,ip1,Ttmp) &
        !$OMP reduction(max:mad)
        do j=1,Gr%Nt ! use open BC for lat.
          do i=1,Gr%Np ! use periodic BC for lon.
            ! Do not smooth anchor points.
            if(.not. isAnchor(i,j)) then
              jm1 = j-1
              jp1 = j+1
              if (j == 1)    jm1 = 1
              if (j == Gr%Nt) jp1 = Gr%Nt

              im1 = i-1
              ip1 = i+1
              if (i == 1)    im1 = Gr%Np
              if (i == Gr%Np) ip1 = 1

              Ttmp =(temp(im1,jm1)+temp(im1,j)+temp(im1,jp1) &
                   + temp(i  ,jm1)+temp(i  ,j)+temp(i  ,jp1) &
                   + temp(ip1,jm1)+temp(ip1,j)+temp(ip1,jp1))/9.D0
              mad  = max(abs(Q(i,j)-Ttmp),mad)
              Q(i,j) = Ttmp
            endif
          enddo
        enddo
        call FixPole(Gr,Q)
        if(mad<thres) exit
      enddo
    end subroutine conductance_smooth

    !Enforce pole condition that all values at the same point are equal
    subroutine FixPole(Gr,Q)
      type(mixGrid_T), intent(in)    :: Gr
      real(rp)       , intent(inout) :: Q(Gr%Np,Gr%Nt)

      real(rp) :: Qavg

      Qavg = sum(Q(:,2))/Gr%Np
      Q(:,1) = Qavg
    end subroutine FixPole

    subroutine conductance_margin(G,array,arraymar)
      ! Conductance boundary processing.
      type(mixGrid_T), intent(in) :: G      
      real(rp),dimension(G%Np,G%Nt), intent(in) :: array
      real(rp),dimension(G%Np+4,G%Nt+4), intent(out) :: arraymar
      arraymar(3:G%Np+2,3:G%Nt+2) = array
      ! reflective BC for lat
      arraymar(3:G%Np+2,1) = array(:,3)
      arraymar(3:G%Np+2,2) = array(:,2)
      arraymar(3:G%Np+2,G%Nt+3) = array(:,G%Nt-1)
      arraymar(3:G%Np+2,G%Nt+4) = array(:,G%Nt-2)
      ! periodic BC for lon
      arraymar(1:2,:) = arraymar(G%Np+1:G%Np+2,:)
      arraymar(G%Np+3:G%Np+4,:) = arraymar(3:4,:)
    end subroutine conductance_margin

    subroutine SetMIXgamma(gamma)
      !Routine to change MIX-gamma on the fly if necessary
      real(rp), intent(in) :: gamma
      MIXgamma = gamma
    end subroutine SetMIXgamma

  end module mixconductance
