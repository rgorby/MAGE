module mixconductance
  use mixdefs
  use mixtypes
  use earthhelper
  use gcmtypes
  use gcminterp
  use math
  
  implicit none

  real(rp), dimension(:,:), allocatable, private :: tmpD,tmpC ! used for chilling in Fedder95. Declare it here so we can allocate in init.
  real(rp), dimension(:,:), allocatable, private :: JF0,RM,RRdi ! used for zhang15
  real(rp), dimension(:,:), allocatable, private :: tmpE,tmpF ! used for smoothing precipitation avg_eng and num_flux
  real(rp), dimension(:,:), allocatable, private :: Kc ! used for multi-reflection modification

  !Replacing some hard-coded inline values (bad) w/ module private values (slightly less bad)
  real(rp), parameter, private :: maxDrop = 20.0 !Hard-coded max potential drop [kV]
  real(rp), private :: RinMHD = 0.0 !Rin of MHD grid (0 if not running w/ MHD)
  logical , private :: doRobKap = .true. !Use Kaeppler+ 15 correction to SigH/SigP from Robinson

  contains
    subroutine conductance_init(conductance,Params,G)
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

      if (.not. allocated(tmpD)) allocate(tmpD(G%Np,G%Nt))
      if (.not. allocated(tmpC)) allocate(tmpC(G%Np,G%Nt))  

      if (.not. allocated(JF0))  allocate(JF0 (G%Np,G%Nt))      
      if (.not. allocated(RM))   allocate(RM  (G%Np,G%Nt))      
      if (.not. allocated(RRdi)) allocate(RRdi(G%Np,G%Nt))      
      if (.not. allocated(tmpE)) allocate(tmpE(G%Np+4,G%Nt+4)) ! for boundary processing.
      if (.not. allocated(tmpF)) allocate(tmpF(G%Np+4,G%Nt+4))
      if (.not. allocated(Kc))   allocate(Kc  (G%Np,G%Nt))

      RinMHD = Params%RinMHD

    end subroutine conductance_init

    subroutine conductance_euv(conductance,G,St)
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

      tmpC = 0.D0
      tmpD = 0.D0

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
      conductance%deltaE = min(maxDrop,conductance%deltaE)

      ! floor on total energy
      St%Vars(:,:,AVG_ENG) = max(conductance%E0 + conductance%deltaE,1.D-8)

      where  ( conductance%deltaE > 0. )
         St%Vars(:,:,NUM_FLUX) = conductance%phi0*(8.D0-7.D0*exp(-conductance%deltaE/7.D0/conductance%E0))
      elsewhere 
         St%Vars(:,:,NUM_FLUX) = conductance%phi0*exp(conductance%deltaE/conductance%E0)
      end where

    end subroutine conductance_fedder95

    subroutine conductance_zhang15(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      
      real(rp) :: signOfY, signOfJ

      tmpC = 0.D0
      tmpD = 0.D0
      JF0 = 0.D0
      
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

      if (conductance%doChill) then
         ! MHD density replaced with gallagher where it's lower      
         ! and temperature changed correspondingly
          tmpD = max(G%D0*Mp_cgs,St%Vars(:,:,DENSITY))
          tmpC = St%Vars(:,:,SOUND_SPEED)*sqrt(St%Vars(:,:,DENSITY)/tmpD)
       else
          tmpD = St%Vars(:,:,DENSITY)
          tmpC = St%Vars(:,:,SOUND_SPEED)
      end if

      call conductance_auroralmask(conductance,G,signOfY)

!Flag_up
      conductance%E0 = conductance%alpha*Mp_cgs*heFrac*erg2kev*tmpC**2
      conductance%phi0 = sqrt(kev2erg)/(heFrac*Mp_cgs)**1.5D0*conductance%beta*sqrt(heFrac*1836.152674)*0.39894228*tmpD*sqrt(conductance%E0)
      ! conductance%phi0 = sqrt(kev2erg)/(heFrac*Mp_cgs)**1.5D0*conductance%beta*tmpD*sqrt(conductance%E0)
      
      JF0 = min( 1.D-4*signOfJ*(St%Vars(:,:,FAC)*1.e-6)/eCharge/(conductance%phi0), RM*0.99 )

      !NOTE: conductance%drift should be turned off when using RCM for diffuse
      where ( JF0 > 1. )
      ! limit the max potential energy drop to 20 [keV]
         conductance%deltaE = min( 0.5*conductance%E0*(RM - 1.D0)*dlog((RM-1.D0)/(RM-JF0)),maxDrop)
         St%Vars(:,:,Z_NFLUX) = JF0*conductance%phi0
      elsewhere
         conductance%deltaE = 0.
         St%Vars(:,:,Z_NFLUX) = conductance%phi0*conductance%drift
      end where

      ! floor on total energy
      St%Vars(:,:,Z_EAVG) = max(2.0*conductance%E0 + conductance%deltaE,1.D-8)
      St%Vars(:,:,Z_NFLUX) = St%Vars(:,:,Z_NFLUX)!*conductance%AuroraMask

      ! Apply Zhang15 precipitation to main arrays
      St%Vars(:,:,AVG_ENG)  = St%Vars(:,:,Z_EAVG) ! [keV]
      St%Vars(:,:,NUM_FLUX) = St%Vars(:,:,Z_NFLUX)! [#/cm^2/s]
      
    end subroutine conductance_zhang15

    subroutine conductance_rcmono(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      call conductance_zhang15(conductance,G,St)
      ! If there is no potential drop OR mono eflux is too low, use RCM precipitation instead.
      where(conductance%deltaE<=0.0)
         St%Vars(:,:,AVG_ENG)  = max(St%Vars(:,:,IM_EAVG),1.D-8) ! [keV]
         St%Vars(:,:,NUM_FLUX) = St%Vars(:,:,IM_EFLUX)/(St%Vars(:,:,AVG_ENG)*kev2erg) ! [ergs/cm^2/s]
         St%Vars(:,:,Z_NFLUX)  = -1.0 ! for diagnostic purposes since full Z15 does not currently work.
      end where

    end subroutine conductance_rcmono

    subroutine conductance_rcmfed(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      real(rp) :: E2Th, E2Tl
      E2Th = 2.0
      E2Tl = 1.0

      ! Use totally Fedder precipitation if deltaE > 2Te.
      call conductance_fedder95(conductance,G,St)
      tmpC = conductance%deltaE/conductance%E0
      where(tmpC<=E2Th.and.tmpC>=E2Tl) ! Linearly combine RCM and Fedder where 1<=deltaE/Te<=2.
         St%Vars(:,:,AVG_ENG)  = max(((E2Th-tmpC)*St%Vars(:,:,IM_EAVG)+(tmpC-E2Tl)*St%Vars(:,:,AVG_ENG))/(E2Th-E2Tl),1.D-8) ! [keV]
         St%Vars(:,:,NUM_FLUX) = ((E2Th-tmpC)*St%Vars(:,:,IM_EFLUX)/(St%Vars(:,:,AVG_ENG)*kev2erg)+(tmpC-E2Tl)*St%Vars(:,:,NUM_FLUX))/(E2Th-E2Tl) ! [ergs/cm^2/s]/[ergs]
         St%Vars(:,:,Z_NFLUX)  = -0.5 ! for diagnostic purposes since full Z15 does not currently work.
      elsewhere(tmpC<E2Tl) ! Use totally RCM precipitation if deltaE<Te.
         St%Vars(:,:,AVG_ENG)  = max(St%Vars(:,:,IM_EAVG),1.D-8) ! [keV]
         St%Vars(:,:,NUM_FLUX) = St%Vars(:,:,IM_EFLUX)/(St%Vars(:,:,AVG_ENG)*kev2erg) ! [ergs/cm^2/s]
         St%Vars(:,:,Z_NFLUX)  = -1.0 ! for diagnostic purposes since full Z15 does not currently work.
      end where

    end subroutine conductance_rcmfed

    subroutine conductance_aurora(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      ! note, this assumes that fedder has been called prior
      conductance%engFlux = kev2erg*St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)  ! Energy flux in ergs/cm^2/s
      conductance%deltaSigmaP = SigmaP_Robinson(St%Vars(:,:,AVG_ENG),conductance%engFlux)
      conductance%deltaSigmaH = SigmaH_Robinson(St%Vars(:,:,AVG_ENG),conductance%engFlux)


    end subroutine conductance_aurora

    ! George Khazanov's multiple reflection(MR) corrections
    subroutine conductance_mr(conductance,St)
      type(mixConductance_T), intent(in) :: conductance
      type(mixState_T), intent(inout) :: St

      ! Modify the diffuse precipitation energy and mean energy based on equation (5) in Khazanov et al. [2019JA026589]
      ! Kc = 3.36-exp(0.597-0.37*Eavg+0.00794*Eavg^2)
      ! Eavg_c = 0.073+0.933*Eavg-0.0092*Eavg^2
      ! Nflx_c = Eflx_c/Eavg_c = Eflx*Kc/Eavg_c = Nflx*Eavg*Kc/Eavg_c
      ! where Eavg is the mean energy in [keV] before MR, Eavg_c is the modified mean energy. Kc is the ratio between modified enflux and unmodified enflux.
      ! print *, "doMR for diffuse electron precipitation."
      Kc = 1.D0
      where(conductance%deltaE<=0.0.and.St%Vars(:,:,AVG_ENG)<=30.and.St%Vars(:,:,AVG_ENG)>=0.5)
      ! The formula was derived from E_avg between 0.5 and 30 keV. Kc becomes negative when E_avg > 47-48 keV.
         Kc = 3.36 - exp(0.597-0.37*St%Vars(:,:,AVG_ENG)+0.00794*St%Vars(:,:,AVG_ENG)**2)
         St%Vars(:,:,NUM_FLUX) = St%Vars(:,:,NUM_FLUX)*St%Vars(:,:,AVG_ENG) ! store Eflux
         St%Vars(:,:,AVG_ENG)  = 0.073+0.933*St%Vars(:,:,AVG_ENG)-0.0092*St%Vars(:,:,AVG_ENG)**2 ! [keV]
         St%Vars(:,:,NUM_FLUX) = Kc*St%Vars(:,:,NUM_FLUX)/St%Vars(:,:,AVG_ENG)
      end where
!      conductance%deltaSigmaP = (2.16-0.87*exp(-0.16*St%Vars(:,:,AVG_ENG)))*conductance%deltaSigmaP
!      conductance%deltaSigmaH = (1.87-0.54*exp(-0.16*St%Vars(:,:,AVG_ENG)))*conductance%deltaSigmaH
    end subroutine conductance_mr

    !Cleaned up routine to calculate MR-augmented eavg and nflx
    subroutine AugmentMR(eavg,nflx,eavgMR,nflxMR)
      real(rp), intent(in)  :: eavg  ,nflx
      real(rp), intent(out) :: eavgMR,nflxMR

      real(rp) :: eflux,Kc

      if ( (eavg<=30.0) .and. (eavg>=0.5) ) then
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

    !Calculate mirror ratio array
    !NOTE: Leaving this to be done every time at every lat/lon to accomodate improved model later
   subroutine GenMirrorRatio(G)
      type(mixGrid_T), intent(in) :: G

      real(rp) :: mlat,mlon
      integer  :: i,j

      if (RinMHD > 0) then
         !Calculate actual mirror ratio
         !NOTE: Should replace this w/ actual inner boundary field strength
         do j=1,G%Nt
            do i=1,G%Np
               mlat = PI/2 - G%t(i,j)
               mlon = G%p(i,j)

               RM(i,j) = MirrorRatio(mlat,RinMHD)
            enddo
         enddo
      else
         !Set mirror ratio everywhere no matter the inner boundary to 10
         !Note: This may have been okay for simulating magnetospheres 40 years ago, but it's 2021 now
         RM = 10.0
      endif
   end subroutine GenMirrorRatio

    subroutine conductance_total(conductance,G,St,gcm,h)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      type(gcm_T),optional,intent(in) :: gcm
      integer,optional,intent(in) :: h

      call GenMirrorRatio(G)

      ! always call fedder to fill in AVG_ENERGY and NUM_FLUX
      ! even if const_sigma, we still have the precip info that way

      ! compute EUV though because it's used in fedder
      call conductance_euv(conductance,G,St)
      select case ( conductance%aurora_model_type )
         case (FEDDER)
            call conductance_fedder95(conductance,G,St)
         case (ZHANG)
            call conductance_zhang15(conductance,G,St)
         case (RCMONO)
            call conductance_rcmono(conductance,G,St)
         case (RCMFED)
            call conductance_rcmfed(conductance,G,St)
         case default
            stop "The aurora precipitation model type entered is not supported."
      end select

      ! correct for multiple reflections if you're so inclined
      if (conductance%doMR) call conductance_mr(conductance,St)

      ! Smooth precipitation energy and flux before calculating conductance.
      if (conductance%doAuroralSmooth) call conductance_auroralsmooth(St,G,conductance)
      
      if (present(gcm)) then
         !write(*,*) 'going to apply!'
         call apply_gcm2mix(gcm,St,h)
         St%Vars(:,:,SIGMAP) = max(conductance%pedmin,St%Vars(:,:,SIGMAP))
         St%Vars(:,:,SIGMAH) = max(conductance%hallmin,St%Vars(:,:,SIGMAH))
         !St%Vars(:,:,SIGMAH) = min(max(conductance%hallmin,St%Vars(:,:,SIGMAH)),&
         !     St%Vars(:,:,SIGMAP)*conductance%sigma_ratio)
      else if (conductance%const_sigma) then
         !write(*,*) "conductance: const_sigma"
         St%Vars(:,:,SIGMAP) = conductance%ped0
         St%Vars(:,:,SIGMAH) = 0.D0
      else
         !write(*,*) "conductance: aurora"
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

    !Returns Robinson's SigP from eavg [kEv] and eflux [ergs/cm^2]
   elemental function SigmaP_Robinson(eavg,eflux) result(SigP)
      real(rp), intent(in) :: eavg,eflux
      real(rp) :: SigP
      SigP = 40.0*eavg*sqrt(eflux)/(16.0 + eavg**2.0)
   end function SigmaP_Robinson

   !Returns Robinson's SigH from eavg [kEv] and eflux [ergs/cm^2]
   !NOTE: Extra correction from Fedder
   elemental function SigmaH_Robinson(eavg,eflux) result(SigH)
      real(rp), intent(in) :: eavg,eflux
      real(rp) :: SigH
      real(rp) :: SigP

      SigP = SigmaP_Robinson(eavg,eflux)
      if (doRobKap) then
         !Kaeppler+ 2015
         SigH = 0.57*SigP*(eavg**0.53)
      else
         SigH = 0.45*SigP*(eavg**0.85)/(1.0 + 0.0025*eavg**2.0) !Includes extra Fedder correction to curb values for high eavg
      endif
   end function SigmaH_Robinson

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

      ! do smoothing
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

    subroutine conductance_margin(G,array,arraymar)
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

  end module mixconductance
