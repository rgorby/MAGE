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

!  real(rp), dimension(:,:), allocatable, private :: beta_RCM,alpha_RCM,gtype_RCM ! two-dimensional beta based on RCM fluxes.

  !Replacing some hard-coded inline values (bad) w/ module private values (slightly less bad)
!  real(rp), parameter, private :: maxDrop = 20.0 !Hard-coded max potential drop [kV]
!  real(rp), parameter, private :: eTINY = mixeTINY ! Floor of average energy [keV]
!  real(rp), parameter, private :: Ne_floor = 0.03e6 ! minimum Ne in [/m^3] when evaluating the linearized FL relation.
!  real(rp), parameter, private :: Ne_psp = 10.0e6 ! Ne threshold for the plasmasphere in [/m^3].
!  real(rp), parameter, private :: GuABNF = 1.e7 ! Gussenhoven+[1983] Auroral Boundary Number Flux in [#/cm^s/s].
!  real(rp), private :: RinMHD = 0.0 !Rin of MHD grid (0 if not running w/ MHD)
!  real(rp), private :: MIXgamma
!  real(rp), private :: beta_inp

  contains

    subroutine conductance_init(conductance,Params,G)
      ! Initialize conductance variables. Called in mixmain.F90/subroutine init_mix.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixParams_T)     , intent(in)    :: Params
      type(mixGrid_T)       , intent(in)    :: G

      ! define these module-wide variables so we don't have to pass the Params object to all the conductance functions
      ! this is similar to how it was done in MIX
      conductance%euv_model_type    = Params%euv_model_type
      conductance%et_model_type     = Params%et_model_type
!      conductance%alpha             = Params%alpha
!      conductance%beta              = Params%beta
!      conductance%R                 = Params%R
      conductance%F107              = Params%F107
      conductance%pedmin            = Params%pedmin
      conductance%hallmin           = Params%hallmin
      conductance%sigma_ratio       = Params%sigma_ratio
      conductance%ped0              = Params%ped0
      conductance%const_sigma       = Params%const_sigma
!      conductance%doRamp            = Params%doRamp
!      conductance%doChill           = Params%doChill
      conductance%doStarlight       = Params%doStarlight
!      conductance%doMR              = Params%doMR
!      conductance%doAuroralSmooth   = Params%doAuroralSmooth      
      conductance%apply_cap         = Params%apply_cap
!      conductance%aurora_model_type = Params%aurora_model_type
      conductance%doEMA             = Params%doEMA

      if (.not. allocated(conductance%zenith)) allocate(conductance%zenith(G%Np,G%Nt))
      if (.not. allocated(conductance%coszen)) allocate(conductance%coszen(G%Np,G%Nt))
      if (.not. allocated(conductance%euvSigmaP)) allocate(conductance%euvSigmaP(G%Np,G%Nt))
      if (.not. allocated(conductance%euvSigmaH)) allocate(conductance%euvSigmaH(G%Np,G%Nt))
      if (.not. allocated(conductance%deltaSigmaP)) allocate(conductance%deltaSigmaP(G%Np,G%Nt))
      if (.not. allocated(conductance%deltaSigmaH)) allocate(conductance%deltaSigmaH(G%Np,G%Nt))

!      if (.not. allocated(conductance%rampFactor)) allocate(conductance%rampFactor(G%Np,G%Nt))
!      if (.not. allocated(conductance%ares)) allocate(conductance%ares(G%Np,G%Nt))
!      if (.not. allocated(conductance%deltaE)) allocate(conductance%deltaE(G%Np,G%Nt))
!      if (.not. allocated(conductance%E0)) allocate(conductance%E0(G%Np,G%Nt))
!      if (.not. allocated(conductance%phi0)) allocate(conductance%phi0(G%Np,G%Nt))
!      if (.not. allocated(conductance%engFlux)) allocate(conductance%engFlux(G%Np,G%Nt))

!      if (.not. allocated(conductance%avgEng)) allocate(conductance%avgEng(G%Np,G%Nt))
!      if (.not. allocated(conductance%drift)) allocate(conductance%drift(G%Np,G%Nt))      
!      if (.not. allocated(conductance%AuroraMask)) allocate(conductance%AuroraMask(G%Np,G%Nt))      
!      if (.not. allocated(conductance%PrecipMask)) allocate(conductance%PrecipMask(G%Np,G%Nt))    

      ! these arrays are global and should not be! reallocate them
!      if(allocated(beta_RCM)) deallocate(beta_RCM)
!      if(allocated(alpha_RCM)) deallocate(alpha_RCM)
!      if(allocated(gtype_RCM)) deallocate(gtype_RCM)

!      allocate(beta_RCM(G%Np,G%Nt))
!      allocate(alpha_RCM(G%Np,G%Nt))
!      allocate(gtype_RCM(G%Np,G%Nt))

!      call SetMIXgamma(Params%gamma)
      ! MHD inner boundary, used to calc mirror ratio.
!      RinMHD = Params%RinMHD
      ! Te/Tmhd
!      alpha_RCM = 1.0/(tiote_RCM+1.0)
      ! Loss cone rate
!      beta_RCM  = conductance%beta
!      beta_inp  = conductance%beta
      ! RCM grid weight: 1. Totally on closed RCM; 0. Totally outside RCM.
      ! if conductance_IM_GTYPE is not called, gtype_RCM has all zero, MHD values have a weight of 1.
!      gtype_RCM = 0.0

    end subroutine conductance_init

    subroutine conductance_total(conductance,G,St,gcm,h)
      ! The main subroutine in module mixconductance.
      ! It derives auroral precipitation and conductance caused by both solar EUV and aurora.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      type(gcm_T),optional,intent(in) :: gcm
      integer,optional,intent(in) :: h

      real(rp), dimension(:,:), allocatable :: SigH0,SigP0 !Old values of SigH/SigP
      real(rp) :: dT,upTau,dnTau,wAvgU,wAvgD
      integer :: i,j

      !Save old Sigs
      allocate(SigH0(G%Np,G%Nt))
      allocate(SigP0(G%Np,G%Nt))
      SigP0 = St%Vars(:,:,SIGMAP)
      SigH0 = St%Vars(:,:,SIGMAH)

!      ! Compute EUV though because it's used in fedder
!      call conductance_euv(conductance,G,St)

!      ! Compute auroral precipitation flux
!      call dragonking_total(conductance,G,St)
      
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

      !Before applying cap, do optional exponential smoothing
      if ( (.not. conductance%const_sigma) .and. conductance%doEMA ) then
         !Want 95% of weight to come from last tau seconds
         !Lazily setting values here
         dT = 5.D0

         dnTau = 30.D0 ![s], fall timescale (eg recombination)
         upTau = 5.D0  ![s], increase timescale (eg state averaging)

         wAvgD = 1.D0 - exp(-3.D0*dT/max(dT,dnTau))
         wAvgU = 1.D0 - exp(-3.D0*dT/max(dT,upTau))
         !Throttle how fast conductance drops (ie lazy recombination timescale)
         !$OMP PARALLEL DO default(shared) &
         !$OMP private(i,j)
         do j=1,G%Nt
            do i=1,G%Np
               if(St%Vars(i,j,SIGMAP) < SigP0(i,j)) then
                  !Local conductance dropping
                  St%Vars(i,j,SIGMAP) = wAvgD*St%Vars(i,j,SIGMAP) + (1-wAvgD)*SigP0(i,j)
               else
                  !Local conductance increasing
                  St%Vars(i,j,SIGMAP) = wAvgU*St%Vars(i,j,SIGMAP) + (1-wAvgU)*SigP0(i,j)
               endif
               if(St%Vars(i,j,SIGMAH) < SigH0(i,j)) then
                  !Local conductance dropping
                  St%Vars(i,j,SIGMAH) = wAvgD*St%Vars(i,j,SIGMAH) + (1-wAvgD)*SigH0(i,j)
               else
                  !Local conductance increasing
                  St%Vars(i,j,SIGMAH) = wAvgU*St%Vars(i,j,SIGMAH) + (1-wAvgU)*SigH0(i,j)
               endif
            enddo
         enddo

         deallocate(SigP0)
         deallocate(SigH0)
      endif

      ! Apply cap
      if ((conductance%apply_cap).and.(.not. conductance%const_sigma).and.(.not. present(gcm))) then
         St%Vars(:,:,SIGMAP) = max(conductance%pedmin,St%Vars(:,:,SIGMAP))
         St%Vars(:,:,SIGMAH) = min(max(conductance%hallmin,St%Vars(:,:,SIGMAH)),&
              St%Vars(:,:,SIGMAP)*conductance%sigma_ratio)
      endif

    end subroutine conductance_total

    subroutine conductance_euv(conductance,G,St)
      ! Derive solar EUV conductance: euvSigmaP, euvSigmaH.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      conductance%coszen = G%x*cos(St%tilt)+sqrt(1.-G%x**2-G%y**2)*sin(St%tilt) ! as it should be
      conductance%zenith = acos(conductance%coszen)

      select case ( conductance%euv_model_type )
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

    subroutine conductance_aurora(conductance,G,St)
      ! Use Robinson formula to get Pedersen and Hall conductance from electron precipitation.
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      real(rp), dimension(:,:), allocatable :: engFlux
      allocate(engFlux(G%Np,G%Nt))

      engFlux = kev2erg*St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)  ! Energy flux in ergs/cm^2/s
      conductance%deltaSigmaP = SigmaP_Robinson(St%Vars(:,:,AVG_ENG),engFlux)
      conductance%deltaSigmaH = SigmaH_Robinson(St%Vars(:,:,AVG_ENG),engFlux)
    end subroutine conductance_aurora

  end module mixconductance
