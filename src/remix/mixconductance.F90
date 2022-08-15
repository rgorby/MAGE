module mixconductance
  use mixdefs
  use mixtypes
  use earthhelper
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
  real(rp), private :: RinMHD = 0.0 !Rin of MHD grid (0 if not running w/ MHD)
  real(rp), private :: MIXgamma
  logical , private :: doDrift = .false. !Whether to add drift term from Zhang

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
      if (.not. allocated(beta_RCM)) allocate(beta_RCM(G%Np,G%Nt))
      if (.not. allocated(alpha_RCM)) allocate(alpha_RCM(G%Np,G%Nt))
      if (.not. allocated(gtype_RCM)) allocate(gtype_RCM(G%Np,G%Nt))

      call SetMIXgamma(Params%gamma)
      RinMHD = Params%RinMHD
      ! alpha_RCM and alpha_beta replace conductance_alpha/beta in zhang15.
      ! Use default xml input if not using rcmhd.
      alpha_RCM = 1.0/(tiote_RCM+1.0)
      beta_RCM  = conductance%beta
      gtype_RCM = 0.0 ! if conductance_IM_GTYPE is not called, gtype_RCM has all zero, MHD values have a weight of 1.

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

    subroutine conductance_zhang15(conductance,G,St,dorcmO)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      logical, optional, intent(in) :: dorcmO
      
      real(rp) :: signOfY, signOfJ
      logical :: dorcm
      real(rp), dimension(G%Np,G%Nt) :: Pe_MHD, Ne_MHD, Pe_RMD, Ne_RMD
      
      if(present(dorcmO)) then
         dorcm = dorcmO
      else ! default is NOT use RCM thermal flux for mono derivation.
         dorcm = .false. 
      endif

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

      Pe_MHD = 0.1/MIXgamma*alpha_RCM*tmpD*tmpC**2 ! electron pressure from MHD side in [Pa].
      Ne_MHD = tmpD/(Mp_cgs*heFrac)*1.0D6      ! electron number density from MHD side in [/m^3].
      if(.not.dorcm) then ! default Zhang15 using MHD thermal flux only.
         conductance%E0   = 2.0/kev2J*Pe_MHD/Ne_MHD     ! Mean energy from MHD electron fluxes in [keV]. E_avg = 2*kT for Maxwellian.
         conductance%phi0 = beta_RCM*sqrt(Pe_MHD*Ne_MHD/(2.0D-3*pi*Me_cgs))*1.0D-4     ! Thermal number flux from MHD electron fluxes in [#/cm^2/s].
      else ! Trigger this part by including dorcm=.true. when calling zhang15.
         ! similarly, E0 is a ratio and should NOT be merged. 
         ! Derive it from merged pressure and density instead.
         ! See kaiju wiki for derivations of the coefficients.
         ! https://bitbucket.org/aplkaiju/kaiju/wiki/userGuide/derivation_of_precipitation
         Pe_RMD = gtype_RCM*St%Vars(:,:,IM_EPRE) + (1.0-gtype_RCM)*Pe_MHD ! Merged electron pressure in [Pa].
         Ne_RMD = gtype_RCM*St%Vars(:,:,IM_EDEN) + (1.0-gtype_RCM)*Ne_MHD ! Merged electron number density in [/m^3].
         conductance%E0   = 2.0/kev2J*Pe_RMD/Ne_RMD     ! Mean energy from merged electron fluxes in [keV].
         conductance%phi0 = beta_RCM*sqrt(Pe_RMD*Ne_RMD/(2.0D-3*pi*Me_cgs))*1.0D-4     ! Thermal number flux from merged electron fluxes in [#/cm^2/s].
      endif
      
      JF0 = min( 1.D-4*signOfJ*(St%Vars(:,:,FAC)*1.e-6)/eCharge/(conductance%phi0), RM*0.99 )

      !NOTE: conductance%drift should be turned off when using RCM for diffuse
      if (.not. doDrift) then
         conductance%drift = 1.0 !Remove dep.
      endif

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

    subroutine conductance_alpha_beta(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      real(rp), dimension(G%Np,G%Nt) :: phi0_rcmz
      integer :: i,j

      ! In MHD, use the same alpha/beta with RCM.
      ! Calculate beta from RCM fluxes.
      ! Default values are from xml in conductance_init.
      ! It's still necessary to initialize alpha/beta_RCM here because they may be using an old value at some points
      ! if IM_EAVG or IM_EPRE or EM_EDEN no longer satisfies the if criteria. Better to use default background beta.
      alpha_RCM = 1.0/(tiote_RCM+1.0)
      beta_RCM = conductance%beta
      phi0_rcmz = sqrt(St%Vars(:,:,IM_EPRE)*St%Vars(:,:,IM_EDEN)/(Me_cgs*1e-3*2*pi))*1.0e-4 ! sqrt([Pa]*[#/m^3]/[kg]) = sqrt([#/m^4/s^2]) = 1e-4*[#/cm^2/s]
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j)
      do j=1,G%Nt
         do i=1,G%Np
            if(phi0_rcmz(i,j)>TINY) then
               beta_RCM(i,j) = St%Vars(i,j,IM_ENFLX)/phi0_rcmz(i,j)
            elseif(St%Vars(i,j,IM_GTYPE) > 0.5) then
               beta_RCM(i,j) = 0.0
            endif
         enddo
      enddo
      St%Vars(:,:,IM_BETA) = min(beta_RCM,1.0)

    end subroutine conductance_alpha_beta

    subroutine conductance_IM_GTYPE(G,St)
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(in) :: St
      real(rp), dimension(3,3) :: A33
      logical, dimension(3,3) :: isG33 !isG33 = (A3>0.0 .and. A3<1.0)
      integer :: i,j,it,MaxIter,im1,ip1,jm1,jp1
      real(rp) :: mad,Ttmp
      real(rp) :: temp(G%Np,G%Nt)

      MaxIter = 5
      
      ! Iterative diffusion algorithm to smooth out IM_GTYPE: 0/1 are boundary cells. Evolve with nine-cell mean. 
      ! Otherwise it only has three values, 0, 0.5, and 1.0, 
      ! which causes discontinuities in the merged precipitation.
      gtype_RCM = St%Vars(:,:,IM_GTYPE) ! supposed to be between 0 and 1.
      do it=1,MaxIter
         mad = 0.D0 ! max abs difference from last iteration.
         temp = gtype_RCM
         do j=1,G%Nt ! use open BC for lat.
            if(j==1) then
               jm1 = 1
            else
               jm1 = j-1
            endif
            if(j==G%Nt) then
               jp1 = G%Nt
            else
               jp1 = j+1
            endif
            do i=1,G%Np ! use periodic BC for lon.
               if(i==1) then
                  im1 = G%Np
               else
                  im1 = i-1
               endif
               if(i==G%Np) then
                  ip1 = 1
               else
                  ip1 = i+1
               endif
               if(St%Vars(i,j,IM_GTYPE)>0.01 .and. St%Vars(i,j,IM_GTYPE)<0.99) then
               ! Use 0.01/0.99 as the boundary for this numerical diffusion because 
               ! the interpolation from RCM to REMIX would slightly modify the 0/1 boundary.
                  Ttmp =(temp(im1,jm1)+temp(im1,j)+temp(im1,jp1) &
                       + temp(i  ,jm1)+temp(i  ,j)+temp(i  ,jp1) &
                       + temp(ip1,jm1)+temp(ip1,j)+temp(ip1,jp1))/9.D0
                  mad  = max(abs(gtype_RCM(i,j)-Ttmp),mad)
                  gtype_RCM(i,j) = Ttmp
               endif
            enddo
         enddo
         if(mad<0.05) exit
      enddo
      gtype_RCM = min(gtype_RCM,1.0)

    end subroutine conductance_IM_GTYPE

    subroutine conductance_rcmhd(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      integer :: i,j
      logical :: isMono
      real(rp) :: mhd_eavg, mhd_nflx, mhd_eflx, rmd_nflx, rmd_eflx, rmd_eavg
      real(rp) :: rmd_eavg_fin, rmd_nflx_fin, rmd_SigP, mhd_SigP
      logical :: dorcm = .true.

      ! derive RCM grid weighting based on that passed from RCM and smooth it with five iterations of numerical diffusion.
      call conductance_IM_GTYPE(G,St)

      ! derive spatially varying beta using RCM precipitation and thermal fluxes. Need IM_GTYPE.
      call conductance_alpha_beta(conductance,G,St)

      ! derive MHD/mono precipitation with zhang15 but include RCM thermal flux to the source by using dorcm=.true.
      call conductance_zhang15(conductance,G,St,dorcm)

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,isMono) &
      !$OMP private(mhd_eavg, mhd_nflx, mhd_eflx, rmd_nflx, rmd_eflx, rmd_eavg)&
      !$OMP private(rmd_eavg_fin, rmd_nflx_fin, rmd_SigP, mhd_SigP)
      do j=1,G%Nt
         do i=1,G%Np
            isMono = conductance%deltaE(i,j) > 0.0    !Potential drop

            mhd_eavg = St%Vars(i,j,Z_EAVG)
            mhd_nflx = St%Vars(i,j,Z_NFLUX)
            mhd_eflx = St%Vars(i,j,Z_EAVG)*St%Vars(i,j,Z_NFLUX)*kev2erg

            rmd_nflx = St%Vars(i,j,IM_ENFLX)*gtype_RCM(i,j)+mhd_nflx*(1.0-gtype_RCM(i,j))
            rmd_eflx = St%Vars(i,j,IM_EFLUX)*gtype_RCM(i,j)+mhd_eflx*(1.0-gtype_RCM(i,j))
            if(rmd_nflx>TINY) then
               rmd_eavg = max(rmd_eflx/(rmd_nflx*kev2erg),1.0e-8)
            else
               rmd_eavg = 0.0
            endif

            if (.not. isMono) then
               !No potential drop, just use merged precipitation.
               St%Vars(i,j,AVG_ENG ) = rmd_eavg
               St%Vars(i,j,NUM_FLUX) = rmd_nflx
               St%Vars(i,j,AUR_TYPE) = AT_RMnoE
               cycle
            endif

            !If still here, we have a potential drop
            !Decide between the two by taking one that gives highest Sig-P
            if (conductance%doMR) then ! be careful here is using nflux.
               call AugmentMR(rmd_eavg,rmd_nflx,rmd_eavg_fin,rmd_nflx_fin) !Correct for MR
            else
               !No corrections
               rmd_eavg_fin = rmd_eavg
               rmd_nflx_fin = rmd_nflx
            endif

            rmd_SigP = SigmaP_Robinson(rmd_eavg_fin,kev2erg*rmd_eavg_fin*rmd_nflx_fin)
            mhd_SigP = SigmaP_Robinson(mhd_eavg    ,kev2erg*mhd_eavg    *mhd_nflx    )

            if (mhd_SigP>rmd_SigP) then
               St%Vars(i,j,AVG_ENG ) = mhd_eavg
               St%Vars(i,j,NUM_FLUX) = mhd_nflx
               St%Vars(i,j,AUR_TYPE) = AT_RMono
            else
               !RMD diffuse is still better than MHD + puny potential drop
               St%Vars(i,j,AVG_ENG ) = rmd_eavg !Use un-augmented value since MR gets called later
               St%Vars(i,j,NUM_FLUX) = rmd_nflx
               conductance%deltaE(i,j) = 0.0 !Wipe out potential drop since it don't matter (otherwise MR won't happen if desired)
               St%Vars(i,j,AUR_TYPE) = AT_RMfnE
            endif
            
         enddo
      enddo

    end subroutine conductance_rcmhd

    subroutine conductance_rcmonoK(conductance,G,St)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      integer :: i,j
      logical :: isRCM,isMono,isMHD
      real(rp) :: rcm_eavg,rcm_nflx,mhd_eavg,mhd_nflx,rcm_eavg_fin,rcm_nflx_fin,rcm_SigP,mhd_SigP,rcm_eavg0,rcm_nflx0
      real(rp) :: pK = 0.5 !Cut-off for "interesting" energy for precipitation [keV]

      ! If using rcmono, beta will be uniform and specified in xml or default.
      ! RCM thermal flux will NOT be used for mono by setting gtype_RCM = 0.0
      ! to be safe, make beta_RCM=specified beta although conductance_alpha_beta won't be called in rcmono mode.
      beta_RCM  = conductance%beta
      call conductance_zhang15(conductance,G,St)

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,isRCM,isMono,isMHD,rcm_SigP,mhd_SigP) &
      !$OMP private(rcm_eavg,rcm_nflx,mhd_eavg,mhd_nflx,rcm_eavg_fin,rcm_nflx_fin)
      do j=1,G%Nt
         do i=1,G%Np
            isMono = conductance%deltaE(i,j) > 0.0    !Potential drop
            isMHD  = conductance%E0(i,j)     > pK     !Have "hot" MHD information, ie worth putting into Robinson
            ! IM_GTYPE>0.5 covers low lat RCM grid. IM_EAVG>1.0e-8 covers buffre region with meaningful RCM precipitation.
            isRCM  = St%Vars(i,j,IM_GTYPE) > 0.5 .or. St%Vars(i,j,IM_EAVG)>1.0e-8

            !Cases: isMono/isRCM/isMHD = 8 cases
            !- No RCM: 
            if (.not. isRCM) then
               !*/F/* = 4 cases
               !If we don't have RCM info then MHD is the only game in town, use Zhang
               St%Vars(i,j,AVG_ENG ) = St%Vars(i,j,Z_EAVG) ! [keV]
               St%Vars(i,j,NUM_FLUX) = St%Vars(i,j,Z_NFLUX)! [#/cm^2/s]
               St%Vars(i,j,AUR_TYPE) = AT_MHD ! AT_MHD=1,AT_RCM,AT_RMnoE,AT_RMfnE,AT_RMono
               cycle
               !NOTE: Should we handle ~isRCM and ~isMHD case separately?
            endif
            
            !If still here we have RCM information
            ! Rcmono is obsolete here because it assumes rcm passes eflux and eavg,
            !     and it does not involve any merging based on rcm grid.
            !     It only chooses whichever is present and whichever gives higher SigP when both are present.
            !     Keep it here as a historical record. Use rcmhd where nflux is merged.
            rcm_eavg = St%Vars(i,j,IM_EAVG)
            if(St%Vars(i,j,IM_EAVG)>TINY) then
               rcm_nflx = St%Vars(i,j,IM_EFLUX)/(St%Vars(i,j,IM_EAVG)*kev2erg)
            else
               rcm_nflx = 0.0
            endif
            mhd_eavg = St%Vars(i,j,Z_EAVG)
            mhd_nflx = St%Vars(i,j,Z_NFLUX)
            
            !- No (hot) MHD, but RCM info (ie low-lat plasmasphere)
            if (.not. isMHD) then !T/T/F & F/T/F
               !NOTE: Split this case into isMono and apply drop to RCM?
               !For now just use RCM
               St%Vars(i,j,AVG_ENG ) = rcm_eavg
               St%Vars(i,j,NUM_FLUX) = rcm_nflx
               St%Vars(i,j,AUR_TYPE) = AT_RCM
               cycle               
            endif

            !Remaining, have both RCM and (hot) MHD. May or may not have pot drop
            !T/T/T and F/T/T

            if (.not. isMono) then
               !F/T/T
               !Have RCM info and no drop, just use RCM
               St%Vars(i,j,AVG_ENG ) = rcm_eavg
               St%Vars(i,j,NUM_FLUX) = rcm_nflx
               St%Vars(i,j,AUR_TYPE) = AT_RMnoE
               cycle
            endif

            !If still here, we have both RCM info and a potential drop
            !Decide between the two by taking one that gives highest Sig-P
            if (conductance%doMR) then
               call AugmentMR(rcm_eavg,rcm_nflx,rcm_eavg_fin,rcm_nflx_fin) !Correct for MR
            else
               !No corrections
               rcm_eavg_fin = rcm_eavg 
               rcm_nflx_fin = rcm_nflx
            endif

            rcm_SigP = SigmaP_Robinson(rcm_eavg_fin,kev2erg*rcm_eavg_fin*rcm_nflx_fin)
            mhd_SigP = SigmaP_Robinson(mhd_eavg    ,kev2erg*mhd_eavg    *mhd_nflx    )

            if (mhd_SigP>rcm_SigP) then
               St%Vars(i,j,AVG_ENG ) = mhd_eavg
               St%Vars(i,j,NUM_FLUX) = mhd_nflx
               St%Vars(i,j,AUR_TYPE) = AT_RMono
            else
               !RCM diffuse is still better than MHD + puny potential drop
               St%Vars(i,j,AVG_ENG ) = rcm_eavg !Use un-augmented value since MR gets called later
               St%Vars(i,j,NUM_FLUX) = rcm_nflx
               conductance%deltaE(i,j) = 0.0 !Wipe out potential drop since it don't matter (otherwise MR won't happen if desired)
               St%Vars(i,j,AUR_TYPE) = AT_RMfnE
            endif

         enddo
      enddo

    end subroutine conductance_rcmonoK

    subroutine conductance_rcmono(conductance,G,St)
      ! Keep a record of even older rcmono in master as of 20220303. 
      ! Slight changes to AUR_TYPE assignment. Note alpha/beta are effectively diff by 2/gamma and 1/sqrt(2).
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St

      call conductance_zhang15(conductance,G,St)
      ! If there is no potential drop OR mono eflux is too low, use RCM precipitation instead.
      St%Vars(:,:,AUR_TYPE) = AT_MHD
      where(conductance%deltaE<=0.0)
         St%Vars(:,:,AVG_ENG)  = max(St%Vars(:,:,IM_EAVG),1.D-8) ! [keV]
         St%Vars(:,:,NUM_FLUX) = St%Vars(:,:,IM_EFLUX)/(St%Vars(:,:,AVG_ENG)*kev2erg) ! [ergs/cm^2/s]
         St%Vars(:,:,AUR_TYPE) = AT_RCM
         !St%Vars(:,:,Z_NFLUX)  = -1.0 ! for diagnostic purposes since full Z15 does not currently work.
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

      ! note, this assumes that fedder95/zhang15/rcmhd/rcmono/rcmfed has been called prior
      ! **********************************************************************************
      ! ********** Diffuse precipitation from RCM has been divided by 2 in ***************
      ! ********** rcm_subs.F90/subroutine kdiffPrecip for each hemisphere. **************
      ! **********************************************************************************
      conductance%engFlux = kev2erg*St%Vars(:,:,AVG_ENG)*St%Vars(:,:,NUM_FLUX)  ! Energy flux in ergs/cm^2/s
      conductance%deltaSigmaP = SigmaP_Robinson(St%Vars(:,:,AVG_ENG),conductance%engFlux)
      conductance%deltaSigmaH = SigmaH_Robinson(St%Vars(:,:,AVG_ENG),conductance%engFlux)

    end subroutine conductance_aurora

    ! George Khazanov's multiple reflection(MR) corrections
    subroutine conductance_mr(conductance,G,St)
      type(mixConductance_T), intent(in) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      real(rp), dimension(G%Np,G%Nt) :: Kc

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
   subroutine GenMirrorRatio(G,St,doIGRFO)
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

   subroutine conductance_total(conductance,G,St,gcm,h)
      type(mixConductance_T), intent(inout) :: conductance
      type(mixGrid_T), intent(in) :: G
      type(mixState_T), intent(inout) :: St
      type(gcm_T),optional,intent(in) :: gcm
      integer,optional,intent(in) :: h

      call GenMirrorRatio(G,St)

      ! always call fedder to fill in AVG_ENERGY and NUM_FLUX
      ! even if const_sigma, we still have the precip info that way

      ! compute EUV though because it's used in fedder
      call conductance_euv(conductance,G,St)
      select case ( conductance%aurora_model_type )
         case (FEDDER)
            call conductance_fedder95(conductance,G,St)
         case (ZHANG)
            doDrift = .true.
            call conductance_zhang15(conductance,G,St)
         case (RCMONO)
            doDrift = .false.
            call conductance_rcmono(conductance,G,St)
         case (RCMHD)
            doDrift = .false.
            call conductance_rcmhd(conductance,G,St)
         case (RCMFED)
            call conductance_rcmfed(conductance,G,St)
         case default
            stop "The aurora precipitation model type entered is not supported."
      end select

      ! correct for multiple reflections if you're so inclined
      if (conductance%doMR) call conductance_mr(conductance,G,St)

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

    !Routine to change MIX-gamma on the fly if necessary
    subroutine SetMIXgamma(gamma)
      real(rp), intent(in) :: gamma
      MIXgamma = gamma
    end subroutine SetMIXgamma

  end module mixconductance
