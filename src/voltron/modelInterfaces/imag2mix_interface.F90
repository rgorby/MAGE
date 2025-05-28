! module to convert imag data to mhd data
module imag2mix_interface
    use mixdefs
    use mixgeom
    use volttypes
    use cmiutils
    use planethelper
    use gridutils
    use shellInterp

    implicit none

    type(mixGrid_T), private :: rai2mixG
    integer, private :: Np_mix, Nt_mix, Np_rai, Nt_rai, Npc_rai, Ntc_rai
    logical, private :: isInit = .true.

    contains

    subroutine init_raiju_mix(imagApp,remixApp)
      ! called by subroutine initializeFromGamera in module voltapp
!      type(voltApp_T), intent(inout) :: vApp
      class(imagCoupler_T), intent(in) :: imagApp
      type(mixApp_T), intent(inout) :: remixApp

      real(rp), dimension(:,:), allocatable :: raijup, raijut ! Np x Nt, remix-style 2-D arrays to hold the RAIJU grid
      integer :: i, j

      !      associate(remixApp=>vApp%remixApp, imagApp=>vApp%imagApp )
      Nt_mix = remixApp%ion(NORTH)%G%Nt
      Np_mix = remixApp%ion(NORTH)%G%Np

      SELECT TYPE (imagA=>imagApp)
         TYPE IS (raijuCoupler_T)
         ! in here you can treat imagType as if it is type raijuCoupler_T, and it points to vApp%imagApp
         Np_rai = imagA%raiApp%Grid%shGrid%Np
         Nt_rai = imagA%raiApp%Grid%shGrid%Nt

         ! Np x Nt, transposed relative to mix grid.
         if (.not.allocated(raijup)) allocate(raijup(Np_rai, Nt_rai))
         if (.not.allocated(raijut)) allocate(raijut(Np_rai, Nt_rai))

         ! construct the 2-D grid
         !! thc/phc: (Nt or Np) [radians] grid centers
         !! th (theta) is colatitude and runs from north pole toward south
         do j=1,Np_rai
           raijut(j,:) = imagA%raiApp%Grid%shGrid%thc
         enddo

         !! Phi is longitude, with zero/2pi at 12 MLT
         do i=1,Nt_rai
           raijup(:,i) = imagA%raiApp%Grid%shGrid%phc
         enddo

         ! call remix grid constructor
         call init_grid_fromTP(rai2mixG,raijut(1:Np_rai-1,:),raijup(1:Np_rai-1,:),isSolverGrid=.false.)
         Npc_rai = rai2mixG%Np ! Np_rai-1
         Ntc_rai = rai2mixG%Nt
      CLASS DEFAULT
         WRITE (*,*) "Imag Coupler is an unsupported type"
         stop
      END SELECT
!      end associate

    end subroutine init_raiju_mix

   subroutine CoupleIMagToMix(vApp)
      class(voltApp_T), intent(inout) :: vApp
      integer :: i,j
      real(rp) :: x1,x2,x1c,x2c
      logical :: isTasty
      real(rp) :: imP(nVars_imag2mix)
      integer, dimension(9) :: varNames = (/IM_EAVG,IM_ENFLX,IM_EFLUX,IM_GTYPE,IM_EDEN,IM_EPRE,IM_NPSP,IM_THCON,IM_PHCON/)

      ! Zero out before conditional updates
      do i=1,size(varNames)
         vApp%remixApp%ion(NORTH)%St%Vars(:,:,varNames(i)) = 0.0
         vApp%remixApp%ion(SOUTH)%St%Vars(:,:,varNames(i)) = 0.0
      enddo
      ! G%t(G%Np,G%Nt), G%p = modulo((atan2(G%y,G%x)+2*pi),(2*pi)) 
      ! G%x = sin(G%t)*cos(G%p)
      ! G%y = sin(G%t)*sin(G%p)
      associate(G=>vApp%remixApp%ion(NORTH)%G, ion=>vApp%remixApp%ion)
         !Loop over and get imag data
         !$OMP PARALLEL DO default(shared) collapse(2) &
         !$OMP private(i,j,x1,x2,x1c,x2c,imP,isTasty)
            do j=1,G%Np
               do i=1,G%Nt
                  ! NH mapping
                  x1 = G%t(j,i) ! G%t is colatitude.
                  x2 = G%p(j,i) ! need to be [0,2pi]
                  imP = 0.0_rp
                  ! Here we use the global voltron grid to first find the conjugate point in the SH.
                  call InterpShellVar_TSC_pnt(vApp%shGrid, vApp%State%tubeShell%latc, x1, x2, x1c)
                  call InterpShellVar_TSC_pnt(vApp%shGrid, vApp%State%tubeShell%lonc, x1, x2, x2c)
                  ion(NORTH)%St%Vars(j,i,IM_THCON) = PI/2.0_rp - x1c ! conjugate co-lat in radians, 0-pi ! State%thcon(i0,j0)
                  ion(NORTH)%St%Vars(j,i,IM_PHCON) = x2c ! conjugate long in radians, 0-2pi ! State%phcon(i0,j0)
                  !call vApp%imagApp%getMomentsPrecip(x1,x2,imP,isTasty)
                  call superSampleMomentsPrecip(G, i, j, 5, imP, isTasty)
                  ! gtype is 1 when isTasty is true, but can be 0 or 0.5 when isTasty is false.
                  ion(NORTH)%St%Vars(j,i,IM_GTYPE) = imP(RAI_GTYPE)      ! [0~1]
                  if (isTasty) then
                     ion(NORTH)%St%Vars(j,i,IM_EAVG ) = imP(RAI_EAVG )      ! [keV]
                     ion(NORTH)%St%Vars(j,i,IM_ENFLX) = imP(RAI_ENFLX)      ! [#/cm^2/s]
                     ion(NORTH)%St%Vars(j,i,IM_EFLUX) = imP(RAI_EFLUX)      ! [ergs/cm^2/s]
                     ion(NORTH)%St%Vars(j,i,IM_EDEN ) = imP(RAI_EDEN )      ! [#/m^3]
                     ion(NORTH)%St%Vars(j,i,IM_EPRE ) = imP(RAI_EPRE )      ! [Pa]
                     ion(NORTH)%St%Vars(j,i,IM_NPSP ) = imP(RAI_NPSP )      ! [#/m^3]
                  endif

                  ! SH mapping
                  ! Earlier we used a simple mirro mapping, i.e.:
                  ! ion(SOUTH)%St%Vars(j,i,IM_EAVG ) = ion(NORTH)%St%Vars(G%Np:1:-1,j,i,IM_EAVG )
                  x1 = PI-G%t(G%Np+1-j,i)     ! G%t is colatitude. SH colat should be 0.75pi-pi if remix llb is 0.25pi (45 deg).
                  x2 = G%p(G%Np+1-j,i)        ! need to be [0,2pi]
                  imP = 0.0_rp
                  ! Here we use the global voltron grid to first find the conjugate point in the NH.
                  call InterpShellVar_TSC_pnt(vApp%shGrid, vApp%State%tubeShell%latc, x1, x2, x1c)
                  call InterpShellVar_TSC_pnt(vApp%shGrid, vApp%State%tubeShell%lonc, x1, x2, x2c)
                  ion(SOUTH)%St%Vars(j,i,IM_THCON) = PI/2.0_rp - x1c ! conjugate co-lat in radians, 0-pi ! State%thcon(i0,j0)
                  ion(SOUTH)%St%Vars(j,i,IM_PHCON) = x2c ! conjugate long in radians, 0-2pi ! State%phcon(i0,j0)

                  ! Note x1c is the conjugate lat in the NH (rather than colat). 
                  ! x1c is returned as 0 for open field lines. 
                  ! Only interpolate if x1c is within the remix grid.
                  if (x1c>=PI/2.0_rp - G%t(1,G%Nt)) then
                     !call vApp%imagApp%getMomentsPrecip(PI/2.0_rp-x1c,x2c,imP,isTasty)
                     call superSampleMomentsPrecip(G, i, j, 5, imP, isTasty, thetaConO=PI/2.0_rp-x1c, phiConO=x2c)
                     ! gtype is 1 when isTasty is true, but can be 0 or 0.5 when isTasty is false.
                     ion(SOUTH)%St%Vars(j,i,IM_GTYPE) = imP(RAI_GTYPE)      ! [0~1]
                     if (isTasty) then
                        ion(SOUTH)%St%Vars(j,i,IM_EAVG ) = imP(RAI_EAVG )      ! [keV]
                        ion(SOUTH)%St%Vars(j,i,IM_ENFLX) = imP(RAI_ENFLX)      ! [#/cm^2/s]
                        ion(SOUTH)%St%Vars(j,i,IM_EFLUX) = imP(RAI_EFLUX)      ! [ergs/cm^2/s]
                        ion(SOUTH)%St%Vars(j,i,IM_EDEN ) = imP(RAI_EDEN )      ! [#/m^3]
                        ion(SOUTH)%St%Vars(j,i,IM_EPRE ) = imP(RAI_EPRE )      ! [Pa]
                        ion(SOUTH)%St%Vars(j,i,IM_NPSP ) = imP(RAI_NPSP )      ! [#/m^3]
                     endif
                  endif
               enddo !i
            enddo !j
      end associate


      contains

      subroutine superSampleMomentsPrecip(G, i, j, nPnts, imP_avg, isTasty, thetaConO, phiConO)
         !! Super sample a remix cell to account for disparate resolutions between remix and imag model
         !! Should ultimately be replaced with a flux-conserving interpolation method eventually
         !! Note: For conjugate points, we are assuming the same cell size as the remix point's origin, which can have some issues
         !!   In order to do better cell area, need to calculate all conjugate points first and then pass all neighbors, do fancy stuff
         type(mixGrid_T), intent(in) :: G
         !class(voltApp_T), intent(in) :: vApp  ! Some sort of intent complaint, now we just use it from parent scope
         integer, intent(in) :: i, j
         integer, intent(in) :: nPnts
         real(rp), intent(inout) :: imP_avg(nVars_imag2mix)
         logical, intent(inout) :: isTasty
         real(rp), intent(in), optional :: thetaConO, phiConO
            !! If wanting conjugate points, provide them here. Otherwise will use G%t and G%p at point i,j

         logical :: doSouth
         integer :: n, nGood
         integer :: iPnt, jPnt
         logical :: isTastyPnt
         real(rp) :: thetaIn, phiIn
            !! Theta and phi vals we are evaluating (either G%t or thetaConO, and same for phi)
         real(rp) :: thMin, thMax, phMin, phMax
         real(rp), dimension(:), allocatable :: theta_subdiv, phi_subdiv
         real(rp) :: imP_pnt(nVars_imag2mix)

         imP_avg = 0.0
         nGood = 0
         if (nPnts < 2) then
            write(*,*)"ERROR in imag2mix_interface.F90:superSampleMomentsPrecip"
            write(*,*)"nPnts must be at least 2, preferrably higher"
            stop
         endif

         if (present(thetaConO) .and. present(phiConO)) then
            thetaIn = thetaConO
            phiIn = phiConO
         else
            thetaIn = G%t(j,i)
            phiIn = G%p(j,i)
         endif

         allocate(theta_subdiv(nPnts))
         allocate(phi_subdiv  (nPnts))

         ! Calculate remix cell bounds for index i,j
         ! Lazy handling of i bounds
         if (i == 1 .or. i == G%Nt) then
            theta_subdiv = thetaIn
         else
            thMin = thetaIn - 0.5_rp*G%dt(j,i-1)
            thMax = thetaIn + 0.5_rp*G%dt(j,i  )
            theta_subdiv = genThetaSubdiv(thMin, thMax, nPnts)
         endif
         ! Gen phi vals
         if (j == 1) then
            phMin = phiIn - 0.5_rp*G%dp(G%Np,i)
         else
            phMin = phiIn - 0.5_rp*G%dp(j-1,i)
         endif
         phMax = phiIn + 0.5_rp*G%dp(j  ,i)
         do n=1,nPnts
            phi_subdiv = phMin + (n-1)*1.0_rp/(nPnts-1)*(phMax - phMin)  ! Don't forget to cast something to real so we don't do integer division
         enddo

         ! Now that we have our points we can start getting our values
         do iPnt=1,nPnts
            do jPnt=1,nPnts
               call vApp%imagApp%getMomentsPrecip(theta_subdiv(iPnt), phi_subdiv(jPnt), imP_pnt, isTastyPnt)
               if (.not. isTastyPnt) then
                  cycle
               endif
               ! Otherwise we got good data
               nGood = nGood + 1
               imP_avg = imP_avg + imP_pnt
            enddo
         enddo
         ! Finish average and make final isTasty decision
         if (nGood > 0) then
            isTasty = .true.

            ! For some variables, don't divide by nGood, divide by # of all points sampled
            ! e.g. in case where 1/nPnts^2 were good it gets applied to whole cell when really the rest are telling you there should be no precip there
            ! So you actually to wanna include some effect from bad points
            !IM_EAVG,IM_ENFLX,IM_EFLUX,IM_GTYPE,IM_EDEN,IM_EPRE,IM_NPSP,IM_THCON,IM_PHCON
            imP_avg(RAI_ENFLX) = imP_avg(RAI_ENFLX)/nPnts**2
            imP_avg(RAI_EFLUX) = imP_avg(RAI_EFLUX)/nPnts**2
            imP_avg(RAI_EDEN)  = imP_avg(RAI_EDEN )/nPnts**2
            imP_avg(RAI_EPRE)  = imP_avg(RAI_EPRE )/nPnts**2
            imP_avg(RAI_NPSP)  = imP_avg(RAI_NPSP )/nPnts**2
            imP_avg(RAI_EAVG)  = imP_avg(RAI_EFLUX) / imP_avg(RAI_ENFLX)
            imP_avg(RAI_GTYPE) = imP_avg(RAI_GTYPE)/nGood
            imP_avg(RAI_THCON) = imP_avg(RAI_THCON)/nGood
            imP_avg(RAI_PHCON) = imP_avg(RAI_PHCON)/nGood
            
         else
            isTasty = .false.
         endif

      end subroutine superSampleMomentsPrecip

      function genThetaSubdiv(theta0, theta1, nPnts) result(theta_subdiv)
         !! Return theta array between theta0 and theta1 where areas of cells will be equal
         !! Expecting Theta0 < Theta1
         real(rp), intent(in) :: theta0, theta1
         integer, intent(in) :: nPnts

         integer :: i
         real(rp) :: dCos, arg
         real(rp), dimension(nPnts) :: theta_subdiv

         dCos = cos(theta0) - cos(theta1)

         theta_subdiv(1) = theta0
         theta_subdiv(nPnts) = theta1
         do i=2,nPnts-1
            arg = cos(theta0) - (i-1)*1.0_rp/(nPnts-1)*dCos  ! Don't forget to cast something to real so we don't do integer division
            theta_subdiv(i) = acos(arg)
         enddo

      end function genThetaSubdiv

   end subroutine CoupleIMagToMix

    subroutine mapRaijuToRemix(vApp)
      type(voltApp_T), intent(inout) :: vApp
!        type(raijuCoupler_T), intent(inout) :: imagApp
!        type(mixApp_T), intent(inout) :: remixApp

        real(rp), dimension(:,:,:), allocatable :: rai_fluxes, mix_fluxes
        real(rp), dimension(:,:), allocatable :: mix_flux
        real(rp), dimension(:,:), allocatable :: mixt, mixp
        real(rp), dimension(:,:), allocatable :: raijup, raijut ! Np x Nt, remix-style 2-D arrays to hold the RAIJU grid
        real(rp), dimension(:), allocatable :: phc, thc
        integer :: Nf = nVars_imag2mix
        integer :: i,j
        type(Map_T) :: raiMap

        ! collect raiju fluxes.
        ! in getMomentsPrecip: allocate(rai_fluxes (is:ie,js:je,nVars_imag2mix)), (Nt_rai, Np_rai, Nf)
!        call vApp%imagApp%getMomentsPrecip(rai_fluxes)

        associate(remixApp=>vApp%remixApp ) !, imagApp=>vApp%imagApp
        allocate(mix_fluxes(Np_mix,Nt_mix,Nf))
        allocate(mix_flux(Np_mix,Nt_mix))
        mix_fluxes = 0.0_rp
        mix_flux   = 0.0_rp

        mixt = remixApp%ion(NORTH)%G%t ! G%t(G%Np,G%Nt)
        mixp = remixApp%ion(NORTH)%G%p
        
        ! NH mapping to remix
        call mix_set_map(rai2mixG,remixApp%ion(NORTH)%G,raiMap)
        do i=1,Nf
          call mix_map_grids(raiMap,transpose(rai_fluxes(:,1:Npc_rai,i)), mix_flux)
          mix_fluxes(:,:,i) = mix_flux
        enddo

        remixApp%ion(NORTH)%St%Vars(:,:,IM_EAVG ) = mix_fluxes(:,:,RAI_EAVG )      ! [keV]
        remixApp%ion(NORTH)%St%Vars(:,:,IM_ENFLX) = mix_fluxes(:,:,RAI_ENFLX)      ! [#/cm^2/s]
        remixApp%ion(NORTH)%St%Vars(:,:,IM_EFLUX) = mix_fluxes(:,:,RAI_EFLUX)      ! [ergs/cm^2/s]
        remixApp%ion(NORTH)%St%Vars(:,:,IM_GTYPE) = mix_fluxes(:,:,RAI_GTYPE)      ! [0~1]
        remixApp%ion(NORTH)%St%Vars(:,:,IM_EDEN ) = mix_fluxes(:,:,RAI_EDEN )      ! [#/m^3]
        remixApp%ion(NORTH)%St%Vars(:,:,IM_EPRE ) = mix_fluxes(:,:,RAI_EPRE )      ! [Pa]
        remixApp%ion(NORTH)%St%Vars(:,:,IM_NPSP ) = mix_fluxes(:,:,RAI_NPSP )      ! [#/m^3]

        ! SH mapping
        mix_fluxes = 0.0_rp
        call mapRaijuSToRemix(rai_fluxes,mixt,mixp,mix_fluxes)
        remixApp%ion(SOUTH)%St%Vars(:,:,IM_EAVG ) = mix_fluxes(Np_mix:1:-1,:,RAI_EAVG )      ! [keV]
        remixApp%ion(SOUTH)%St%Vars(:,:,IM_ENFLX) = mix_fluxes(Np_mix:1:-1,:,RAI_ENFLX)      ! [#/cm^2/s]
        remixApp%ion(SOUTH)%St%Vars(:,:,IM_EFLUX) = mix_fluxes(Np_mix:1:-1,:,RAI_EFLUX)      ! [ergs/cm^2/s]
        remixApp%ion(SOUTH)%St%Vars(:,:,IM_GTYPE) = mix_fluxes(Np_mix:1:-1,:,RAI_GTYPE)      ! [0~1]
        remixApp%ion(SOUTH)%St%Vars(:,:,IM_EDEN ) = mix_fluxes(Np_mix:1:-1,:,RAI_EDEN )      ! [#/m^3]
        remixApp%ion(SOUTH)%St%Vars(:,:,IM_EPRE ) = mix_fluxes(Np_mix:1:-1,:,RAI_EPRE )      ! [Pa]
        remixApp%ion(SOUTH)%St%Vars(:,:,IM_NPSP ) = mix_fluxes(Np_mix:1:-1,:,RAI_NPSP )      ! [#/m^3]
        end associate
    end subroutine mapRaijuToRemix

    subroutine mapRaijuSToRemix(rai_fluxes,mixt,mixp,mix_fluxes)
      ! Directly map from irregular raiju SH grid to ReMIX.
      real(rp), dimension(Nt_rai, Np_rai,nVars_imag2mix), intent(in)  :: rai_fluxes
      real(rp), dimension(Np_mix, Nt_mix), intent(in) :: mixt, mixp
      real(rp), dimension(Np_mix, Nt_mix,nVars_imag2mix), intent(out) :: mix_fluxes
      
      real(rp), dimension(Nt_rai, Np_rai) :: colatc, glongc
      real(rp), dimension(:,:), allocatable :: mixtE, mixpE
      real(rp), dimension(Np_mix, Nt_mix) :: Ainvdwgt2
      real(rp) :: dlat, delt, delp, invdwgt
      integer :: i, j, i0, j0, il, iu, jp, dj
  
      ! Source grid: latc is negative. colatc is positive from ~15 to 75 deg. Note latc=0 for open field lines.
      colatc = PI-rai_fluxes(:,:,RAI_THCON) ! RAI_THCON is conjugate co-lat in radians pi/2-pi, PI - RAI_THCON is -> 0-pi/2
      glongc = rai_fluxes(:,:,RAI_PHCON) ! conjugate long in radians, 0-2pi, need to double check if they are consistent at lon=0
  
      ! Destination grid: remix Grid.
      dlat = mixt(1,2)-mixt(1,1)
      ! dj is the ratio of rcm dlon to remix dlon, i.e. min number of rcm/raiju cells to collect.
      ! now with 360 raiju/rcm longitudinal cells, dj is 1 with quad res remix grid.
      dj = nint(dble(Np_mix)/dble(Np_rai))
      ! make an extended mixt/mixp grid for easier periodic boundary processing.
      allocate(mixtE(1-dj:Np_mix+dj,1:Nt_mix))
      allocate(mixpE(1-dj:Np_mix+dj,1:Nt_mix))
      mixtE(1:Np_mix,:) = mixt
      mixtE(1-dj:0,:) = mixt(Np_mix+1-dj:Np_mix,:)
      mixtE(Np_mix+1:Np_mix+dj,:) = mixt(1:dj,:)
      mixpE(1:Np_mix,:) = mixp
      mixpE(1-dj:0,:) = mixp(Np_mix+1-dj:Np_mix,:)
      mixpE(Np_mix+1:Np_mix+dj,:) = mixp(1:dj,:)
  
      ! Mapping: remix dlat is ~10x of rcm, dlon is ~1/3.6 of rcm. Remix lat is from 0-45 deg. RCM is from 15-75 deg.
      ! For each rcm SH point, find the nearest remix lat. If it's not too far away (within dlat) then
      ! find the nearest remix lon. Assign rcm contribution to the nearest lat shell within 2 rcm dlon.
      ! The difference is due to remix dlat is larger while dlon is smaller. Need to make sure all remix grids have some contribution from rcm.
      ! Lastly, normalize the contribution by total IDW.
      Ainvdwgt2 = 0.0_rp
      mix_fluxes = 0.0_rp
      !$OMP PARALLEL DO default(shared) collapse(2) &
      !$OMP private(i,j,i0,il,iu,j0,jp,delt,delp,invdwgt) &
      !$OMP reduction(+:Ainvdwgt2,mix_fluxes)
      do j=1,Np_rai
         do i=1,Nt_rai
            i0 = minloc(abs(mixt(1,:)-colatc(i,j)),1) ! Find the nearest remix colat index for rcm colatc(i,j)
            if(mixt(1,i0)<=colatc(i,j)) then ! If the nearest remix colat is < rcm colatc, only collect rcm to this colat and its next grid.
               il=i0
               iu=min(i0+1,Nt_mix)
            else ! Otherwise, collect from this point and its neighbor lat.
               il=max(i0-1,1)
               iu=i0
            endif
            do i0=il,iu 
               ! For any remix grid, interpolate if rcm lat is within dlat away
               if(abs(mixt(1,i0)-colatc(i,j))<dlat) then 
                  jp = minloc(abs(mixp(:,1)-glongc(i,j)),1)
                  ! 1 <= jp <= Np_mix
                  ! jp-dj>= 1-dj; jp+dj<= Np_mix+dj
                  ! mixtE/mixpE is from 1-dj:Np_mix+dj
                  do j0=jp-dj,jp+dj
                     delt = abs(mixtE(j0,i0)-colatc(i,j))
                     delp = abs((mixpE(j0,i0)-glongc(i,j)))*sin(mixtE(j0,i0))
                     invdwgt = 1.0_rp/sqrt(delt**2+delp**2)
                     mix_fluxes(j0,i0,:) = mix_fluxes(j0,i0,:) + rai_fluxes(i,j,:)*invdwgt
                     Ainvdwgt2(j0,i0) = Ainvdwgt2(j0,i0) + invdwgt
                  enddo
               endif
            enddo
  !          endif
         enddo
      enddo
      !$OMP PARALLEL DO default(shared) collapse(2) &
      !$OMP private(i0,j0)
      do j0=1,Np_mix
         do i0=1,Nt_mix
            if(Ainvdwgt2(j0,i0)>0.0_rp) then
               mix_fluxes(j0,i0,:) = mix_fluxes(j0,i0,:)/Ainvdwgt2(j0,i0)
            endif
         enddo
      enddo
    end subroutine mapRaijuSToRemix

end module

