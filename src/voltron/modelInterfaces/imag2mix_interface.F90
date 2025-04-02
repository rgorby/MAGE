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
         print *,'Np_rai/Nt_rai: ',Np_rai,Nt_rai
         Nt_rai = imagA%raiApp%Grid%shGrid%ie - imagA%raiApp%Grid%shGrid%is+1
         Np_rai = imagA%raiApp%Grid%shGrid%je - imagA%raiApp%Grid%shGrid%js+1
         print *,'is-ie+1/js-je+1: ',Np_rai,Nt_rai

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
        integer :: i1,j1,j1p1
        real(rp), dimension(4) :: F
        real(rp), dimension(:,:), allocatable :: F1
        logical :: isbad = .false.

        ! collect raiju fluxes.
        ! in getMomentsPrecip: allocate(rai_fluxes (is:ie,js:je,nVars_imag2mix)), (Nt_rai, Np_rai, Nf)
        call vApp%imagApp%getMomentsPrecip(rai_fluxes)

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

        ! debug
        allocate(F1(Npc_rai,Nt_rai))
        F1 = transpose(rai_fluxes(:,1:Npc_rai,RAI_EFLUX))
        do i=1,Np_mix
         do j=1,Nt_mix
            if(isnan(mix_fluxes(i,j,RAI_EFLUX))) then
               isbad = .true.
               i1 = raiMap%I1(i,j)
               j1 = raiMap%J1(i,j)
               if (j1.eq.size(F1,1)) then
                  j1p1=1
               else 
                  j1p1=j1+1
               end if
               if (i1 == size(F1,1)) then
                 F = (/ F1(j1,i1), F1(j1p1,i1), F1(j1,i1), F1(j1p1,i1)/)
               elseif (i1 ==0) then
                 F = (/ F1(j1,i1+1), F1(j1p1,i1+1), F1(j1,i1+1), F1(j1p1,i1+1)/)
               else
                 F = (/ F1(j1,i1), F1(j1p1,i1), F1(j1,i1+1), F1(j1p1,i1+1)/)
               endif
               print *,'ldong_20250325 mix_fluxes nan: i=',i,' j=',j,' mix_fluxes=',mix_fluxes(i,j,RAI_EFLUX), &
               ' i1/j1=',i1,j1,' Map=',raiMap%M(i,j,:),' F=',F
            endif
         enddo
        enddo
        if(isbad) then
            print *,'got NaNs in mapRaijuToRemix'
            stop
        endif

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

