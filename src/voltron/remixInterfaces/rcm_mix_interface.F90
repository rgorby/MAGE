module rcm_mix_interface
  use mixdefs
  use mixgeom
  use volttypes
  use rcm_mhd_interfaces
  
  implicit none

  type(mixGrid_T), private :: rcmG, rcmG_mixstyle, rcmGS
  real(rp), dimension(:,:), allocatable :: rcmTopod

contains

  subroutine init_rcm_mix(rcmApp,imag2mix)
    type(rcm_mhd_T),intent(in) :: rcmApp
    type(imag2Mix_T), intent(inout) :: imag2mix
    real(rp), dimension(:,:), allocatable :: rcmp, rcmt ! remix-style 2-D arrays to hold the RCM grid
    integer :: i, j, Np, Nt

    Np = size(rcmApp%glong)
    Nt = size(rcmApp%gcolat)

  !Create imag2mix object
    !Note, imag2mix objects are in RCM order (gcolat,lon)
    allocate(imag2mix%gcolat(Nt))
    allocate(imag2mix%glong (Np))
    imag2mix%gcolat = rcmApp%gcolat
    imag2mix%glong  = rcmApp%glong
    allocate(imag2mix%eflux(Nt,Np))
    allocate(imag2mix%iflux(Nt,Np))
    allocate(imag2mix%eavg (Nt,Np))
    allocate(imag2mix%iavg (Nt,Np))
    allocate(imag2mix%latc (Nt,Np))
    allocate(imag2mix%lonc (Nt,Np))
    allocate(imag2mix%fac  (Nt,Np))
    allocate(imag2mix%inIMag(Nt,Np))
    allocate(imag2mix%inIMagActive(Nt,Np))
    allocate(imag2mix%inIMagBuffer(Nt,Np))
    imag2mix%inIMag(:,:) = .false.
    imag2mix%inIMagActive(:,:) = .false.
    imag2mix%inIMagBuffer(:,:) = .false.
    imag2mix%isInit = .true.

  !Now do remix mapping
    if (.not.allocated(rcmp)) allocate(rcmp(Np,Nt))
    if (.not.allocated(rcmt)) allocate(rcmt(Np,Nt))

    ! construct the 2-D grid
    do j=1,Np
       rcmt(j,:) = rcmApp%gcolat
    enddo

    do i=1,Nt
       rcmp(:,i) = rcmApp%glong
    enddo

    !write(*,*) "===================================",rcmApp%glong(1:3),rcmApp%glong(Np-2:Np),rcmApp%gcolat(1:3),rcmApp%gcolat(Nt-2:Nt)

    ! call remix grid constructor
    call init_grid_fromTP(rcmG,rcmt,rcmp,isSolverGrid=.false.)
    call init_grid_fromTP(rcmG_mixstyle,rcmt(1:Np-1,:),rcmp(1:Np-1,:),isSolverGrid=.false.)

  end subroutine init_rcm_mix

  subroutine map_rcm_mix(voltApp,rcmPot)
    type(voltApp_T), intent(in) :: voltApp 
    real(rp), dimension(:,:), allocatable, intent(inout) :: rcmPot
    type(Map_T) :: rcmMap  

    ! do mapping here since in geo the RCM grid will be moving
    call mix_set_map(voltApp%remixApp%ion(NORTH)%G,rcmG,rcmMap)
    call mix_map_grids(rcmMap,voltApp%remixApp%ion(NORTH)%St%Vars(:,:,POT),rcmPot)
    !Convert from kV to V
    rcmPot = (1.0e+3)*rcmPot

  end subroutine map_rcm_mix

  !Take fluxes from RCM and use for conductance
  subroutine mapIMagToRemix(imag2mix,remixApp)
    type(imag2Mix_T), intent(inout) :: imag2mix
    type(mixApp_T), intent(inout) :: remixApp
    type(Map_T) :: rcmMap, rcmMapS
    real(rp),dimension(:,:),allocatable :: rcmEflux_mix,rcmEavg_mix,rcmTopod_mix
    real(rp), dimension(:,:), allocatable :: efluxS, eavgS, rcmTopodS ! for SH mapping. will add ifluxS and iavgS later.
    integer :: ii, jj, kk, Nt, Np
    integer :: SHmaptype 
    ! # of steps for mapping RCM SH precipitation (may make it an option in XML later): 
    ! 0. direct mirror mapping using NH results; 
    ! 1. map from irregular RCM SH to remix; 
    ! 2. map from irregular RCM SH to a regular equivalent then map to remix as the NH does.
    SHmaptype=1

    if ( (.not. imag2mix%isInit) .or. (.not. imag2mix%isFresh) ) return

    !Pull info and do cool stuff here
    ! do mapping here since in geo the RCM grid will be moving
    ! FIXME: if we do RCM in SM, though, this is not necessary (can set map in the init routine above)

    Nt = size(imag2mix%latc,1)
    Np = size(imag2mix%latc,2) ! imag2mix%latc (Nt,Np)
    if (.not.allocated(rcmTopod)) allocate(rcmTopod(Nt,Np))
    rcmTopod = 0.D0
    where(imag2mix%inIMagActive)
       rcmTopod = 1.0
    elsewhere(imag2mix%inIMagBuffer)
       rcmTopod = 0.5
    endwhere
    print *,'ldong_20211105 rcmTopod',minval(rcmTopod),maxval(rcmTopod)

    call mix_set_map(rcmG_mixstyle,remixApp%ion(NORTH)%G,rcmMap)
    associate(rcmNt=>rcmG_mixstyle%Nt,rcmNp=>rcmG_mixstyle%Np)
    call mix_map_grids(rcmMap,transpose(imag2mix%eflux(:,1:rcmNp)),rcmEflux_mix)
    call mix_map_grids(rcmMap,transpose(imag2mix%eavg(:,1:rcmNp)),rcmEavg_mix)
    call mix_map_grids(rcmMap,transpose(rcmTopod),rcmTopod_mix)
    end associate

    remixApp%ion(NORTH)%St%Vars(:,:,IM_EAVG)  = rcmEavg_mix*1e-3 ! [eV -> keV]
    remixApp%ion(NORTH)%St%Vars(:,:,IM_EFLUX) = rcmEflux_mix     ! [ergs/cm^2/s]
    remixApp%ion(NORTH)%St%Vars(:,:,IM_TOPOD) = rcmTopod_mix

    ! Southern Hemisphere Mapping
    if(SHmaptype==1) then
       call mapIMagSToRemix(imag2mix,remixApp,efluxS,eavgS,rcmTopodS)
       rcmEavg_mix  = transpose(eavgS)
       rcmEflux_mix = transpose(efluxS)
       rcmTopod_mix = transpose(rcmTopodS)
    elseif(SHmaptype==2) then
       call mapIMagSToIMag(imag2mix,efluxS,eavgS) ! need updates to deal with inIMagActive and inIMagBuffer. But SHmaptype=2 is never used.
       call mix_set_map(rcmGS,remixApp%ion(NORTH)%G,rcmMapS)
       associate(rcmNt=>rcmGS%Nt,rcmNp=>rcmGS%Np)
       call mix_map_grids(rcmMapS,transpose(efluxS(:,1:rcmNp-1)),rcmEflux_mix)
       call mix_map_grids(rcmMapS,transpose(eavgS(:,1:rcmNp-1)), rcmEavg_mix)
       end associate
    endif

    associate(Nt=>remixApp%ion(SOUTH)%G%Nt,Np=>remixApp%ion(SOUTH)%G%Np)
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_EAVG)  = rcmEavg_mix(Np:1:-1,:)*1e-3 ! [eV -> keV]
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_EFLUX) = rcmEflux_mix(Np:1:-1,:)
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_TOPOD) = rcmTopod_mix(Np:1:-1,:)
    end associate

! For proton precipitation (all zero for now)
    rcmEflux_mix=0.0
    rcmEavg_mix=0.0
    associate(rcmNt=>rcmG_mixstyle%Nt,rcmNp=>rcmG_mixstyle%Np)
    call mix_map_grids(rcmMap,transpose(imag2mix%iflux(:,1:rcmNp)),rcmEflux_mix)
    call mix_map_grids(rcmMap,transpose(imag2mix%iavg(:,1:rcmNp)),rcmEavg_mix)
    end associate
    remixApp%ion(NORTH)%St%Vars(:,:,IM_IAVG)  = rcmEavg_mix*1e-3 ! [eV -> keV]
    remixApp%ion(NORTH)%St%Vars(:,:,IM_IFLUX) = rcmEflux_mix
    associate(Nt=>remixApp%ion(SOUTH)%G%Nt,Np=>remixApp%ion(SOUTH)%G%Np)
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_IAVG)  = rcmEavg_mix(Np:1:-1,:)*1e-3 ! [eV -> keV]
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_IFLUX) = rcmEflux_mix(Np:1:-1,:)
    end associate

    !Set toggle and ignore it until isFresh toggled back
    imag2mix%isFresh = .false.
  end subroutine mapIMagToRemix

  subroutine mapIMagSToRemix(imag2mix,remixApp,efluxS,eavgS,rcmTopodS)
  ! Directly map from irregular RCM SH grid to ReMIX.
    type(imag2Mix_T), intent(in) :: imag2mix
    type(mixApp_T), intent(inout) :: remixApp
    real(rp), dimension(:,:), allocatable, intent(inout) :: efluxS, eavgS, rcmTopodS
    real(rp), dimension(:,:), allocatable :: colatc, glongc, rcmt, rcmp, Ainvdwgt2
    real(rp) :: dlat, delt, delp, invdwgt
    integer :: i, j, Np, Nt, i0, j0, NpS, NtS, jl, ju, il, iu, jp, dj

    Nt = size(imag2mix%latc,1)
    Np = size(imag2mix%latc,2) ! imag2mix%latc (Nt,Np)
    if (.not.allocated(colatc)) allocate(colatc(Nt,Np))
    if (.not.allocated(glongc)) allocate(glongc(Nt,Np))
    ! Source grid: latc is negative. colatc is positive from ~15 to 75 deg. Note latc=0 for open field lines.
    colatc = PI/2 + imag2mix%latc 
    glongc = imag2mix%lonc

    ! Destination grid: remix Grid.
    rcmt = remixApp%ion(NORTH)%G%t
    rcmp = remixApp%ion(NORTH)%G%p
    NpS  = size(rcmt,1)
    NtS  = size(rcmt,2)
    dlat = rcmt(1,2)-rcmt(1,1)
    dj = nint(dble(NpS)/dble(Np)) ! ratio of rcm dlon to remix dlon.

    ! Mapping: remix dlat is ~10x of rcm, dlon is ~1/3.6 of rcm. Remix lat is from 0-45 deg. RCM is from 15-75 deg.
    ! For each rcm SH point, find the nearest remix lat. If it's not too far away (within dlat) then
    ! find the nearest remix lon. Assign rcm contribution to the nearest lat shell within 2 rcm dlon.
    ! The difference is due to remix dlat is larger while dlon is smaller. Need to make sure all remix grids have some contribution from rcm.
    ! Lastly, normalize the contribution by total IDW.
    if (.not.allocated(efluxS)) allocate(efluxS(NtS,NpS))
    if (.not.allocated(eavgS))  allocate(eavgS(NtS,NpS))
    if (.not.allocated(rcmTopodS))  allocate(rcmTopodS(NtS,NpS))
    if (.not.allocated(Ainvdwgt2))  allocate(Ainvdwgt2(NtS,NpS))
    efluxS = 0.0
    eavgS = 0.0
    rcmTopodS = 0.0
    Ainvdwgt2 = 0.0
    !$OMP PARALLEL DO default(shared) collapse(2) &
    !$OMP private(i,j,i0,il,iu,j0,jl,ju,jp,delt,delp,invdwgt) &
    !$OMP reduction(+:efluxS,eavgS,Ainvdwgt2,rcmTopodS)
    do j=1,Np
       do i=1,Nt
!          if(imag2mix%eflux(i,j)>0.0) then
             i0 = minloc(abs(rcmt(1,:)-colatc(i,j)),1)
             if(rcmt(1,i0)<=colatc(i,j)) then
                il=i0
                iu=min(i0+1,NtS)
             else
                il=max(i0-1,1)
                iu=i0
             endif
             do i0=il,iu 
                if(abs(rcmt(1,i0)-colatc(i,j))<dlat) then
                   jp = minloc(abs(rcmp(:,1)-glongc(i,j)),1)
                   jl = max(jp-dj,1) ! dj used to 2.
                   ju = min(jp+dj,NpS)
                   if(jp<=dj) then  ! The code here may be optimized to be more concise.
                     do j0=NpS-(dj-jp),NpS
                       delt = abs(rcmt(j0,i0)-colatc(i,j))
                       delp = abs((rcmp(j0,i0)-glongc(i,j)))*sin(rcmt(j0,i0))
                       invdwgt = 1./sqrt(delt**2+delp**2)
                       efluxS(i0,j0) = efluxS(i0,j0) + imag2mix%eflux(i,j)*invdwgt
                       eavgS(i0,j0)  = eavgS(i0,j0)  + imag2mix%eavg(i,j)*invdwgt
                       rcmTopodS(i0,j0)  = rcmTopodS(i0,j0)  + rcmTopod(i,j)*invdwgt
                       Ainvdwgt2(i0,j0)  = Ainvdwgt2(i0,j0)  + invdwgt
                     enddo
                   elseif(jp>NpS-dj) then
                     do j0=1,dj-(NpS-jp)
                       delt = abs(rcmt(j0,i0)-colatc(i,j))
                       delp = abs((rcmp(j0,i0)-glongc(i,j)))*sin(rcmt(j0,i0))
                       invdwgt = 1./sqrt(delt**2+delp**2)
                       efluxS(i0,j0) = efluxS(i0,j0) + imag2mix%eflux(i,j)*invdwgt
                       eavgS(i0,j0)  = eavgS(i0,j0)  + imag2mix%eavg(i,j)*invdwgt
                       rcmTopodS(i0,j0)  = rcmTopodS(i0,j0)  + rcmTopod(i,j)*invdwgt
                       Ainvdwgt2(i0,j0)  = Ainvdwgt2(i0,j0)  + invdwgt
                     enddo
                   endif
                   do j0=jl,ju
                      delt = abs(rcmt(j0,i0)-colatc(i,j))
                      delp = abs((rcmp(j0,i0)-glongc(i,j)))*sin(rcmt(j0,i0))
                      invdwgt = 1./sqrt(delt**2+delp**2)
                      efluxS(i0,j0) = efluxS(i0,j0) + imag2mix%eflux(i,j)*invdwgt
                      eavgS(i0,j0)  = eavgS(i0,j0)  + imag2mix%eavg(i,j)*invdwgt
                      rcmTopodS(i0,j0)  = rcmTopodS(i0,j0)  + rcmTopod(i,j)*invdwgt
                      Ainvdwgt2(i0,j0)  = Ainvdwgt2(i0,j0)  + invdwgt
                   enddo
                endif
             enddo
!          endif
       end do
    end do
    !$OMP PARALLEL DO default(shared) collapse(2) &
    !$OMP private(i0,j0)
    do j0=1,NpS
       do i0=1,NtS
          if(Ainvdwgt2(i0,j0)>0.0) then
             efluxS(i0,j0) = efluxS(i0,j0)/Ainvdwgt2(i0,j0)
             eavgS(i0,j0) = eavgS(i0,j0)/Ainvdwgt2(i0,j0)
             rcmTopodS(i0,j0) = rcmTopodS(i0,j0)/Ainvdwgt2(i0,j0)
          endif
       end do
    end do
  end subroutine mapIMagSToRemix

  subroutine mapIMagSToIMag(imag2mix,efluxS,eavgS)
    type(imag2Mix_T), intent(in) :: imag2mix
    real(rp), dimension(:,:), allocatable, intent(inout) :: efluxS, eavgS
    real(rp), dimension(:,:), allocatable :: colatc, glongc, rcmt, rcmp, Ainvdwgt2
    real(rp) :: colatmin, dlat, dlon, delt, delp, invdwgt, Ainvdwgt
    integer :: i, j, Np, Nt, i0, j0, NpS, NtS, il, iu, jl, ju

! RCM NH grid is regular but not uniform. The lat spacing increases toward low latitude (high colat) 
! from 0.015 deg to 0.6 deg, on average (75-15)/200 = 0.3 deg. The lon spacing is uniform, 360/100 = 3.6 deg.
! When constructing a regular grid for SH, may keep the highest lat resolution in the high lat end. The low
! lat end would be no different.

    Nt = size(imag2mix%latc,1)
    Np = size(imag2mix%latc,2) ! imag2mix%latc (Nt,Np)
    if (.not.allocated(colatc)) allocate(colatc(Nt,Np))
    if (.not.allocated(glongc)) allocate(glongc(Nt,Np))
    ! latc is negative. colatc is positive from ~15 to 75 deg. Note latc=0 for open field lines.
    colatc = PI/2 + imag2mix%latc 
    glongc = imag2mix%lonc
     
    ! Step 1. Determine the highest lat/lowest colat in SH conjugate grid.
    colatmin = PI/2
    do j=1,Np
       do i=1,Nt
          if(colatc(i,j)<colatmin) then
             colatmin = colatc(i,j)
          end if
       end do
    end do

    ! Step 2. Create a regular grid that covers all latc/lonc. 
    !     Keep the same size for the SH grid if SH polar cap is larger.
    !     Need to expand for the SH conjugate grid otherwise. 
    !     Note NH colat spacing is 0.07, 0.015, 0.018 for the first four shells.
    !     Use 0.07 deg as dlat for grid expansion. Value of dlat will change later
    dlat = imag2mix%gcolat(2)-imag2mix%gcolat(1) 
    dlon = imag2mix%glong(3)-imag2mix%glong(2)
    NtS = Nt+ceiling(max((imag2mix%gcolat(1)-colatmin),0.0)/dlat)
    NpS = Np
    if (.not.allocated(rcmt)) allocate(rcmt(NpS,NtS))
    if (.not.allocated(rcmp)) allocate(rcmp(NpS,NtS))

    do j=1,NpS
       do i=1,NtS-Nt
          rcmt(j,i)=colatmin+(i-1)*dlat ! Expanded grid.
          rcmp(j,i)=imag2mix%glong(j)
       end do
       do i=NtS-Nt+1,NtS
          rcmt(j,i)=imag2mix%gcolat(i-NtS+Nt)
          rcmp(j,i)=imag2mix%glong(j)
       end do
    end do
    call init_grid_fromTP(rcmGS,rcmt,rcmp,isSolverGrid=.false.)

    ! Step 3. Return efluxS and eavgS on that expanded regular grid using nearest neighbors with inverse distance weighting (IDW).
    !     Neighbors within two grid spacings will be used for IDW average.
    if (.not.allocated(efluxS)) allocate(efluxS(NtS,NpS))
    if (.not.allocated(eavgS))  allocate(eavgS(NtS,NpS))
    if (.not.allocated(Ainvdwgt2))  allocate(Ainvdwgt2(NtS,NpS))
    efluxS = 0.0
    eavgS = 0.0
    Ainvdwgt2 = 0.0
    ! Traverse the irregular grid. Cells in the destination regular grid can be determined which are within 2 spacings.
    !$OMP PARALLEL DO default(shared) collapse(2) &
    !$OMP private(i,j,i0,j0,il,iu,jl,ju,delt,delp,invdwgt) &
    !$OMP reduction(+:efluxS,eavgS,Ainvdwgt2)
    do j=1,Np
       do i=1,Nt
          if(imag2mix%eflux(i,j)>0.0) then
             i0 = minloc(abs(rcmt(1,:)-colatc(i,j)),1)
             il = max(i0-2,1)
             iu = min(i0+2,NtS)
             j0 = minloc(abs(rcmp(:,1)-glongc(i,j)),1)
             jl = max(j0-2,1)
             ju = min(j0+2,NpS)
             do i0=il,iu ! Warning: this range in lat is too broad and would result in latitudinal expansion.
                do j0=jl,ju
                   delt = abs(rcmt(j0,i0)-colatc(i,j))
                   delp = abs((rcmp(j0,i0)-glongc(i,j)))*sin(rcmt(j0,i0))
                   invdwgt = 1./sqrt(delt**2+delp**2)
                   efluxS(i0,j0) = efluxS(i0,j0) + imag2mix%eflux(i,j)*invdwgt
                   eavgS(i0,j0)  = eavgS(i0,j0)  + imag2mix%eavg(i,j)*invdwgt
                   Ainvdwgt2(i0,j0)  = Ainvdwgt2(i0,j0)  + invdwgt
                enddo
             end do
          endif
       end do
    end do
    !$OMP PARALLEL DO default(shared) collapse(2) &
    !$OMP private(i0,j0)
    do j0=1,NpS
       do i0=1,NtS
          if(Ainvdwgt2(i0,j0)>0.0) then
             efluxS(i0,j0) = efluxS(i0,j0)/Ainvdwgt2(i0,j0)
             eavgS(i0,j0) = eavgS(i0,j0)/Ainvdwgt2(i0,j0)
          endif
       end do
    end do
  end subroutine mapIMagSToIMag

end module  rcm_mix_interface
