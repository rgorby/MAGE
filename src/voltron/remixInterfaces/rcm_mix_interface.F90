module rcm_mix_interface
  use mixdefs
  use mixgeom
  use volttypes
  use rcm_mhd_interfaces
  
  implicit none

  type(mixGrid_T), private :: rcmG, rcmG_mixstyle, rcmGS

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
    allocate(imag2mix%isClosed(Nt,Np))
    imag2mix%isClosed(:,:) = .false.
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
    call init_grid_fromTP(rcmG,rcmt,rcmp,SOLVER_GRID=.false.)
    call init_grid_fromTP(rcmG_mixstyle,rcmt(1:Np-1,:),rcmp(1:Np-1,:),SOLVER_GRID=.false.)

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
    real(rp),dimension(:,:),allocatable :: rcmEflux_mix,rcmEavg_mix
    real(rp), dimension(:,:), allocatable :: efluxS, eavgS ! for SH mapping. will add ifluxS and iavgS later.
    integer :: ii, jj, kk, Nt, Np

    if ( (.not. imag2mix%isInit) .or. (.not. imag2mix%isFresh) ) return

    !Pull info and do cool stuff here
    ! do mapping here since in geo the RCM grid will be moving
    ! FIXME: if we do RCM in SM, though, this is not necessary (can set map in the init routine above)

    associate(rcmNt=>rcmG_mixstyle%Nt,rcmNp=>rcmG_mixstyle%Np)

    call mix_set_map(rcmG_mixstyle,remixApp%ion(NORTH)%G,rcmMap)
    call mix_map_grids(rcmMap,transpose(imag2mix%eflux(:,1:rcmNp)),rcmEflux_mix)
    call mix_map_grids(rcmMap,transpose(imag2mix%eavg(:,1:rcmNp)),rcmEavg_mix)
    rcmEflux_mix = max(rcmEflux_mix,1.D-10)

    end associate

    remixApp%ion(NORTH)%St%Vars(:,:,IM_EAVG)  = rcmEavg_mix*1e-3 ! [eV -> keV]
    remixApp%ion(NORTH)%St%Vars(:,:,IM_EFLUX) = rcmEflux_mix

    ! Southern Hemisphere Mapping
    ! 1. map eflux and eavg to rcmG* based on NN. Inputs: imag2mix%eflux, eavg, latc, lonc. Outputs: efluxS, eavgS
    call mapIMagToIMagS(imag2mix,efluxS,eavgS)

    ! 2. Use the same interpolation procedure for NH.
    associate(rcmNt=>rcmGS%Nt,rcmNp=>rcmGS%Np)
    call mix_set_map(rcmGS,remixApp%ion(NORTH)%G,rcmMapS)
    call mix_map_grids(rcmMapS,transpose(efluxS(:,1:rcmNp-1)),rcmEflux_mix)
    call mix_map_grids(rcmMapS,transpose(eavgS(:,1:rcmNp-1)), rcmEavg_mix)
    rcmEflux_mix = max(rcmEflux_mix,1.D-10)
    end associate

    ! 3. Like before: set grid size
    associate(Nt=>remixApp%ion(SOUTH)%G%Nt,Np=>remixApp%ion(SOUTH)%G%Np)
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_EAVG)  = rcmEavg_mix(Np:1:-1,:)*1e-3 ! [eV -> keV]
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_EFLUX) = rcmEflux_mix(Np:1:-1,:)

    end associate

! For proton precipitation
    associate(rcmNt=>rcmG_mixstyle%Nt,rcmNp=>rcmG_mixstyle%Np)
    rcmEflux_mix=0.0
    rcmEavg_mix=0.0
    call mix_map_grids(rcmMap,transpose(imag2mix%iflux(:,1:rcmNp)),rcmEflux_mix)
    call mix_map_grids(rcmMap,transpose(imag2mix%iavg(:,1:rcmNp)),rcmEavg_mix)
    end associate
    rcmEflux_mix = max(rcmEflux_mix,1.D-10)
    remixApp%ion(NORTH)%St%Vars(:,:,IM_IAVG)  = rcmEavg_mix*1e-3 ! [eV -> keV]
    remixApp%ion(NORTH)%St%Vars(:,:,IM_IFLUX) = rcmEflux_mix
    associate(Nt=>remixApp%ion(SOUTH)%G%Nt,Np=>remixApp%ion(SOUTH)%G%Np)
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_IAVG)  = rcmEavg_mix(Np:1:-1,:)*1e-3 ! [eV -> keV]
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_IFLUX) = rcmEflux_mix(Np:1:-1,:)
    end associate

    !Set toggle and ignore it until isFresh toggled back
    imag2mix%isFresh = .false.
  end subroutine mapIMagToRemix

  subroutine mapIMagToIMagS(imag2mix,efluxS,eavgS)
    type(imag2Mix_T), intent(in) :: imag2mix
    real(rp), dimension(:,:), allocatable, intent(inout) :: efluxS, eavgS
    real(rp), dimension(:,:), allocatable :: colatc, glongc, xc, yc, xS, yS, rcmt, rcmp
    real(rp) :: dist, distmin, colatmin, colatmax, dlat, rlat, rlon
    integer :: i, j, Np, Nt, i0, j0, NpS, NtS, im, jm

    Nt = size(imag2mix%latc,1)
    Np = size(imag2mix%latc,2) ! imag2mix%latc (Nt,Np)
    if (.not.allocated(colatc)) allocate(colatc(Nt,Np))
    if (.not.allocated(glongc)) allocate(glongc(Nt,Np))
    if (.not.allocated(xc)) allocate(xc(Nt,Np))
    if (.not.allocated(yc)) allocate(yc(Nt,Np))
    colatc = PI/2 + imag2mix%latc
    glongc = imag2mix%lonc
    dlat = (imag2mix%gcolat(Nt)-imag2mix%gcolat(1))/(Nt-1)
    
    ! Step 1. Create a regular grid that covers all latc/lonc. 
! Need to first determine the size in colat by assuming the same lat resolution and using the same lon grid.
    colatmin = PI/2
    colatmax = 0.0
    do j=1,Np
       do i=1,Nt
          if(colatc(i,j)<PI/2) then
             colatmin = min(colatmin,colatc(i,j))
             colatmax = max(colatmax,colatc(i,j))
          end if
          xc(i,j)=sin(colatc(i,j))*cos(glongc(i,j))
          yc(i,j)=sin(colatc(i,j))*sin(glongc(i,j))
       end do
    end do
    NpS = Np
    NtS = ceiling((colatmax-colatmin)/dlat)
    if (.not.allocated(rcmt)) allocate(rcmt(NtS,NpS))
    if (.not.allocated(rcmp)) allocate(rcmp(NtS,NpS))
    if (.not.allocated(xS)) allocate(xS(NtS,NpS))
    if (.not.allocated(yS)) allocate(yS(NtS,NpS))

    do j=1,NpS
       do i=1,NtS
          rcmt(i,j)=colatmin+(i-1)*dlat ! rlat is actually co-lat.
          rcmp(i,j)=imag2mix%glong(j)
       end do
    end do
    call init_grid_fromTP(rcmGS,transpose(rcmt),transpose(rcmp),SOLVER_GRID=.false.)
    xS=transpose(rcmGS%x)
    yS=transpose(rcmGS%y)

    ! Step 3. Return efluxS and eavgS on that expanded regular grid using nearest neighbor.
    ! Find nearest xc/yc for each xS/yS. Can increase precision by decreasing dlat.
    if (.not.allocated(efluxS)) allocate(efluxS(NtS,NpS))
    if (.not.allocated(eavgS))  allocate(eavgS(NtS,NpS))
    do j0=1,NpS
       do i0=1,NtS
          distmin = 1.0D5
          do j=1,Np
             do i=1,Nt
                dist=(xS(i0,j0)-xc(i,j))**2+(yS(i0,j0)-yc(i,j))**2
                if(dist < distmin) then
                   distmin = dist
                   im = i
                   jm = j
                end if
             end do
          end do
          if(distmin<0.025) then !lat cell size: 60deg/200/180*pi=0.005. lon cell size: 360deg/100/180*pi*sin(15-75)=0.016
             efluxS(i0,j0) = imag2mix%eflux(im,jm)
             eavgS(i0,j0)  = imag2mix%eavg(im,jm)
          else
             efluxS(i0,j0) = 0.0
             eavgS(i0,j0)  = 0.0
          endif
       end do
    end do
  end subroutine mapIMagToIMagS


end module  rcm_mix_interface
