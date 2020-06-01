module rcm_mix_interface
  use mixdefs
  use mixgeom
  use volttypes
  use rcm_mhd_interfaces
  
  implicit none

  type(mixGrid_T), private :: rcmG
  type(mixGrid_T), private :: rcmG_mixstyle

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

    write(*,*) "===================================",rcmApp%glong(1:3),rcmApp%glong(Np-2:Np),rcmApp%gcolat(1:3),rcmApp%gcolat(Nt-2:Nt)

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
    type(Map_T) :: rcmMap  
    real(rp),dimension(:,:),allocatable :: rcmEflux_mix,rcmEavg_mix
    integer :: ii, jj, kk

    if ( (.not. imag2mix%isInit) .or. (.not. imag2mix%isFresh) ) return

    !Pull info and do cool stuff here
    ! do mapping here since in geo the RCM grid will be moving
    ! FIXME: if we do RCM in SM, though, this is not necessary (can set map in the init routine above)

!    print *, 'imag2mix%eflux min/max: ',minval(imag2mix%eflux),maxval(imag2mix%eflux)
!    print *, 'imag2mix%eavg  min/max: ',minval(imag2mix%eavg),maxval(imag2mix%eavg)
!    print *, 'imag2mix%iflux min/max: ',minval(imag2mix%iflux),maxval(imag2mix%iflux)
!    print *, 'imag2mix%iavg  min/max: ',minval(imag2mix%iavg),maxval(imag2mix%iavg)

    associate(rcmNt=>rcmG_mixstyle%Nt,rcmNp=>rcmG_mixstyle%Np)

    call mix_set_map(rcmG_mixstyle,remixApp%ion(NORTH)%G,rcmMap)
!    imag2mix%eflux = max(imag2mix%eflux, 1.D-10) ! floor eflux before interpolating.
    call mix_map_grids(rcmMap,transpose(imag2mix%eflux(:,1:rcmNp)),rcmEflux_mix)
    call mix_map_grids(rcmMap,transpose(imag2mix%eavg(:,1:rcmNp)),rcmEavg_mix)
!    where(abs(rcmEflux_mix)<1e-10)
!       rcmEflux_mix = 0.D0
!    end where
!    kk = 0
!    do ii=1,remixApp%ion(NORTH)%G%Np
!      do jj=1,remixApp%ion(NORTH)%G%Nt
!        if(rcmEflux_mix(ii,jj)<0.D0) then
!          write(*,"(a30,1x,i4,i4,e15.7,e15.7)") 'Negative rcmEavg_mix: ',ii,jj,rcmEavg_mix(ii,jj)*1e-3,rcmEflux_mix(ii,jj)
!          kk=kk+1
!        endif
!      enddo
!    enddo
!    print *, 'Number of negative rcmEflux_mix: ', kk
!    print '(a30,e15.7,e15.7,a25,e15.7,e15.7)', 'rcmEflux_mix min/max: ',minval(rcmEflux_mix),maxval(rcmEflux_mix),' rcmEavg_mix min/max: ',minval(rcmEavg_mix)*1e-3,maxval(rcmEavg_mix)*1e-3
!    print '(a30,e15.7,e15.7,a25,e15.7,e15.7)', 'imag2mix%eflux min/max: ',minval(imag2mix%eflux),maxval(imag2mix%eflux),' imag2mix%eavg min/max: ',minval(imag2mix%eavg),maxval(imag2mix%eavg)
!    print *, rcmG_mixstyle%Nt, rcmG_mixstyle%Np, rcmG%Nt, rcmG%Np, shape(imag2mix%eflux), shape(rcmEflux_mix)
! rcmG_mixstyle%Nt=200, rcmG_mixstyle%Np=98,rcmG%Nt=200,rcmG%Np=99, shape(imag2mix%eflux)=200,99, shape(rcmEflux_mix)=360,45.
    rcmEflux_mix = max(rcmEflux_mix,1.D-10)

    end associate

    remixApp%ion(NORTH)%St%Vars(:,:,IM_EAVG)  = rcmEavg_mix*1e-3 ! [eV -> keV] ldong_20200424
    remixApp%ion(NORTH)%St%Vars(:,:,IM_EFLUX) = rcmEflux_mix

    ! FIXME: eventually should map to south along B-field

    ! set grid size
    associate(Nt=>remixApp%ion(SOUTH)%G%Nt,Np=>remixApp%ion(SOUTH)%G%Np)

    remixApp%ion(SOUTH)%St%Vars(:,:,IM_EAVG)  = rcmEavg_mix(Np:1:-1,:)*1e-3 ! [eV -> keV] ldong_20200424
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
    remixApp%ion(NORTH)%St%Vars(:,:,IM_IAVG)  = rcmEavg_mix*1e-3 ! [eV -> keV] ldong_20200424
    remixApp%ion(NORTH)%St%Vars(:,:,IM_IFLUX) = rcmEflux_mix
    associate(Nt=>remixApp%ion(SOUTH)%G%Nt,Np=>remixApp%ion(SOUTH)%G%Np)
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_IAVG)  = rcmEavg_mix(Np:1:-1,:)*1e-3 ! [eV -> keV] ldong_20200424
    remixApp%ion(SOUTH)%St%Vars(:,:,IM_IFLUX) = rcmEflux_mix(Np:1:-1,:)
    end associate

    !Set toggle and ignore it until isFresh toggled back
    imag2mix%isFresh = .false.
  end subroutine mapIMagToRemix

end module  rcm_mix_interface
