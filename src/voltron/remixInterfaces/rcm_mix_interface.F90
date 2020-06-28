module rcm_mix_interface
  use mixdefs
  use mixgeom
  use volttypes
  use rcm_mhd_interfaces
  
  implicit none

  type(mixGrid_T), private :: rcmG
  type(Map_T), private :: rcmMap  

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

    ! call remix grid constructor
    call init_grid_fromTP(rcmG,rcmt,rcmp,SOLVER_GRID=.false.)
  end subroutine init_rcm_mix

  subroutine map_rcm_mix(voltApp,rcmPot)
    type(voltApp_T), intent(in) :: voltApp 
    real(rp), dimension(:,:), allocatable, intent(inout) :: rcmPot

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

    if ( (.not. imag2mix%isInit) .or. (.not. imag2mix%isFresh) ) return

    !Pull info and do cool stuff here

    !Set toggle and ignore it until isFresh toggled back
    imag2mix%isFresh = .false.
  end subroutine mapIMagToRemix

end module  rcm_mix_interface
