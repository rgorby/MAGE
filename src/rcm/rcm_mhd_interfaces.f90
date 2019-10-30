!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> RCM Ionosphere grid definition, coupled model variables and helper
!! functions.
module rcm_mhd_interfaces
  USE rcm_precision
  USE Rcm_mod_subs, ONLY : isize, jsize, jwrap, pi, colat, aloct
  implicit none
  type rcm_mhd_T
  integer(iprec) :: nLat_ion 
  integer(iprec) :: nLon_ion 
  real(rprec),allocatable :: gcolat(:) !> RCM Latitude grid points
  real(rprec),allocatable :: glong(:)  !> RCM Longitude grid points

  real(rprec),allocatable :: pot(:,:)     !> Potential; received from MHD [Volts]
  real(rprec),allocatable :: eng_avg(:,:) !> Average Energy (sent to MIX Coupler/Solver)
  real(rprec),allocatable :: flux(:,:)    !> Energy Flux (sent to MIX Coupler/Solver)
  real(rprec),allocatable :: fac(:,:)     !> Total FAC density (sent to MIX Coupler/Solver)A
  real(rprec),allocatable :: Pave(:,:)    ! MHD supplied average pressure on Pa
  real(rprec),allocatable :: Nave(:,:)   ! MHD supplied average density in #/m^3
  real(rprec),allocatable :: Vol(:,:)     ! MHD supplied flux tube volume, -ve => open fieldline - [(Pa^1/gamma)*m/T]
  real(rprec),allocatable :: X_bmin(:,:,:)! MHD supplied location of Bmin surface, x,y,z in meters
  real(rprec),allocatable :: Bmin(:,:)    ! MHD supplied  bmin strenght in T
  real(rprec),allocatable :: beta_average(:,:)    ! MHD field line averaged plasma beta (\int 2mu0P/B^3ds/B/\int ds/B)
  integer(iprec),allocatable :: iopen(:,:) ! MHD supplied mask open/closed field line (-1: closed; 1: open; 1: else)
  
  real(rprec),allocatable :: Prcm(:,:)    ! RCM supplied pressure on Pa
  real(rprec),allocatable :: Nrcm(:,:)    ! RCM supplied density in #/m^3

  real(rprec),allocatable :: sigmap(:,:)
  real(rprec),allocatable :: sigmah(:,:)

  end type rcm_mhd_T
end module rcm_mhd_interfaces
