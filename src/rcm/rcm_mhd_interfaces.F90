!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> RCM Ionosphere grid definition, coupled model variables and helper
!! functions.
module rcm_mhd_interfaces
    USE kdefs, ONLY : strLen
    USE rcm_precision
    USE Rcm_mod_subs, ONLY : isize, jsize, kcsize,jwrap, doRCMVerbose
    USE rcmdefs, ONLY : RCMTOPCLOSED,RCMTOPNULL,RCMTOPOPEN
    implicit none
    integer(iprec), parameter :: RCMINIT=0,RCMADVANCE=1,RCMRESTART=2,RCMWRITERESTART=-2,RCMWRITEOUTPUT=-3,RCMWRITETIMING=-1
    integer(iprec), parameter :: RCMCOLDSTART=10
    logical :: doColdstart =.true.

    !Scaling parameters
    real(rprec), parameter :: rcmPScl = 1.0e+9 !Convert Pa->nPa
    real(rprec), parameter :: rcmNScl = 1.0e-6 !Convert #/m3 => #/cc
    real(rprec), parameter :: rcmBScl = 1.0e+9 !Convert T to nT

    type rcm_mhd_T
        real(rprec) :: llBC !MHD low-latitude boundary (radians)
        real(rprec) :: dtCpl !Current coupling timescale (can change), [s]
        integer(iprec) :: nLat_ion 
        integer(iprec) :: nLon_ion
        real(rprec) :: planet_radius ! m
        real(rprec) :: iono_radius ! m 
        real(rprec),allocatable :: gcolat(:) !> RCM Latitude grid points
        real(rprec),allocatable :: glong(:)  !> RCM Longitude grid points
        real(rprec),allocatable :: pot(:,:)     !> Potential; received from MHD [Volts]
        real(rprec),allocatable :: eng_avg(:,:,:) !> Average Energy (sent to MIX Coupler/Solver)
        real(rprec),allocatable :: flux(:,:,:)    !> Energy Flux (sent to MIX Coupler/Solver)
        real(rprec),allocatable :: nflx(:,:,:)    !> Number Flux (sent to MIX Coupler/Solver)
        real(rprec),allocatable :: fac(:,:)     !> Total FAC density (sent to MIX Coupler/Solver)A
        real(rprec),allocatable :: Pave(:,:)    ! MHD supplied average pressure on Pa
        real(rprec),allocatable :: Nave(:,:)    ! MHD supplied average density in #/m^3
        real(rprec),allocatable :: Vol(:,:)     ! MHD supplied flux tube volume, -ve => open fieldline - [Re/T]
        real(rprec),allocatable :: X_bmin(:,:,:)! MHD supplied location of Bmin surface, x,y,z in meters
        real(rprec),allocatable :: Bmin(:,:)    ! MHD supplied  bmin strenght in T
        real(rprec),allocatable :: beta_average(:,:)    ! MHD field line averaged plasma beta (\int 2mu0P/B^3ds/B/\int ds/B)
        integer(iprec),allocatable :: iopen(:,:) ! MHD supplied mask open/closed field line (-1: closed; 1: open; 0: else)

        real(rprec),allocatable :: Prcm(:,:)    ! RCM supplied pressure in Pa
        real(rprec),allocatable :: Nrcm(:,:)    ! RCM supplied density in #/m^3
        real(rprec),allocatable :: Npsph(:,:)   ! RCM supplied plasmasphere density in #/m^3
        real(rprec),allocatable :: sigmap(:,:)
        real(rprec),allocatable :: sigmah(:,:)
        real(rprec),allocatable :: oxyfrac(:,:)   ! O+ fraction of MHD number density
        real(rprec),allocatable :: Percm(:,:) ! RCM electron (only) pressure in Pa
        !Conjugate mapping, lat/lon of conjugate point mapped
        real(rprec),allocatable :: latc(:,:)
        real(rprec),allocatable :: lonc(:,:)

        !Field line arc length [Re]
        real(rprec),allocatable :: Lb(:,:)
        !Alfven Bounce timescale [s]
        real(rprec),allocatable :: Tb(:,:)
        !Loss cone size [rad]
        real(rprec),allocatable :: losscone(:,:)
        !Curvature radius [Rp]
        real(rprec),allocatable :: radcurv(:,:)
        !Information about MHD ingestion
        logical, allocatable :: toMHD(:,:)
        !RCM confidence weight, [0,1]
        real(rprec),allocatable :: wIMAG(:,:)

        integer(iprec),allocatable :: nTrc(:,:) !Number of steps on this flux-tube

        !Arrays to hold error in D,P => eta => D',P'. Storing X'/X
        real(rprec), allocatable,dimension(:,:) :: errD,errP
        
        !Information to sync restarts w/ MHD
        integer(iprec) :: rcm_nOut,rcm_nRes !Indices for output/restart
        character(len=strLen) :: rcm_runid
  
        !Some simple quantities for keeping track of RCM energy channels
        real(rprec) :: MaxAlam = 0.0
        
        !Current pressure floor from MHD [nPa]
        real(rprec) :: pFloor = 0.0

        integer(iprec) :: NkT = kcsize !Current number of used channels
    end type rcm_mhd_T

    contains
        !Copy A (RCM/MHD-sized) into B (RCM-sized) and wrap (fill periodic)
        subroutine EmbiggenWrap(rmA,rcmA)
          REAL(rprec), intent(in)    :: rmA (isize,jsize-jwrap+1)
          REAL(rprec), intent(inout) :: rcmA(isize,jsize)

          INTEGER(iprec) :: j

          rcmA(:,jwrap:jsize) = rmA(:,:)
          do j=1,jwrap-1
            rcmA(:,j) = rcmA(:,jsize-jwrap+j)
          enddo

        end subroutine EmbiggenWrap

        !Same as above, but for int
        subroutine EmbiggenWrapI(rmA,rcmA)
          INTEGER(iprec), intent(in)    :: rmA (isize,jsize-jwrap+1)
          INTEGER(iprec), intent(inout) :: rcmA(isize,jsize)

          INTEGER(iprec) :: j

          rcmA(:,jwrap:jsize) = rmA(:,:)
          do j=1,jwrap-1
            rcmA(:,j) = rcmA(:,jsize-jwrap+j)
          enddo

        end subroutine EmbiggenWrapI
        
        !Copy rcmA (RCM-sized) into rmA (RCM/MHD-sized)
        subroutine Unbiggen(rcmA,rmA)
          REAL(rprec), intent(in)       :: rcmA(isize,jsize)
          REAL(rprec), intent(inout)    :: rmA (isize,jsize-jwrap+1)
          
          rmA(:,:) = rcmA(:,jwrap:jsize)

      end subroutine Unbiggen

end module rcm_mhd_interfaces
