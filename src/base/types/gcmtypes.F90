module gcmtypes
  use mixtypes
  
  implicit none
  
  integer,parameter :: gcm2mix_nvar = 2 ! SigmaP, SigmaH, Neutral Winds(eventually)
  integer,parameter :: mix2gcm_nvar = 3 ! pot, eng, flux

  integer,parameter :: GCMSIGMAP = 1,GCMSIGMAH = 2, GCMNORTH=NORTH, GCMSOUTH=SOUTH
  integer,parameter :: GCMhemispheres=2

  
  type var_T
      real(rp), dimension(:,:),allocatable :: var
  end type var_T

  type gcm_T
      !type(mixGrid_T) :: rGEOGrid,tGEOGrid,tSMGrid ! GEO&SM Grid for tgcm, remix GEO grid
      !type(mixGrid_T) :: mixGrid,gcmGrid
      type(Map_T), allocatable, dimension(:) :: r2tMaps, t2rMaps
      type(mixGrid_T), allocatable, dimension(:) :: SM,GEO
      type(var_T), dimension(GCMhemispheres,gcm2mix_nvar) :: gcmInput,mixInput
      type(var_T), dimension(GCMhemispheres,mix2gcm_nvar) :: gcmOutput,mixOutput
      !real(rp), dimension(:,:,:,:), allocatable :: mixInput, mixOutput
      !real(rp), dimension(:,:,:,:), allocatable :: gcmInput, gcmOutput
      integer :: nlat,nlon,nlev,nhlat,ntime,order,nhemi,lonshift
      integer, dimension(:), allocatable :: outlist
      integer, dimension(:), allocatable :: t2N, t2S
      real(rp), dimension(:), allocatable :: time,lev,lon,clat,lat
      real(rp),dimension(:,:),allocatable :: gx,gy
      real(rp),dimension(:,:,:),allocatable :: glon,gclat
      integer :: cplStep = 1
      !character(len=strlen) :: mix2gcmH5,gcm2mixH5,mix2gcmLock,gcm2mixLock
      character(len=strlen) :: mix2gcmH5 = "mix4gcm.h5"
      character(len=strlen) :: mix2gcmLock = "mix2gcm.txt"
      character(len=strlen) :: gcm2mixH5 = "gcm4mix.h5"
      character(len=strlen) :: gcm2mixLock = "gcm2mix.txt"
      logical :: isRestart
  end type gcm_T

  
end module gcmtypes
