module gcminterp
  use mixdefs
  use mixtypes
  use gcmtypes
  use mixgeom
  use math
  use ioh5
  use mixio
  use kai2geo
  use clocks

  implicit none

  type(mixGrid_T),allocatable,private :: gcmG
  character(len=strLen), dimension(nVars), private :: mixVarNames
  character(len=strLen), dimension(nVars), private :: mixUnitNames  

  contains

    subroutine init_gcm(gcm,ion,isRestart)
      type(gcm_T),intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      logical, intent(in) :: isRestart

      write(*,*) "start init_gcm"

      if (.not. allocated(gcm%outlist)) allocate(gcm%outlist(mix2gcm_nvar))
      gcm%outlist(1) = POT
      gcm%outlist(2) = AVG_ENG
      gcm%outlist(3) = NUM_FLUX
      call initGCMNames()

      gcm%isRestart = isRestart

      call init_gcm_mix(gcm,ion)
      
      
    end subroutine init_gcm

  subroutine initGCMNames()
    ! NOTE: these have to be in the same order as the
    ! variable enumerator in mixdefs
    ! LAZY COPY PASTE HERE
    mixVarNames(POT)           = "Potential"
    mixUnitNames(POT)          = "kV"
    mixVarNames(FAC)           = "Field-aligned current"
    mixUnitNames(FAC)          = "muA/m**2"
    mixVarNames(SIGMAP)        = "Pedersen conductance"
    mixUnitNames(SIGMAP)       = "S"
    mixVarNames(SIGMAH)        = "Hall conductance"
    mixUnitNames(SIGMAH)       = "S"
    mixVarNames(SOUND_SPEED)   = "Sound speed"
    mixUnitNames(SOUND_SPEED)  = "cm/s"
    mixVarNames(DENSITY)       = "Density"    
    mixUnitNames(DENSITY)      = "g/cm^3"
    mixVarNames(AVG_ENG)       = "Average energy"
    mixUnitNames(AVG_ENG)      = "keV"
    mixVarNames(NUM_FLUX)      = "Number flux"
    mixUnitNames(NUM_FLUX)     = "1/cm^2 s"
    mixVarNames(NEUTRAL_WIND)  = "Neutral wind"
    mixUnitNames(NEUTRAL_WIND) = "cm/s"
    mixVarNames(EFIELD)        = "Electric field"
    mixUnitNames(EFIELD)       = "mV/m" 
    mixVarNames(IM_EFLUX)      = "IM Energy flux"
    mixUnitNames(IM_EFLUX)     = "ergs/cm^2 s"
    mixVarNames(IM_EAVG)       = "IM average energy"
    mixUnitNames(IM_EAVG)      = "keV" ! add *1e-3 in rcm_mix_interface.F90
    mixVarNames(IM_IFLUX)      = "IM Energy flux proton"
    mixUnitNames(IM_IFLUX)     = "ergs/cm^2 s"
    mixVarNames(IM_IAVG)       = "IM average energy proton"
    mixUnitNames(IM_IAVG)      = "keV" ! add *1e-3 in rcm_mix_interface.F90
    mixVarNames(Z_NFLUX)       = "Zhang number flux"
    mixUnitNames(Z_NFLUX)      = "1/cm^2 s"
    mixVarNames(Z_EAVG)        = "Zhang average energy"
    mixUnitNames(Z_EAVG)       = "keV"
    mixVarNames(CRPOT)         = "Corotation Potential"
    mixUnitNames(CRPOT)        = "kV"
    mixVarNames(TPOT)          = "Total Potential"
    mixUnitNames(TPOT)         = "kV"
    mixVarNames(IM_GTYPE)      = "RCM grid type"
    mixUnitNames(IM_GTYPE)     = "0-1"
    mixVarNames(AUR_TYPE)      = "Auroral model type"
    mixUnitNames(AUR_TYPE)     = "Zhang Fedder RCM RCMZ"
    mixVarNames(IM_BETA)       = "RCM beta"
    mixUnitNames(IM_BETA)      = "0-1"
    mixVarNames(IM_EDEN)       = "RCM electron density"
    mixUnitNames(IM_EDEN)      = "#/m^3"
    mixVarNames(IM_EPRE)       = "RCM electron pressure"
    mixUnitNames(IM_EPRE)      = "Pa"
    mixVarNames(IM_ENFLX)      = "IM Number flux"
    mixUnitNames(IM_ENFLX)     = "1/cm^2 s"
    mixVarNames(IM_INFLX)      = "IM Number flux proton"
    mixUnitNames(IM_INFLX)     = "1/cm^2 s"
    mixVarNames(DELTAE)        = "Mono potential drop"
    mixUnitNames(DELTAE)       = "kV"
  end subroutine initGCMNames

    subroutine init_gcm_mix(gcm,ion)!,remixApp
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      !type(mixApp_T), intent(in) :: remixApp
      real(rp), dimension(:,:), allocatable :: gcmp,gcmt
      integer :: i,j,h,v,Np,Nt,Nh,N,nlon,nlat,nhlatN,nhlatS
      integer, allocatable :: positive_values(:)
      integer, dimension(:), allocatable :: indicesN,indicesS
      character (len=strLen) :: H5file, lockStr
      logical :: fExist = .false., doSP = .false.
      type(IOVAR_T), dimension(MAXMIXIOVAR) :: IOVars

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Prepare for reading
      H5file  = "init_"//trim(gcm%gcm2mixH5)
      lockStr = "init_"//trim(gcm%gcm2mixLock)
      write(*,*) "mix: waiting for ",trim(lockStr)," to be appear"
      do while (.not.fExist)
        inquire(file=lockStr,exist=fExist)
        call sleep(1)
      end do
      write(*,*) "mix: ",trim(lockStr)," is here!"

      write(*,*) 'mix: Reading in gcm file: ', trim(H5file)," step: init"
      call CheckFileOrDie(H5file,"Couldn't find input h5 file. Exiting...")
      call ClearIO(IOVars) !Reset IO chain

      !write(*,*) 'Start scalar read'
      !List variables to read
      
      !Scalars (need to specify integers), order doesn't matter
      
      !Arrays

      !1D Arrays
      call AddInVar(IOVars,"Lat")
      call AddInVar(IOVars,"Lon")

      !Now do actual reading
      call ReadVars(IOVars,doSP,H5File)

      !Find grid dimension
      N = FindIO(IOVars,"Lon")
      nlon = IOVars(N)%dims(1)
      N = FindIO(IOVars,"Lat")
      nlat = IOVars(N)%dims(1)
      Nh = GCMhemispheres
      
      if (.not.allocated(gcm%lat)) allocate(gcm%lat(nlat))
      if (.not.allocated(gcm%lon)) allocate(gcm%lon(nlon))

      call IOArray1DFill(IOVars,"Lat",gcm%lat)
      call IOArray1DFill(IOVars,"Lon",gcm%lon)

      !call CheckAndKill(lockStr)
      
      write(*,*) "mix: Done init GCM read"

      !Save dimension information
      gcm%nlon  = nlon
      gcm%nlat  = nlat
      gcm%nhemi = Nh
      
      ! Let's try something
      ! Build a list of index that says these are N and these are S
      ! Then use that to map to N/S in gcm
      !indicesN = merge( 0, [ ( i, i = 1, size( gcm%lat ) ) ], gcm%lat < 0 )
      !indicesS = merge( 0, [ ( i, i = 1, size( gcm%lat ) ) ], gcm%lat > 0 )
      !nhlatN = count(indicesN)
      !nhlatS = count(indicesS)
      !if ( .not. allocated(gcm%t2N) ) allocate(gcm%t2N(nhlatN))
      !if ( .not. allocated(gcm%t2S) ) allocate(gcm%t2S(nhlatS))
      !gcm%t2N = pack( indicesN, indicesN /= 0 )
      !gcm%t2S = pack( indicesS, indicesS /= 0 )
      gcm%t2N = pack([(i,i=1,nlat)], gcm%lat > 10)
      gcm%t2S = pack([(i,i=1,nlat)], gcm%lat < -10)
      gcm%lonshift = findloc(gcm%lon(:nlon-1),0.,1)-1

      ! Going to assume the hemispheres are symetrically sized for now
      if ( size(gcm%t2N) /= size(gcm%t2S) ) write(*,*) 'WE GOT ASYMMETRY'
      gcm%nhlat = size(gcm%t2N)

      if (.not.allocated(gcm%gclat)) allocate(gcm%gclat(gcm%nlon,gcm%nhlat,Nh))
      if (.not.allocated(gcm%glon)) allocate(gcm%glon(gcm%nlon,gcm%nhlat,Nh))

      ! Convert things from latitude to colatitude (and funky colat for south)
      do i=1,gcm%nlon
        gcm%gclat(i,:,GCMNORTH) = 90.-gcm%lat(gcm%t2N)
        gcm%gclat(i,:,GCMNORTH) = gcm%gclat(i,gcm%nhlat:1:-1,GCMNORTH) ! reverse order. Ascending.
        gcm%gclat(i,:,GCMSOUTH) = 90.+gcm%lat(gcm%t2S) !For mapping to work, we're going to use remix's funky southern colat
      enddo
        
      ! This complicated looking thing converts [-180,180] longitude into [0,360] and also shifts it so that longitude starts at 0
      do j=1,gcm%nhlat
        gcm%glon(:gcm%nlon-1,j,GCMNORTH) = modulo(cshift(gcm%lon(:nlon-1),gcm%lonshift)+360.,360.)!modulo((atan2(G2%y,G2%x)+2*pi),(2*pi))
        gcm%glon(:gcm%nlon-1,j,GCMSOUTH) = modulo(cshift(gcm%lon(:nlon-1),gcm%lonshift)+360.,360.)
        ! cover periodic longitude point
        gcm%glon(gcm%nlon,j,GCMNORTH) = 360. ! no matter what, last point is periodic point
        gcm%glon(gcm%nlon,j,GCMSOUTH) = 360. ! hard coding the value of last point for now
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Np = gcm%nlon
      Nt = gcm%nhlat
      Nh = gcm%nhemi

      !Now do remix mapping
      if (.not.allocated(gcm%GEO)) allocate(gcm%GEO(Nh))
      if (.not.allocated(gcm%SM)) allocate(gcm%SM(Nh))
      if (.not.allocated(gcm%r2tMaps)) allocate(gcm%r2tMaps(Nh))
      if (.not.allocated(gcm%t2rMaps)) allocate(gcm%t2rMaps(Nh))

      call init_grid_fromTP(gcm%GEO(GCMNORTH),gcm%gclat(:,:,GCMNORTH)*deg2rad,gcm%glon(:,:,GCMNORTH)*deg2rad,isSolverGrid=.true.)
      call init_grid_fromTP(gcm%GEO(GCMSOUTH),gcm%gclat(:,:,GCMSOUTH)*deg2rad,gcm%glon(:,:,GCMSOUTH)*deg2rad,isSolverGrid=.true.)

      do h=1,Nh
        do v=1,gcm2mix_nvar
          if (.not.allocated(gcm%mixInput(h,v)%var)) allocate(gcm%mixInput(h,v)%var(ion(h)%G%Np,ion(h)%G%Nt))
          if (.not.allocated(gcm%gcmInput(h,v)%var)) allocate(gcm%gcmInput(h,v)%var(gcm%nlon,gcm%nhlat))
        end do
        do v=1,mix2gcm_nvar
          if (.not. allocated(gcm%gcmOutput(h,v)%var)) allocate(gcm%gcmOutput(h,v)%var(gcm%nlon,gcm%nhlat))
        end do
      enddo
    end subroutine init_gcm_mix

    subroutine coupleGCM2MIX(gcm,ion,do2way,mjd,time)
      type(mixIon_T),dimension(:),intent(inout) :: ion
      type(gcm_T), intent(inout) :: gcm
      logical,optional,intent(in) :: do2way
      real(rp), optional, intent(in) :: time, mjd

      ! maybe create the SM and GEO list of points here
      ! since the destination grid does not need to be structured
      ! can do a simple loop over all grid points and transform
      ! transform all remix points to GEO
      ! transform all gcm points to SM


      !Must MIX export first.  TIEGCM will also import first.
      call Tic("Export")
      call exportgcm(ion,gcm,mjd,time)
      call Toc("Export")

      !Import gcm data
      call Tic("Import")
      call importgcm(gcm, ion)
      call Toc("Import")

      if (gcm%isRestart) gcm%isRestart=.false.
      ! We have a separate internal counter for coupling here. 
      ! This may be used later on for WACCM-X coupling which is desync from remix coupling time
      ! TIEGCM coupling time is also 5s while WACCM-X will couple at 1 min default
      gcm%cplStep = gcm%cplStep + 1
        
    end subroutine coupleGCM2MIX

    subroutine importgcm(gcm,ion)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      integer :: h,ymod1

      call Tic("Read")
      call ReadH5gcm(gcm)
      call Toc("Read")      

      ! The weird ymod here is to undo the funky southern hemisphere colat issue.
      do h=1,size(ion)
        if (h == GCMSOUTH) then
          ymod1 = -1
        else
          ymod1 = 1
        endif
        !ymod1 = 1
        !write(*,*) "SM -> GEO START:",h,ymod1
        call transform_grid(ion(h)%G,ion(h)%Ggeo,iSMtoGEO,h,ym1=ymod1)
        !write(*,*) "SM -> GEO END: ",h
      end do
      
      !Map the data to MIX grid
      call mapGCM2MIX(gcm,ion)

    end subroutine

    subroutine exportgcm(ion,gcm,mjd,time)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T), dimension(:),intent(inout) :: ion
      real(rp), optional, intent(in) :: time, mjd
      integer :: h,ymod2

      ! The weird ymod here is to undo the funky southern hemisphere colat issue.
      do h=1,gcm%nhemi
        if (h == GCMSOUTH) then
          ymod2 = -1
        else
          ymod2 = 1
        endif
        !ymod2 = 1
        !write(*,*) "GEO -> SM START:",h,ymod2
        call transform_grid(gcm%GEO(h),gcm%SM(h),iGEOtoSM,h,ym2=ymod2)
        !write(*,*) "GEO -> SM END: ",h
      end do

      ! Map from mix grid to gcm grid
      call mapMIX2GCM(ion,gcm)

      ! write the coupling file
      call writeMIX2GCM(ion,gcm,mjd,time)

    end subroutine

    subroutine ReadH5gcm(gcm)
      type(gcm_T), intent(inout) :: gcm
      type(IOVAR_T), dimension(MAXMIXIOVAR) :: IOVars
      character (len=strLen) :: H5file, lockStr
      integer :: ntime,nlon,nlat,nlev,nlathemi,nhemi
      integer :: i,j,t,v,ordering
      integer :: N,Nr,Ndim4
      real(rp), dimension(:), allocatable :: time,lev,etac
      real(rp), dimension(:,:),allocatable :: x,y
      real(rp), dimension(:,:,:), allocatable :: var2d

      logical :: doSP = .false. !Restarts are always double precision !doSP = Do single precision
      logical :: fExist = .false.

      !Prepare for reading
      call Tic("Wait4File")
      H5file = gcm%gcm2mixH5
      lockStr = gcm%gcm2mixLock    
      write(*,*) "mix: waiting for ",trim(lockStr)," to be appear"
      do while (.not.fExist)
        inquire(file=lockStr,exist=fExist)
        call sleep(1)
      end do
      write(*,*) "mix: ",trim(lockStr)," is here!"
      call Toc("Wait4File")

      call Tic("Load_Arrays")
      write(*,*) 'mix: Reading in gcm file: ', trim(H5file)," step: ",gcm%cplStep
      call CheckFileOrDie(H5file,"Couldn't find input h5 file. Exiting...")
      call ClearIO(IOVars) !Reset IO chain


      !write(*,*) 'Start scalar read'
      !List variables to read
      
      !Scalars (need to specify integers), order doesn't matter
      
      !Arrays

      !1D Arrays
      call AddInVar(IOVars,"Lat")
      call AddInVar(IOVars,"Lon")

      !2D Arrays
      call AddInVar(IOVars,"Pedersen Conductance")
      call AddInVar(IOVars,"Hall Conductance")
      
      !3D Arrays
      
      !4D Arrays
      
      !Now do actual reading
      call ReadVars(IOVars,doSP,H5File)
      
      !The arrays are read in in the opposite order as listed by netcdf.
      !The time array is the last dimension
      if (.not. allocated(var2d)) allocate(var2d(gcm%nlat,gcm%nlon,gcm2mix_nvar))

      !Pull Scalars

      !Pull 1D Arrays

      !Pull 2D Arrays
      call IOArray2DFill(IOVars,"Pedersen Conductance",var2d(:,:,GCMSIGMAP))
      call IOArray2DFill(IOVars,"Hall Conductance",var2d(:,:,GCMSIGMAH))
      
      !Pull 3D Arrays
      
      !Pull 4D Arrays

      ! Split global GCM array into two hemispheres
      ! then shift the array so that longitude starts at 0
      do v=1,gcm2mix_nvar
        do i=1,gcm%nhlat
          gcm%gcmInput(GCMNORTH,v)%var(:gcm%nlon-1,i) = cshift(var2d(gcm%t2N(gcm%nhlat-i+1),:gcm%nlon-1,v),gcm%lonshift)
          gcm%gcmInput(GCMSOUTH,v)%var(:gcm%nlon-1,i) = cshift(var2d(gcm%t2S(i),:gcm%nlon-1,v),gcm%lonshift)
          gcm%gcmInput(GCMNORTH,v)%var(gcm%nlon,i) = gcm%gcmInput(GCMNORTH,v)%var(1,i)
          gcm%gcmInput(GCMSOUTH,v)%var(gcm%nlon,i) = gcm%gcmInput(GCMSOUTH,v)%var(1,i)
        enddo
        !write(*,*) "SHAPES: ",shape(var2d),shape(gcm%gcmInput(GCMNORTH,v)%var)
        !write(*,*) "var2d: ",maxval(var2d(:,:,v)),minval(var2d(:,:,v)),v
        !write(*,*) "GCMINPUT NORTH: ",maxval(gcm%gcmInput(GCMNORTH,v)%var),minval(gcm%gcmInput(GCMNORTH,v)%var),v
        !write(*,*) "GCMINPUT SOUTH: ",maxval(gcm%gcmInput(GCMSOUTH,v)%var),minval(gcm%gcmInput(GCMSOUTH,v)%var),v
      end do

      if (allocated(var2d)) deallocate(var2d)
      call CheckAndKill(lockStr)
      call Toc("Load_Arrays")

      !Noisy Diagnostics below
      !write(*,*) 'Time array1: ',gcm%time
      !write(*,*) 'Lon array1: ',gcm%glon(:,1)!modulo(deg2rad*gcm%glon+2*pi,(2*pi))
      !write(*,*) 'Lat array1: ',gcm%glat(1,:)
      !write(*,*) 'Lat array2: ',gcm%glat(1,::-1)
      !write(*,*) 'maxLat: ',maxval(gcm%glat)
      !write(*,*) 'minLat: ',minval(gcm%glat)
      !write(*,*) 'Colat array1: ',gcm%gcolat
      !write(*,*) 'input shape1: ',shape(gcm%gcmInput(GCMNORTH,GCMSIGMAP)%var)
      !write(*,*) 'input shape2: ',shape(var2d)
      !write(*,*) 'SigmaP maxval: ',maxval(gcm%gcmInput(GCMNORTH,GCMSIGMAP)%var(:,:)),maxval(gcm%mixInput(GCMNORTH,GCMSIGMAP)%var(:,:))
      !write(*,*) 'SigmaH maxval: ',maxval(gcm%gcmInput(GCMNORTH,GCMSIGMAH)%var(:,:)),maxval(gcm%mixInput(GCMNORTH,GCMSIGMAH)%var(:,:))

    end subroutine ReadH5gcm

    subroutine mapGCM2MIX(gcm,ion)
      type(mixIon_T),dimension(:),intent(inout) :: ion
      type(gcm_T), intent(inout) :: gcm
      type(Map_T) :: Map
      
      real(rp), dimension(:,:), allocatable :: F
      integer :: h,v

      do h=1,gcm%nhemi
        call mix_set_map(gcm%GEO(h),ion(h)%Ggeo,gcm%t2rMaps(h))
        do v=1,gcm2mix_nvar
          call mix_map_grids(gcm%t2rMaps(h),gcm%gcmInput(h,v)%var(:,:),F)
          select case (v)
          case (GCMSIGMAP)
            gcm%mixInput(h,v)%var(:,:) = F
          case (GCMSIGMAH)
            gcm%mixInput(h,v)%var(:,:) = F
          end select
        end do
      end do
    end subroutine mapGCM2MIX
    
    subroutine mapMIX2GCM(ion,gcm)
      type(mixIon_T),dimension(:),intent(inout) :: ion
      type(gcm_T), intent(inout) :: gcm
      type(Map_T) :: Map
      
      real(rp), dimension(:,:), allocatable :: F
      integer :: h,v

      do h=1,gcm%nhemi
        call mix_set_map(ion(h)%G,gcm%SM(h),gcm%r2tMaps(h))
        do v=1,mix2gcm_nvar
          call mix_map_grids(gcm%r2tMaps(h),ion(h)%St%Vars(:,:,gcm%outlist(v)),F)
          gcm%gcmOutput(h,v)%var(:,:) = F
        end do
      end do
    end subroutine mapMIX2GCM
    
    subroutine apply_gcm2mix(gcm,St,h)
      type(gcm_T),intent(in) :: gcm
      type(mixState_T),intent(inout) :: St
      integer :: v,h
      
      do v=1,gcm2mix_nvar
        select case (v)
        case (GCMSIGMAP)
          St%Vars(:,:,SIGMAP) = gcm%mixInput(h,v)%var(:,:)
        case (GCMSIGMAH)
          St%Vars(:,:,SIGMAH) = gcm%mixInput(h,v)%var(:,:)
        end select
      end do
      
    end subroutine apply_gcm2mix

    subroutine writeMIX2GCM(I,gcm,mjd,time)
      use mixio
      type(mixIon_T),dimension(:),intent(in) :: I
      type(gcm_T),intent(in) :: gcm
      real(rp), optional, intent(in) :: time, mjd
      character(len=strLen) :: vStr

      type(IOVAR_T), dimension(MAXMIXIOVAR) :: IOVars

      integer :: v,h,n0,cplStep
      character(len=strLen) :: gStr,uStr,hStr
      logical :: doDump = .true.,fExist=.true.
      real(rp) :: cpcp = 0.0
      real(rp), dimension(:,:),allocatable :: xc,yc
      real(rp), dimension(:,:,:),allocatable :: var2d
      
      character(len=strLen) :: cplStr,lockStr


      if (.not. allocated(var2d)) allocate(var2d(gcm%nlat,gcm%nlon,mix2gcm_nvar))
      var2d = 0.
      ! Before we start, we collapse to 1 globe instead of 2 hemispheres
      do v=1,mix2gcm_nvar
        var2d(gcm%t2N(gcm%nhlat:1:-1),:gcm%nlon-1,v) = transpose(cshift(gcm%gcmOutput(GCMNORTH,v)%var(:gcm%nlon-1,:),-1*gcm%lonshift))
        var2d(gcm%t2S,:gcm%nlon-1,v) = transpose(cshift(gcm%gcmOutput(GCMSOUTH,v)%var(:gcm%nlon-1,:),-1*gcm%lonshift))
        var2d(gcm%t2N,gcm%nlon,v) = var2d(gcm%t2N,1,v)
        var2d(gcm%t2S,gcm%nlon,v) = var2d(gcm%t2S,1,v)
      end do

      cplStr = gcm%mix2gcmH5
      lockStr = gcm%mix2gcmLock
      cplStep = gcm%cplStep

      inquire(file=lockStr,exist=fExist)

      !write(*,*) "waiting for ",trim(lockStr)," to be disappear"
      do while (fExist)
         inquire(file=lockStr,exist=fExist)
         call sleep(1)
      end do
      write(*,*) trim(lockStr)," is gone so creating ",trim(cplStr),' start coupling @ ',mjd
      call CheckAndKill(cplStr)
      
      !Reset IO chain
      call ClearIO(IOVars)
              
      do v=1,mix2gcm_nvar
        ! NOTE: assuming initMIXNames got called before
        vStr = trim(mixVarNames(gcm%outlist(v)))
        uStr = trim(mixUnitNames(gcm%outlist(v)))
        
        call AddOutVar(IOVars,vStr,var2d(:,:,v),uStr)
      enddo

      ! now add time
      if (present(time)) call AddOutVar(IOVars,"time",time)
      if (present(mjd))  call AddOutVar(IOVars,"MJD",mjd)
      ! also add tilt
      call AddOutVar(IOVars,"tilt",I(NORTH)%St%tilt)

      ! add grid info
      call AddOutVar(IOVars,"lat",gcm%lat)
      call AddOutVar(IOVars,"lon",gcm%lon)

      ! add cpcp
      call AddOutVar(IOVars,"nCPCP",maxval(I(NORTH)%St%Vars(:,:,POT))-minval(I(NORTH)%St%Vars(:,:,POT)))
      call AddOutVar(IOVars,"sCPCP",maxval(I(SOUTH)%St%Vars(:,:,POT))-minval(I(SOUTH)%St%Vars(:,:,POT)))    
      
      !Write out the chain (to root)
      write(gStr,'(A,I0)') "Step#", cplStep
      call WriteVars(IOVars,.false.,cplStr,gStr)
     
      write(*,*) "Done making ",trim(cplStr)," so locking"
      open(303,file=trim(lockStr))
        write(303,*) mjd
      close(303)
      write(*,*) trim(lockStr)," done, so go ahead"

      if (allocated(var2d)) deallocate(var2d)

    end subroutine writeMIX2GCM
end module gcminterp

