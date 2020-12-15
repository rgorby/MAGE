module gcminterp
  use mixdefs
  use mixtypes
  use gcmtypes
  use mixgeom
  use math
  use ioh5
  use mixio
  !use geopack

  implicit none

  type(mixGrid_T), dimension(:),allocatable,private :: gcmG

  contains

    subroutine init_gcm(gcm,isRestart)
      type(gcm_T),intent(inout) :: gcm
      logical, intent(in) :: isRestart

      gcm%isRestart = isRestart
      !write(*,*) "init_gcm: ",gcm%isRestart
      if (gcm%cplStep == 1) then
        call CheckAndKill(gcm%mix2gcmLock)
        call CheckAndKill(gcm%mix2gcmH5)
      endif
    end subroutine init_gcm

    subroutine init_gcm_mix(gcmApp,ion)!,remixApp
      type(gcm_T), intent(inout) :: gcmApp
      type(mixIon_T),dimension(:),intent(in) :: ion
      !type(mixApp_T), intent(in) :: remixApp
      real(rp), dimension(:,:), allocatable :: gcmp,gcmt
      integer :: i,j,h,Np,Nt,Nh

    !Create gcm object
      !call read_from_netcdf(gcmApp)
      ! let's assume gcm grid in the form of a mix grid for now
      ! not sure if true
      !gcmApp%remixGrid = remixApp%
      !call mix_interpolant(gcmApp%remixGrid)
      !write(*,*) 'Hello from init_gcm_mix'
      call readH5gcm(gcmApp)

      Np = gcmApp%nlon
      Nt = gcmApp%nlat
      Nh = gcmApp%nhemi
      if (.not.allocated(gcmG)) allocate(gcmG(Nh))
    !Now do remix mapping

      if (.not.allocated(gcmp)) allocate(gcmp(Np,Nt))
      if (.not.allocated(gcmt)) allocate(gcmt(Np,Nt))

      if (.not.allocated(gcmApp%r2tMaps)) allocate(gcmApp%r2tMaps(Nh))
      if (.not.allocated(gcmApp%t2rMaps)) allocate(gcmApp%t2rMaps(Nh))

      ! construct the 2-D grid
      !do j=1,Np
      !   gcmt(j,:) = deg2rad*gcmApp%gcolat
      !enddo

      !do i=1,Nt
      !   gcmp(:,i) = deg2rad*(gcmApp%glon)+pi
      !   write(*,*) 'gcolat check: ',gcmt(9,i)
      !enddo

      !do j=1,Np
      !  do i=1,Nt
      !    gcmt(j,i) = !deg2rad*(gcmApp%gcolat(i)-gcmApp%gcolat(1)) !Shift pole
      !    gcmp(j,i) = modulo(deg2rad*gcmApp%glon(j)+2*pi,(2*pi))  !Positive only.
      !  end do
      !end do

      !do h=1,Nh
      !call init_grid_fromTP(gcmG(GCMNORTH),gcmApp%glat,gcmApp%glon,isSolverGrid=.true.)
      !call init_grid_fromTP(gcmG(GCMSOUTH),gcmApp%glat,gcmApp%glon,isSolverGrid=.true.)
      call init_grid_fromXY(gcmG(GCMNORTH),gcmApp%gx,gcmApp%gy,isSolverGrid=.true.)
      call init_grid_fromXY(gcmG(GCMSOUTH),gcmApp%gx,gcmApp%gy,isSolverGrid=.true.)
      !write(*,*) "Array Check: ",gcmG(GCMSOUTH)%Np,gcmG(GCMSOUTH)%Nt,shape(gcmApp%gx)
        ! call remix grid constructor
        !write(*,*) 'do interpolant ',h
        !call mix_interpolant(gcmG(h))
      do h=1,Nh
        do i=1,gcm2mix_nvar
          if (.not.allocated(gcmApp%mixInput(h,i)%var)) allocate(gcmApp%mixInput(h,i)%var(ion(h)%G%Np,ion(h)%G%Nt))
        end do
      enddo
      !write(*,*) 'Grid check 1: ',gcmG(h)%
    end subroutine init_gcm_mix

    subroutine coupleGCM2MIX(gcm,ion,do2way,mjd,time)
      type(mixIon_T),dimension(:),intent(inout) :: ion
      type(gcm_T), intent(inout) :: gcm
      logical,optional,intent(in) :: do2way
      real(rp), optional, intent(in) :: time, mjd
      
      !Must MIX export first.  TIEGCM will also import first.
      call writeMIX2GCM(ion,gcm%mix2gcmH5,gcm%mix2gcmLock,gcm%cplStep,mjd,time)

      !Import gcm data
      if (gcm%cplStep == 1) then
        call init_gcm_mix(gcm,ion)
      else
        call ReadH5gcm(gcm)
      end if
      !Map the data to MIX grid (redundant for now)
      call mapGCM2MIX(gcm,ion)

      if (gcm%isRestart) gcm%isRestart=.false.
      gcm%cplStep = gcm%cplStep + 1
        
    end subroutine coupleGCM2MIX

    subroutine ReadH5gcm(gcm)
      type(gcm_T), intent(inout) :: gcm
      type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
      character (len=strLen) :: H5file, lockStr
      integer :: ntime,nlon,nlat,nlev,nlathemi,nhemi
      integer :: i,j,t,v,ordering
      integer :: N,Nr,Ndim4
      real(rp), dimension(:), allocatable :: time,lev,etac
      real(rp), dimension(:,:),allocatable :: x,y
      real(rp), dimension(:,:,:,:), allocatable :: var2d

      logical :: doSP = .false. !Restarts are always double precision!Do single precision
      logical :: fExist = .false.

      !Prepare for reading
      !H5file = "LITd-p32_mar26_2003_80F_no_tsoft_sech_tie_2003-03-28T12-01-00_2003-03-28T12-30-00.nc4"
      H5file = gcm%gcm2mixH5
      lockStr = gcm%gcm2mixLock    
      write(*,*) "mix: waiting for ",trim(lockStr)," to be appear"
      do while (.not.fExist)
        inquire(file=lockStr,exist=fExist)
        call sleep(1)
      end do
      write(*,*) "mix: ",trim(lockStr)," is here!"

      write(*,*) 'mix: Reading in gcm file: ', trim(H5file)," step: ",gcm%cplStep
      call CheckFileOrDie(H5file,"Couldn't find input h5 file. Exiting...")
      call ClearIO(IOVars) !Reset IO chain

      !write(*,*) 'Start scalar read'
      !List variables to read
      
      !Scalars (need to specify integers), order doesn't matter
      
      !Arrays

      !1D Arrays

      !2D Arrays
      call AddInVar(IOVars,"Grid X")
      call AddInVar(IOVars,"Grid Y")
      call AddInVar(IOVars,"Pedersen Conductance NORTH")
      call AddInVar(IOVars,"Pedersen Conductance SOUTH")
      call AddInVar(IOVars,"Hall Conductance NORTH")
      call AddInVar(IOVars,"Hall Conductance SOUTH")
      
      !3D Arrays
      
      !4D Arrays
      
      !write(*,*) 'ReadVars'
      !Now do actual reading
      call ReadVars(IOVars,doSP,H5File)
      
      !Find grid dimension
      N = FindIO(IOVars,"Grid X")
      nlon = IOVars(N)%dims(2)
      nlat = IOVars(N)%dims(1)
      !nlev = IOVars(N)%dims(3)
      !ntime = IOVars(N)%dims(4)
      nhemi = GCMhemispheres
      
      !Save dimension information
      gcm%nlon  = nlon
      gcm%nlat  = nlat
      gcm%nhemi = nhemi

      !allocate(etac(lon))
      !write(*,*) 'Allocate arrays: ',nlon,nlat
      !The arrays are read in in the opposite order as listed by netcdf.
      !The time array is the last dimension if (.not. allocated(var2d)) 
      allocate(var2d(nlat,nlon,nhemi,gcm2mix_nvar))
      if (.not. allocated(x)) allocate(x(nlat,nlon))
      if (.not. allocated(y)) allocate(y(nlat,nlon))
      if (.not. allocated(gcm%gx)) allocate(gcm%gx(nlon,nlat))
      if (.not. allocated(gcm%gy)) allocate(gcm%gy(nlon,nlat))
      if (.not. allocated(gcm%gcmInput(GCMNORTH,GCMSIGMAP)%var)) allocate(gcm%gcmInput(GCMNORTH,GCMSIGMAP)%var(nlon,nlat))
      if (.not. allocated(gcm%gcmInput(GCMSOUTH,GCMSIGMAP)%var)) allocate(gcm%gcmInput(GCMSOUTH,GCMSIGMAP)%var(nlon,nlat))
      if (.not. allocated(gcm%gcmInput(GCMNORTH,GCMSIGMAH)%var)) allocate(gcm%gcmInput(GCMNORTH,GCMSIGMAH)%var(nlon,nlat))
      if (.not. allocated(gcm%gcmInput(GCMSOUTH,GCMSIGMAH)%var)) allocate(gcm%gcmInput(GCMSOUTH,GCMSIGMAH)%var(nlon,nlat))

      !Pull Scalars

      !write(*,*) 'Pull 1D Arrays'
      !Pull 1D Arrays

      !Pull 2D Arrays
      call IOArray2DFill(IOVars,"Grid X",x)
      call IOArray2DFill(IOVars,"Grid Y",y)
      call IOArray2DFill(IOVars,"Pedersen Conductance NORTH",var2d(:,:,GCMNORTH,GCMSIGMAP))
      call IOArray2DFill(IOVars,"Pedersen Conductance SOUTH",var2d(:,:,GCMSOUTH,GCMSIGMAP))
      call IOArray2DFill(IOVars,"Hall Conductance NORTH",var2d(:,:,GCMNORTH,GCMSIGMAH))
      call IOArray2DFill(IOVars,"Hall Conductance SOUTH",var2d(:,:,GCMSOUTH,GCMSIGMAH))
      
      !Pull 3D Arrays
      
      !Pull 4D Arrays
      
      do i = 1,nlon
        do j = 1,nlat
          do v = 1,gcm2mix_nvar
            gcm%gx(i,j) = x(j,i)
            gcm%gy(i,j) = y(j,i)
            gcm%gcmInput(GCMNORTH,v)%var(i,j) = var2d(j,i,GCMNORTH,v)
            gcm%gcmInput(GCMSOUTH,v)%var(i,j) = var2d(j,i,GCMSOUTH,v)
          end do
        end do
      end do

      if (allocated(var2d)) deallocate(var2d)
      call CheckAndKill(lockStr)

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
      !write(*,*) 'SigmaP maxval: ',maxval(gcm%gcmInput(GCMNORTH,GCMSIGMAP)%var(:,:)),maxval(gcmApp%mixInput(GCMNORTH,GCMSIGMAP)%var(:,:))
      !write(*,*) 'SigmaH maxval: ',maxval(gcm%gcmInput(GCMNORTH,GCMSIGMAH)%var(:,:)),maxval(gcmApp%mixInput(GCMNORTH,GCMSIGMAH)%var(:,:))

    end subroutine ReadH5gcm

    subroutine mapGCM2MIX(gcmApp,ion)
      type(mixIon_T),dimension(:),intent(inout) :: ion
      type(gcm_T), intent(inout) :: gcmApp
      !type(mixState_T), intent(inout) :: St
      type(Map_T) :: Map
      
      real(rp), dimension(:,:), allocatable :: F
      integer :: h,v

      ! GEOPACK transforms require setting UT time
      !call RECALC(utime(1),utime(2),utime(3),utime(4),utime(5))

      !write(*,*) "mapGCM2MIX CPL STEP: ", gcmApp%cplStep

      do h=1,2
        if (gcmApp%cplStep == 1) then
          ! Plan A, write a modified mix_set_map as gcm_set_map
          !write(*,*) 'Mapping r2t'
          call gcm_set_map(ion(h)%G,gcmG(h),gcmApp%r2tMaps(h),iGEOtoSM)
          !write(*,*) 'Mapping t2r'
          call gcm_set_map(gcmG(h),ion(h)%G,gcmApp%t2rMaps(h),iSMtoGEO)
        endif
        do v=1,gcm2mix_nvar
          call mix_map_grids(gcmApp%t2rMaps(h),gcmApp%gcmInput(h,v)%var(:,:),F)
          select case (v)
          case (GCMSIGMAP)
            gcmApp%mixInput(h,v)%var(:,:) = F
          case (GCMSIGMAH)
            gcmApp%mixInput(h,v)%var(:,:) = F
          end select
        end do
      end do
      !write(*,*) 'SigmaP maxval: ',maxval(gcmApp%gcmInput(GCMNORTH,GCMSIGMAP)%var(:,:)),maxval(gcmApp%mixInput(GCMNORTH,GCMSIGMAP)%var(:,:))
      !write(*,*) 'SigmaH maxval: ',maxval(gcmApp%gcmInput(GCMNORTH,GCMSIGMAH)%var(:,:)),maxval(gcmApp%mixInput(GCMNORTH,GCMSIGMAH)%var(:,:))
    end subroutine mapGCM2MIX
    
    subroutine apply_gcm2mix(gcm,St,h)
      type(gcm_T),intent(in) :: gcm
      type(mixState_T),intent(inout) :: St
      integer :: v,h
      
      do v=1,gcm2mix_nvar
        select case (v)
        case (GCMSIGMAP)
          !write(*,*) 'SigmaP size: ',shape(St%Vars(:,:,SIGMAP)),shape(gcm%mixInput(h,v)%var(:,:))
          St%Vars(:,:,SIGMAP) = gcm%mixInput(h,v)%var(:,:)
        case (GCMSIGMAH)
          !write(*,*) 'SigmaH size: ',shape(St%Vars(:,:,SIGMAH)),shape(gcm%mixInput(h,v)%var(:,:))
          St%Vars(:,:,SIGMAH) = gcm%mixInput(h,v)%var(:,:)
          !write(*,*) 'SigmaH maxval: ',maxval(gcm%mixInput(h,v)%var(:,:))
        end select
      end do
      
      !write(*,*) 'Done applying GCM!'
    end subroutine apply_gcm2mix

    subroutine gcm_set_map(G1,G2,Map,transform)
      ! Alternative to creating and holding onto an entire mapping for intermediate grid
      ! This is a slightly modified version of the mix_set_map that includes a grid transformation
      ! Prior to calling gcm_set_map, must define recalc
      type(mixGrid_T), intent(in) :: G1,G2
      type(Map_T), intent(inout) :: Map
      integer, intent(in) :: transform
  
      integer :: i1,j1,i2,j2,j,i
      real(rp), dimension(4) :: x
      real(rp), dimension(3) :: gsm,gnew
      real(rp) :: pnew,tnew,zin
      
      !write(*,*) 'Np: ',G2%Np,' Nt: ',G2%Nt
      if (.not. allocated(Map%M)) allocate(Map%M(G2%Np,G2%Nt,4))
      if (.not. allocated(Map%I1)) allocate(Map%I1(G2%Np,G2%Nt))
      if (.not. allocated(Map%J1)) allocate(Map%J1(G2%Np,G2%Nt))

      ! loop over points onto which we're interpolating
      do i2=1,G2%Nt
         do j2=1,G2%Np
            zin = sin(G2%t(j2,i2))
            
            select case (transform)
            case (iGEOtoSM)
                !call GEOGSW_08(G2%x(j2,i2),G2%y(j2,i2),zin,gsm(1),gsm(2),gsm(3),iGEOtoGSM)
                !call SMGSW_08(gnew(1),gnew(2),gnew(3),xggsm(1)sm,gsm(2),gsm(3),iGSMtoSM)
            case (iSMtoGEO)
                !call SMGSW_08(G2%x(j2,i2),G2%y(j2,i2),zin,gsm(1),gsm(2),gsm(3),iSMtoGSM)
                !call GEOGSW_08(gnew(1),gnew(2),gnew(3),gsm(1),gsm(2),gsm(3),iGSMtoGEO)
            end select
            !!!!! TEMPORARY UNTIL GEOPACK WORKS !!!!!
            gnew(1) = G2%x(j2,i2)
            gnew(2) = G2%y(j2,i2)
            gnew(3) = zin
            !write(*,*) 'gnew: ',gnew
            !!!!! END OF TEMP CODE !!!!!
            
            tnew = asin(sqrt(gnew(1)**2+gnew(2)**2))
            pnew = modulo((atan2(gnew(2),gnew(1))+2*pi),(2*pi)) ! any pole problems??
            
            call mix_search(G1,pnew,tnew,j1,i1)

            if ((i1 == 0) ) then
              !write(*,*) '??? ',pnew,tnew,j1,i1,gnew,j2,i2,G2%t(j2,i2),G2%p(j2,i2)
              i1=1
            end if

            Map%I1(j2,i2) = i1
            Map%J1(j2,i2) = j1
            x = (/1._rp, pnew,tnew, pnew*tnew /)

            ! treat the poles
            if ((i1.eq.1).and.(G1%t(1,1)).lt.TINY) then
               ! just average over the three vertices of the triangle
               ! could make it more fancy to interpolate using a0+a1*t+a2*t*p interpolant, but why bother???
               Map%M(j2,i2,:) = 1._rp/3._rp*(/0.5_rp, 0.5_rp, 1._rp, 1._rp/)
            else if ((i1.eq.G1%Nt-1).and.(pi-G1%t(1,G1%Nt-1)).lt.TINY) then
               Map%M(j2,i2,:) = 1._rp/3._rp*(/1._rp, 1._rp, 0.5_rp, 0.5_rp/)               
            else 
               Map%M(j2,i2,:) = matmul(G1%Interpolant(j1,i1,:,:),x)
            end if
         end do
      end do
    end subroutine gcm_set_map
end module gcminterp
