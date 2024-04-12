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

  contains

    subroutine init_gcm_grid(gcm,ion)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(in) :: ion
      
      integer :: Np, Nt, Nh, i,j,h,v
      Np = gcm%nlon
      Nt = gcm%nhlat
      Nh = gcm%nhemi

      ! Let's try something
      ! Build a list of index that says these are N and these are S
      ! Then use that to map to N/S in gcm
      gcm%t2N = pack([(i,i=1,gcm%nlat)], gcm%lat > 10)
      gcm%t2S = pack([(i,i=1,gcm%nlat)], gcm%lat < -10)
      gcm%lonshift = findloc(gcm%lon(:gcm%nlon-1),0.,1)-1

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
        gcm%glon(:gcm%nlon-1,j,GCMNORTH) = modulo(cshift(gcm%lon(:gcm%nlon-1),gcm%lonshift)+360.,360.)!modulo((atan2(G2%y,G2%x)+2*pi),(2*pi))
        gcm%glon(:gcm%nlon-1,j,GCMSOUTH) = modulo(cshift(gcm%lon(:gcm%nlon-1),gcm%lonshift)+360.,360.)
        ! cover periodic longitude point
        gcm%glon(gcm%nlon,j,GCMNORTH) = 360. ! no matter what, last point is periodic point
        gcm%glon(gcm%nlon,j,GCMSOUTH) = 360. ! hard coding the value of last point for now
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
    end subroutine init_gcm_grid

    subroutine process_gcmimport(gcm,ion)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      
      integer :: v, i, h, ymod1

      call Tic("Processing")
      call Tic("Shifting")
      ! Split global GCM array into two hemispheres
      ! then shift the array so that longitude starts at 0
      do v=1,gcm2mix_nvar
        do i=1,gcm%nhlat
          gcm%gcmInput(GCMNORTH,v)%var(:gcm%nlon-1,i) = cshift(gcm%invar2d(gcm%t2N(gcm%nhlat-i+1),:gcm%nlon-1,v),gcm%lonshift)
          gcm%gcmInput(GCMSOUTH,v)%var(:gcm%nlon-1,i) = cshift(gcm%invar2d(gcm%t2S(i),:gcm%nlon-1,v),gcm%lonshift)
          gcm%gcmInput(GCMNORTH,v)%var(gcm%nlon,i) = gcm%gcmInput(GCMNORTH,v)%var(1,i)
          gcm%gcmInput(GCMSOUTH,v)%var(gcm%nlon,i) = gcm%gcmInput(GCMSOUTH,v)%var(1,i)
        enddo
      end do
      call Toc("Shifting")

      call Tic("Mapping")
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
      call Toc("Mapping")
      call Toc("Processing")

    end subroutine process_gcmimport

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

end module gcminterp

