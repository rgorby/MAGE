module gcm_mpi
  use gcminterp
  use mpi_f08
  use mpidefs

  implicit none

  interface exportgcm
    module procedure exportgcmmpi
  end interface
  interface importgcm
    module procedure importgcmmpi
  end interface

contains

    subroutine init_gcm_mpi(gcm,ion,isRestart)
      type(gcm_T),intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      logical, intent(in) :: isRestart

      write(*,*) "start init_gcm"

      if (.not. allocated(gcm%outlist)) allocate(gcm%outlist(mix2gcm_nvar))
      gcm%outlist(1) = POT
      gcm%outlist(2) = AVG_ENG
      gcm%outlist(3) = NUM_FLUX

      gcm%isRestart = isRestart

      call init_gcm_grid(gcm,ion)
      
    end subroutine init_gcm_mpi

    subroutine init_gcm_mix_mpi(gcm,gcmCplComm,gcmCplRank)
      type(gcm_T), intent(inout) :: gcm
      type(MPI_Comm) :: gcmCplComm
      integer, intent(in) :: gcmCplRank

      integer :: nlon,nlat,Nh,ierr

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Find grid dimension
 
      write(*,*) "MIXCPL Waiting for GCM Grid"

      call mpi_recv(nlat, 1, MPI_INTEGER, gcmCplRank, (tgcmId+voltId)*100+1, gcmCplComm, MPI_STATUS_IGNORE,ierr)
      write(*,*) "MIXCPL GOT NLAT: ",nlat
      
      call mpi_recv(nlon, 1, MPI_INTEGER, gcmCplRank, (tgcmId+voltId)*100+2, gcmCplComm, MPI_STATUS_IGNORE,ierr)
      write(*,*) "MIXCPL GOT NLON: ",nlon

      Nh = GCMhemispheres
      
      if (.not.allocated(gcm%lat)) allocate(gcm%lat(nlat))
      if (.not.allocated(gcm%lon)) allocate(gcm%lon(nlon))


      call mpi_recv(gcm%lat, nlat, MPI_DOUBLE_PRECISION, gcmCplRank, (tgcmId+voltId)*100+3, gcmCplComm, MPI_STATUS_IGNORE,ierr)
      write(*,*) "MIXCPL GOT GLAT: ",gcm%lat
      call mpi_recv(gcm%lon, nlon, MPI_DOUBLE_PRECISION, gcmCplRank, (tgcmId+voltId)*100+4, gcmCplComm, MPI_STATUS_IGNORE,ierr)
      write(*,*) "MIXCPL GOT GLON: ",gcm%lon
      
      write(*,*) "MIXCPL GOT GRID INFO: ",nlat, nlon

      !Save dimension information
      gcm%nlon  = nlon
      gcm%nlat  = nlat
      gcm%nhemi = Nh

    end subroutine init_gcm_mix_mpi

    subroutine coupleGCM2MIX(gcm,ion,mjd,time,gcmCplComm,myRank)
      type(mixIon_T),dimension(:),intent(inout) :: ion
      type(gcm_T), intent(inout) :: gcm
      real(rp), intent(in) :: time, mjd
      type(MPI_Comm), optional :: gcmCplComm
      integer, optional, intent(in) :: myRank

      ! maybe create the SM and GEO list of points here
      ! since the destination grid does not need to be structured
      ! can do a simple loop over all grid points and transform
      ! transform all remix points to GEO
      ! transform all gcm points to SM

      !Skip MPI exchange on first restart
      if (present(gcmCplComm) .and. gcm%isRestart) then
        gcm%isRestart = .False.
        return
      endif
      !Must MIX export first.  TIEGCM will also import first.
      call Tic("Export")
      if (present(gcmCplComm)) then
        if (gcm%isRestart) then
          gcm%isRestart = .False.
          return
        endif
        call exportgcm(ion,gcm,mjd,time,gcmCplComm,myRank)
      else
        write(*,*) "Are we trying to Export to Couple GCM?"
      endif
      call Toc("Export")

      !Import gcm data
      call Tic("Import")
      if (present(gcmCplComm)) then
        call importgcm(gcm, ion, gcmCplComm,myRank)
      else
        write(*,*) "Are we trying to Import to Couple GCM?"
      endif
      call process_gcmimport(gcm,ion)
      call Toc("Import")

      if (gcm%isRestart) gcm%isRestart=.false.
      ! We have a separate internal counter for coupling here. 
      ! This may be used later on for WACCM-X coupling which is desync from remix coupling time
      ! TIEGCM coupling time is also 5s while WACCM-X will couple at 1 min default
      gcm%cplStep = gcm%cplStep + 1
        
    end subroutine coupleGCM2MIX


    subroutine importgcmmpi(gcm,ion, gcmCplComm,gcmCplRank)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T),dimension(:),intent(inout) :: ion
      type(MPI_Comm) :: gcmCplComm
      integer, intent(in) :: gcmCplRank

      integer :: v,ierr,length
      character( len = MPI_MAX_ERROR_STRING) :: message

      call Tic("MpiExchange")
      write(*,*) "MIX: IMPORT GCM"

      if (gcmCplComm /= MPI_COMM_NULL) then
        call Tic("Passing")
        if (.not. allocated(gcm%invar2d)) allocate(gcm%invar2d(gcm%nlat,gcm%nlon,gcm2mix_nvar))
        do v=1,gcm2mix_nvar
          call mpi_recv(gcm%invar2d(:,:,v), gcm%nlat*gcm%nlon, MPI_DOUBLE_PRECISION, gcmCplRank, (tgcmId+voltId)*100+v, gcmCplComm, MPI_STATUS_IGNORE,ierr)
          if(ierr /= MPI_Success) then
              call MPI_Error_string( ierr, message, length, ierr)
              print *,message(1:length)
              call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
          end if
        enddo
        call Toc("Passing")

      endif
      call Toc("MpiExchange")      

    end subroutine importgcmmpi

    subroutine exportgcmmpi(ion,gcm,mjd,time,gcmCplComm,gcmCplRank)
      type(gcm_T), intent(inout) :: gcm
      type(mixIon_T), dimension(:),intent(inout) :: ion
      real(rp), intent(in) :: time, mjd
      type(MPI_Comm) :: gcmCplComm
      integer, intent(in) :: gcmCplRank
      
      integer :: h,ymod2,v,ierr,length
      character( len = MPI_MAX_ERROR_STRING) :: message

      write(*,*) "MIX: EXPORT GCM"

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

      if (gcmCplComm /= MPI_COMM_NULL) then
      
        ! Prepare the export data
        if (.not. allocated(gcm%outvar2d)) allocate(gcm%outvar2d(gcm%nlat,gcm%nlon,mix2gcm_nvar))
        gcm%outvar2d = 0.
        ! Before we start, we collapse to 1 globe instead of 2 hemispheres
        do v=1,mix2gcm_nvar
          gcm%outvar2d(gcm%t2N(gcm%nhlat:1:-1),:gcm%nlon-1,v) = transpose(cshift(gcm%gcmOutput(GCMNORTH,v)%var(:gcm%nlon-1,:),-1*gcm%lonshift))
          gcm%outvar2d(gcm%t2S,:gcm%nlon-1,v) = transpose(cshift(gcm%gcmOutput(GCMSOUTH,v)%var(:gcm%nlon-1,:),-1*gcm%lonshift))
          gcm%outvar2d(gcm%t2N,gcm%nlon,v) = gcm%outvar2d(gcm%t2N,1,v)
          gcm%outvar2d(gcm%t2S,gcm%nlon,v) = gcm%outvar2d(gcm%t2S,1,v)
        end do

        ! Send the coupling data
        do v=1,mix2gcm_nvar
        write(*,*) " MIXCPL: ", gcmCplRank,(tgcmId+voltId)*100+v,gcmCplComm,gcm%nlat,gcm%nlon
          call mpi_send(gcm%outvar2d(:,:,v), gcm%nlat*gcm%nlon, MPI_DOUBLE_PRECISION, gcmCplRank, (tgcmId+voltId)*100+v, gcmCplComm, ierr)
          if(ierr /= MPI_Success) then
              call MPI_Error_string( ierr, message, length, ierr)
              print *,message(1:length)
              call mpi_Abort(MPI_COMM_WORLD, 1, ierr)
          end if
        enddo

      endif
      
    end subroutine exportgcmmpi
end module gcm_mpi
