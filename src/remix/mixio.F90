! Input/output for mix

module mixio
  use hdf5
  use mixdefs  
  use mixtypes

  implicit none

  !Necessary for HDF5 routines
  integer :: herror

  contains

    ! get UT: int array (Year, Month, Day, Hr, Min, Sec)
    subroutine getUT(fname,simtime)
      character(len=*),intent(in) :: fname
      integer, dimension(6), target, intent(out) :: simtime
      
      integer(HID_T) :: h5fId,attrId
#ifdef NO2003_HDF5
      integer(HSIZE_T), dimension(1) :: dims
#else
      TYPE(C_PTR) :: f_ptr
#endif

      call checkFile(fname)
      !Open file
      call h5fopen_f(fname, H5F_ACC_RDONLY_F, h5fId, herror)

      call h5aopen_f(h5fId,"SimTime",attrId,herror)
#ifdef NO2003_HDF5
      call h5aread_f(attrId,H5T_NATIVE_INTEGER,simtime,dims,herror)
#else
      f_ptr = C_LOC(simtime)
      call h5aread_f(attrId,H5T_NATIVE_INTEGER,f_ptr,herror)
#endif
      
      call h5aclose_f(attrId,herror)
      call h5fclose_f(h5fId, herror)
    end subroutine getUT

    subroutine readVar(fname,varname,var)
      character(len=*),intent(in) :: fname
      character(len=*),intent(in) :: varname
      
      integer(HID_T) :: h5fId,dsId,dspaceId
      integer(HSIZE_T), dimension(2) :: dims,maxdims

#ifndef NO2003_HDF5
      TYPE(C_PTR) :: f_ptr
#endif

      real, dimension(:,:), target, allocatable :: var

      call checkFile(fname)
      !Open file
      call h5fopen_f(fname, H5F_ACC_RDONLY_F, h5fId, herror)

      call h5dopen_f(h5fId,varname,dsId,herror)
      call h5dget_space_f(dsId,dspaceId,herror);
      call h5sget_simple_extent_dims_f(dspaceId, dims, maxdims, herror) 

      if (.not. allocated(var)) then 
         allocate(var(1:dims(1),1:dims(2)))
      end if

#ifdef NO2003_HDF5
      call h5dread_f(dsId,H5T_NATIVE_REAL,var,dims,herror)
#else
      f_ptr = C_LOC(var)
      call h5dread_f(dsId,H5T_NATIVE_REAL,f_ptr,herror)
#endif
      
      call h5dclose_f(dsId,herror)
      call h5fclose_f(h5fId, herror)
    end subroutine readVar

    subroutine writeMIX(fname,I)
      type(mixIon_T),dimension(:),intent(in) :: I
      character(len=*), intent(in) :: fname
      integer(HID_T) :: h5fId,h5gId
      integer :: h

      call h5open_f(herror) !Setup H5 Fortran interface            
      call h5fcreate_f(fname,H5F_ACC_TRUNC_F, h5fId, herror) 
      call h5fopen_f(fname,H5F_ACC_RDWR_F, h5fId, herror)
      do h=1,size(I)
         if (h.eq.NORTH) then
            call h5gcreate_f (h5fId,"NORTH",h5gId,herror)
         else if (h.eq.SOUTH) then
            call h5gcreate_f (h5fId,"SOUTH",h5gId,herror)
         else
            write(*,*) "writeMIX: Wrong hemisphere identifier. Stopping..."
            stop
         end if 
         
         call writeState(h5gId,I(h)%St)
         call writeGrid(h5gId,I(h)%G)

         call h5gclose_f (h5gId,herror)
      end do

      call h5fclose_f(h5fId, herror)
      call h5close_f(herror)  ! Close H5 Fortran interface
    end subroutine writeMIX


    subroutine writeState(h5gId,St)
      type(mixState_T), intent(in) :: St
      integer(HID_T), intent(in) :: h5gId
      integer :: v

      do v=1,nVars
         select case (v)
            case (POT)
               call writeVar(h5gId,St%Vars(:,:,v),"Potential","kV")
            case (FAC)
               call writeVar(h5gId,St%Vars(:,:,v),"Field-aligned current","muA/m**2")
            case (SIGMAP)
               call writeVar(h5gId,St%Vars(:,:,v),"Pedersen conductance","S")
            case (SIGMAH)
               call writeVar(h5gId,St%Vars(:,:,v),"Hall conductance","S")
            case (AVG_ENG)
               call writeVar(h5gId,St%Vars(:,:,v),"Average energy","keV")
            case (NUM_FLUX)
               call writeVar(h5gId,St%Vars(:,:,v),"Number flux","1/cm^2 s")
         end select
      enddo

    end subroutine writeState

    subroutine writeGrid(h5gId,G)
      type(mixGrid_T), intent(in) :: G
      integer(HID_T) :: h5gId

      call writeVar(h5gId,G%x,'X',"Ri")
      call writeVar(h5gId,G%y,'Y',"Ri")
    end subroutine writeGrid

    subroutine writeVar(gId,var,varName,units)
      integer(HID_T), intent(in) :: gId
      real(rp), dimension(:,:), intent(in) :: var
      character(len=*),intent(in) :: varName
      character(len=*),intent(in) :: units
      
      integer(HID_T) :: dsId,dspaceId,aspaceId,attrId,atypeId
      integer(HSIZE_T), dimension(2) :: dims
      integer(HSIZE_T), dimension(1) :: adims = [1]
      integer(SIZE_T) :: attrlen
      integer :: rank = 2, arank = 1

      attrlen = len(units)      

      dims = shape(var)
      ! create data space
      call h5screate_simple_f(rank, dims,dspaceId,herror)
      !Create data set and write data
      call h5dcreate_f(gId,varName,H5T_NATIVE_REAL, dspaceId,dsId, herror)

      ! data set attributes
      call h5screate_simple_f(arank, adims, aspaceId,herror)
      call h5tcopy_f(H5T_NATIVE_CHARACTER, atypeId, herror)
      call h5tset_size_f(atypeId,attrlen, herror)
      call h5acreate_f (dsId,"Units",atypeId,aspaceId,attrId,herror)
      call h5awrite_f (attrId,atypeId,trim(units),adims,herror) 
      call h5aclose_f (attrId,herror)

      ! note type conversion to io_real (single, typically)
      call h5dwrite_f(dsId, H5T_NATIVE_REAL,real(var,iop),dims, herror)

      !Close up shop
      call h5dclose_f(dsId,herror)
      call h5sclose_f(dspaceId,herror)
    end subroutine writeVar

    ! check H5 file existence
    subroutine checkFile(fname)
      character(len=*),intent(in) :: fname
      logical :: fExist

      inquire(file=fname,exist=fExist)
      if (.not.(fExist)) then
         write(*,"(a,a)") "Cannot read file ", fname
         stop
      endif
    end subroutine checkFile

    subroutine readMHDVar(fname,varname,var)
      character(len=*),intent(in) :: fname
      character(len=*),intent(in) :: varname
      
      integer(HID_T) :: h5fId,dsId,dspaceId
      integer(HSIZE_T), dimension(3) :: dims,maxdims

#ifndef NO2003_HDF5
      TYPE(C_PTR) :: f_ptr
#endif

      real, dimension(:,:,:), target, allocatable :: var

      call checkFile(fname)
      !Open file
      call h5fopen_f(fname, H5F_ACC_RDONLY_F, h5fId, herror)

      call h5dopen_f(h5fId,varname,dsId,herror)
      call h5dget_space_f(dsId,dspaceId,herror);
      call h5sget_simple_extent_dims_f(dspaceId, dims, maxdims, herror) 

      if (.not. allocated(var)) then 
         allocate(var(1:dims(1),1:dims(2),1:dims(3)))
      end if

#ifdef NO2003_HDF5
      call h5dread_f(dsId,H5T_NATIVE_REAL,var,dims,herror)
#else
      f_ptr = C_LOC(var)
      call h5dread_f(dsId,H5T_NATIVE_REAL,f_ptr,herror)
#endif
      
      call h5dclose_f(dsId,herror)
      call h5fclose_f(h5fId, herror)
    end subroutine readMHDVar

    ! Version that reads variables stored under group
    ! E.g., gamera stores X, Y, Z as higher level datasets
    ! while physical vars are stored under groups "Step#0" etc.
    subroutine readGMHDVar(fname,groupname,varname,var)
      character(len=*),intent(in) :: fname
      character(len=*),intent(in) :: groupname
      character(len=*),intent(in) :: varname

      
      integer(HID_T) :: h5fId,dsId,dspaceId,gId
      integer(HSIZE_T), dimension(3) :: dims,maxdims

#ifndef NO2003_HDF5
      TYPE(C_PTR) :: f_ptr
#endif

      real, dimension(:,:,:), target, allocatable :: var

      call checkFile(fname)

      !Open file
      call h5fopen_f(fname, H5F_ACC_RDONLY_F, h5fId, herror)

      ! Open group
      call h5gopen_f(h5fId,trim(groupname),gId,herror)


      call h5dopen_f(gId,varname,dsId,herror)
      call h5dget_space_f(dsId,dspaceId,herror);
      call h5sget_simple_extent_dims_f(dspaceId, dims, maxdims, herror) 

      if (.not. allocated(var)) then 
         allocate(var(1:dims(1),1:dims(2),1:dims(3)))
      end if

#ifdef NO2003_HDF5
      call h5dread_f(dsId,H5T_NATIVE_REAL,var,dims,herror)
#else
      f_ptr = C_LOC(var)
      call h5dread_f(dsId,H5T_NATIVE_REAL,f_ptr,herror)
#endif
      
      call h5dclose_f(dsId,herror)
      call h5gclose_f(gId,herror)
      call h5fclose_f(h5fId, herror)
    end subroutine readGMHDVar

    subroutine potMinMax(I)
      type(mixIon_T),dimension(:),intent(in) :: I

      integer :: h ! hemisphere counter
      do h=1,size(I)
         write(*,*) 'Hemisphere:',h,'; Min/Max potential',minval(I(h)%St%Vars(:,:,POT)),maxval(I(h)%St%Vars(:,:,POT))
      end do
      
    end subroutine potMinMax
    
end module mixio
