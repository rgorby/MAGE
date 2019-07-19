module mix2mhd_interface
  use mixdefs
  use mixtypes
  use mixgeom
  use mixinterp
  use mixio
  use mixmain
  use ioH5
  use mixapp
  use gamapp
  
  implicit none

contains

     subroutine convertRemixToGamera(gameraApp, remixApp)
        type(mixApp_T), intent(inout) :: remixApp
        type(gamApp_T), intent(inout) :: gameraApp

        ! convert the "remixOutputs" variable to inEijk and inExyz, which are in
        ! Gamera coordinates
        integer :: i

         ! populate potential on gamera grid
         remixApp%gPsi = 0.0
         do i=1,remixApp%PsiShells+1
            remixApp%gPsi(i,:,gameraApp%Grid%ks:gameraApp%Grid%ke/2+1)   = remixApp%mixOutput(i,:,:,MHDPSI,NORTH)
            remixApp%gPsi(i,:,gameraApp%Grid%ke/2+1:gameraApp%Grid%ke+1) = remixApp%mixOutput(i,:,:,MHDPSI,SOUTH)
         enddo

        ! add corotation
        call CorotationPot(gameraApp%Model, gameraApp%Grid, remixApp%gPsi)

        ! find the remix BC to write data into
        SELECT type(iiBC=>gameraApp%Grid%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                call Ion2MHD(gameraApp%Model,gameraApp%Grid,remixApp%gPsi,iiBC%inEijk,iiBC%inExyz,remixApp%rm2g)
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in remix IC'
                stop
        END SELECT

    end subroutine convertRemixToGamera

  subroutine init_mix_mhd_interface(remixApp,mhdJGrid,mhdPsiGrid,optFilename)
    type(mixApp_T), intent(inout) :: remixApp
    real(rp), dimension(:,:,:,:,:), intent(in) :: mhdJGrid,mhdPsiGrid ! (i,j,k,x-z,hemisphere)
    character(len=*), optional, intent(in) :: optFilename    

    integer :: l,h ! h for hemisphere
    type(mixGrid_T) :: mhdGfpd,mhdG
    type(Map_T) :: Map
    real(rp), dimension(:,:), allocatable :: mhdt, mhdp, mhdtFpd, mhdpFpd
    real(rp) :: mhd_Rin  ! the radius of the shell given by the MHD grid

    if(present(optFilename)) then
        ! read from the prescribed file
        call init_mix(remixApp%ion,hmsphrs,remixApp%conductance,optFilename)
    else
        call init_mix(remixApp%ion,hmsphrs,remixApp%conductance)
    endif

    remixApp%PsiShells = size(mhdPsiGrid,1)
    remixApp%JShells = size(mhdJGrid,1)
    allocate(remixApp%PsiMaps(remixApp%PsiShells))
    allocate(remixApp%JMaps(remixApp%JShells))
    do h=1,size(remixApp%ion)
       ! set up interpolation map(s) for mhd2mix
       do l=1,remixApp%JShells
          call mix_mhd_grid(mhdJGrid(l,:,:,:,h),mhdt,mhdp,mhdtFpd,mhdpFpd,mhd_Rin)
          call init_grid_fromTP(mhdGfpd,mhdtFpd,mhdpFpd,.false.)  
          call flip_grid(remixApp%ion(h)%G,remixApp%mixGfpd,mhd_Rin) ! storing flipped
          ! grid for MIX only to
          ! use in mhd2mix
          ! below for zeroing
          ! out equatorward of
          ! MHD boundary, if
          ! necessary
          call mix_set_map(mhdGfpd,remixApp%mixGfpd,Map) 
          remixApp%JMaps(l) = Map
       end do

       ! set up interpolation map(s) for mix2mhd
       do l=1,remixApp%PsiShells
          call mix_mhd_grid(mhdPsiGrid(l,:,:,:,h),mhdt,mhdp,mhdtFpd,mhdpFpd,mhd_Rin)
          call init_grid_fromTP(mhdG,mhdt,mhdp,.false.)
          call mix_set_map(remixApp%ion(h)%G,mhdG,Map)
          remixApp%PsiMaps(l) = Map
       enddo
    enddo

  end subroutine init_mix_mhd_interface

  subroutine mix_mhd_output(ion,mhdvarsin,hmsphrs,time)
    type(mixIon_T),dimension(:),intent(inout) :: ion ! I for ionosphere (is an array of 1 or 2 elements for north and south)
    real(rp), dimension(:,:,:,:,:),intent(in) :: mhdvarsin
    integer, dimension(:), intent(in) :: hmsphrs ! array of integers marking hemispheres for the I object array.
    real(rp), intent(in) :: time

    character(strLen) :: fnstr,fname,vID
    real(rp), save :: next_t = 0.0
    logical :: mixOut,isThere,fExist
    integer, save :: step = 0
    real(rp) :: cpcp(2) = 0.0 

    type(IOVAR_T), dimension(MAXMIXIOVAR) :: IOVars

    mixOut = any(ion(:)%P%dtOut > 0.0)

    !Save CPCP for diagnostics
    cpcp(NORTH) = maxval(mhdvarsin(1,:,:,MHDPSI,NORTH))-minval(mhdvarsin(1,:,:,MHDPSI,NORTH))
    cpcp(SOUTH) = maxval(mhdvarsin(1,:,:,MHDPSI,SOUTH))-minval(mhdvarsin(1,:,:,MHDPSI,SOUTH))

    if (time >= next_t) then ! not in use jet
        write(*,*) '----- CMI -----'
        write(*,'(a,2f8.3)') 'N/S CPCP [kV] = ', cpcp(NORTH), cpcp(SOUTH)

        if (mixOut)then
          write(fnstr,'(I0.6)') floor(time/minval(ion(:)%P%dtOut)) !step
          fname = 'mixtest'//trim(fnstr)//'.h5'

          inquire(file=trim(fname),exist=fExist)
          if (.not. fExist) then
            !If the file doesn't exist
            call writeMIX(fname,ion,hmsphrs)

            !Add extra attribute information to output
            vID = "t"
            isThere = ioExist(fname,trim(vID))
            if (.not. isThere) then
              call ClearIO(IOVars)
              call AddOutVar(IOVars,"t"   ,time)
              call AddOutVar(IOVars,"ts"  ,step)
              call AddOutVar(IOVars,"nCPCP",cpcp(NORTH))
              call AddOutVar(IOVars,"sCPCP",cpcp(SOUTH))
              call WriteVars(IOVars,.true.,fname)
            endif !isThere
          endif !File exists
          
        end if

        next_t = next_t + minval(ion(:)%P%dtOut)

    end if

    step = step +1

  end subroutine mix_mhd_output

  ! assume what's coming here is mhdvars(i,j,k,var,hemisphere)
  ! thus the transposes below
  subroutine mhd2mix(remixApp)
    type(mixApp_T), intent(inout) :: remixApp

    real(rp), dimension(:,:), allocatable :: F
    integer :: l,h ! hemisphere
    integer :: v ! mhd var

    if (size(remixApp%mixInput,5).ne.size(remixApp%ion)) then
       write(*,*) 'The number of hemispheres in mhdvars is different from the size of the MIX ionosphere object. I am stopping.'
       stop
    end if

    do h=1,size(remixApp%ion)
       do v=1,size(remixApp%mixInput,4)
          do l=1,remixApp%JShells ! here we loop over Jshells but always use the last one (F)
             ! note the transpose to conform to the MIX layout (phi,theta)
             call mix_map_grids(remixApp%JMaps(l),transpose(remixApp%mixInput(l,:,:,v,h)),F)  

             ! note, cleaning MHD vars equatorward of the MHD boundary
             ! if the MIX boundary is equatorward of MHD boundary
             select case (v)
             case (MHDJ)
                ! zero out the current
                where (remixApp%mixGfpd%mask.eq.-1) F=0._rp
             case (MHDD, MHDC)
                ! set density and sound speed to min values
                ! this helps with conductdance calculation
                where (remixApp%mixGfpd%mask.eq.-1) F=minval(remixApp%mixInput(l,:,:,v,h))
             end select
          end do

          select case (v)
          case (MHDJ)
             remixApp%ion(h)%St%Vars(:,:,FAC) = F
          case (MHDD)
             remixApp%ion(h)%St%Vars(:,:,DENSITY) = F
          case (MHDC)
             remixApp%ion(h)%St%Vars(:,:,SOUND_SPEED) = F
          end select
       end do
    end do

    call run_mix(remixApp%ion,remixApp%tilt,remixApp%conductance)
  end subroutine mhd2mix

  subroutine mix2mhd(remixApp)
    type(mixApp_T), intent(inout) :: remixApp
    integer :: l,h ! hemisphere
    integer :: v ! mhd var

    real(rp), dimension(:,:), allocatable :: gPsi_tmp  ! gamera potential

    ! map to potential shells
    do h=1,size(remixApp%ion)
       do v=1,size(remixApp%mixOutput,4)  ! mirroring mix2mhd here, but for now only one variable 
          do l=1,remixApp%PsiShells
             call mix_map_grids(remixApp%PsiMaps(l),remixApp%ion(h)%St%Vars(:,:,POT),gPsi_tmp)  
             ! this is going back to gamera
             select case (v)
             case (MHDPSI)
                remixApp%mixOutput(l,:,:,v,h) = gPsi_tmp   ! north
             end select
          enddo
       enddo
    enddo
  end subroutine mix2mhd

  ! This function gets the MHD grid from MHD directly
  subroutine mix_mhd_grid(mhdg,t,p,tFpd,pFpd,Rinner)
    real(rp), dimension(:,:,:), intent(in) :: mhdg ! MHD grid cell coords
    real(rp), dimension(:,:), allocatable, intent(out) :: t,p,tFpd,pFpd
    real(rp), dimension(:,:), allocatable :: r
    real(rp), intent(out) :: Rinner
    real(rp), dimension(:,:), allocatable :: dipRatio
    integer :: nj, nk2

    nj = size(mhdg,1); nk2 = size(mhdg,2) ! nk2+1 for Psi shells but keep the notation for brevity

    ! Also define t, p variables in the gamera coordinate space (x-axis is the spherical axis)
    ! "Fpd" suffix stands for "flipped"
    if (.not.allocated(t)) allocate(t(nj,nk2))      
    if (.not.allocated(p)) allocate(p(nj,nk2))      
    if (.not.allocated(r)) allocate(r(nj,nk2))      
    if (.not.allocated(dipRatio)) allocate(dipRatio(nj,nk2))      

    if (.not.allocated(tFpd)) allocate(tFpd(nk2,nj))
    if (.not.allocated(pFpd)) allocate(pFpd(nk2,nj))      

    ! NOTE: this is already mapped to the ionosphere
    r = sqrt(mhdg(:,:,1)**2+mhdg(:,:,2)**2+mhdg(:,:,3)**2)

    ! capture for the rare (but possible!) case where the psiShell
    ! we're mapping to is below the ionosphere
    ! this may happen, e.g., for the bottom-most ghost cell inside the Gamera inner boundary.
    ! In this case, we simply set theta=pi/2 so mix_map_grids (in mixinterp.F90) will simply 
    ! extrapolate equatorward from the MIX low lat boundary
    dipRatio = sqrt(mhdg(:,:,1)**2+mhdg(:,:,2)**2)/r**1.5
    where( dipRatio < 1. )
       t = asin(dipRatio)
    elsewhere
       t = PI/2.
    end where

    p = modulo(atan2(mhdg(:,:,2),mhdg(:,:,1)),2*pi)
    Rinner = sum(r)/size(r)

    ! NOTE, transposition is necessary to make sure the phi coordinate goes first
    ! because all the codes, particulalry, interpolation (mixinterp) assume that
    tFpd = transpose(acos(mhdg(:,:,1)/r))
    pFpd = transpose(modulo(atan2(mhdg(:,:,3),mhdg(:,:,2)),2*pi))

    ! Fixes for south
    if (mhdg(nj/2,nk2/2,2).lt.0) then ! pick the pole and see if
       ! z<0. this works even if nk2
       ! is nk2+1 (for PsiShells)
       ! For south, phi obtained this way goes from pi to 2pi.  so we
       ! subtract pi because we always assume one hemisphere when
       ! mapping
       pFpd = pFpd-pi ! this is phi in the lfm/gamera space

       ! FIXME: make sure this actually works
       p = modulo(atan2(-mhdg(:,:,2),mhdg(:,:,1)),2*pi)
    end if
  end subroutine mix_mhd_grid

  ! This is an ancillary function to get data from gamera dump
  ! In production, the mhd grid would be coming from gamera directly
  subroutine mhd_fromFile(fname,shells,mhdg) 
    character(len=*),intent(in) :: fname
    real(rp), dimension(:,:,:,:,:), allocatable, intent(out) :: mhdg ! MHD grid cell coords (i,j,k,x-z,hemisphere)
    integer :: shells ! number of shells to extract
    real(iop), dimension(:,:,:), allocatable :: gx, gy, gz
    integer :: nip1,njp1,nkp1,nj,nk2,h

    call readMHDVar(fname,"X",gx)
    call readMHDVar(fname,"Y",gy)
    call readMHDVar(fname,"Z",gz)

    nip1 = size(gx,1); njp1 = size(gx,2); nkp1 = size(gx,3)
    nj = njp1-1; nk2 = (nkp1-1)/2

    if (.not.allocated(mhdg)) allocate(mhdg(shells,nj,nk2,3,2))

    ! Loop over hemispheres
    do h=1,2
       ! offset of (h-1)*nk2 is 0 for north (h=1) and nk/2 for south (h=2)
       mhdg(1:shells,:,:,1,h) = 0.125*(&
            gx(1:shells,1:nj,    1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gx(1:shells,1:nj,    2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gx(1:shells,2:njp1,  1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gx(1:shells,2:njp1,  2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gx(2:shells+1,1:nj,  1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gx(2:shells+1,1:nj,  2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gx(2:shells+1,2:njp1,1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gx(2:shells+1,2:njp1,2+(h-1)*nk2:nk2+1+(h-1)*nk2))

       mhdg(1:shells,:,:,2,h) = 0.125*(&
            gy(1:shells,1:nj,    1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gy(1:shells,1:nj,    2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gy(1:shells,2:njp1,  1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gy(1:shells,2:njp1,  2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gy(2:shells+1,1:nj,  1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gy(2:shells+1,1:nj,  2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gy(2:shells+1,2:njp1,1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gy(2:shells+1,2:njp1,2+(h-1)*nk2:nk2+1+(h-1)*nk2))

       mhdg(1:shells,:,:,3,h) = 0.125*(&
            gz(1:shells,1:nj,    1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gz(1:shells,1:nj,    2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gz(1:shells,2:njp1,  1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gz(1:shells,2:njp1,  2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gz(2:shells+1,1:nj,  1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gz(2:shells+1,1:nj,  2+(h-1)*nk2:nk2+1+(h-1)*nk2)+&
            gz(2:shells+1,2:njp1,1+(h-1)*nk2:nk2  +(h-1)*nk2)+&
            gz(2:shells+1,2:njp1,2+(h-1)*nk2:nk2+1+(h-1)*nk2))
    enddo
  end subroutine mhd_fromFile

end module mix2mhd_interface
