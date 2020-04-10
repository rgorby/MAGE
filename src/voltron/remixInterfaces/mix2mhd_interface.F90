module mix2mhd_interface
  use volttypes
  use mixdefs
  use mixtypes
  use mixgeom
  use mixinterp
  use mixio
  use mixmain
  use ioH5
  use gamapp
  use msphutils
  use mixinterfaceutils
  use uservoltic ! required to have IonInnerBC_T defined
  
  implicit none

  integer, parameter :: mix2mhd_varn = 1  ! for now just the potential is sent back

contains

     subroutine init_mix2Mhd(mix2mhd, remixApp, gameraApp)
        type(mix2Mhd_T), intent(inout) :: mix2mhd
        type(mixApp_T), intent(inout) :: remixApp
        type(gamApp_T), intent(inout) :: gameraApp

        integer :: i,j,k,iG,h,l
        real(rp) :: mhd_Rin
        real(rp), allocatable, dimension(:,:,:,:,:) :: mhdPsiGrid
        real(rp), allocatable, dimension(:,:) :: mhdt, mhdp, mhdtFpd, mhdpFpd
        type(mixGrid_T) :: mhdG
        type(Map_T) :: Map

        mix2Mhd%rm2g = gB0*gV0*gx0*1.0e-12 !Scaling factor for remix potential [kV]

        ! allocate remix arrays
        allocate(mix2mhd%gPsi(1:mix2mhd%PsiShells+1,gameraApp%Grid%js:gameraApp%Grid%je+1,gameraApp%Grid%ks:gameraApp%Grid%ke+1))
        allocate(mhdPsiGrid(1:mix2mhd%PsiShells+1, gameraApp%Grid%js:gameraApp%Grid%je+1, gameraApp%Grid%ks:gameraApp%Grid%ke/2+1, 1:3, 1:2))
        allocate(mix2mhd%mixOutput(1:mix2mhd%PsiShells+1, gameraApp%Grid%js:gameraApp%Grid%je+1, gameraApp%Grid%ks:gameraApp%Grid%ke/2+1, 1:mix2mhd_varn, 1:2))
        allocate(mix2mhd%PsiMaps(mix2mhd%PsiShells))

        ! get those grid coordinates (corner centers for Psi)
        do k=gameraApp%Grid%ks,gameraApp%Grid%ke+1
            do j=gameraApp%Grid%js,gameraApp%Grid%je+1
                ! note, PsiShells give shell numbers based on cell centers per our
                ! convenion
                ! thus, no -1 below
                do i=1,mix2mhd%PsiShells+1
                    iG = mix2mhd%PsiStart+i-1

                    ! note conversion to Rion units which are expected on the remix
                    ! side
                    if (k<=gameraApp%Grid%ke/2+1) then
                        mhdPsiGrid(i,j,k,:,NORTH) = gameraApp%Grid%xyz(iG,j,k,:)/Rion
                    endif
                    if (k>=gameraApp%Grid%ke/2+1) then
                        mhdPsiGrid(i,j,k-gameraApp%Grid%ke/2,:,SOUTH) = gameraApp%Grid%xyz(iG,j,k,:)/Rion
                    endif
                enddo
            enddo
        enddo

        do h=1,size(remixApp%ion)
            ! set up interpolation map(s) for mix2mhd
            do l=1,mix2mhd%PsiShells
                call mix_mhd_grid(mhdPsiGrid(l,:,:,:,h),mhdt,mhdp,mhdtFpd,mhdpFpd,mhd_Rin)
                call init_grid_fromTP(mhdG,mhdt,mhdp,.false.)
                call mix_set_map(remixApp%ion(h)%G,mhdG,Map)
                mix2mhd%PsiMaps(l) = Map
            enddo
        enddo
 
        deallocate(mhdPsiGrid)

     end subroutine init_mix2Mhd

     subroutine convertRemixToGamera(mix2mhd, remixApp, gameraApp, doCorotO)
        type(mix2Mhd_T), intent(inout) :: mix2mhd
        type(mixApp_T), intent(inout) :: remixApp
        type(gamApp_T), intent(inout) :: gameraApp
        logical, intent(in), optional :: doCorotO

        ! convert the "remixOutputs" variable to inEijk and inExyz, which are in
        ! Gamera coordinates
        integer :: i
        logical :: doCorot

         ! populate potential on gamera grid
         mix2mhd%gPsi = 0.0
         do i=1,mix2mhd%PsiShells+1
            mix2mhd%gPsi(i,:,gameraApp%Grid%ks:gameraApp%Grid%ke/2+1)   = mix2mhd%mixOutput(i,:,:,MHDPSI,NORTH)
            mix2mhd%gPsi(i,:,gameraApp%Grid%ke/2+1:gameraApp%Grid%ke+1) = mix2mhd%mixOutput(i,:,:,MHDPSI,SOUTH)
         enddo

        if (present(doCorotO)) then
          doCorot = doCorotO
        else
          doCorot = .true.
        endif

        ! add corotation
        if (doCorot) call CorotationPot(gameraApp%Model, gameraApp%Grid, mix2mhd%gPsi)

        ! find the remix BC to write data into
        SELECT type(iiBC=>gameraApp%Grid%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                call Ion2MHD(gameraApp%Model,gameraApp%Grid,mix2mhd%gPsi,iiBC%inEijk,iiBC%inExyz,mix2mhd%rm2g)
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in remix IC'
                stop
        END SELECT

    end subroutine convertRemixToGamera

  subroutine mapRemixToGamera(mix2mhd, remixApp)
    type(mix2Mhd_T), intent(inout) :: mix2mhd
    type(mixApp_T), intent(inout) :: remixApp
    integer :: l,h ! hemisphere
    integer :: v ! mhd var

    real(rp), dimension(:,:), allocatable :: gPsi_tmp  ! gamera potential

    ! map to potential shells
    do h=1,size(remixApp%ion)
       do v=1,size(mix2mhd%mixOutput,4)  ! mirroring mix2mhd here, but for now only one variable 
          do l=1,mix2mhd%PsiShells
             call mix_map_grids(mix2mhd%PsiMaps(l),remixApp%ion(h)%St%Vars(:,:,POT),gPsi_tmp)  
             ! this is going back to gamera
             select case (v)
             case (MHDPSI)
                mix2mhd%mixOutput(l,:,:,v,h) = gPsi_tmp   ! north
             end select
          enddo
       enddo
    enddo
  end subroutine mapRemixToGamera

end module mix2mhd_interface
