!Various routines to handle coupling between Gamera magnetosphere and remix ionosphere
!General coupling consists of XXX

module mhd2mix_interface
    use types
    use mixtypes
    use gamutils
    use math
    use clocks
    use xml_input
    use gridutils
    use msphutils !Get scaling values
    !Debugging stuff
    use ioH5
    use hdf5
    use h5lt
    use gamapp
    
    implicit none

    ! how many variables are we sending (should be consistent with the enumerator in mixdefs.F90)
    integer, parameter :: mhd2mix_varn = 3

    type mhd2Mix_T

        ! data for gamera -> remix conversion
        real(rp), dimension(:,:,:,:,:), allocatable :: mixInput
        real(rp), dimension(:,:,:,:), allocatable :: gJ
        type(Map_T), allocatable, dimension(:) :: Jmaps
        integer :: JStart = 2, JShells = 1

    end type mhd2mix_T

    contains

    subroutine init_mhd2Mix(mhd2mix, gameraApp)
        type(mhd2Mix_T), intent(inout) :: mhd2mix
        type(gamApp_T), intent(inout) :: gameraApp

        integer :: i,j,k,iG
        real(rp) :: xc,yc,zc
        real(rp), allocatable, dimension(:,:,:,:,:) :: mhdJGrid

        mhd2Mix%rm2g = gB0*gV0*gx0*1.0e-12 !Scaling factor for remix potential [kV]

        ! allocate remix arrays
        allocate(mhd2Mix%gJ(1:mhd2Mix%JShells, gameraApp%Grid%js:gameraApp%Grid%je, gameraApp%Grid%ks:GameraApp%Grid%ke, 1:NDIM))
        allocate(mhdJGrid(1:mhd2Mix%JShells, gameraApp%Grid%js:gameraApp%Grid%je, gameraApp%Grid%ks:gameraApp%Grid%ke/2, 1:3, 1:2))
        allocate(mhd2Mix%mixInput(1:mhd2Mix%JShells, gameraApp%Grid%js:gameraApp%Grid%je, gameraApp%Grid%ks:gameraApp%Grid%ke/2, 1:mhd2mix_varn, 1:2))
        allocate(mhd2Mix%JMaps(mhd2Mix%JShells))

        ! get those grid coordinates (cell centers for Jp)
        do k=gameraApp%Grid%ks,gameraApp%Grid%ke
            do j=gameraApp%Grid%js,gameraApp%Grid%je
                do i=1,remixApp%JShells
                    iG = remixApp%JStart+i-1
                    call cellCenter(gameraApp%Grid,iG,j,k,xc,yc,zc)

                    ! note conversion to Rion units which are expected on the remix
                    ! side
                    if (k<=gameraApp%Grid%ke/2) then
                        mhdJGrid(i,j,k,:,NORTH) = [xc,yc,zc]/Rion
                    else
                        mhdJGrid(i,j,k-ke/2,:,SOUTH) = [xc,yc,zc]/Rion
                    endif
                enddo
            enddo
        enddo

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
        end do

        deallocate(mhdJGrid)

    end subroutine init_mhd2Mix

    subroutine convertGameraToRemix(gameraApp, remixApp)
        type(mixApp_T), intent(inout) :: remixApp
        type(gamApp_T), intent(in) :: gameraApp

        ! convert incoming gamera data to the "remixInputs" variable
        real(rp) ::  B0mag,Bi2m
        integer :: i,j,k,iG
        real(rp) :: xc,yc,zc

        real(rp) :: Con(NVAR)
        real(rp) :: Cs


        !write(*,*) 'Prepping remix data at T = ',Model%t

        !Only working on Bxyz from perturbation
        !B0 in inner region (where we care for remix)
        !Assumd to be current-free
        call GetShellJ(gameraApp%Model, gameraApp%Grid, gameraApp%State%Bxyz, remixApp%gJ)

        !Now loop over only shells and populate mhdvars to be sent to mix
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,iG,j,k,B0mag,Bi2m,xc,yc,zc,Con,Cs)
        do k=gameraApp%Grid%ks,gameraApp%Grid%ke
            do j=gameraApp%Grid%js,gameraApp%Grid%je
                do i=1,remixApp%JShells
                    iG = remixApp%JStart+i-1
                    B0mag = sqrt(dot_product(gameraApp%Grid%B0(iG,j,k,:),gameraApp%Grid%B0(iG,j,k,:)))
                    ! get cell centers
                    call cellCenter(gameraApp%Grid,iG,j,k,xc,yc,zc)
                    !Get ion2mag scaling factor for field
                    Bi2m = BIon2Mag(xc,yc,zc)

                    if (k<=gameraApp%Grid%ke/2) then
                        !!! NOTE: assuming gB0 in nT and gx0 in m and gv0 in m/s

                        ! note conversion to microA/m^2
                        remixApp%mixInput(i,j,k,MHDJ,NORTH) = dot_product(remixApp%gJ(i,j,k,:),gameraApp%Grid%B0(iG,j,k,:)/B0mag)*Bi2m*(gB0/gx0*1.e4/4/PI)
                        ! note conversion to g/cm^3
                        remixApp%mixInput(i,j,k,MHDD,NORTH) = gameraApp%State%Gas(iG,j,k,DEN,BLK)*(gB0*1.e-7/gv0)**2/4/pi
                        ! get sound speed first
                        Con = gameraApp%State%Gas(iG,j,k,:,BLK)
                        call CellPress2Cs(gameraApp%Model,Con,Cs)
                        remixApp%mixInput(i,j,k,MHDC,NORTH) = Cs*gv0*1.e2
                    else
                        remixApp%mixInput(i,j,k-gameraApp%Grid%ke/2,MHDJ,SOUTH) = dot_product(remixApp%gJ(i,j,k,:),gameraApp%Grid%B0(iG,j,k,:)/B0mag)*Bi2m*(gB0/gx0*1.e4/4/PI)
                        remixApp%mixInput(i,j,k-gameraApp%Grid%ke/2,MHDD,SOUTH) = gameraApp%State%Gas(iG,j,k,DEN,BLK)*(gB0*1.e-7/gv0)**2/4/pi
                        ! get sound speed first
                        Con = gameraApp%State%Gas(iG,j,k,:,BLK)
                        call CellPress2Cs(gameraApp%Model,Con,Cs)
                        remixApp%mixInput(i,j,k-ke/2,MHDC,SOUTH) = Cs*gv0*1.e2
                    endif
                enddo
            enddo
        enddo

    end subroutine convertGameraToRemix

    !Calculate Ion->Mag scaling factor, xyz are in units of Re (like Gamera)
    ! note for Bi2m to be correct, xc,yc,zc,r need to be in
    ! units of Rion where we solve for the potential
    ! Rion is defined in msphutils in code units

    function BIon2Mag(x,y,z) result(Bi2m)
        real(rp), intent(in) :: x,y,z
        real(rp) :: Bi2m

        real(rp) :: rScl,zScl,zor2

        rScl = norm2([x,y,z])/Rion
        zScl = z/Rion
        zor2 = (zScl/rScl)**2.0

        Bi2m = (rScl**3.0)*sqrt(1+3*(1-zor2)/rScl)/sqrt(1+3*zor2)
        ! bion2bmag in LFM-para/src/interfaces/MHDBoundaryInterface.C
        !Bi2m  = r**3*sqrt(1+3*(1-(1-(zc/r)**2)/r))/sqrt(1+3*(zc/r)**2) 


    end function BIon2Mag

    function need_remix_BC(Model,Grid)
        type(Model_T) ,intent(in) :: Model
        type(Grid_T)  ,intent(in) :: Grid
        logical need_remix_BC
     
        need_remix_BC = Grid%ijkShift(1) == 0

    end function need_remix_BC

end module mhd2mix_interface
