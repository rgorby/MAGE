!Various routines to handle coupling between Gamera magnetosphere and remix ionosphere
!General coupling consists of XXX

module mhd_mix_interface
    use types
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
    
    !Remix modules
    use mix_mhd_interface, ONLY: mixIon_T, NORTH,SOUTH, &
           MHDC, MHDD, MHDJ, MHDPSI

    implicit none

    !!!!!!!!!!! CMI exchange variables !!!!!!!!!!!!!!!!!!!
    !
    ! Gamera normalization
    ! Scaling factor for remix potential [kV]
    real(rp) :: rm2g
    
    integer, parameter :: mhd2mix_varn = 3  ! how many variables are we sending (should be consistent with the enumerator in mixdefs.F90)
    integer, parameter :: mix2mhd_varn = 1  ! for now just the potential is sent back
    
    ! Electric potential on Gamera inner shell (node-centered)
    real(rp), allocatable, dimension(:,:,:) :: gPsi
    ! field aligned currrent
    real(rp), allocatable, dimension(:,:,:,:) :: gJ

    ! MHD grids to pass to mix
    real(rp), allocatable, dimension(:,:,:,:,:) :: mhdJGrid,mhdPsiGrid
    
    real(rp) :: tilt = -9999._rp

    contains

    !Initialize everything necessary for CMI
    subroutine InitCMI(Model,Grid)
        type(Model_T),intent(in) :: Model
        type(Grid_T) ,intent(in) :: Grid

        integer :: i,j,k,iG
        integer :: h ! hemisphere
        real(rp) :: xc,yc,zc

         rm2g = gB0*gV0*gx0*1.0e-12 !Scaling factor for remix potential [kV]
        
        if (JpSh > 1) then
           ! support is there on the MIX side (although not tested),
           ! but it will simply choose the last shell always, becasue
           ! it doesn't really know what to do with multiple shells.
           write(*,*) 'You are entering a world of pain, this will not work!'
           stop
        endif

        !Allocate outgoing grid
        allocate(mhdJGrid(1:JpSh,Grid%js:Grid%je,Grid%ks:Grid%ke/2,3,2))      ! cell centered
        allocate(mhdPsiGrid(1:PsiSh+1,Grid%js:Grid%je+1,Grid%ks:Grid%ke/2+1,3,2)) ! corner centered

        !i = local index (1-...)
        !iG = global index, ie to map into Grid array
        ! get those grid coordinates (cell centers for Jp)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=1,JpSh
                    iG = JpSt+i-1
                    call cellCenter(Grid,iG,j,k,xc,yc,zc)

                    ! note conversion to Rion units which are expected on the remix side
                    if (k<=Grid%ke/2) then
                        mhdJGrid(i,j,k,:,NORTH) = [xc,yc,zc]/Rion
                    else
                        mhdJGrid(i,j,k-Grid%ke/2,:,SOUTH) = [xc,yc,zc]/Rion
                    endif
                enddo
            enddo
        enddo

        ! get those grid coordinates (corner centers for Psi)
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                ! note, PsiShells give shell numbers based on cell centers per our convenion
                ! thus, no -1 below
                do i=1,PsiSh+1
                    iG = PsiSt+i-1

                    ! note conversion to Rion units which are expected on the remix side
                    if (k<=Grid%ke/2+1) then
                        mhdPsiGrid(i,j,k,:,NORTH) = [Grid%x(iG,j,k),Grid%y(iG,j,k),Grid%z(iG,j,k)]/Rion
                    endif
                    if (k>=Grid%ke/2+1) then
                        mhdPsiGrid(i,j,k-Grid%ke/2,:,SOUTH) = [Grid%x(iG,j,k),Grid%y(iG,j,k),Grid%z(iG,j,k)]/Rion
                    endif
                enddo
            enddo
        enddo

        !Setup holders for Gamera calculation of fac, density, Cs and potential
        !call allocGridVec(Model,Grid,gJ,.false.,NDIM)
        allocate(gJ(1:JpSh,Grid%js:Grid%je,Grid%ks:Grid%ke,1:NDIM))
        allocate(mhdvarsout(1:JpSh,    Grid%js:Grid%je,  Grid%ks:Grid%ke/2,  mhd2mix_varn,2))  
        allocate(mhdvarsin(1:PsiSh+1,Grid%js:Grid%je+1,Grid%ks:Grid%ke/2+1,mix2mhd_varn,2))
        allocate(gPsi(PsiSh+1,Grid%js:Grid%je+1,Grid%ks:Grid%ke+1)) 

        ! Init mix (both hemispheres) and interpolation maps
        !call init_mix_mhd_interface(ion,hmsphrs,mhdJGrid,mhdPsiGrid)

    end subroutine InitCMI

    !Update for CMI, called every timestep
    subroutine UpdateCMI(Model,Grid,inEijk,inExyz)
        type(Model_T),intent(in) :: Model
        type(Grid_T) ,intent(in) :: Grid
        real(rp), intent(inout) :: inEijk(1:PsiSh+1,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), intent(inout) :: inExyz(1:PsiSh,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)

        real(rp) ::  B0mag,r,Bi2m
        integer :: i,j,k,iG
        real(rp) :: xc,yc,zc
        real(rp) :: Con(NVAR)
        real(rp) :: Cs

         ! populate potential on gamera grid
         do i=1,PsiSh+1
            gPsi(i,:,Grid%ks:Grid%ke/2+1)   = mhdvarsin(i,:,:,MHDPSI,NORTH)
            gPsi(i,:,Grid%ke/2+1:Grid%ke+1) = mhdvarsin(i,:,:,MHDPSI,SOUTH)
         enddo

         ! add corotation
         call CorotationPot(Model,Grid,gPsi)

         !Use potential to set E field values
         inEijk = 0.0
         inExyz = 0.0
         ! FIXME: rename to avoid confusion with mix2mhd
         call Ion2MHD(Model,Grid,gPsi,inEijk,inExyz,rm2g)
        
        
    end subroutine UpdateCMI

    !Update for CMI, called every timestep
    subroutine PrepRemixData(Model,Grid,State)
        type(Model_T) ,intent(in) :: Model
        type(Grid_T)  ,intent(in) :: Grid
        type(State_T) ,intent(in) :: State

        real(rp) ::  B0mag,Bi2m
        integer :: i,j,k,iG
        real(rp) :: xc,yc,zc
        
        real(rp) :: Con(NVAR)
        real(rp) :: Cs

        
        !write(*,*) 'Prepping remix data at T = ',Model%t

        !Only working on Bxyz from perturbation 
        !B0 in inner region (where we care for remix)
        !Assumd to be current-free
        call GetShellJ(Model,Grid,State%Bxyz,gJ)

        !Now loop over only shells and populate mhdvars to be sent to mix
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,iG,j,k,B0mag,Bi2m,xc,yc,zc,Con,Cs)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do i=1,JpSh
                    iG = JpSt+i-1
                    B0mag = sqrt(dot_product(Grid%B0(iG,j,k,:),Grid%B0(iG,j,k,:)))
                    ! get cell centers
                    call cellCenter(Grid,iG,j,k,xc,yc,zc)
                    !Get ion2mag scaling factor for field
                    Bi2m = BIon2Mag(xc,yc,zc)

                    if (k<=Grid%ke/2) then
                        !!! NOTE: assuming gB0 in nT and gx0 in m and gv0 in m/s

                        ! note conversion to microA/m^2
                        mhdvarsout(i,j,k,MHDJ,NORTH) = dot_product(gJ(i,j,k,:),Grid%B0(iG,j,k,:)/B0mag)*Bi2m*(gB0/gx0*1.e4/4/PI)
                        ! note conversion to g/cm^3
                        mhdvarsout(i,j,k,MHDD,NORTH) = State%Gas(iG,j,k,DEN,BLK)*(gB0*1.e-7/gv0)**2/4/pi
                        ! get sound speed first
                        Con = State%Gas(iG,j,k,:,BLK)
                        call CellPress2Cs(Model,Con,Cs)
                        mhdvarsout(i,j,k,MHDC,NORTH) = Cs*gv0*1.e2
                    else
                        mhdvarsout(i,j,k-Grid%ke/2,MHDJ,SOUTH) = dot_product(gJ(i,j,k,:),Grid%B0(iG,j,k,:)/B0mag)*Bi2m*(gB0/gx0*1.e4/4/PI)
                        mhdvarsout(i,j,k-Grid%ke/2,MHDD,SOUTH) = State%Gas(iG,j,k,DEN,BLK)*(gB0*1.e-7/gv0)**2/4/pi
                        ! get sound speed first
                        Con = State%Gas(iG,j,k,:,BLK)
                        call CellPress2Cs(Model,Con,Cs)
                        mhdvarsout(i,j,k-Grid%ke/2,MHDC,SOUTH) = Cs*gv0*1.e2
                    endif
                enddo
            enddo
        enddo
        
    end subroutine PrepRemixData

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

end module mhd_mix_interface
