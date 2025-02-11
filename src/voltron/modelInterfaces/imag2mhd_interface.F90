! module to convert imag data to mhd data
module imag2mhd_interface
    use volttypes
    use cmiutils
    use planethelper
    use gridutils
    use shellInterp

    implicit none

    contains

    ! have this short interface function in case we want to alter how we fill source data, or use something other than chimp and imagApp
    subroutine CoupleSourceToMhd(vApp)
        class(voltApp_T), intent(inout) :: vApp

        real(rp), dimension(:,:,:)  , allocatable :: gPsi
        real(rp), dimension(:,:,:,:), allocatable :: Eijk
        logical , dimension(:,:,:)  , allocatable :: isGoodCC

        integer :: n,i,j,k
        real(rp) :: x1,x2,x1c,di,dj,dk,Rion,t,eScl
        real(rp), dimension(NDIM) :: xyz
        real(rp) :: Qs(8) !Squish all the corners of a single cell
        real(rp) :: imW(IM_D_RING:IM_TSCL) !All imag variables
        logical :: isTasty,isG


        !Add projection to ion instead of invlat
        Rion = vApp%shGrid%radius !Voltron grid radius in Rx
        eScl = vApp%gApp%rm2g !Scaling for remix kV

        associate(Gr=>vApp%gApp%Grid,gModel=>vApp%gApp%Model)
            Gr%Gas0 = 0.0

            allocate(gPsi(Gr%isg:Gr%ieg+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1)) !All nodes
            gPsi = 0.0

            !Loop over "real" nodes and get total potential from voltron grid
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,x1,x2,x1c,xyz,isG)
            do k=Gr%ks,Gr%ke+1
                do j=Gr%js,Gr%je+1
                    do i=Gr%isg,Gr%ie+1

                        if (i < Gr%is) then
                            !Use dipole projection
                            xyz = Gr%xyz(i,j,k,:) !Gamera grid corner
                            call Proj2Rad(xyz,Rion,x1,x2)
                            isG = .true.
                        else
                            x1  = vApp%chmp2mhd%xyzSquish(i,j,k,1)
                            x2  = vApp%chmp2mhd%xyzSquish(i,j,k,2)
                            isG = vApp%chmp2mhd%isGood(i,j,k)
                        endif
                        x1c = PI/2 - x1
                        if (isG) then
                            !Get potential in kV
                            call InterpShellVar_TSC_pnt(vApp%shGrid,vApp%State%potential_total,&
                                                        x1c,x2,gPsi(i,j,k) )
                        endif
                    enddo
                enddo !j
            enddo !k

            !Create edge (ijk) E and then convert to cell-centered Exyz
            allocate(Eijk(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM))
            Eijk = 0.0

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,di,dj,dk)
            do k=Gr%ksg,Gr%keg
                do j=Gr%jsg,Gr%jeg
                    do i=Gr%isg,Gr%ieg
                        di = Gr%edge(i,j,k,IDIR)
                        dj = Gr%edge(i,j,k,JDIR)
                        dk = Gr%edge(i,j,k,KDIR)
                        Eijk(i,j,k,IDIR) = -( gPsi(i+1,j,k) - gPsi(i,j,k) )/di
                        Eijk(i,j,k,JDIR) = -( gPsi(i,j+1,k) - gPsi(i,j,k) )/dj
                        Eijk(i,j,k,KDIR) = -( gPsi(i,j,k+1) - gPsi(i,j,k) )/dk

                        Eijk(i,j,k,:) = Eijk(i,j,k,:)/eScl
                        
                    enddo
                enddo
            enddo !k

            !Convert edge E fields to xyz and store
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k)
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%isg,Gr%ieg
                        Gr%Gas0(i,j,k,IONEX:IONEZ) = CellExyz(gModel,Gr,Eijk,i,j,k)
                    enddo
                enddo !j
            enddo !k

            do n=IONEX,IONEZ
                call FillGhostsCC(gModel,Gr,Gr%Gas0(:,:,:,n))
            enddo

            allocate(isGoodCC(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg))
            isGoodCC = .false.
            !Get cell-centered projections

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,xyz,x1,x2,Qs)
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%isg,Gr%ie
                        
                        if (i<Gr%is) then
                            !Below inner boundary, do dipole projection
                            isGoodCC(i,j,k) = .true.
                            xyz = Gr%xyzcc(i,j,k,:) !Gamera grid center
                            call Proj2Rad(xyz,Rion,x1,x2)
                            Gr%Gas0(i,j,k,PROJLAT) = x1
                            Gr%Gas0(i,j,k,PROJLON) = x2

                        else
                            !Get value from xyzsquish

                            if ( all(vApp%chmp2mhd%isGood(i:i+1,j:j+1,k:k+1)) ) then
                                !All values are good, so just do this thing
                                call SquishCorners(vApp%chmp2mhd%xyzSquish(i:i+1,j:j+1,k:k+1,1),Qs)
                                Gr%Gas0(i,j,k,PROJLAT) = ArithMean(Qs)
                                call SquishCorners(vApp%chmp2mhd%xyzSquish(i:i+1,j:j+1,k:k+1,2),Qs)
                                Gr%Gas0(i,j,k,PROJLON) = CircMean(Qs)
                                isGoodCC(i,j,k) = .true.
                            else
                                Gr%Gas0(i,j,k,PROJLAT) = 0.0
                                Gr%Gas0(i,j,k,PROJLON) = 0.0
                                isGoodCC(i,j,k) = .false.
                            endif !corner projection => center
                        endif !inner-i vs. not

                    enddo
                enddo !j
            enddo !k

            call FillGhostsCC(gModel,Gr,Gr%Gas0(:,:,:,PROJLAT))
            call FillGhostsCC(gModel,Gr,Gr%Gas0(:,:,:,PROJLON))

            t = (vApp%gApp%Model%t)*(vApp%gApp%Model%Units%gT0)

            !Loop over and get imag data
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP private(i,j,k,x1c,x2,imW,isTasty)
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    do i=Gr%isg,Gr%ie
                        Gr%Gas0(i,j,k,IM_D_RING:IM_P_RING) = 0.0
                        Gr%Gas0(i,j,k,IM_D_COLD:IM_P_COLD) = 0.0
                        Gr%Gas0(i,j,k,IM_TSCL) = 0.0

                        if (isGoodCC(i,j,k)) then
                            !Get imag value here
                            imW = 0.0
                            x1c = PI/2 - Gr%Gas0(i,j,k,PROJLAT)
                            x2  = Gr%Gas0(i,j,k,PROJLON)

                            call vApp%imagApp%getMoments(x1c,x2,t,imW,isTasty)
                            if (isTasty) then
                                !Density/pressure coming back in #/cc and nPa
                                !Ingestion timescale coming back in seconds
                                Gr%Gas0(i,j,k,IM_D_RING) = imW(IM_D_RING)
                                Gr%Gas0(i,j,k,IM_D_COLD) = imW(IM_D_COLD)
                                Gr%Gas0(i,j,k,IM_P_RING) = imW(IM_P_RING)/vApp%gApp%Model%Units%gP0
                                Gr%Gas0(i,j,k,IM_P_COLD) = imW(IM_P_COLD)/vApp%gApp%Model%Units%gP0
                                Gr%Gas0(i,j,k,IM_TSCL  ) = imW(IM_TSCL  )/vApp%gApp%Model%Units%gT0
                            endif
                            
                        endif
                    enddo !i
                enddo
            enddo !k

            do n=IM_D_RING,IM_TSCL
                call FillGhostsCC(gModel,Gr,Gr%Gas0(:,:,:,n))
            enddo
        end associate



        !---
        contains

            !Take a full-sized cell-centered grid w/ active values defined and fill ghosts
            subroutine FillGhostsCC(Model,Gr,Q)
                type(Model_T), intent(in) :: Model
                type(Grid_T) , intent(in) :: Gr
                real(rp), intent(inout) :: Q(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg)

                integer :: i,j,k,ip,jp,kp
                logical :: isActive

                !$OMP PARALLEL DO default(shared) collapse(2) &
                !$OMP private(i,j,k,isActive,ip,jp,kp)
                do k=Gr%ksg,Gr%keg
                    do j=Gr%jsg,Gr%jeg
                        do i=Gr%isg,Gr%ieg
                            isActive = (j >= Gr%js) .and. (j <= Gr%je) .and. &
                                       (k >= Gr%ks) .and. (k <= Gr%ks)
                            if (isActive) cycle
                            !If still here map this ghost to active and set value based on active
                            call lfmIJKcc(Model,Gr,i,j,k,ip,jp,kp)
                            Q(i,j,k) = Q(ip,jp,kp)
                        enddo
                    enddo !j
                enddo !k

            end subroutine FillGhostsCC

            !Project xyz along dipole to R0 and return lat (x1) and lon (x2)
            subroutine Proj2Rad(xyz,R0,x1,x2)
                real(rp), intent(in ) :: xyz(NDIM), R0
                real(rp), intent(out) :: x1,x2

                real(rp), dimension(NDIM) :: xyz0

                xyz0 = DipoleShift(xyz,R0)
                x1 = asin(xyz0(ZDIR)/R0) !Lat
                x2 = katan2(xyz0(YDIR),xyz0(XDIR)) !katan => [0,2pi] instead of [-pi,pi]

            end subroutine Proj2Rad 

    end subroutine

    ! subroutine convertImagToGamera(gApp, vApp)
    !     class(gamCoupler_T), intent(inout) :: gApp
    !     type(voltApp_T), intent(inout) :: vApp

    !     integer :: i,j,k,Nk
    !     real(rp) :: x1,x2,t
    !     !real(rp) :: imW(NVARIMAG),Qs(8),xyz(NDIM)
    !     logical :: isTasty
    !     !NOTE: Using isgTasy to handle isg:is region (not covered by chimp)
    !     logical , dimension(:,:,:), allocatable :: isgTasty !Embiggen-ed logical array


    !     !TODO: Think about what time to evaluate at
    !     t = gApp%Model%t*gApp%Model%Units%gT0

    !     gApp%Grid%Gas0 = 0.0 !Just set it all to zero

    ! ! !Proceed in two steps
    ! ! ! 1) Get ingestion values at each node (cell corner)
    ! ! ! 2) Loop over cells and average from corners to cell centers

    ! !     if(size(gApp%SrcNC,1) .ne. (vApp%chmp2mhd%iMax+1)) then
    ! !         deallocate(gApp%SrcNC)
    ! !         !NOTE: Embiggening SrcNC to include isg:is region
    ! !         allocate(gApp%SrcNC(gApp%Grid%isg:vApp%chmp2mhd%iMax+1,gApp%Grid%js:gApp%Grid%je+1,gApp%Grid%ks:gApp%Grid%ke+1,1:NVARIMAG))
    ! !     endif

    ! !     associate(Gr=>gApp%Grid,chmp2mhd=>vApp%chmp2mhd,SrcNC=>gApp%SrcNC)

    ! !     allocate(isgTasty(Gr%isg:chmp2mhd%iMax+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1))
    ! !     !Create local storage for cell corner imW's
    ! !     SrcNC = 0.0
    ! !     isgTasty = .false.
    ! !     chmp2mhd%isEdible = .false.

    ! ! ! 1) Cell corner ingestion
    ! !     !$OMP PARALLEL DO default(shared) collapse(2) &
    ! !     !$OMP schedule(dynamic) &
    ! !     !$OMP private(i,j,k,x1,x2,imW,isTasty)
    ! !     do k=Gr%ks,Gr%ke+1
    ! !         do j=Gr%js,Gr%je+1
    ! !             do i=Gr%is,chmp2mhd%iMax+1

    ! !                 if (chmp2mhd%isGood(i,j,k)) then
    ! !                     !Good projection, let's get some values
    ! !                     x1 = PI/2 - chmp2mhd%xyzSquish(i,j,k,1)  ! Colat
    ! !                     x2 = chmp2mhd%xyzSquish(i,j,k,2)
    ! !                     call vApp%imagApp%getMoments(x1,x2,t,imW,isTasty)
    ! !                 else
    ! !                     !Projection wasn't good, nothing good to eat
    ! !                     imW = 0.0
    ! !                     isTasty = .false.

    ! !                 endif !isGood
    ! !                 SrcNC(i,j,k,:) = imW
    ! !                 chmp2mhd%isEdible(i,j,k) = isTasty
    ! !                 isgTasty(i,j,k) = isTasty !Save in embiggen-ed array
    ! !             enddo !i loop
    ! !         enddo
    ! !     enddo

    ! !     !Embiggen to include i-ghosts for certain cases
    ! !     if ( (vApp%isEarth) .and. (vApp%prType == LLPROJ) ) then
    ! !         !Do lat-lon (dipole) projections for gamera ghost cells
    ! !         !$OMP PARALLEL DO default(shared) collapse(2) &
    ! !         !$OMP schedule(dynamic) &
    ! !         !$OMP private(i,j,k,x1,x2,xyz,imW,isTasty)
    ! !         do k=Gr%ks,Gr%ke+1
    ! !             do j=Gr%js,Gr%je+1
    ! !                 do i=Gr%isg,Gr%is-1
    ! !                     !Dipole project
    ! !                     xyz = Gr%xyz(i,j,k,:) !Gamera grid corner
    ! !                     x1 = PI/2 - InvLatitude(xyz)  ! Colat
    ! !                     x2 = katan2(xyz(YDIR),xyz(XDIR)) !katan => [0,2pi] instead of [-pi,pi]
    ! !                     call vApp%imagApp%getMoments(x1,x2,t,imW,isTasty)
    ! !                     SrcNC(i,j,k,:) = imW
    ! !                     isgTasty(i,j,k) = isTasty
    ! !                 enddo
    ! !             enddo
    ! !         enddo
    ! !     endif


    ! ! ! 2) Corners => Centers
    ! !     Gr%Gas0 = 0.0

    ! !     !$OMP PARALLEL DO default(shared) collapse(2) &
    ! !     !$OMP schedule(dynamic) &
    ! !     !$OMP private(i,j,k,imW,Qs)
    ! !     do k=Gr%ks,Gr%ke
    ! !         do j=Gr%js,Gr%je
    ! !             do i=Gr%isg,chmp2mhd%iMax


    ! !                 if ( all(isgTasty(i:i+1,j:j+1,k:k+1)) ) then
    ! !                 !Density and pressure
    ! !                     call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMDEN),Qs)
    ! !                     imW(IMDEN) = ArithMean(Qs)
    ! !                     call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMPR) ,Qs)
    ! !                     imW(IMPR)  = ArithMean(Qs)
    ! !                 !x1 and x2
    ! !                     call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMX1),Qs)
    ! !                     imW(IMX1) = ArithMean(Qs)
    ! !                     call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMX2),Qs)
    ! !                     imW(IMX2) = CircMeanDeg(Qs)

    ! !                 !Timescale
    ! !                     call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMTSCL),Qs)
    ! !                     if ( all(Qs>TINY) ) then
    ! !                         imW(IMTSCL) = ArithMean(Qs)
    ! !                     else if (any(Qs>TINY)) then
    ! !                         imW(IMTSCL) = sum(Qs,mask=(Qs>TINY))/count(Qs>TINY)
    ! !                     else
    ! !                         imW(IMTSCL) = vApp%DeepDT
    ! !                     endif
    ! !                 else
    ! !                     !Not good to eat
    ! !                     imW = 0.0
    ! !                 endif

    ! !                 imW(IMTSCL) = max(imW(IMTSCL),vApp%DeepDT)

    ! !                 !Do scaling and store
    ! !                 !density/pressure coming back in #/cc and nPa
    ! !                 !ingestion timescale coming back in seconds
    ! !                 Gr%Gas0(i,j,k,IMDEN ) = imW(IMDEN)
    ! !                 Gr%Gas0(i,j,k,IMPR  ) = imW(IMPR)/gApp%Model%Units%gP0
    ! !                 Gr%Gas0(i,j,k,IMX1  ) = imW(IMX1)
    ! !                 Gr%Gas0(i,j,k,IMX2  ) = imW(IMX2)
    ! !                 Gr%Gas0(i,j,k,IMTSCL) = imW(IMTSCL)/gApp%Model%Units%gT0

    ! !             enddo !i loop
    ! !         enddo
    ! !     enddo

    ! ! !Now do some touch up at the axis and get outta here

    ! !     !Do averaging for first cell next to singularity
    ! !     !Do for +/- X pole and density/pressure
    ! !     Nk = Gr%ke-Gr%ks+1
    ! !     !$OMP PARALLEL DO default(shared) &
    ! !     !$OMP private(i,imW)
    ! !     do i=Gr%is,chmp2mhd%iMax
    ! !         !+X pole
    ! !         imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN),Nk)
    ! !         imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ),Nk)
    ! !         Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN) = imW(IMDEN)
    ! !         Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ) = imW(IMPR )

    ! !         !-X pole
    ! !         imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN),Nk)
    ! !         imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ),Nk)
    ! !         Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN) = imW(IMDEN)
    ! !         Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ) = imW(IMPR )

    ! !         !-X pole
    ! !         imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN),Nk)
    ! !         imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMPR ),Nk)
    ! !         Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN) = imW(IMDEN)
    ! !         Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMPR ) = imW(IMPR )
    ! !     enddo

    ! !     end associate

    !     contains
    !         function AvgOverGood(Q,Nk) result(Qavg)
    !             real(rp), intent(in), dimension(Nk) :: Q
    !             integer , intent(in) :: Nk

    !             real(rp) :: Qavg
    !             integer :: Nkg

    !             if ( any(Q>TINY) ) then
    !                 Nkg = count(Q>TINY)
    !                 Qavg = sum(Q,mask=(Q>TINY))/Nkg
    !             else
    !                 Qavg = 0.0
    !             endif

    !         end function AvgOverGood

    ! end subroutine

end module

