! module to convert imag data to mhd data
module imag2mhd_interface
    use volttypes
    use cmiutils
    use planethelper
    
    implicit none

    contains

    ! have this short interface function in case we want to alter how we fill source data, or use something other than chimp and imagApp
    subroutine CoupleSourceToMhd(vApp)
        class(voltApp_T), intent(inout) :: vApp

        call convertImagToGamera(vApp%gApp, vApp)

    end subroutine

    subroutine convertImagToGamera(gApp, vApp)
        class(gamCoupler_T), intent(inout) :: gApp
        type(voltApp_T), intent(inout) :: vApp

        integer :: i,j,k,Nk
        real(rp) :: x1,x2,t
        real(rp) :: imW(NVARIMAG),Qs(8),xyz(NDIM)
        logical :: isTasty
        !NOTE: Using isgTasy to handle isg:is region (not covered by chimp)
        logical , dimension(:,:,:), allocatable :: isgTasty !Embiggen-ed logical array

    !Proceed in two steps
    ! 1) Get ingestion values at each node (cell corner)
    ! 2) Loop over cells and average from corners to cell centers

        !TODO: Think about what time to evaluate at
        t = gApp%Model%t*gApp%Model%Units%gT0

        if(size(gApp%SrcNC,1) .ne. (vApp%chmp2mhd%iMax+1)) then
            deallocate(gApp%SrcNC)
            !NOTE: Embiggening SrcNC to include isg:is region
            allocate(gApp%SrcNC(gApp%Grid%isg:vApp%chmp2mhd%iMax+1,gApp%Grid%js:gApp%Grid%je+1,gApp%Grid%ks:gApp%Grid%ke+1,1:NVARIMAG))
        endif

        associate(Gr=>gApp%Grid,chmp2mhd=>vApp%chmp2mhd,SrcNC=>gApp%SrcNC)

        allocate(isgTasty(Gr%isg:chmp2mhd%iMax+1,Gr%js:Gr%je+1,Gr%ks:Gr%ke+1))
        !Create local storage for cell corner imW's
        SrcNC = 0.0
        isgTasty = .false.
        chmp2mhd%isEdible = .false.

    ! 1) Cell corner ingestion
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,x1,x2,imW,isTasty)
        do k=Gr%ks,Gr%ke+1
            do j=Gr%js,Gr%je+1
                do i=Gr%is,chmp2mhd%iMax+1

                    if (chmp2mhd%isGood(i,j,k)) then
                        !Good projection, let's get some values
                        x1 = PI/2 - chmp2mhd%xyzSquish(i,j,k,1)  ! Colat
                        x2 = chmp2mhd%xyzSquish(i,j,k,2)
                        call vApp%imagApp%getMoments(x1,x2,t,imW,isTasty)
                    else
                        !Projection wasn't good, nothing good to eat
                        imW = 0.0
                        isTasty = .false.

                    endif !isGood
                    SrcNC(i,j,k,:) = imW
                    chmp2mhd%isEdible(i,j,k) = isTasty
                    isgTasty(i,j,k) = isTasty !Save in embiggen-ed array
                enddo !i loop
            enddo
        enddo

        !Embiggen to include i-ghosts for certain cases
        if ( (vApp%isEarth) .and. (vApp%prType == LLPROJ) ) then
            !Do lat-lon (dipole) projections for gamera ghost cells
            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,k,x1,x2,xyz,imW,isTasty)
            do k=Gr%ks,Gr%ke+1
                do j=Gr%js,Gr%je+1
                    do i=Gr%isg,Gr%is-1
                        !Dipole project
                        xyz = Gr%xyz(i,j,k,:) !Gamera grid corner
                        x1 = PI/2 - InvLatitude(xyz)  ! Colat
                        x2 = katan2(xyz(YDIR),xyz(XDIR)) !katan => [0,2pi] instead of [-pi,pi]
                        call vApp%imagApp%getMoments(x1,x2,t,imW,isTasty)
                        SrcNC(i,j,k,:) = imW
                        isgTasty(i,j,k) = isTasty
                    enddo
                enddo
            enddo
        endif


    ! 2) Corners => Centers
        Gr%Gas0 = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k,imW,Qs)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%isg,chmp2mhd%iMax


                    if ( all(isgTasty(i:i+1,j:j+1,k:k+1)) ) then
                    !Density and pressure
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMDEN),Qs)
                        imW(IMDEN) = ArithMean(Qs)
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMPR) ,Qs)
                        imW(IMPR)  = ArithMean(Qs)
                    !x1 and x2
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMX1),Qs)
                        imW(IMX1) = ArithMean(Qs)
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMX2),Qs)
                        imW(IMX2) = CircMeanDeg(Qs)

                    !Timescale
                        call SquishCorners(SrcNC(i:i+1,j:j+1,k:k+1,IMTSCL),Qs)
                        if ( all(Qs>TINY) ) then
                            imW(IMTSCL) = ArithMean(Qs)
                        else if (any(Qs>TINY)) then
                            imW(IMTSCL) = sum(Qs,mask=(Qs>TINY))/count(Qs>TINY)
                        else
                            imW(IMTSCL) = vApp%DeepDT
                        endif
                    else
                        !Not good to eat
                        imW = 0.0
                    endif

                    imW(IMTSCL) = max(imW(IMTSCL),vApp%DeepDT)

                    !Do scaling and store
                    !density/pressure coming back in #/cc and nPa
                    !ingestion timescale coming back in seconds
                    Gr%Gas0(i,j,k,IMDEN ) = imW(IMDEN)
                    Gr%Gas0(i,j,k,IMPR  ) = imW(IMPR)/gApp%Model%Units%gP0
                    Gr%Gas0(i,j,k,IMX1  ) = imW(IMX1)
                    Gr%Gas0(i,j,k,IMX2  ) = imW(IMX2)
                    Gr%Gas0(i,j,k,IMTSCL) = imW(IMTSCL)/gApp%Model%Units%gT0

                enddo !i loop
            enddo
        enddo

    !Now do some touch up at the axis and get outta here

        !Do averaging for first cell next to singularity
        !Do for +/- X pole and density/pressure
        Nk = Gr%ke-Gr%ks+1
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,imW)
        do i=Gr%is,chmp2mhd%iMax
            !+X pole
            imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN),Nk)
            imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ),Nk)
            Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN) = imW(IMDEN)
            Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ) = imW(IMPR )

            !-X pole
            imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN),Nk)
            imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ),Nk)
            Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMDEN) = imW(IMDEN)
            Gr%Gas0(i,Gr%js,Gr%ks:Gr%ke,IMPR ) = imW(IMPR )

            !-X pole
            imW(IMDEN) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN),Nk)
            imW(IMPR ) = AvgOverGood(Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMPR ),Nk)
            Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMDEN) = imW(IMDEN)
            Gr%Gas0(i,Gr%je,Gr%ks:Gr%ke,IMPR ) = imW(IMPR )
        enddo

        end associate

        contains
            function AvgOverGood(Q,Nk) result(Qavg)
                real(rp), intent(in), dimension(Nk) :: Q
                integer , intent(in) :: Nk

                real(rp) :: Qavg
                integer :: Nkg

                if ( any(Q>TINY) ) then
                    Nkg = count(Q>TINY)
                    Qavg = sum(Q,mask=(Q>TINY))/Nkg
                else
                    Qavg = 0.0
                endif

            end function AvgOverGood

    end subroutine

end module

