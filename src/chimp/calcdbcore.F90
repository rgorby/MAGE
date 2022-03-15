module calcdbcore
	use chmpdefs
	use chmpunits
    use ebtypes
    use calcdbutils
    use clocks
    use geopack
    
	implicit none

    !Bounds for supermag indices [geomagnetic latitude]
    real(rp), parameter, private :: SMLowLat = +40.0
    real(rp), parameter, private :: SMHiLat  = +80.0
    real(rp), parameter, private :: SMRLat   =  50.0

	contains

    !Calculate contribution from BSGrid to ground
    subroutine BS2Gr(Model,t,ebState,magBS,ionBS,facBS,gGr)
        type(chmpModel_T), intent(in) :: Model
        real(rp)         , intent(in) :: t
        type(ebState_T)  , intent(in) :: ebState
        type(BSGrid_T), intent(inout) :: magBS,ionBS,facBS
        type(grGrid_T), intent(inout) :: gGr

        type(BSGrid_T) :: magltBS !Squashed magBS
        real(rp) :: mjd

        mjd = ioTabMJD(ebState%ebTab,t)
        call MJDRecalc(mjd) !Setup geopack for this time

        !Remove far away points in magnetospheric grid and make remaining contiguous
        if (magBS%isActive) then
            call Tic("Compactify")
            call Compactify(magBS,gGr%rMax,magltBS)
            call Toc("Compactify")

            call Tic("BSMag")
            call BSIntegral(magltBS,gGr,gGr%dbMAG_xyz)
            call Toc("BSMag")
        else
            gGr%dbMAG_xyz = 0.0
        endif

        if (ionBS%isActive) then
            call Tic("BSIon")
            call BSIntegral(ionBS,gGr,gGr%dbION_xyz)
            call Toc("BSIon")
        else
            gGr%dbION_xyz = 0.0
        endif

        if (facBS%isActive) then
            call Tic("BSFac")
            call BSIntegral(facBS,gGr,gGr%dbFAC_xyz)
            call Toc("BSFac")
        else
            gGr%dbFAC_xyz = 0.0
        endif

    !We've done all the work to get dB-XYZ (SM), now get some SM coordinates/indices

        !Calculate magnetic coordinates/northward deflection
        call CalcMagCoordinates(Model,t,gGr)

        !Calculate synthetic supermag indices
        call CalcSuperMAGIndices(Model,t,gGr)

    !Map dB-XYZ (SM) to dB-XYZ (GEO) if desired
        if (gGr%doGEO) then
            call Tic("BSRemap")
            call BSRemap(gGr,gGr%dbMAG_xyz)
            call BSRemap(gGr,gGr%dbION_xyz)
            call BSRemap(gGr,gGr%dbFAC_xyz)

            call Toc("BSRemap")
        endif
    end subroutine BS2Gr

!------
    subroutine CalcMagCoordinates(Model,t,gGr)
        type(chmpModel_T), intent(in) :: Model
        real(rp)         , intent(in) :: t
        type(grGrid_T), intent(inout) :: gGr

        real(rp), dimension(NDIM) :: x0,nhat,dbXYZ
        real(rp) :: rad,theta,phi
        integer  :: i,j,k

        !Calculate smlat/lon and dBn

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,x0,rad,theta,phi,nhat,dbXYZ)
        do k=1,gGr%Nz
            do j=1,gGr%NLon
                do i=1,gGr%NLat
                    x0 = gGr%SMxyzC(i,j,k,XDIR:ZDIR) !Cell center of ground grid
                    rad = norm2(x0)
                    gGr%smlat(i,j,k) = asin(x0(ZDIR)/rad)       *180.0/PI !geomagnetic latitude
                    gGr%smlon(i,j,k) = katan2(x0(YDIR),x0(XDIR))*180.0/PI !geomagnetic longitude
                    theta = acos (x0(ZDIR)/rad)
                    phi   = atan2(x0(YDIR),x0(XDIR))
                    nhat = -[cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)] !- thetahat from SM
                    dbXYZ = gGr%dbMAG_xyz(i,j,k,:) + gGr%dbION_xyz(i,j,k,:) + gGr%dbFAC_xyz(i,j,k,:)
                    gGr%dBn(i,j,k) = dot_product(dBxyz,nhat)
                enddo
            enddo
        enddo

    end subroutine CalcMagCoordinates

    !Remove far away points and make contiguous
    subroutine Compactify(xBS,rMax,xsBS)
        type(BSGrid_T), intent(in)    :: xBS
        real(rp)      , intent(in)    :: rMax
        type(BSGrid_T), intent(inout) :: xsBS

        logical, dimension(:), allocatable :: isG
        integer :: n,np

        !Identify good points
        allocate(isG(xBS%NumP))
        isG = .false.

        do n=1,xBS%NumP
            isG(n) = ( norm2(xBS%XYZcc(n,:)) <= rMax )
        enddo

        !Create sub-grid for only good points
        xsBS%NumP = count(isG)
        call BSSubInit(xsBS,xsBS%NumP)

        !Scrape values into new sub-grid
        xsBS%jScl = xBS%jScl

        np = 1
        do n=1,xBS%NumP !Loop over original grid
            if (isG(n)) then
                !Good point, copy this
                xsBS%XYZcc(np,:) = xBS%XYZcc(n,:)
                xsBS%Jxyz (np,:) = xBS%Jxyz (n,:)
                xsBS%dV   (np)   = xBS%dV   (n)
                np = np+1
            endif
        enddo

    end subroutine Compactify

    !Remap ground grid vectors to GEO
    subroutine BSRemap(gGr,dbXYZ)
        type(grGrid_T), intent(in) :: gGr
        real(rp), intent(inout) :: dbXYZ(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)

        integer :: i,j,k
        real(rp), dimension(NDIM) :: sm,geo

        if (.not. gGr%doGEO) return
        
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,sm,geo)
        do k=1,gGr%Nz
            do j=1,gGr%NLon
                do i=1,gGr%NLat
                    sm = dbXYZ(i,j,k,:) !SM ground dB
                    call SM2GEO(sm,geo)
                    dbXYZ(i,j,k,:) = geo
                enddo
            enddo
        enddo !k

    end subroutine BSRemap
    
    !Do individual BS integral using BSGr
    subroutine BSIntegral(xBS,gGr,dbXYZ)
        type(BSGrid_T), intent(in) :: xBS
        type(grGrid_T), intent(in) :: gGr
        real(rp), intent(inout) :: dbXYZ(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)
        
        integer :: iG,jG,kG,nS
        real(rp), dimension(NDIM) :: x0,xCC,ddB,J,R
        real(rp) :: dV,r3

        dbXYZ = 0.0
        
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(iG,jG,kG,nS) &
        !$OMP private(dV,r3,x0,xCC,ddB,J,R) 
        do kG=1,gGr%Nz
            do jG=1,gGr%NLon
                do iG=1,gGr%NLat
                    x0 = gGr%SMxyzC(iG,jG,kG,:) !Cell center of ground grid

                    do nS=1,xBS%NumP
                        xCC = xBS%XYZcc(nS,:) !Location of source contribution
                        R = x0-xCC !Vector pointing from source to destination/station
                        
                        r3 = norm2(R)**3.0
                        dV = xBS%dV(nS)
                        J = xBS%Jxyz(nS,:) !Current contribution

                        !Avoid array temporary
                        !NOTE: This current vector is in SM coordinates
                        ddB(XDIR) = ( J(YDIR)*R(ZDIR) - J(ZDIR)*R(YDIR) )/r3
                        ddB(YDIR) = ( J(ZDIR)*R(XDIR) - J(XDIR)*R(ZDIR) )/r3
                        ddB(ZDIR) = ( J(XDIR)*R(YDIR) - J(YDIR)*R(XDIR) )/r3

                        !dbXYZ(iG,jG,kG,:) = dbXYZ(iG,jG,kG,:) + xBS%jScl*dV*ddB
                        dbXYZ(iG,jG,kG,:) = dbXYZ(iG,jG,kG,:) + dV*ddB !Pull out overall scaling
                    enddo !nS
                    !Do scaling factor here just once
                    dbXYZ(iG,jG,kG,:) = xBS%jScl*dbXYZ(iG,jG,kG,:)

                enddo !iG
            enddo !jG
        enddo !kG

    end subroutine BSIntegral

    !Calculate supermag indices
    subroutine CalcSuperMAGIndices(Model,t,gGr)
        type(chmpModel_T), intent(in) :: Model
        real(rp)         , intent(in) :: t
        type(grGrid_T), intent(inout) :: gGr

        integer :: i,j,k
        integer, dimension(2) :: ijC
        real(rp), dimension(:,:,:), allocatable :: Bncorr
        real(rp), dimension(:,:), allocatable :: mlatIJ,mlonIJ,BnIJ,BncIJ
        logical , dimension(:,:), allocatable :: I_00,I_06,I_12,I_18
        logical , dimension(:,:), allocatable :: I_UL,I_R,IRLT

    !Allocate arrays
        allocate(Bncorr(gGr%NLat,gGr%NLon,gGr%Nz))
        !Sliced down arrays
        allocate(mlatIJ  (gGr%NLat,gGr%NLon))
        allocate(mlonIJ  (gGr%NLat,gGr%NLon))
        allocate(BnIJ    (gGr%NLat,gGr%NLon))
        allocate(BncIJ   (gGr%NLat,gGr%NLon))

        !Logical masks for different regions (lazy, but meh)
        allocate(I_UL(gGr%NLat,gGr%NLon))
        allocate(I_R (gGr%NLat,gGr%NLon))
        allocate(I_00(gGr%NLat,gGr%NLon))
        allocate(I_06(gGr%NLat,gGr%NLon))
        allocate(I_12(gGr%NLat,gGr%NLon))
        allocate(I_18(gGr%NLat,gGr%NLon))

        allocate(IRLT(gGr%NLat,gGr%NLon))

        !-Calculate corrected northward
        Bncorr = gGr%dBn/cos(gGr%smlat*PI/180.0)

    !Slice down to specific level
        k = 1 !Do calculation at lowest level

        mlatIJ = gGr%smlat(:,:,k)
        mlonIJ = gGr%smlon(:,:,k)
        BnIJ   = gGr%dBn  (:,:,k)
        BncIJ  = Bncorr   (:,:,k)

        !Create mask arrays
        I_UL = (mlatIJ <= SMHiLat) .and. (mlatIJ >= SMLowLat) !Treat the southern hemisphere like garbage
        I_R  = (mlatIJ <= +SMRLat) .and. (mlatIJ >=  -SMRLat)

        I_00 = (mlonIJ >= 135) .and. (mlonIJ <= 225)
        I_06 = (mlonIJ >= 225) .and. (mlonIJ <= 315)
        I_12 = (mlonIJ >= 315) .or.  (mlonIJ <=  45) !Straddling mlon=0
        I_18 = (mlonIJ >=  45) .and. (mlonIJ <= 135)

    !Get indices
    !AL/AU
        ijC = minloc(BnIJ,mask=I_UL)
        gGr%SML      = BnIJ  (ijC(1),ijC(2))
        gGr%SML_MLat = mlatIJ(ijC(1),ijC(2))
        gGr%SML_MLon = mlonIJ(ijC(1),ijC(2))
        
        ijC = maxloc(BnIJ,mask=I_UL)
        gGr%SMU      = BnIJ  (ijC(1),ijC(2))
        gGr%SMU_MLat = mlatIJ(ijC(1),ijC(2))
        gGr%SMU_MLon = mlonIJ(ijC(1),ijC(2))

    !AL/U-LTs (use min/max of raw Bn)
        !NOTE: These are technically quadrant values so not precisely the supermag octant SML/SMU-LTs
        IRLT = I_UL .and. I_00
        gGr%SML_00 = minval(BnIJ,mask=IRLT)
        gGr%SMU_00 = maxval(BnIJ,mask=IRLT)

        IRLT = I_UL .and. I_06
        gGr%SML_06 = minval(BnIJ,mask=IRLT)
        gGr%SMU_06 = maxval(BnIJ,mask=IRLT)

        IRLT = I_UL .and. I_12
        gGr%SML_12 = minval(BnIJ,mask=IRLT)
        gGr%SMU_12 = maxval(BnIJ,mask=IRLT)

        IRLT = I_UL .and. I_18
        gGr%SML_18 = minval(BnIJ,mask=IRLT)
        gGr%SMU_18 = maxval(BnIJ,mask=IRLT)

    !SMR-LTs (use avg lat-corrected Bn)
        IRLT = I_R .and. I_00
        gGr%SMR_00 = sum(BncIJ,mask=IRLT)/count(IRLT)

        IRLT = I_R .and. I_06
        gGr%SMR_06 = sum(BncIJ,mask=IRLT)/count(IRLT)

        IRLT = I_R .and. I_12
        gGr%SMR_12 = sum(BncIJ,mask=IRLT)/count(IRLT)

        IRLT = I_R .and. I_18
        gGr%SMR_18 = sum(BncIJ,mask=IRLT)/count(IRLT)

    !Calculate derived indices
        gGr%SME =  gGr%SMU - gGr%SML
        gGr%SMO = (gGr%SMU + gGr%SML)/2
        gGr%SMR = 0.25*(gGr%SMR_00+gGr%SMR_06+gGr%SMR_12+gGr%SMR_18)

    end subroutine CalcSuperMAGIndices

end module calcdbcore
