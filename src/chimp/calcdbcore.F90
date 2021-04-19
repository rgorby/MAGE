module calcdbcore
	use chmpdefs
	use chmpunits
    use ebtypes
    use calcdbutils
    use clocks
    use geopack
    use ebtabutils
    
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

        !Remove far away points in magnetospheric grid and make remaining contiguous
        call Tic("Compactify")
        call Compactify(magBS,gGr%rMax,magltBS)
        call Toc("Compactify")

        call Tic("BSMag")
        call BSIntegral(magltBS,gGr,gGr%dbMAG_xyz)
        call Toc("BSMag")

        call Tic("BSIon")
        call BSIntegral(ionBS,gGr,gGr%dbION_xyz)
        call Toc("BSIon")

        call Tic("BSFac")
        call BSIntegral(facBS,gGr,gGr%dbFAC_xyz)
        call Toc("BSFac")

        !We've done all the work to get dB-XYZ (SM)
        !Before doing anything else calculate auroral indices
        call CalcAuroralIndices(Model,t,gGr)

        !Map dB-XYZ (SM) to dB-XYZ (GEO) if desired
        if (gGr%doGEO) then
            mjd = MJDAt(ebState%ebTab,t)
            call MJDRecalc(mjd) !Setup geopack for this time
            
            call Tic("BSRemap")
            call BSRemap(gGr,gGr%dbMAG_xyz)
            call BSRemap(gGr,gGr%dbION_xyz)
            call BSRemap(gGr,gGr%dbFAC_xyz)

            call Toc("BSRemap")
        endif
    end subroutine BS2Gr

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
        
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,sm,geo)
        do k=1,gGr%Nz
            do j=1,gGr%NLon
                do i=1,gGr%NLat
                    sm = dbXYZ(i,j,k,:) !SM ground dB
                    call SM2GEO(sm(XDIR),sm(YDIR),sm(ZDIR),geo(XDIR),geo(YDIR),geo(ZDIR))
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

    !Calculate auroral indices
    subroutine CalcAuroralIndices(Model,t,gGr)
        type(chmpModel_T), intent(in) :: Model
        real(rp)         , intent(in) :: t
        type(grGrid_T), intent(inout) :: gGr

        integer :: i,j,k
        integer, dimension(2) :: ijC
        real(rp), dimension(NDIM) :: x0,dBxyz,dBrtp
        real(rp) :: Bn
        real(rp), dimension(:,:), allocatable :: mlatIJ,mlonIJ,BnIJ,BncorrIJ
        logical , dimension(:,:), allocatable :: I_UL,I_R,I_00,I_06,I_12,I_18,IRLT

        k = 1 !Do calculation at lowest level

        !Allocate arrays
        allocate(mlatIJ  (gGr%NLat,gGr%NLon))
        allocate(mlonIJ  (gGr%NLat,gGr%NLon))
        allocate(BnIJ    (gGr%NLat,gGr%NLon))
        allocate(BncorrIJ(gGr%NLat,gGr%NLon))

        !Logical masks for different regions
        allocate(I_UL(gGr%NLat,gGr%NLon))
        allocate(I_R (gGr%NLat,gGr%NLon))
        allocate(I_00(gGr%NLat,gGr%NLon))
        allocate(I_06(gGr%NLat,gGr%NLon))
        allocate(I_12(gGr%NLat,gGr%NLon))
        allocate(I_18(gGr%NLat,gGr%NLon))
        allocate(IRLT(gGr%NLat,gGr%NLon))

        !Get array values (to easily do minloc/maxloc)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,x0,dBxyz,dBrtp,Bn)
        do j=1,gGr%NLon
            do i=1,gGr%NLat
                x0 = gGr%SMxyzC(i,j,k,:) !Cell center of ground grid
                mlatIJ(i,j) = asin(x0(ZDIR)/norm2(x0)) *180.0/PI !geomagnetic latitude
                mlonIJ(i,j) = katan2(x0(YDIR),x0(XDIR))*180.0/PI !geomagnetic longitude
                !Get total SM-XYZ deflection
                dBxyz = gGr%dbMAG_xyz(i,j,k,:) + gGr%dbION_xyz(i,j,k,:) + gGr%dbFAC_xyz(i,j,k,:)
                !Convert to spherical coordinates
                dBrtp = xyz2rtp(x0,dBxyz)
                !Bn = -dB_theta-SM
                Bn = -dBrtp(2) !Deflection in direction of geomagnetic north
                BnIJ(i,j) = Bn
                BncorrIJ(i,j) = Bn/cos(mlatIJ(i,j))

            enddo
        enddo

        !Create mask arrays
        I_UL = (mlatIJ <= SMHiLat) .and. (mlatIJ >= SMLowLat)
        I_R  = (mlatIJ <= +SMRLat) .and. (mlatIJ >=  -SMRLat)

        I_00 = (mlonIJ >= 135) .and. (mlonIJ <= 225)
        I_06 = (mlonIJ >= 225) .and. (mlonIJ <= 315)
        I_12 = (mlonIJ >= 315) .or.  (mlonIJ <=  45) !Straddling mlon=0
        I_18 = (mlonIJ >=  45) .and. (mlonIJ <= 135)

    !Get indices
        !AL
        ijC = minloc(BnIJ,mask=I_UL)
        gGr%SML      = BnIJ  (ijC(1),ijC(2))
        gGr%SML_MLat = mlatIJ(ijC(1),ijC(2))
        gGr%SML_MLon = mlonIJ(ijC(1),ijC(2))
        !AU
        ijC = maxloc(BnIJ,mask=I_UL)
        gGr%SMU      = BnIJ  (ijC(1),ijC(2))
        gGr%SMU_MLat = mlatIJ(ijC(1),ijC(2))
        gGr%SMU_MLon = mlonIJ(ijC(1),ijC(2))

        !SMR-LTs
        IRLT = I_R .and. I_00
        gGr%SMR_00 = sum(BnIJ,mask=IRLT)/count(IRLT)

        IRLT = I_R .and. I_06
        gGr%SMR_06 = sum(BnIJ,mask=IRLT)/count(IRLT)

        IRLT = I_R .and. I_12
        gGr%SMR_12 = sum(BnIJ,mask=IRLT)/count(IRLT)

        IRLT = I_R .and. I_18
        gGr%SMR_18 = sum(BnIJ,mask=IRLT)/count(IRLT)

    !Calculate derived indices
        gGr%SME =  gGr%SMU - gGr%SML
        gGr%SMO = (gGr%SMU + gGr%SML)/2
        gGr%SMR = 0.25*(gGr%SMR_00+gGr%SMR_06+gGr%SMR_12+gGr%SMR_18)

    end subroutine CalcAuroralIndices

end module calcdbcore
