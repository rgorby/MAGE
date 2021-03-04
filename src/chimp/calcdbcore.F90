module calcdbcore
	use chmpdefs
	use chmpunits
    use ebtypes
    use calcdbutils
    use clocks
    use geopack

	implicit none

	contains

    !Calculate contribution from BSGrid to ground
    subroutine BS2Gr(Model,magBS,ionBS,facBS,gGr)
        type(chmpModel_T), intent(in) :: Model
        type(BSGrid_T), intent(inout) :: magBS,ionBS,facBS
        type(grGrid_T), intent(inout) :: gGr

        type(BSGrid_T) :: magltBS !Squashed magBS

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

        if (gGr%doGEO) then
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
                        
                        R = xCC-x0 !R = x_src - x_station
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


end module calcdbcore
