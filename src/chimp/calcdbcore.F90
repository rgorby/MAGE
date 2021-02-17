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
        type(BSGrid_T), intent(in) :: magBS,ionBS,facBS
        type(grGrid_T), intent(inout) :: gGr

        call Tic("BSMag")
        call BSIntegral(magBS,gGr,gGr%dbMAG_xyz)
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
                        if ( norm2(xCC) > gGr%rMax ) cycle
                        
                        R = xCC-x0 !R = x_src - x_station
                        r3 = norm2(R)**3.0
                        dV = xBS%dV(nS)
                        J = xBS%Jxyz(nS,:) !Current contribution

                        !Avoid array temporary
                        ddB(XDIR) = ( J(YDIR)*R(ZDIR) - J(ZDIR)*R(YDIR) )/r3
                        ddB(YDIR) = ( J(ZDIR)*R(XDIR) - J(XDIR)*R(ZDIR) )/r3
                        ddB(ZDIR) = ( J(XDIR)*R(YDIR) - J(YDIR)*R(XDIR) )/r3

                        !NOTE: This current vector is in SM coordinates

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
