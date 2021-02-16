module calcdbcore
	use chmpdefs
	use chmpunits
    use ebtypes
    use calcdbutils
    use clocks
    
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

    end subroutine BS2Gr

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
                    x0 = gGr%xyzC(iG,jG,kG,:) !Cell center of ground grid

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

                        dbXYZ(iG,jG,kG,:) = dbXYZ(iG,jG,kG,:) + xBS%jScl*dV*ddB
                    enddo !nS

                enddo !iG
            enddo !jG
        enddo !kG

    end subroutine BSIntegral


end module calcdbcore
