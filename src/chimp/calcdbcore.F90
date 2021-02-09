module calcdbcore
	use chmpdefs
	use chmpunits
    use ebtypes
    use calcdbutils
    use ebinterp
    
	implicit none

	contains

	!Calculate magnetospheric DB contribution
	subroutine CalcMagDB(Model,ebState,gGr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
		type(sphGrid_T), intent(inout) :: gGr

        real(rp), dimension(:,:,:,:), allocatable :: Jxyz
        real(rp) :: w1,w2,dV,B0,r3
        real(rp), dimension(NDIM) :: x0,xCC,ddB
        integer :: iG,jG,kG,iM,jM,kM
        
        B0 = bScale()
        call GetTWgts(Model,ebState,Model%t,w1,w2)

    	allocate(Jxyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM))
        !$OMP PARALLEL WORKSHARE
    	Jxyz = w1*ebState%eb1%Jxyz + w2*ebState%eb2%Jxyz
        !$OMP END PARALLEL WORKSHARE
	
		gGr%dbMAG_xyz = 0.0

    	associate( ebGr=>ebState%ebGr )
    	!Big-ass loop, loop over grid cells we want dB at and then loop over MHD grid contribution
    	!iG,jG,kG = ground grid
    	!iM,jM,kM = MHD grid

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(iG,jG,kG,iM,jM,kM) &
        !$OMP private(dV,r3,x0,xCC,ddB)
    	do kG=1,gGr%Nz
    		do jG=1,gGr%NLon
    			do iG=1,gGr%NLat
    				x0 = gGr%xyzC(iG,jG,kG,:) !Cell center of ground grid

    				do kM=ebGr%ks,ebGr%ke
    					do jM=ebGr%js,ebGr%je
    						do iM=ebGr%is+gGr%i0-1,ebGr%ie
    							xCC = ebGr%xyzcc(iM,jM,kM,:) !MHD grid cell center
    							if ( norm2(xCC) > gGr%rMax ) cycle
                                
                                r3 = norm2(x0-xCC)**3.0
    							dV = ebGr%dV(iM,jM,kM)
    							!Get differential contribution
                                !Avoid array temporary
                                ddB(XDIR) = -( Jxyz(iM,jM,kM,YDIR)*xCC(ZDIR) - Jxyz(iM,jM,kM,ZDIR)*xCC(YDIR) )/r3
                                ddB(YDIR) = -( Jxyz(iM,jM,kM,ZDIR)*xCC(XDIR) - Jxyz(iM,jM,kM,XDIR)*xCC(ZDIR) )/r3
                                ddB(ZDIR) = -( Jxyz(iM,jM,kM,XDIR)*xCC(YDIR) - Jxyz(iM,jM,kM,YDIR)*xCC(XDIR) )/r3

                                !Pulling out overall scaling
    							
                                gGr%dbMAG_xyz(iG,jG,kG,:) = gGr%dbMAG_xyz(iG,jG,kG,:) + dV*ddB

    						enddo !iM
    					enddo !jM
    				enddo !kM

    			enddo !iG
    		enddo !jG
    	enddo !kG

        !Moving overall scaling factor to secondary loop to only apply once
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(iG,jG,kG)
        do kG=1,gGr%Nz
            do jG=1,gGr%NLon
                do iG=1,gGr%NLat
                    gGr%dbMAG_xyz(iG,jG,kG,:) = B0*gGr%dbMAG_xyz(iG,jG,kG,:)/(4.0*PI)
                enddo !iG
            enddo !jG
        enddo !kG

    	end associate
    							
	end subroutine CalcMagDB


    !Lazy function to return scaling, should be replaced
    function bScale() result(B0)
        real(rp) :: B0
        real(rp) :: mu0,Mp,x0,u0,t0,d0,p0
        mu0 = 4*PI*1e-7
        Mp = 1.67e-27 ![kg]
        x0 = 1*6.38e6 ![m]   - RE
        u0 = 100e3    ![m/s] - 100 km/s
        t0 = x0/u0 ![s]   -
        d0 = Mp*1e6 ! [kg/m^3] - 1 particle/cc
        p0 = d0*u0*u0 ![N/m^2]
        B0 = sqrt(mu0*d0*u0*u0)*1e9 ! [nT]
    end function bScale

    !Calculate ionospheric DB contribution
    subroutine CalcIonDB(Model,ebState,gGr,ionGr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(sphGrid_T), intent(inout) :: gGr
        type(ionGrid_T), intent(inout) :: ionGr

        gGr%dbION_xyz = 0.0
        !CALCDB-TODO: Add the BS integral
    end subroutine CalcIonDB

    !Calculate ionospheric DB contribution
    subroutine CalcFacDB(Model,ebState,gGr,facGr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(sphGrid_T), intent(inout) :: gGr
        type(facGrid_T), intent(inout) :: facGr

        gGr%dbFAC_xyz = 0.0
        !CALCDB-TODO: Add the BS integral

    end subroutine CalcFacDB

end module calcdbcore
