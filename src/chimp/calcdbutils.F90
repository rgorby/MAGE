module calcdbutils
	use kdefs
	use chmpdefs
	use ebtypes
    use calcdbtypes
    
	implicit none

    real(rp), parameter :: RIon = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880
	real(rp) :: dzGG = 60.0 !Default height spacing [km]

	contains

	subroutine facGridInit(Model,ebState,rmState,facGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(facGrid_T)  , intent(inout) :: facGrid

        integer :: Np,Nth

        Np = rmState%Np
        Nth = rmState%Nth

        !Allocate and zero out facGrid
        allocate(facGrid%XYZcc(Np,Nth,rSegs,2,NDIM))
        allocate(facGrid%Jxyz (Np,Nth,rSegs,2,NDIM))
        allocate(facGrid%dV   (Np,Nth,rSegs,2))
        facGrid%XYZcc = 0.0
        facGrid%Jxyz  = 0.0
        facGrid%dV    = 0.0
        facGrid%Np = Np
        facGrid%Nth = Nth
        facGrid%rSegs = rSegs

        !CALCDB-TODO: Need to add code here to define facGrid%XYZcc and dV
        !NOTE: XYZ from remix is in units of Rion, not Re

	end subroutine facGridInit

    subroutine ionGridInit(Model,ebState,rmState,ionGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(ionGrid_T)  , intent(inout) :: ionGrid

        integer :: Np,Nth

        Np = rmState%Np
        Nth = rmState%Nth

        allocate(ionGrid%XYZcc(Np,Nth,2,NDIM))
        allocate(ionGrid%Jxyz (Np,Nth,2,NDIM))
        allocate(ionGrid%dS   (Np,Nth,2)) !Surface area per patch

        ionGrid%XYZcc = 0.0
        ionGrid%Jxyz  = 0.0
        ionGrid%dS    = 0.0
        ionGrid%Np = Np
        ionGrid%Nth = Nth
        !CALCDB-TODO: Need to add code here to define grid (XYZcc and dS)
        !NOTE: XYZ from remix is in units of Rion, not Re

    end subroutine ionGridInit

    !Initialize holders for Bios-Savart contributions
    subroutine BSGridInit(Model,ebState,rmState,magBS,ionBS,facBS)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(BSGrid_T), intent(inout) :: magBS,ionBS,facBS

        integer :: NMag,NIon,NFac
        real(rp) :: B0
    !Mag BS grid
        B0 = bScale()
        !Number of cells
        NMag = (ebState%ebGr%Nip)*(ebState%ebGr%Njp)*(ebState%ebGr%Nkp)
        call BSSubInit(magBS,NMag)
        magBS%jScl = B0/(4.0*PI) !Scaling factor

        !CALCDB-TODO: Add ionBS%jScl and facBS%jScl values
    !Ion BS grid
        !Number of cells
        NIon = (rmState%Np)*(rmState%Nth)*(2) !Include N/S hemispheres
        call BSSubInit(ionBS,NIon)
        ionBS%jScl = 1.0
    !FAC BS grid
        !Number of cells
        NFac = NIon*rSegs
        call BSSubInit(facBS,NFac)
        facBS%jScl = 1.0

        contains

        !Initialize BS grid w/ N points
        subroutine BSSubInit(xBS,N)
            type(BSGrid_T), intent(inout) :: xBS
            integer, intent(in) :: N
            xBS%NumP = N
            allocate(xBS%XYZcc(N,NDIM))
            allocate(xBS%Jxyz (N,NDIM))
            allocate(xBS%dV(N))
            xBS%XYZcc = 0.0
            xBS%Jxyz = 0.0
            xBS%dV = 0.0
        end subroutine BSSubInit

    end subroutine BSGridInit

	!Using a rmState (remix data), fill facGrid Jxyz
	subroutine facGridUpdate(Model,ebState,rmState,facGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(facGrid_T)  , intent(inout) :: facGrid

        !CALCDB-TODO: Write this
        facGrid%Jxyz = 0.0

    end subroutine facGridUpdate

    !Using a rmState (remix data), fill ionGrid Jxyz
    subroutine ionGridUpdate(Model,ebState,rmState,ionGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(ionGrid_T)  , intent(inout) :: ionGrid

        !CALCDB-TODO: Write this
        !Need to calculate E field from potential in both N/S hemispheres
        !Then get J from E and SigP/SigH

        ionGrid%Jxyz = 0.0

    end subroutine ionGridUpdate

    !Set rmState given properly set 4 hemispheres and temporal weights
    subroutine hemi2rm(rmState,w1,w2)
        type(rmState_T)  , intent(inout) :: rmState
        real(rp), intent(in) :: w1,w2

        rmState%nFac  = w1*rmState%rmN1%xFac  + w2*rmState%rmN2%xFac 
        rmState%nPot  = w1*rmState%rmN1%xPot  + w2*rmState%rmN2%xPot 
        rmState%nSigP = w1*rmState%rmN1%xSigP + w2*rmState%rmN2%xSigP
        rmState%nSigH = w1*rmState%rmN1%xSigH + w2*rmState%rmN2%xSigH
        rmState%sFac  = w1*rmState%rmS1%xFac  + w2*rmState%rmS2%xFac 
        rmState%sPot  = w1*rmState%rmS1%xPot  + w2*rmState%rmS2%xPot 
        rmState%sSigP = w1*rmState%rmS1%xSigP + w2*rmState%rmS2%xSigP
        rmState%sSigH = w1*rmState%rmS1%xSigH + w2*rmState%rmS2%xSigH


    end subroutine hemi2rm

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

end module calcdbutils
