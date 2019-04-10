!Routines to handle interpolation on grid
!Mapping to reference element, weight calculation for interp and derivatives

module gridinterp
    use chmpdefs
    use ebtypes
    use gridloc
    
    implicit none

    !Weight type
    abstract interface
        function Wgt1D_T(eta) result(wgts)
            Import :: rp
            real(rp), intent(in) :: eta
            real(rp) :: wgts(-1:1)
        end function Wgt1D_T
    end interface

    !Default weight functions
    procedure(Wgt1D_T), pointer :: Wgt1D, Wgt1Dp

    logical, parameter :: doPole=.false.
    logical, parameter :: doCornerCut = .false.
    integer, parameter :: Nw=3 !Interpolation stencil size in each dimension

    contains

    !Given point xyz in cell i,j,k find eta,zeta,psi values
    !Force bounds on output ezp, |x| <= 1/2
    function Map2ezp(xyz,ijk,Model,ebGr) result(ezp)
        real(rp), intent(in) :: xyz(NDIM)
        integer, intent(in) :: ijk(NDIM)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in)   :: ebGr

        real(rp), dimension(NDIM) :: ezp, delx,lfmC
        integer :: n,i0,j0,k0
        integer :: ns,ne
        
        i0 = ijk(IDIR) ; j0 = ijk(JDIR) ; k0 = ijk(KDIR)
        
        delx = xyz-ebGr%xyzcc(i0,j0,k0,:) !xyz displacement from cell center
        ns = IDIR
        ne = KDIR

        !Do some fixes for geometry
        if (ebGr%GrID == EGGGRID .or. ebGr%GrID == LFMGRID) then
            lfmC = lfmCoords(xyz) !Get r,phi,thx coords
            !do K for either lfm/egg
            !Set angular coordinate manually (using rotation about x)

            ezp(KDIR) = (lfmC(KDIR)/locAux%dTh) + 0.5-1.0*k0
            if (ebGr%GrID == EGGGRID) then
                !Finish i/j and get out of here for speed
                ezp(JDIR) = (lfmC(JDIR)/locAux%dPhi) + 0.5-1.0*j0
                ezp(IDIR) = dot_product(ebGr%Tix(i0,j0,k0,IDIR,:),delx)
                !Enforce bounds and get out
                do n=1,NDIM
                    ezp(n) = min(ezp(n), 0.5)
                    ezp(n) = max(ezp(n),-0.5)
                enddo

                return

            else
                !LFM grid
                ne = JDIR
            endif
            
        endif

        !Do what's left (possibly everything) generically
        do n=ns,ne
            !Use Tix mapping to go to ezp coords
            ezp(n) = dot_product(ebGr%Tix(i0,j0,k0,n,:),delx)
        enddo

        !Enforce bounds here
        do n=1,NDIM
            ezp(n) = min(ezp(n), 0.5)
            ezp(n) = max(ezp(n),-0.5)
        enddo

    end function Map2ezp

    !-----------------------
    !Main 3D weight calculation
    !Given point xyz in cell ijk, calculate TSC weights
    !Map to reference element, calculate weights in continous index variables
    subroutine GetWeights(xyz,ijk,W,Model,ebGr)
        real(rp), intent(in) :: xyz(NDIM)
        integer, intent(in) :: ijk(NDIM)
        real(rp), intent(out) :: W(Nw,Nw,Nw)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in)   :: ebGr

        real(rp) :: ezp(NDIM) !Eta,zeta,psi coordinates
        real(rp), dimension(Nw) :: wE,wZ,wP !1D weights
        integer :: i,j,k

        !W = 0.0
        !W = 1.0/(Nw**3.0)
        !Do mapping to ezp coordinates
        ezp = Map2ezp(xyz,ijk,Model,ebGr)
        
        !Get 1D weights
        wE = Wgt1D(ezp(IDIR))
        wZ = Wgt1D(ezp(JDIR))
        wP = Wgt1D(ezp(KDIR))

        !Create multi-D weights
        do k=1,Nw
            do j=1,Nw
                do i=1,Nw
                    W(i,j,k) = wE(i)*wZ(j)*wP(k)
                enddo
            enddo
        enddo

        if (doCornerCut) call CutCorners(ijk,W,Model,ebGr)
        if (doPole) call wgtPole(ijk,W,Model,ebGr)
    end subroutine GetWeights

    !Zero out central weight
    subroutine wgtPole(ijk,W,Model,ebGr)
        integer, intent(in) :: ijk(NDIM)
        real(rp), intent(inout) :: W(-1:1,-1:1,-1:1)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in)   :: ebGr

        if ( (ijk(JDIR) == 1) .or. (ijk(JDIR) == ebGr%je) )then
            W(0,0,0) = 0.0
            W = W/sum(W)
        endif
    end subroutine wgtPole

    !Zero out weights of corner ghosts (side ghosts are fine)
    subroutine CutCorners(ijk,W,Model,ebGr)
        integer, intent(in) :: ijk(NDIM)
        real(rp), intent(inout) :: W(-1:1,-1:1,-1:1)
        type(chmpModel_T), intent(in) :: Model
        type(ebGrid_T), intent(in)   :: ebGr

        integer :: i,j,k,d1,d2,d3
        integer, dimension(-1:1,-1:1,-1:1) :: iW,jW,kW,gW

        i = ijk(IDIR) ; j = ijk(JDIR) ; k = ijk(KDIR)

        !TODO: Bail early

        iW = 0
        jW = 0
        kW = 0
        !Set xW to 1 where reliant on x ghost

        !I
        if (i == 1) then
            iW(-1,:,:) = 1
        else if (i == ebGr%Nip) then
            iW( 1,:,:) = 1
        endif

        !J
        if (j == 1) then
            jW(:,-1,:) = 1
        else if (j == ebGr%Njp) then
            jW(:, 1,:) = 1
        endif

        !K
        if (k == 1) then
            kW(:,:,-1) = 1
        else if (k == ebGr%Nkp) then
            kW(:,:, 1) = 1
        endif
        gW = iW+jW+kW

        if (maxval(gW)>=2) then
            !Only use cross
            do d1=-1,1
                do d2=-1,1
                    do d3=-1,1
                        if (abs(d1)+abs(d2)+abs(d3) > 1) then
                            W(d1,d2,d3) = 0.0
                        endif
                    enddo
                enddo
            enddo
            !Renormalize
            W = W/sum(W)

        endif


    end subroutine CutCorners

!-----------------------
!Weight functions
    
    !1D linear interpolation weights
    function lin1D(eta) result(wtsc)
        real(rp), intent(in) :: eta
        real(rp) :: wtsc(-1:1)
        real(rp) :: ceta
        ceta = eta
        ceta = max(-0.5,ceta)
        ceta = min(0.5,ceta)

        wtsc = 0.0
        if (eta >= 0) then
            !Right-sided
            wtsc(0) = 1-eta
            wtsc(1) = eta
        else
            !Left-sided
            wtsc( 0) = 1+eta
            wtsc(-1) = -eta
        endif
    end function lin1D

    !1D triangular shaped cloud weights
    !1D weights for triangular shaped cloud interpolation
    !Assuming on -1,1 reference element, dx=1
    !Check for degenerate cases ( |eta| > 0.5 )
    function tsc1D(eta) result(wtsc)
        real(rp), intent(in) :: eta
        real(rp) :: wtsc(-1:1)

        wtsc(-1) = 0.5*(0.5-eta)**2.0
        wtsc( 1) = 0.5*(0.5+eta)**2.0
        wtsc( 0) = 0.75 - eta**2.0

    end function tsc1D

    !1D weights for derivatives of TSC weights
    function tsc1Dp(eta) result(wp)
        real(rp), intent(in) :: eta
        real(rp) :: wp(-1:1)

        wp(-1) = -(0.5 - eta)
        wp( 1) =   0.5 + eta
        wp( 0) = -2.0*eta
    end function tsc1Dp

    !Parabolic interpolation
    function quad1D(eta) result(wtsc)
        real(rp), intent(in) :: eta
        real(rp) :: wtsc(-1:1)
        real(rp) :: ceta
        ceta = eta
        ceta = max(-0.5,ceta)
        ceta = min(0.5,ceta)

        wtsc = 0.0
        wtsc(-1) = 0.5*ceta*(ceta-1)
        wtsc( 0) = -1.0*(ceta-1)*(ceta+1)
        wtsc(+1) = 0.5*ceta*(ceta+1)

    end function quad1D
end module gridinterp