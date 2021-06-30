module streamutils
    use chmpdefs
    use ebtypes
    use math
    use ebinterp
    
    implicit none

    integer , parameter, private :: NumRKF45 = 6
    real(rp), parameter, private :: TolRKF45 = 1.0e-3 !Units of planetary radius
    real(rp), dimension(NumRKF45), parameter, private :: & !Coefficients for RKF45
                    RKF45_LO = [25.0/216,0.0,1408.0/2565 ,2197.0/4104  ,-1.0/5 ,0.0   ], &
                    RKF45_HO = [16.0/216,0.0,6656.0/12825,28561.0/56430,-9.0/50,2.0/55]

    !Holder for data defining point on grid
    type GridPoint_T
        real(rp) :: xyz(NDIM)
        real(rp) :: t,dl !Time/lengthscale of current cell
        integer :: ijkG(NDIM) !Guess for cell location
    end type GridPoint_T

    !Generic one-step streamline routine
    abstract interface
        subroutine OneStep_T(gpt,Model,ebState,eps,h,dx)
            import :: rp,NDIM,GridPoint_T,chmpModel_T,ebState_T
            type(GridPoint_T), intent(inout) :: gpt
            type(chmpModel_T), intent(in)    :: Model
            type(ebState_T)  , intent(in)    :: ebState
            real(rp), intent(in) :: eps
            real(rp), intent(inout) :: h,dx(NDIM)
        end subroutine OneStep_T
    end interface

    procedure(OneStep_T), pointer :: StreamStep=>Step_RKF45

    contains

    subroutine Step_FE(gpt,Model,ebState,eps,h,dx)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)

        real(rp), dimension(NDIM) :: x0
        real(rp) :: t
        integer, dimension(NDIM) :: ijkG
        logical :: isGood

        x0 = gpt%xyz
        ijkG = gpt%ijkG
        t = gpt%t

        dx = h*MagHat(x0,t,Model,ebState,isGood,gpt%ijkG)
    end subroutine Step_FE

    subroutine Step_RKF45(gpt,Model,ebState,eps,h,dx)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)

        real(rp) :: eMax = 1.5 !Largest step size (units of local cell)
        real(rp), dimension(NDIM) :: x0,x2,x3,x4,x5,x6,k1,k2,k3,k4,k5,k6
        real(rp), dimension(NDIM) :: dxLO,dxHO
        real(rp) :: t,ddx,sScl,absh
        integer, dimension(NDIM) :: ijkG
        logical :: isGoods(NumRKF45)

        isGoods = .false.
        x0 = gpt%xyz
        ijkG = gpt%ijkG
        t = gpt%t

        k1 = h*MagHat(x0,t,Model,ebState,isGoods(1),ijkG)
        
        if (all(isGoods(1:1))) then
        	x2 = x0 + k1/4.0
        	k2 = h*MagHat(x2,t,Model,ebState,isGoods(2),ijkG)
        endif
        if (all(isGoods(1:2))) then
        	x3 = x0 + (3.0/32)*k1 + (9.0/32)*k2
        	k3 = h*MagHat(x3,t,Model,ebState,isGoods(3),ijkG)
        endif
        if (all(isGoods(1:3))) then
        	x4 = x0 + (1932.0/2197)*k1 - (7200.0/2197)*k2 + (7296.0/2197)*k3
        	k4 = h*MagHat(x4,t,Model,ebState,isGoods(4),ijkG)
        endif
        if (all(isGoods(1:4))) then
        	x5 = x0 + (439.0/216)*k1 - 8.0*k2 + (3680.0/513)*k3 - (845.0/4104)*k4
        	k5 = h*MagHat(x5,t,Model,ebState,isGoods(5),ijkG)
        endif
        if (all(isGoods(1:5))) then
        	x6 = x0 - (8.0/27)*k1 + 2.0*k2 - (3544.0/2565)*k3 + (1859.0/4104)*k4 - (11.0/40)*k5
        	k6 = h*MagHat(x6,t,Model,ebState,isGoods(6),ijkG)
        endif

        if (.not. all(isGoods)) then
            !At least one step was bad
            dx = 0.0
            h = 0.0
            return
        endif

        !If we're still here, the step was good so let's calculate
        dxLO = RKF45_LO(1)*k1 + RKF45_LO(2)*k2 + RKF45_LO(3)*k3 + RKF45_LO(4)*k4 + RKF45_LO(5)*k5 + RKF45_LO(6)*k6
        dxHO = RKF45_HO(1)*k1 + RKF45_HO(2)*k2 + RKF45_HO(3)*k3 + RKF45_HO(4)*k4 + RKF45_HO(5)*k5 + RKF45_HO(6)*k6

        ddx = max( norm2(dxLO-dxHO),TINY )
        sScl = 0.84*( (TolRKF45/ddx)**(0.25) )

        !Now calculate new step
        absh = abs(h)*sScl !Optimal value according to math
        !Clamp min/max step based on fraction of cell size
        call ClampValue(absh,eps*gpt%dl,eMax*gpt%dl)
        h = sign(absh,h)

        !Using HO step (but see some comments about using LO for stiff problems)
        dx = dxHO

    end subroutine Step_RKF45

    !Return magnetic field unit vector at given point/time
    function MagHat(xyz,t,Model,ebState,isIn,ijkG) result(bhat)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out) :: isIn
        integer, intent(in) :: ijkG(NDIM)
        real(rp) :: bhat(NDIM)

        real(rp), dimension(NDIM) :: B

        B = fldInterp(xyz,t,Model,ebState,BFLD,isIn,ijkG)
        bhat = normVec(B)

    end function MagHat

    function getDiag(ebGr,ijk) result (dl)
        type(ebGrid_T), intent(in)   :: ebGr
        integer, intent(in) :: ijk(NDIM)
        real(rp) :: dl
        integer :: i,j,k
        i = ijk(IDIR) ; j = ijk(JDIR) ; k = ijk(KDIR)

        dl = norm2( ebGr%xyz(i+1,j+1,k+1,:)-ebGr%xyz(i,j,k,:) )

    end function getDiag

    subroutine cleanStream(fL)
        type(fLine_T), intent(inout) :: fL

        integer :: i
        if (allocated(fL%xyz)) deallocate(fL%xyz)
        if (allocated(fL%ijk)) deallocate(fL%ijk)
        do i=0,NumVFL
            if (allocated(fL%lnVars(i)%V)) deallocate(fL%lnVars(i)%V)
        enddo
        fL%x0 = 0.0
        fL%Nm = 0
        fL%Np = 0
        fL%isGood = .false.
    end subroutine cleanStream                
end module streamutils