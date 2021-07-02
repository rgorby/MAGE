module streamutils
    use chmpdefs
    use ebtypes
    use math
    use ebinterp
    
    implicit none

    integer , parameter, private :: NumRKF45 = 6
    integer , parameter, private :: NumBS23 = 4
    real(rp), parameter, private :: StreamTol = 1.0e-3 !Units of planetary radius
    real(rp), dimension(NumRKF45), parameter, private :: & !Coefficients for RKF45
                    RKF45_LO = [25.0/216,0.0,1408.0/2565 ,2197.0/4104  ,-1.0/5 ,0.0   ], &
                    RKF45_HO = [16.0/216,0.0,6656.0/12825,28561.0/56430,-9.0/50,2.0/55]
    real(rp), parameter, private :: eMax = 8.0 !Units of local cell

    !Holder for data defining point on grid
    type GridPoint_T
        real(rp) :: xyz(NDIM)
        real(rp) :: t,dl !Time/lengthscale of current cell
        integer :: ijkG(NDIM) !Guess for cell location
    end type GridPoint_T

    !Generic one-step streamline routine
    !NOTE: TO handle FSAL methods, possible we get first B and output B at final position
    abstract interface
        subroutine OneStep_T(gpt,Model,ebState,eps,h,dx,iB,oB)
            import :: rp,NDIM,GridPoint_T,chmpModel_T,ebState_T
            type(GridPoint_T), intent(inout) :: gpt
            type(chmpModel_T), intent(in)    :: Model
            type(ebState_T)  , intent(in)    :: ebState
            real(rp), intent(in) :: eps
            real(rp), intent(inout) :: h,dx(NDIM)
            real(rp), intent(in) , dimension(NDIM), optional :: iB
            real(rp), intent(out), dimension(NDIM), optional :: oB
        end subroutine OneStep_T
    end interface

    procedure(OneStep_T), pointer :: StreamStep=>Step_RKF45
    !procedure(OneStep_T), pointer :: StreamStep=>Step_RK4L

    contains

    subroutine Step_FE(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB

        real(rp), dimension(NDIM) :: x0
        logical :: isGood

        x0 = gpt%xyz
        if (present(iB)) then
            dx = h*normVec(iB)
            isGood = .true.
        else
            dx = h*FastHat(x0,gPt%t,Model,ebState,isGood,gpt%ijkG)
        endif

        if (present(oB)) then
            oB = FastMag(x0+dx,gPt%t,Model,ebState,isGood,gpt%ijkG)
        endif

    end subroutine Step_FE

    subroutine Step_BS23(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB


        real(rp), dimension(NDIM) :: x0,x2,x3,x4,k1,k2,k3,k4
        real(rp), dimension(NDIM) :: dxLO,dxHO,B
        real(rp) :: ddx,sScl,absh
        logical :: isGoods(NumBS23),isG
        real(rp), parameter :: eStpT = 1.0e-3
        real(rp), parameter :: dlMax = 2.0 !Max cell length

        write(*,*) 'This method needs to be debugged!'
        isGoods = .false.
        x0 = gpt%xyz

        if (present(iB)) then
            k1 = h*normVec(iB)
            isGoods(1) = .true.
        else
            k1 = h*FastHat(x0,gPt%t,Model,ebState,isGoods(1),gPt%ijkG)
        endif
        if (all(isGoods(1:1))) then
            x2 = x0 + k1/2.0
            k2 = h*FastHat(x2,gPt%t,Model,ebState,isGoods(2),gPt%ijkG)
        endif
        if (all(isGoods(1:2))) then
            x3 = x0 + (3.0/4)*k2
            k3 = h*FastHat(x3,gPt%t,Model,ebState,isGoods(3),gPt%ijkG)
        endif
        
        if (all(isGoods(1:3))) then
            dxHO = (2*k1+3*k2+4*k3)/9.0
            x4 = x0 + dxHO !Actual solution
            B = FastMag(x4,gPt%t,Model,ebState,isGoods(4),gPt%ijkG)
            if (present(oB)) then
                oB = B
            endif
            k4 = h*normVec(B)
        endif
        !Check for failure
        if (.not. all(isGoods)) then
            !At least one step was bad
            dx = 0.0
            h = 0.0
            if (present(oB)) oB = 0.0
            return
        endif
    !Now finish up
        dxLO = (7*k1+6*k2+8*k3+3*k4)/24.0
        dx = dxHO
        
    ! !Kutta-Merson style
    !     ddx = norm2(dxHO-dxLO)/gpt%dl
    !     if (ddx >= eStpT) then
    !         h = h/2.0 !Reduce step
    !     else if (ddx <= eStpT/64.0) then
    !         h = 2.0*h !Increase step
    !     endif
    !     absh = abs(h)
    !     call ClampValue(absh,eps*gpt%dl,eMax*gpt%dl)
        

    !Embedded opt style
        ! ddx = max( norm2(dxLO-dxHO),TINY )
        ! !sScl = 0.9*0.7*( (StreamTol/ddx)**(0.5) )
        ! !sScl = 0.8*(StreamTol/ddx)**(0.25)

        ! !Now calculate new step
        ! absh = abs(h)*sScl !Optimal value according to math
        ! !Clamp min/max step based on fraction of cell size
        ! call ClampValue(absh,eps*gpt%dl,eMax*gpt%dl)

        h = sign(absh,h)

    end subroutine Step_BS23

    subroutine Step_RK4L(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB

        real(rp), dimension(NDIM) :: Jb,Jb2,Jb3,F1,F2,F3,F4
        real(rp), dimension(NDIM) :: x0,E,B
        real(rp) :: ds,dsmag,MagB,MagJb
        logical :: isGood
        type(gcFields_T) :: gcF

        x0 = gpt%xyz
        !Doing full field calc
        call ebFields(x0,gPt%t,Model,ebState,E,B,gPt%ijkG,gcFields=gcF)

        !Only using sign of h
        MagB = norm2(B)
        MagJb = norm2(gcF%JacB)
        if (MagJb <= TINY) then
            !Field is constant-ish, use local grid size
            dsmag = gpt%dl
        else
            dsmag = MagB/MagJb
        endif
        ds = sign(1.0_rp,h)*eps*min(gpt%dl,dsmag)
        !ds = sign(1.0_rp,h)*min(gpt%dl,eps*dsmag)
        !Save step for next round
        h = ds

        !Convert ds to streamline units
        ds = ds/max(MagB,TINY)

        !Get powers of jacobian
        Jb  = matmul(gcF%JacB,B  )
        Jb2 = matmul(gcF%JacB,Jb )
        Jb3 = matmul(gcF%JacB,Jb2)

        !Calculate steps
        F1 = ds*B
        F2 = F1 + (ds*ds/2)*Jb
        F3 = F2 + (ds*ds*ds/4)*Jb2
        F4 = ds*B + ds*ds*Jb + (ds**3.0/2.0)*Jb2 + (ds**4.0/4.0)*Jb3

        !Advance
        dx = (F1+2*F2+2*F3+F4)/6.0
 
        if (present(oB)) then
            oB = FastMag(x0+dx,gPt%t,Model,ebState,isGood,gpt%ijkG)
        endif

    end subroutine Step_RK4L

    subroutine Step_RKF45(gpt,Model,ebState,eps,h,dx,iB,oB)
        type(GridPoint_T), intent(inout) :: gpt
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)    :: ebState
        real(rp), intent(in) :: eps
        real(rp), intent(inout) :: h,dx(NDIM)
        real(rp), intent(in) , dimension(NDIM), optional :: iB
        real(rp), intent(out), dimension(NDIM), optional :: oB

        real(rp), dimension(NDIM) :: x0,x2,x3,x4,x5,x6,k1,k2,k3,k4,k5,k6
        real(rp), dimension(NDIM) :: dxLO,dxHO
        real(rp) :: ddx,sScl,absh
        logical :: isGoods(NumRKF45),isG

        isGoods = .false.
        x0 = gpt%xyz

        if (present(iB)) then
            k1 = h*normVec(iB)
            isGoods(1) = .true.
        else
            k1 = h*FastHat(x0,gPt%t,Model,ebState,isGoods(1),gPt%ijkG)
        endif

        if (all(isGoods(1:1))) then
        	x2 = x0 + k1/4.0
        	k2 = h*FastHat(x2,gPt%t,Model,ebState,isGoods(2),gPt%ijkG)
        endif
        if (all(isGoods(1:2))) then
        	x3 = x0 + (3.0/32)*k1 + (9.0/32)*k2
        	k3 = h*FastHat(x3,gPt%t,Model,ebState,isGoods(3),gPt%ijkG)
        endif
        if (all(isGoods(1:3))) then
        	x4 = x0 + (1932.0/2197)*k1 - (7200.0/2197)*k2 + (7296.0/2197)*k3
        	k4 = h*FastHat(x4,gPt%t,Model,ebState,isGoods(4),gPt%ijkG)
        endif
        if (all(isGoods(1:4))) then
        	x5 = x0 + (439.0/216)*k1 - 8.0*k2 + (3680.0/513)*k3 - (845.0/4104)*k4
        	k5 = h*FastHat(x5,gPt%t,Model,ebState,isGoods(5),gPt%ijkG)
        endif
        if (all(isGoods(1:5))) then
        	x6 = x0 - (8.0/27)*k1 + 2.0*k2 - (3544.0/2565)*k3 + (1859.0/4104)*k4 - (11.0/40)*k5
        	k6 = h*FastHat(x6,gPt%t,Model,ebState,isGoods(6),gPt%ijkG)
        endif

        if (.not. all(isGoods)) then
            !At least one step was bad
            dx = 0.0
            h = 0.0
            if (present(oB)) oB = 0.0
            return
        endif

        !If we're still here, the step was good so let's calculate
        dxLO = RKF45_LO(1)*k1 + RKF45_LO(2)*k2 + RKF45_LO(3)*k3 + RKF45_LO(4)*k4 + RKF45_LO(5)*k5 + RKF45_LO(6)*k6
        dxHO = RKF45_HO(1)*k1 + RKF45_HO(2)*k2 + RKF45_HO(3)*k3 + RKF45_HO(4)*k4 + RKF45_HO(5)*k5 + RKF45_HO(6)*k6

        ddx = max( norm2(dxLO-dxHO),TINY )
        !sScl = 0.84*( (StreamTol/ddx)**(0.25) )
        sScl = 0.84*( (StreamTol*gpt%dl/ddx)**(0.25) ) !Relative

        !Now calculate new step
        absh = 0.9*abs(h)*sScl !Optimal value according to math
        !Clamp min/max step based on fraction of cell size
        !call ClampValue(absh,eps*gpt%dl,eMax*gpt%dl)
        call ClampValue(absh,eps*gpt%dl,3*abs(h))

        h = sign(absh,h)

        !Using HO step (but see some comments about using LO for stiff problems)
        dx = dxHO
        if (present(oB)) then
            oB = FastMag(x0+dx,gPt%t,Model,ebState,isG,gPt%ijkG)
        endif
    end subroutine Step_RKF45

    !Return magnetic field unit vector at given point/time
    function FastHat(xyz,t,Model,ebState,isIn,ijkG) result(bhat)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out) :: isIn
        integer, intent(inout) :: ijkG(NDIM)
        real(rp) :: bhat(NDIM)

        real(rp), dimension(NDIM) :: B

        B = FastMag(xyz,t,Model,ebState,isIn,ijkG)
        bhat = normVec(B)
    end function FastHat

    !Return magnetic field vector at given point/time
    !NOTE: This reproduces a lot of ebinterp code but is designed to be streamlined/faster
    function FastMag(xyz,t,Model,ebState,isIn,ijkG) result(B)
        real(rp), intent(in) :: xyz(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        logical, intent(out) :: isIn
        integer, intent(inout) :: ijkG(NDIM)
        real(rp), dimension(NDIM) :: B

        real(rp), dimension(NDIM) :: B0
        real(rp), dimension(Nw,Nw,Nw) :: W
        real(rp), dimension(Nw,Nw,Nw,NDIM) :: Qb !Buffer
        integer :: n,i0,j0,k0,i1,i2,i3

        !B = fldInterp(xyz,t,Model,ebState,BFLD,isIn,ijkG)
        !Locate w/ guess
        call locate(xyz,ijkG,Model,ebState%ebGr,isIn,ijkG)
        if (.not. isIn) then
            B = Model%B0(xyz)
            return
        endif

    !If still here, do work. NOTE: Assuming static field to avoid extra work
        B0 = Model%B0(xyz)
        call GetWeights(xyz,ijkG,W,Model,ebState%ebGr)
        i0 = ijkG(IDIR) ; j0 = ijkG(JDIR) ; k0 = ijkG(KDIR)

        Qb = ebState%eb1%dB(i0-1:i0+1,j0-1:j0+1,k0-1:k0+1,XDIR:ZDIR)

        do n=1,NDIM
            !NOTE: Avoiding doing below to avoid array temporary creation
            !B(n) = sum( W(:,:,:)*Qb(:,:,:,n) ) + B0(n)
            B(n) = B0(n)
            do i3=1,Nw
                do i2=1,Nw
                    do i1=1,Nw
                        B(n) = B(n) + W(i1,i2,i3)*Qb(i1,i2,i3,n)
                    enddo
                enddo
            enddo
        enddo

    end function FastMag

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