module gamutils
    use gamtypes

    implicit none
    
    contains

    !Loads a brickette of size (iMax,recLen) 
    !Into a brickette of size (vecLen,recLen)
    !Necessary to reconstruct in direction d, all cells from iS:iS+iMax,j,k
    !Here we recast Q from 3D->1D, and use explicit stride to ensure proper vector instructions
    !NOTE: Assuming that memory being pulled from is aligned
    subroutine LoadBlock(Model,Gr,Qb,Q,iS,j,k,iMax,d)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: iS,j,k,iMax,d
        !Remap Q to 1D from (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), intent(in) :: Q(Gr%Ni*Gr%Nj*Gr%Nk) 
        real(rp), intent(out) :: Qb(vecLen,-Nr2:Nr2-1)
        
        integer :: l,n,l0,nSi,nSj,nSk,nSt
        !DIR$ ASSUME_ALIGNED Qb: ALIGN
        !DIR$ ASSUME_ALIGNED Q: ALIGN
            
        !Mapping of iS,j,k to 1D (l0)
        l0 = ijk2n(Gr,iS,j,k)

        !Strides dep. on direction
        nSi = 1; nSj = Gr%Ni; nSk = Gr%Ni*Gr%Nj

        !Set nSt equal to nSi/nSj/nSk for given direction
        nSt = nSi*ijkD(d,IDIR) + nSj*ijkD(d,JDIR) + nSk*ijkD(d,KDIR)

        do n = -Nr2,Nr2-1
            l = l0 + nSt*n !Jump to next 1:iMax in direction d
            Qb(1:iMax,n) = Q(l:l+iMax-1)

        enddo

        !In comments: long-form loop using Q(isg:,jsg:,ksg:)
        !Does same thing, slower but more readable
        ! do n = -Nr2,Nr2-1
        !     do i=1,iMax
        !         iG = i+iS-1 !Global index into Q
        !         ip = iG + n*ijkD(d,IDIR)
        !         jp = j  + n*ijkD(d,JDIR)
        !         kp = k  + n*ijkD(d,KDIR)
        !         Qb(i,n) = Q(ip,jp,kp)
        !     enddo
        ! enddo   

    end subroutine LoadBlock

    !Same as LoadBlock but for interface-centered variables
    !Not necessarily aligned
    subroutine LoadBlockI(Model,Gr,Qb,Q,iS,j,k,iMax,d)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: iS,j,k,iMax,d
        !Remap Q to 1D from (isg:ieg,jsg:jeg,ksg:keg)
        real(rp), intent(in) :: Q( (Gr%Ni+1)*(Gr%Nj+1)*(Gr%Nk+1) ) 
        real(rp), intent(out) :: Qb(vecLen,-Nr2:Nr2-1)
        
        integer :: l,n,l0,nSi,nSj,nSk,nSt
            
        !Mapping of iS,j,k to 1D (l0)
        l0 = ijk2nI(Gr,iS,j,k)

        !Strides dep. on direction
        nSi = 1; nSj = Gr%Ni+1; nSk = (Gr%Ni+1)*(Gr%Nj+1)

        !Set nSt equal to nSi/nSj/nSk for given direction
        nSt = nSi*ijkD(d,IDIR) + nSj*ijkD(d,JDIR) + nSk*ijkD(d,KDIR)
        Qb = 0.0
        do n = -Nr2,Nr2-1
            l = l0 + nSt*n !Jump to next 1:iMax in direction d
            Qb(1:iMax,n) = Q(l:l+iMax-1)

        enddo

    end subroutine LoadBlockI

    ! subroutine LoadBlockI_Old(Model,Gr,Qb,Q,iS,j,k,iMax,d)
    !     type(Model_T), intent(in) :: Model
    !     type(Grid_T),  intent(in) :: Gr
    !     integer, intent(in) :: iS,j,k,iMax,d
    !     real(rp), intent(in) :: Q(Gr%isg:Gr%ieg+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1)
    !     real(rp), intent(out) :: Qb(vecLen,-Nr2:Nr2-1)

    !     integer :: i,n,iG,ip,jp,kp

    !     Qb = 0.0
    !     do n = -Nr2,Nr2-1
    !         do i=1,iMax
    !             iG = i+iS-1 !Global index into Q
    !             ip = iG + n*ijkD(d,IDIR)
    !             jp = j  + n*ijkD(d,JDIR)
    !             kp = k  + n*ijkD(d,KDIR)

    !             Qb(i,n) = Q(ip,jp,kp)
                
    !         enddo
    !     enddo

    ! end subroutine LoadBlockI_Old

    !Converts i,j,k index (isg:ieg,jsg:jeg,ksg:keg)
    !To 1D index, n= 1:Ni*Nj*Nk)

    function ijk2n(Gr,i,j,k) result (n)
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: i,j,k
        integer :: n

        integer :: ip,jp,kp

        !Convert to ip,jp,kp (1:Ni,j,k)
        ip = (i-Gr%isg)+1
        jp = (j-Gr%jsg)+1
        kp = (k-Gr%ksg)+1

        !Convert to n
        n = ip + (jp-1)*Gr%Ni + (kp-1)*(Gr%Ni*Gr%Nj)

    end function ijk2n

    function ijk2nI(Gr,i,j,k) result (n)
        type(Grid_T),  intent(in) :: Gr
        integer, intent(in) :: i,j,k
        integer :: n

        integer :: ip,jp,kp

        !Convert to ip,jp,kp (1:Ni,j,k)
        ip = (i-Gr%isg)+1
        jp = (j-Gr%jsg)+1
        kp = (k-Gr%ksg)+1

        !Convert to n
        n = ip + (jp-1)*(Gr%Ni+1) + (kp-1)*( (Gr%Ni+1)*(Gr%Nj+1) )

    end function ijk2nI

    subroutine allocState(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        type(State_T), intent(inout) :: State

        if ( .not. allocated(State%Gas) ) then
            allocate( State%Gas(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NVAR,0:Model%nSpc) )
            State%Gas = 0.0
        endif

        if ( .not. allocated(State%magFlux) ) then
            call allocGridVec(Model,Grid,State%magFlux,.true.,NDIM)
            State%magFlux = 0.0
        endif
        
        if ( .not. allocated(State%Efld) ) then
            call allocGridVec(Model,Grid,State%Efld,.false.,NDIM)
            State%Efld = 0.0
        endif

        if ( Model%doResistive .and. .not. allocated(State%Deta) ) then
            call allocGridVec(Model,Grid,State%Deta,.true.,NDIM)
            State%Deta = 0.0
        endif

        if ( .not. allocated(State%Bxyz) ) then
            call allocGridVec(Model,Grid,State%Bxyz,.false.,NDIM)
            State%Bxyz = 0.0
        endif

        State%time = Model%t

    end subroutine allocState

    subroutine deallocState(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        type(State_T), intent(inout) :: State
        
        if ( allocated(State%Gas)     ) deallocate(State%Gas)
        if ( allocated(State%magFlux) ) deallocate(State%magFlux)
        if ( allocated(State%Efld)    ) deallocate(State%Efld)
        if ( allocated(State%Bxyz)    ) deallocate(State%Bxyz)
        if ( allocated(State%Deta)    ) deallocate(State%Deta)

    end subroutine deallocState
    
    !Allocates space for a grid-sized variable (w/ ghosts)
    !If doP1, create all dims+1
    subroutine allocGridVar(Model,Grid,gV,doP1)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp), dimension(:,:,:), allocatable, intent(inout) :: gV
        logical, optional, intent(in) :: doP1

        logical :: doPlus
        if ( .not. allocated(gV) ) then
            !Choose size, either all xlg:xhg or xlg:xhg+1
            if (present(doP1)) then
                doPlus = doP1
            else
                doPlus = .false.
            endif

            if (doPlus) then
                allocate( gV(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1) )
            else
                allocate( gV(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg) )
            endif
            !Go ahead and wipe newly allocated data
            gV = 0.0
            
        endif

    end subroutine allocGridVar

    !Allocates space for a grid-sized vector (w/ ghosts)
    !If doP1, create all dims+1
    !numDim is the number of vector components (ie, for state vector or spatial vector)
    subroutine allocGridVec(Model,Grid,gV,doP1,numDim)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp), dimension(:,:,:,:), allocatable, intent(inout) :: gV
        logical, optional, intent(in) :: doP1
        integer, optional, intent(in) :: numDim

        logical :: doPlus
        integer :: nD

        if ( .not. allocated(gV) ) then
            !Choose size, either all xlg:xhg or xlg:xhg+1
            if (present(doP1)) then
                doPlus = doP1
            else
                doPlus = .false.
            endif
            !Choose number of elements in vector, default to NDIM
            if (present(numDim)) then
                nD = numDim
            else
                nD = NDIM
            endif

            if (doPlus) then
                allocate( gV(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:nD ) )
            else
                allocate( gV(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:nD) )
            endif
            
            !Go ahead and wipe newly allocated data
            gV = 0.0
        endif       

    end subroutine allocGridVec

    !Wipes a gridVec structure
    !If doP1, create all dims+1
    !numDim is the number of vector components (ie, for state vector or spatial vector)
    subroutine wipeGridVec(Model,Grid,gV,doP1,numDim)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp), dimension(:,:,:,:),intent(out) :: gV
        logical, optional, intent(in) :: doP1
        integer, optional, intent(in) :: numDim

        logical :: doPlus
        integer :: nD,n,i,j,k,Ni,Nj,Nk

        if (present(doP1)) then
            doPlus = doP1
        else
            doPlus = .false.
        endif

        !Choose number of elements in vector, default to NDIM
        if (present(numDim)) then
            nD = numDim
        else
            nD = NDIM
        endif

        Ni = Grid%Ni
        Nj = Grid%Nj
        Nk = Grid%Nk

        if (doPlus) then
            Ni = Ni +1
            Nj = Nj +1
            Nk = Nk +1
        endif

        !$OMP PARALLEL DO default(shared) collapse(3)
        !Loop over and wipe array
        do n=1,nD
            do k=1,Nk
                do j=1,Nj
                    do i=1,Ni
                        gV(i,j,k,n) = 0.0
                    enddo
                enddo
            enddo
        enddo

    end subroutine wipeGridVec

    !Converts grid of conserved quantities to primitive
    subroutine Con2Prim(Model,Grid,Con,Prim)
        type(Model_T), intent(in) :: Model
        type(Grid_T),  intent(in) :: Grid
        real(rp), intent(in) :: Con(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NVAR)
        real(rp), intent(out) :: Prim(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NVAR)

        integer :: i,j,k
        real(rp) :: rho,Vx,Vy,Vz,E,KinE

        do k=Grid%ksg, Grid%keg
            do j=Grid%jsg, Grid%jeg
                do i=Grid%isg, Grid%ieg
                    rho = Con(i,j,k,DEN)
                    Vx = Con(i,j,k,MOMX)/rho
                    Vy = Con(i,j,k,MOMY)/rho
                    Vz = Con(i,j,k,MOMZ)/rho
                    E = Con(i,j,k,ENERGY)

                    KinE = 0.5*rho*(Vx**2.0+Vy**2.0+Vz**2.0)
                    Prim(i,j,k,DEN) = rho
                    Prim(i,j,k,VELX) = Vx
                    Prim(i,j,k,VELY) = Vy
                    Prim(i,j,k,VELZ) = Vz
                    Prim(i,j,k,PRESSURE) = (Model%gamma-1)*(E-KinE)
                enddo
            enddo
        enddo

    end subroutine Con2Prim

    subroutine CellC2P(Model,Con,Prim)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in)  :: Con(NVAR)
        real(rp), intent(out) :: Prim(NVAR)

        real(rp) :: rho,Vx,Vy,Vz,E,KinE

        rho = Con(DEN)
        Vx  = Con(MOMX)/rho
        Vy  = Con(MOMY)/rho
        Vz  = Con(MOMZ)/rho
        E   = Con(ENERGY)

        KinE = 0.5*rho*(Vx**2.0+Vy**2.0+Vz**2.0)

        Prim(DEN)      = max(rho,dFloor)
        Prim(VELX)     = Vx
        Prim(VELY)     = Vy
        Prim(VELZ)     = Vz
        Prim(PRESSURE) = max((Model%gamma-1)*(E-KinE),pFloor)
    end subroutine CellC2P

    subroutine CellPress2Cs(Model,Con,Cs)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in)  :: Con(NVAR)
        real(rp), intent(out) :: Cs

        real(rp) :: rho,Vx,Vy,Vz,E,KinE,P

        rho = Con(DEN)
        Vx  = Con(MOMX)/rho
        Vy  = Con(MOMY)/rho
        Vz  = Con(MOMZ)/rho
        E   = Con(ENERGY)

        KinE = 0.5*rho*(Vx**2.0+Vy**2.0+Vz**2.0)

        P = (Model%gamma-1)*(E-KinE)
        Cs = sqrt(Model%gamma*P/rho)
      end subroutine CellPress2Cs

    subroutine CellP2C(Model,Prim,Con)
        type(Model_T), intent(in) :: Model
        real(rp), intent(out) :: Con(NVAR)
        real(rp), intent(in)  :: Prim(NVAR)

        real(rp) :: rho,Vx,Vy,Vz,KinE,IntE

        rho  = max(Prim(DEN),dFloor)
        Vx   = Prim(VELX)
        Vy   = Prim(VELY)
        Vz   = Prim(VELZ)

        IntE = max(Prim(PRESSURE),pFloor)/(Model%gamma-1)
        KinE = 0.5*rho*(Vx**2.0+Vy**2.0+Vz**2.0)

        Con(DEN)    = rho
        Con(MOMX)   = rho*Vx
        Con(MOMY)   = rho*Vy
        Con(MOMZ)   = rho*Vz
        Con(ENERGY) = IntE+KinE
        
    end subroutine CellP2C

end module gamutils
