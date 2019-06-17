!Various utilities to handle solar wind BC

module wind
    use types
    use gamutils
    use math
    use xml_input
    use ioH5
    use gridutils
    use msphutils
    use multifluid

    implicit none

    ! Global Parameters
    integer, parameter :: SWSPC = 1 !SW fluid is always 1st in multifluid
    integer, parameter :: MAXWINDVARS = 20

    !Type for generic solar wind BC (from file or subroutine)
    !Either use discrete tW,Qw(NVAR) series or subroutine
    type, extends(baseBC_T) :: WindBC_T

        character(len=strLen) :: wID !Wind ID string
        logical :: isDiscrete=.false.

        !Holder for subroutine data
        procedure(WindTS_T), pointer, nopass :: getWind => NULL()

        !Holder for discrete data
        real(rp), allocatable :: tW(:), Q(:,:), B(:,:)
        real(rp) :: tMin,tMax
        integer :: NumT

        !Coefficients for solar wind geometry
        real(rp) :: ByC = 0.0, BzC = 0.0

        logical :: doDSW = .false. !Density perturbations in SW
        real(rp) :: dSWAmp = 0.0 !Magnitude of density perturbations (0,1)
        real(rp) :: dtDSW = 0.125 !Period of perturbations
        real(rp) :: tDSW = 0.0 !Next re-perturbation

        integer :: NgW = 4 !Number of solar wind ghost cells
        real(rp), dimension(NDIM) :: xyzW !Position of solar wind time-series (ie L1)

        !Shell values for solar wind fields
        !Size = [NgW,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,NDIM]
        real(rp), dimension(:,:,:,:), allocatable :: BxyzW,ExyzW,VxyzW

        !Shell values for solar wind scalars
        real(rp), dimension(:,:,:), allocatable :: RhoW,PrW,dRhoW

        !Boolean values on outer shell for solar wind influence
        logical, dimension(:,:,:), allocatable :: isWind

        contains

        procedure :: doInit => InitWind
        procedure :: doBC => WindBC

    end type WindBC_T


    !Type for solar wind time-series subroutine
    abstract interface
        subroutine WindTS_T(windBC,Model,t,Rho,Pr,V,B)
            Import :: rp,Model_T,NDIM,WindBC_T
            class(WindBC_T), intent(inout) :: windBC
            type(Model_T), intent(in) :: Model
            real(rp), intent(in) :: t
            real(rp), intent(out) :: Rho,Pr
            real(rp), dimension(NDIM), intent(out) :: V, B
        end subroutine WindTS_T
    end interface

    contains

    subroutine InitWind(bc,Model,Grid,State,xmlInp)
        class(WindBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        !Allocate arrays
        !TODO: Fix the size of these arrays
        call allocGridVec(Model,Grid,bc%BxyzW)
        call allocGridVec(Model,Grid,bc%ExyzW)
        call allocGridVec(Model,Grid,bc%VxyzW)

        call allocGridVar(Model,Grid,bc%RhoW)
        call allocGridVar(Model,Grid,bc%PrW)

        allocate(bc%isWind(bc%NgW,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg))
        bc%isWind = .false.

        !Set time series reference point to Origin by default
        bc%xyzW = 0.0
        call xmlInp%Set_Val(bc%wID,"wind/tsfile","NONE")

        !---------------
        select case (trim(toUpper(bc%wID)))
            case("NONE")
                write(*,*) "No wind TS file specified, relying on user (don't mess this up)"
            case default
                !Set discrete wind function
                bc%isDiscrete = .true.
                bc%getWind => InterpWind
                !Read data into discrete time series
                call readWind(bc,Model,xmlInp)
        end select
        call xmlInp%Set_Val(bc%doDSW,"wind/doDSW",.false.)
        if (bc%doDSW) then
            call xmlInp%Set_Val(bc%dSWAmp,"wind/dSWAmp",0.01_rp)
            call xmlInp%Set_Val(bc%dSWAmp,"wind/dtDSW",bc%dtDSW)
            bc%tDSW = Model%t+bc%dtDSW
            call allocGridVar(Model,Grid,bc%dRhoW)
            bc%dRhoW = 0.0
        endif

        !Zero out initial values
        bc%BxyzW = 0.0
        bc%ExyzW = 0.0
        bc%VxyzW = 0.0
        bc%RhoW = 0.0
        bc%PrW = 0.0

    end subroutine InitWind

    !Set outer shell values for solar wind
    !FIXME
    !For now just doing propagation along x based on constant velocity
    subroutine RefreshWind(windBC,Model,Grid)
        class(WindBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid

        real(rp) :: D,P,DelT
        real(rp), dimension(NDIM) :: nHat,V,B,V0,Vfr,xcc
        real(rp) :: CosF,SinF,CosBp,SinBp
        real(rp) :: bcc,btt, vMag, dSW
        integer :: n,ip,j,k

        !Find current wind info
        call windBC%getWind(windBC,Model,Model%t,D,P,V0,B)

        !Rotate current velocity
        bcc = sqrt(windBC%ByC**2.0 + windBC%BzC**2.0)
        btt = sqrt(1+bcc*2.0)

        CosBp = windBC%ByC/bcc
        SinBp = windBC%BzC/bcc

        CosF = 1.0/btt
        SinF = bcc/btt

        !Now calculate front velocity
        !Simplifying out terms to avoid 0/0
        ! Vfr(XDIR) =  V0(XDIR)*(CosF**2.0)
        ! Vfr(YDIR) = -V0(XDIR)*CosF*SinF*CosBp
        ! Vfr(ZDIR) = -V0(XDIR)*CosF*SinF*SinBp
        Vfr(XDIR) =  V0(XDIR)*(CosF**2.0)
        Vfr(YDIR) = -V0(XDIR)*windBC%ByC/(btt**2.0)
        Vfr(ZDIR) = -V0(XDIR)*windBC%BzC/(btt**2.0)


        vMag = norm2(Vfr)

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,j,k,ip,nHat,xcc,DelT,D,P,V,B,dSW)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                ip = Grid%ie
                nHat = Grid%Tf(ip+1,j,k,NORMX:NORMZ,IDIR) !Outward normal
                do n=1,Model%Ng
                    !Calculate lag time to cell center of this ghost
                    xcc = Grid%xyzcc(ip+n,j,k,:)
                    DelT = dot_product(xcc-windBC%xyzW,Vfr)/vMag**2.0

                    !Set isWind, ie influenced by solar wind or not
                    call getWind(windBC,Model,Model%t-DelT,D,P,V,B)
                    !Note, assuming that Bx is consistent with tilted front
                    !Ie, B(XDIR) = By*ByC + Bz*BzC + Bx0

                    !Check direction
                    if ( dot_product(nHat,V)<=0 ) then
                        windBC%isWind(n,j,k) = .true.
                    else
                        windBC%isWind(n,j,k) = .false.
                    endif
                    if (windBC%doDSW) then
                        dSW = windBC%dRhoW(n,j,k)
                        if (Model%t>=windBC%tDSW) then
                            !Set new perturbation
                            windBC%dRhoW(n,j,k) = genRand(-windBC%dSWAmp,windBC%dSWAmp)
                        endif
                    else
                        dSW = 0.0
                    endif

                    if (windBC%isWind(n,j,k)) then
                        windBC%RhoW(n,j,k) = D*(1.0+dSW)
                        windBC%PrW(n,j,k) = P
                        windBC%VxyzW(n,j,k,:) = V
                        windBC%BxyzW(n,j,k,:) = B
                        windBC%ExyzW(n,j,k,:) = -cross(V,B)
                    else
                        !Do nothing otherwise since values unused
                        !For now just set everything
                        windBC%RhoW(n,j,k) = D*(1.0+dSW)
                        windBC%PrW(n,j,k) = P
                        windBC%VxyzW(n,j,k,:) = V
                        windBC%BxyzW(n,j,k,:) = B
                        windBC%ExyzW(n,j,k,:) = -cross(V,B)

                    endif 
                enddo
            enddo
        enddo
        if (Model%t>=windBC%tDSW) then
            !Set next perturbation time
            windBC%tDSW = windBC%tDSW+windBC%dtDSW
        endif

    end subroutine RefreshWind

    !Do outer I BC for solar wind
    subroutine WindBC(bc,Model,Grid,State)
        class(WindBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: ig,ip,n,j,k,s
        real(rp), dimension(NVAR) :: gW,gCon
        real(rp), dimension(NDIM) :: Bxyz

        !Refresh solar wind shell values
        call RefreshWind(bc,Model,Grid)

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(ig,ip,n,j,k,s) &
        !$OMP private(gW,gCon,Bxyz)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                ip = Grid%ie
                do n=1,Model%Ng
                    ig = Grid%ie+n

                    !Set gCon, conserved variables for SW species
                    if (bc%isWind(n,j,k)) then
                        !Do solar wind values
                        gW(DEN)      = bc%RhoW (n,j,k)
                        gW(PRESSURE) = bc%PrW  (n,j,k)
                        gW(VELX:VELZ)= bc%VxyzW(n,j,k,:)
                        call CellP2C(Model,gW,gCon)
            
                        Bxyz = bc%BxyzW(n,j,k,:)
                    else
                        ! !Use floating BCs from last physical cell
                        ! gCon = State%Gas(ip,j,k,:,BLK)
                        ! Bxyz = CellBxyz(Model,Grid,State%magFlux,ip,j,k)

                        !Testing new option, using SW values for ghost cells but not replacing E field
                        !Do solar wind values
                        gW(DEN)      = bc%RhoW (n,j,k)
                        gW(PRESSURE) = bc%PrW  (n,j,k)
                        gW(VELX:VELZ)= bc%VxyzW(n,j,k,:)
                        call CellP2C(Model,gW,gCon)
            
                        Bxyz = bc%BxyzW(n,j,k,:)

                    endif

                    !FIXME: Properly handle multi-fluid BC for outflow
                    if (Model%doMultiF) then
                        !Use SW or outflow for SW species
                        State%Gas(ig,j,k,:,SWSPC) = gCon
                        do s=SWSPC+1,Model%nSpc
                            !Use zero-grad for all other species
                            State%Gas(ig,j,k,:,s) = State%Gas(ip,j,k,:,s)
                        enddo
                        !Now accumulate
                        call MultiF2Bulk(Model,State%Gas(ig,j,k,:,:))
                    else
                        !Just set bulk
                        State%Gas(ig,j,k,:,BLK) = gCon
                    endif !Multifluid
                    
                    !Set flux/fields based on Bxyz
                    State%Bxyz(ig,j,k,:) = Bxyz
                    ! !NOTE: Using geometry from active grid to ensure smoother stencil
                    State%magFlux(ig+1,j,k,IDIR) = Grid%face(ip+1,j,k,IDIR)*dot_product(Grid%Tf(ip+1,j,k,NORMX:NORMZ,IDIR),Bxyz)
                    State%magFlux(ig  ,j,k,JDIR) = Grid%face(ip  ,j,k,JDIR)*dot_product(Grid%Tf(ip  ,j,k,NORMX:NORMZ,JDIR),Bxyz)
                    State%magFlux(ig  ,j,k,KDIR) = Grid%face(ip  ,j,k,KDIR)*dot_product(Grid%Tf(ip  ,j,k,NORMX:NORMZ,KDIR),Bxyz)

                    ! State%magFlux(ig+1,j,k,IDIR) = Grid%face(ig+1,j,k,IDIR)*dot_product(Grid%Tf(ig+1,j,k,NORMX:NORMZ,IDIR),Bxyz)
                    ! State%magFlux(ig  ,j,k,JDIR) = Grid%face(ig  ,j,k,JDIR)*dot_product(Grid%Tf(ig  ,j,k,NORMX:NORMZ,JDIR),Bxyz)
                    ! State%magFlux(ig  ,j,k,KDIR) = Grid%face(ig  ,j,k,KDIR)*dot_product(Grid%Tf(ig  ,j,k,NORMX:NORMZ,KDIR),Bxyz)

                    
                enddo
            enddo
        enddo

    end subroutine WindBC

    !Fix outer shell electric fields to solar wind values
    subroutine WindEFix(windBC,Model,Grid,State)
        class(windBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k,n
        real(rp), dimension(NDIM) :: Exyz, e1,e2

        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                if (windBC%isWind(1,j,k)) then
                    !Set front-side tangential electric fields to solar wind
                    do i=Grid%ie+1,Grid%ie+1
                        Exyz = windBC%ExyzW(1,j,k,:)
                        do n=JDIR,KDIR
                            !Get edge coordinates
                            call edgeCoords(Model,Grid,i,j,k,n,e1,e2)
                            State%Efld(i,j,k,n) = dot_product(Exyz,e2-e1)
                        enddo

                    enddo

                endif
            enddo
        enddo

    end subroutine WindEFix


    !Read solar wind data from file and initialize WindBC_T (qWind)
    subroutine readWind(windBC,Model,inpXML)
        class(windBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        type(XML_Input_T), intent(in) :: inpXML

        integer :: N
        logical :: fExist
        type(IOVAR_T), dimension(MAXWINDVARS) :: IOVars

        write(*,*) "---------------"
        write(*,*) "Solar wind data"
        write(*,*) "Reading wind data from ", trim(windBC%wID)
        write(*,*) "Assuming input units: t,D,V,P,B = [s],[#/cm3],[m/s],[nPa],[nT]"
        !Check file
        inquire(file=trim(windBC%wID),exist=fExist)
        if (.not. fExist) then
            write(*,*) "Error reading ", trim(windBC%wID), " exiting ..."
            stop
        endif

        !Setup input chain
        call ClearIO(IOVars)
        call AddInVar(IOVars,"T")
        call AddInVar(IOVars,"D")
        call AddInVar(IOVars,"Vx")
        call AddInVar(IOVars,"Vy")
        call AddInVar(IOVars,"Vz")
        call AddInVar(IOVars,"P")
        call AddInVar(IOVars,"Bx")
        call AddInVar(IOVars,"By")
        call AddInVar(IOVars,"Bz")

        !Read data, don't use IO precision
        call ReadVars(IOVars,.false.,windBC%wID)
        N = IOVars(1)%N
        

        !Allocate holders
        windBC%NumT = N
        allocate(windBC%tW(N))
        allocate(windBC%Q(N,NVAR))
        allocate(windBC%B(N,NDIM))

        windBC%tW            = (1/gT0)*IOVars(1)%data
        windBC%Q(:,DEN)      = (1/1.0)*IOVars(2)%data
        windBC%Q(:,VELX)     = (1/gv0)*IOVars(3)%data
        windBC%Q(:,VELY)     = (1/gv0)*IOVars(4)%data
        windBC%Q(:,VELZ)     = (1/gv0)*IOVars(5)%data
        windBC%Q(:,PRESSURE) = (1/gP0)*IOVars(6)%data
        windBC%B(:,XDIR)     = (1/gB0)*IOVars(7)%data
        windBC%B(:,YDIR)     = (1/gB0)*IOVars(8)%data
        windBC%B(:,ZDIR)     = (1/gB0)*IOVars(9)%data

        windBC%tMin = minval(windBC%tW)
        windBC%tMax = maxval(windBC%tW)

        write(*,*) "Finished reading solar wind data"
        write(*,*) "---------------"

    end subroutine readWind

    !Interpolate from qWind data to provide wind BC
    subroutine InterpWind(windBC,Model,t,Rho,Pr,V,B)
        class(WindBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: t
        real(rp), intent(out) :: Rho,Pr
        real(rp), dimension(NDIM), intent(out) :: V, B

        integer :: i0,i1
        real(rp) :: w0,w1,dT
        Rho = 0.0
        Pr = 0.0
        V = 0.0
        B = 0.0

        if (t >= windBC%tMax) then
            i0 = windBC%NumT
            i1 = i0
            w0 = 1.0
            w1 = 0.0

        else if (t <= windBC%tMin) then
            i0 = 1
            i1 = i0
            w0 = 1.0
            w1 = 0.0
        else
            i0 = maxloc(windBC%tW,dim=1,mask=windBC%tW .le. t)
            i1 = i0+1
            dT = windBC%tW(i1)-windBC%tW(i0)
            w0 = (windBC%tW(i1)-t)/dT
            w1 = (t-windBC%tW(i0))/dT
        endif

        Rho = w0*windBC%Q(i0,DEN      ) + w1*windBC%Q(i1,DEN      )
        Pr  = w0*windBC%Q(i0,PRESSURE ) + w1*windBC%Q(i1,PRESSURE )
        V   = w0*windBC%Q(i0,VELX:VELZ) + w1*windBC%Q(i1,VELX:VELZ)
        B   = w0*windBC%B(i0,:        ) + w1*windBC%B(i1,:        )

    end subroutine InterpWind

end module wind
