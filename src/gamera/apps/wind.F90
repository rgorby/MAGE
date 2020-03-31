!Various utilities to handle solar wind BC

module wind
    use gamtypes
    use gamutils
    use math
    use xml_input
    use ioH5
    use gridutils
    use msphutils
    use multifluid
    use bcs

    implicit none

    ! Global Parameters
    integer, parameter :: SWSPC = 1 !SW fluid is always 1st in multifluid

    logical :: doWindInterp = .false.

    !TODO: Remove WindTS_T pointer and call interpwind

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
        real(rp) :: Bx0 = 0.0, ByC = 0.0, BzC = 0.0

        integer :: NgW = 4 !Number of solar wind ghost cells
        real(rp), dimension(NDIM) :: xyzW !Position of solar wind time-series (ie L1)
        real(rp), dimension(NDIM) :: vFr !Front velocity @ current time

        contains

        procedure :: doInit => InitWind
        procedure :: doBC => WindBC

    end type WindBC_T


    !Type for solar wind time-series subroutine
    abstract interface
        subroutine WindTS_T(windBC,Model,t,Rho,Pr,V,B)
            Import :: rp,Model_T,NDIM,WindBC_T
            class(WindBC_T), intent(in) :: windBC
            type(Model_T), intent(in) :: Model
            real(rp), intent(in) :: t
            real(rp), intent(out) :: Rho,Pr
            real(rp), dimension(NDIM), intent(out) :: V, B
        end subroutine WindTS_T
    end interface

    contains

    subroutine InitWind(bc,Model,Grid,State,xmlInp)
        class(WindBC_T), intent(inout) :: bc
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        !Set time series reference point to Origin by default
        bc%xyzW = 0.0
        bc%vFr = [-1.0,0.0,0.0] !Initial front velocity

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

    end subroutine InitWind

    !Set outer shell values for solar wind
    !FIXME
    !For now just doing propagation along x based on constant velocity
    subroutine RefreshWind(windBC,Model,Grid)
        class(WindBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid

        real(rp) :: D,P,DelT
        real(rp), dimension(NDIM) :: nHat,vHat,V,B,V0,xcc
        real(rp) :: CosF,SinF,CosBp,SinBp
        real(rp) :: bcc,btt, vMag
        integer :: n,ip,j,k

        !Find current wind info
        call windBC%getWind(windBC,Model,Model%t,D,P,V0,B)

        !Rotate current velocity
        bcc = sqrt(windBC%ByC**2.0 + windBC%BzC**2.0)
        btt = sqrt(1+bcc**2.0)

        !Guard against div by 0
        if (bcc > TINY) then
            CosBp = windBC%ByC/bcc
            SinBp = windBC%BzC/bcc

            CosF = 1.0/btt
            SinF = bcc/btt

            !Calculate front velocity
            windBC%vFr(XDIR) =  V0(XDIR)*(CosF**2.0)
            windBC%vFr(YDIR) = -V0(XDIR)*CosF*SinF*CosBp
            windBC%vFr(ZDIR) = -V0(XDIR)*CosF*SinF*SinBp
        else
            !ByC and BzC are both nearly zero
            windBC%vFr(XDIR) = V0(XDIR)
            windBC%vFr(YDIR) = 0.0
            windBC%vFr(ZDIR) = 0.0
        endif

        vMag = norm2(windBC%vFr)

    end subroutine RefreshWind

    !Get solar wind at point in space-time
    subroutine GetWindAt(windBC,Model,xyz,t,Rho,Pr,V,B)
        class(WindBC_T), intent(in) :: windBC
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: t
        real(rp), intent(in) :: xyz(NDIM)
        real(rp), intent(out) :: Rho,Pr
        real(rp), dimension(NDIM), intent(out) :: V, B

        real(rp) :: vMag,DelT
        vMag = norm2(windBC%vFr)

        !Calculate lag time to this cell
        DelT = dot_product(xyz-windBC%xyzW,windBC%vFr)/vMag**2.0
        !Rewind by lag time
        call windBC%getWind(windBC,Model,t-DelT,Rho,Pr,V,B)

    end subroutine GetWindAt

    !Given face normal, decide how important solar wind is [0,1]
    !Given point, decide how important solar wind is [0,1]
    function wgtWind(windBC,Model,xyz,t) result(wgt)
        class(WindBC_T), intent(in) :: windBC
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: t
        real(rp), intent(in) :: xyz(NDIM)
        real(rp) :: wgt

        real(rp) :: Rho,P,dTh
        real(rp), dimension(NDIM) :: V,B,vHat,nHat

        !Get radial vector to this point
        nHat = normVec(xyz)
        !Start by getting wind
        call GetWindAt(windBC,Model,xyz,t,Rho,P,V,B)
        vHat = normVec(V) !Normalized velocity direction

        !Get angle between -vHat and nHat
        dTh = acos(dot_product(-vHat,nHat))*180.0/PI !Degrees
        !0 <= dTh <= 180
        !wgt = 1, dTh<=90
        !wgt->0 as dTh=>180
        wgt = RampDown(dTh,90.0_rp,90.0_rp)
        
    end function wgtWind

    !Do outer I BC for solar wind
    subroutine WindBC(bc,Model,Grid,State)
        class(WindBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        real(rp) :: t,D,P,wSW,wMHD,swFlx,inFlx,Cs
        integer :: ip,ig,j,k,n,s
        real(rp), dimension(NDIM) :: xcc,V,B,nHat
        real(rp), dimension(NVAR) :: gW,gW_sw,gW_in,gCon
        
        if (.not. Grid%hasUpperBC(IDIR)) return

        !Refresh solar wind shell values
        call RefreshWind(bc,Model,Grid)
        
        t = Model%t

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(ip,ig,j,k,n,s) &
        !$OMP private(D,P,wSW,wMHD,swFlx,inFlx) &
        !$OMP private(xcc,V,B,nHat,gW,gW_sw,gW_in,gCon,Cs)
        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1
                ip = Grid%ie
                do n=1,Model%Ng
                    ig = Grid%ie+n

                    gW_sw = 0.0
                    gW_in = 0.0

                    !Do cell-centered stuff
                    if (isCellCenterG(Model,Grid,ig,j,k)) then
                    
                        !Get weight from cell center
                        xcc = Grid%xyzcc(ig,j,k,:)
                        wSW = wgtWind(bc,Model,xcc,t)
                        wMHD = 1.0-wSW
                        !Do SW cell centered
                        if (wSW > TINY) then
                            call GetWindAt(bc,Model,xcc,t,D,P,V,B)
                            gW_sw(DEN) = D
                            gW_sw(PRESSURE) = P
                            gW_sw(VELX:VELZ) = V

                        endif !SW inflow

                        !Do MHD outflow cell-centered
                        if (wMHD>TINY) then
                            !Get values from last physical
                            gCon = State%Gas(ip,j,k,:,BLK)
                            call CellC2P(Model,gCon,gW_in)
                            call CellPress2Cs(Model,gCon,Cs)

                            !Do some work on the flow velocity
                            V = gW_in(VELX:VELZ)
                            nHat = Grid%Tf(ip+1,j,k,NORMX:NORMZ,IDIR)

                            !Enforce diode, should be leaving at Mach 1
                            if (dot_product(V,nHat) < Cs) then
                                V = V - Vec2Para(V,nHat) + Cs*nHat
                            endif
                            gW_in(VELX:VELZ) = V
                        endif !MHD outflow

                        !Mix BCs and set final ghost values
                        gW = wSW*gW_sw + wMHD*gW_in
                        call CellP2C(Model,gW,gCon)

                        !FIXME: Handle multi-fluid better
                        if (Model%doMultiF) then
                            !Use calculated values for SW fluid
                            State%Gas(ig,j,k,:,SWSPC) = gCon
                            do s=SWSPC+1,Model%nSpc
                                !Lazy zero grad for rest
                                State%Gas(ig,j,k,:,s) = State%Gas(ip,j,k,:,s)
                            enddo !Fluid loop
                            !Now accumulate
                            call MultiF2Bulk(Model,State%Gas(ig,j,k,:,:))
                        else
                            !Just set bulk
                            State%Gas(ig,j,k,:,BLK) = gCon
                        endif !Multifluid

                    endif !Cell-centered ghosts

                !Now do fluxes
                    !I flux (start at +1)
                    wSW = wgtWind(bc,Model,Grid%xfc(ig+1,j,k,:,IDIR),t)
                    swFlx = WindMagFlux(bc,Model,Grid,t,IDIR,ig+1,j,k)
                    inFlx = State%magFlux(ip+1,j,k,IDIR)
                    State%magFlux(ig+1,j,k,IDIR) = wSW*swFlx + (1-wSW)*inFlx
                    
                    !J flux
                    wSW = wgtWind(bc,Model,Grid%xfc(ig,j,k,:,JDIR),t)
                    swFlx = WindMagFlux(bc,Model,Grid,t,JDIR,ig,j,k)
                    inFlx = State%magFlux(ip,j,k,JDIR)
                    State%magFlux(ig,j,k,JDIR) = wSW*swFlx + (1-wSW)*inFlx
                    

                    !K flux
                    wSW = wgtWind(bc,Model,Grid%xfc(ig,j,k,:,KDIR),t)
                    swFlx = WindMagFlux(bc,Model,Grid,t,KDIR,ig,j,k)
                    inFlx = State%magFlux(ip,j,k,KDIR)
                    State%magFlux(ig,j,k,KDIR) = wSW*swFlx + (1-wSW)*inFlx
                    
                enddo !n
            enddo !j
        enddo !k
        
    end subroutine WindBC


    !Calculate solar wind magnetic flux at ijkdir-face of cell i,j,k
    function WindMagFlux(bc,Model,Grid,t,ijkdir,i,j,k) result(swFlux)
        class(windBC_T), intent(in) :: bc
        type(Model_T)  , intent(in) :: Model
        type(Grid_T)   , intent(in) :: Grid
        real(rp)       , intent(in) :: t
        integer        , intent(in) :: i,j,k,ijkdir

        real(rp) :: swFlux,D,P
        real(rp), dimension(NDIM) :: xyz,V,B
        xyz = Grid%xfc(i,j,k,:,ijkdir)
        call GetWindAt(bc,Model,xyz,t,D,P,V,B)
        swFlux = Project2Face(Model,Grid,B,ijkdir,i,j,k)

    end function WindMagFlux

    !Correct predictor Bxyz
    subroutine WindPredFix(bc,Model,Grid,State) 
        class(windBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: ig,j,k,n,ip
        real(rp), dimension(NDIM) :: xcc,V,swB,inB
        real(rp) :: t,D,P,wSW

        if (.not. Grid%hasUpperBC(IDIR)) return

        !Refresh solar wind shell values
        call RefreshWind(bc,Model,Grid)
        
        t = Model%t

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(j,k,n,ig,ip,xcc,D,P,V,wSW,swB,inB)
        do k=Grid%ksg,Grid%keg
            do j=Grid%js,Grid%je
                ip = Grid%ie
                do n=1,Model%Ng
                    ig = Grid%ie+n
                    xcc = Grid%xyzcc(ig,j,k,:)
                    wSW = wgtWind(bc,Model,xcc,t)

                    call GetWindAt(bc,Model,xcc,t,D,P,V,swB)
                    inB = State%Bxyz(ip,j,k,:)

                    State%Bxyz(ig,j,k,:) = wSW*swB + (1-wSW)*inB
                enddo !n loop
            enddo !j loop
        enddo !k loop

    end subroutine WindPredFix

    !Nudge outer-most physical cell
    subroutine NudgeSW(windBC,Model,Grid,State)
        class(windBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k,ip
        real(rp) :: t,dtSW0,dtSW,wSW,D,P
        real(rp), dimension(NDIM) :: xcc,V,B
        real(rp), dimension(NVAR) :: gCon,gW,gW_sw,gW_mhd

        if (.not. Grid%hasUpperBC(IDIR)) return

        !Refresh solar wind shell values
        call RefreshWind(windBC,Model,Grid)

        dtSW0 = 60.0/Model%Units%gT0 !One minute from SW time series

        ip = Grid%ie !Only doing outer-most cell
        t = Model%t

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(j,k,D,P,wSW,dtSW) &
        !$OMP private(xcc,V,B,gW,gCon,gW_mhd,gW_sw)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je

                xcc = Grid%xyzcc(ip,j,k,:)
                wSW = wgtWind(windBC,Model,xcc,t) !SW weight of last physical cell
                if (wSW>TINY) then
                    !Use nudge timescale (dtSW0) over weight
                    dtSW = max(dtSW0/wSW,Model%dt)

                !Get solar wind state in this cell
                    call GetWindAt(windBC,Model,xcc,t,D,P,V,B)
                    gW_sw(DEN) = D
                    gW_sw(PRESSURE) = P
                    gW_sw(VELX:VELZ) = V
                !Get MHD state in this cell
                    !For multifluid, only nudging SW fluid
                    if (Model%doMultiF) then
                        gCon = State%Gas(ip,j,k,:,SWSPC)
                    else
                        gCon = State%Gas(ip,j,k,:,BLK)
                    endif

                    call CellC2P(Model,gCon,gW_mhd)

                !Nudge primitive MHD state to primitive SW state
                    gW = gW_mhd + (Model%dt/dtSW)*( gW_sw - gW_mhd )
                    !Flip back to conserved and put back in physical cell
                    call CellP2C(Model,gW,gCon)
                    if (Model%doMultiF) then
                        State%Gas(ip,j,k,:,SWSPC) = gCon
                        call MultiF2Bulk(Model,State%Gas(ip,j,k,:,:))
                    else
                        State%Gas(ip,j,k,:,BLK) = gCon
                    endif !Multifluid

                endif !wSW nudge

            enddo
        enddo
        
    end subroutine NudgeSW

    !Fix outer shell electric fields to solar wind values
    subroutine WindEFix(windBC,Model,Grid,State)
        class(windBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k,n
        real(rp), dimension(NDIM) :: xcc,nHat,e1,e2,ecc,Vxyz,Bxyz,swExyz
        real(rp) :: wSW,D,P,mhdExyz,t

        if (.not. Grid%hasUpperBC(IDIR)) return

        !Refresh solar wind shell values
        call RefreshWind(windBC,Model,Grid)

        t = Model%t
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,n,e1,e2,ecc,Vxyz,Bxyz,swExyz,wSW,D,P,mhdExyz)
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=Grid%ie-2,Grid%ie+1

                    do n=IDIR,KDIR
                        !Get edge coordinates
                        call edgeCoords(Model,Grid,i,j,k,n,e1,e2)
                        ecc = 0.5*(e1+e2)

                        !Get weight and wind at edge center
                        wSW = wgtWind(windBC,Model,ecc,t)
                        if (wSW > TINY) then
                            call GetWindAt(windBC,Model,ecc,t,D,P,Vxyz,Bxyz)
                            if (i <= Grid%ie) then
                                !Use full weight only for outer-most shell
                                wSW = wSW/(1.0 + Grid%ie+1 - i)
                            endif
                            !Calculate solar wind E field (not yet projected to edge)
                            swExyz  = -cross(Vxyz,Bxyz)
                            !Get MHD EMF
                            mhdExyz = State%Efld(i,j,k,n)
                            State%Efld(i,j,k,n) = (1.0-wSW)*mhdExyz + wSW*dot_product(swExyz,e2-e1)
                        endif
                    enddo !Edge loop

                    !Add diffusive electric field
                    if (i == Grid%ie+1) then
                        swExyz = DiffuseOuter(Model,Grid,State,i,j,k)
                        State%Efld(i,j,k,:) = State%Efld(i,j,k,:) + swExyz
                    endif

                enddo !i cells
            enddo
        enddo

    end subroutine WindEFix


    !Calculate diffusive electric field
    function DiffuseOuter(Model,Grid,State,i,j,k) result(Ed)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State
        integer, intent(in) :: i,j,k
        real(rp), dimension(NDIM) :: Ed,Jd
        real(rp) :: Vd,db2,db1,dl

        Ed = 0.0
        Jd = 0.0
        !Calculate current
        !Jk = d_i (Bj) - d_j (Bi), db2 - db1 (see fields.F90)
        
        db2 = State%magFlux(i,j,k,JDIR)/Grid%face(i,j,k,JDIR) - State%magFlux(i-1,j,k,JDIR)/Grid%face(i-1,j,k,JDIR)
        db1 = State%magFlux(i,j,k,IDIR)/Grid%face(i,j,k,IDIR) - State%magFlux(i,j-1,k,IDIR)/Grid%face(i,j-1,k,IDIR)
        Jd(KDIR) = db2 - db1

        !Jj = d_k (Bi) - d_i (Bk)
        db2 = State%magFlux(i,j,k,IDIR)/Grid%face(i,j,k,IDIR) - State%magFlux(i,j,k-1,IDIR)/Grid%face(i,j,k-1,IDIR)
        db1 = State%magFlux(i,j,k,KDIR)/Grid%face(i,j,k,KDIR) - State%magFlux(i-1,j,k,KDIR)/Grid%face(i-1,j,k,KDIR)
        Jd(JDIR) = db2 - db1

        Vd = Model%Ca
        dl = Grid%volume(i,j,k)**(1.0/3.0)
        Vd = min(Vd,Model%CFL*dl/Model%dt)

        Ed(IDIR) = 0.0
        Ed(JDIR) = Vd*Jd(JDIR)*Grid%edge(i,j,k,JDIR)
        Ed(KDIR) = Vd*Jd(KDIR)*Grid%edge(i,j,k,KDIR)

    end function DiffuseOuter

    !Read solar wind data from file and initialize WindBC_T (qWind)
    subroutine readWind(windBC,Model,inpXML)
        class(windBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        type(XML_Input_T), intent(in) :: inpXML

        integer :: i,N
        logical :: isByC,isBzC
        real(rp) :: BCoef(3)
        real(rp), parameter :: ergcc2nPa = 1.0e8

        type(IOVAR_T), dimension(MAXWINDVARS) :: IOVars

        if (Model%isLoud) then
            write(*,*) "---------------"
            write(*,*) "Solar wind data"
            write(*,*) "Reading wind data from ", trim(windBC%wID)
            write(*,*) "Assuming input units: t,D,V,T,B = [s],[#/cm3],[m/s],[K],[nT]"
        endif        
        !Make sure file exists
        call CheckFileOrDie(windBC%wID, "Error opening wind file, exiting ...")

        if (.not.(ioExist(trim(windBC%wID),"Temp"))) then
           write(*,*) 'As of 5 October 2019 solar wind temperature, rather than thermal pressure,'
           write(*,*) 'is stored in the solar wind h5 file by the omni2wind script.'
           write(*,*) 'The solar wind file used in this run does not have the "Temp" variable. Quitting...'
           stop
        endif

        !Setup input chain
        call ClearIO(IOVars)
        call AddInVar(IOVars,"T")
        call AddInVar(IOVars,"D")
        call AddInVar(IOVars,"Vx")
        call AddInVar(IOVars,"Vy")
        call AddInVar(IOVars,"Vz")
        call AddInVar(IOVars,"Temp")
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
        ! compute pressure from density and temperature
        ! note, assuming density in /cc and temperature in K
        ! Kbltz is defined in kdefs in erg/K, so convert to nPa
        windBC%Q(:,PRESSURE) = (1/gP0)*windBC%Q(:,DEN)*Kbltz*IOVars(6)%data*ergcc2nPa
        windBC%B(:,XDIR)     = (1/gB0)*IOVars(7)%data
        windBC%B(:,YDIR)     = (1/gB0)*IOVars(8)%data
        windBC%B(:,ZDIR)     = (1/gB0)*IOVars(9)%data

        windBC%tMin = minval(windBC%tW)
        windBC%tMax = maxval(windBC%tW)

    !Now go back in to get coefficients
        !Try to grab all three from file
        call ClearIO(IOVars)
        call AddInVar(IOVars,"ByC",vTypeO=IOREAL)
        call AddInVar(IOVars,"BzC",vTypeO=IOREAL)
        call AddInVar(IOVars,"Bx0",vTypeO=IOREAL)
        !Read data, don't use IO precision
        call ReadVars(IOVars,.false.,windBC%wID)

        BCoef(1:3) = 0.0
        do i=1,3
            if (IOVars(i)%isDone) then
                !This coefficient is present, so grab it
                BCoef(i) = IOVars(i)%data(1)
            endif
        enddo

        !Don't need to scale coefficients, but need to convert Bx0 to code units
        windBC%ByC = BCoef(1)
        windBC%BzC = BCoef(2)
        windBC%Bx0 = BCoef(3)*(1/gB0)

        !Now redo Bx to be consistent with tilted front
        !Bx = Bx0 + ByC*By + BzC*Bz

        windBC%B(:,XDIR) = windBC%Bx0 + windBC%ByC*windBC%B(:,YDIR) + windBC%BzC*windBC%B(:,ZDIR) 

        if (Model%isLoud) then
            write(*,'(a,3f8.3)') ' SW Coefficients (Bx0,ByC,BzC) = ', BCoef(3),BCoef(1),BCoef(2)

            write(*,*) "Finished reading solar wind data"
            write(*,*) "---------------"
        endif        
    end subroutine readWind

    !Interpolate from qWind data to provide wind BC
    subroutine InterpWind(windBC,Model,t,Rho,Pr,V,B)
        class(WindBC_T), intent(in) :: windBC
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

        if (.not. doWindInterp) then
            !Use discrete walls
            w0 = 1.0
            w1 = 0.0
        endif

        Rho = w0*windBC%Q(i0,DEN      ) + w1*windBC%Q(i1,DEN      )
        Pr  = w0*windBC%Q(i0,PRESSURE ) + w1*windBC%Q(i1,PRESSURE )
        V   = w0*windBC%Q(i0,VELX:VELZ) + w1*windBC%Q(i1,VELX:VELZ)
        B   = w0*windBC%B(i0,:        ) + w1*windBC%B(i1,:        )

        !Replace Bx w/ coefficient expansion
        B(XDIR) = windBC%Bx0 + windBC%ByC*B(YDIR) + windBC%BzC*B(ZDIR)
        if (t <= windBC%tMin) then
            B(:) = 0.0
        endif
    end subroutine InterpWind

end module wind
