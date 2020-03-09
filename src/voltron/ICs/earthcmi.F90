!Toy gamera/remix magnetosphere

module uservoltic
    use gamtypes
    use gamdebug
    use gamutils
    use math
    use gridutils
    use xml_input
    use bcs
    use background
    use msphutils
    use wind
    use multifluid

    implicit none
!Earth scaling
!x0 = 6.38e6 m [1 RE]
!rho0 = 1.67e-21 kg/m^3 [1 particle/cc]
!v0 = 10e5 m/s [100 km/s]
!B0 = 4.58 nT
!P0 = 1.67e-11
!t0 = 63.8 second
!M0 = -0.31*1.0e+5/B0

    !Various module variables
    real(rp),private :: Rho0,P0

    !Some knobs for initialization
    logical , private :: doNewIC = .true.
    real(rp), private :: Lc  = 8.0
    real(rp), private :: dLc = 4.0
    real(rp), private :: DInner = 10.0

    ! type for remix BC
    type, extends(baseBC_T) :: IonInnerBC_T

        !Main electric field structures
        real(rp), allocatable, dimension(:,:,:,:) :: inEijk,inExyz

        contains

        procedure :: doInit => InitIonInner
        procedure :: doBC => IonInner

    end type IonInnerBC_T

    contains

    subroutine initUser(Model,Grid,State,inpXML)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: inpXML
        procedure(GasIC_T), pointer :: Wxyz
        procedure(VectorField_T), pointer :: Axyz
        procedure(HackE_T), pointer :: eHack
        procedure(HackStep_T), pointer :: tsHack

        real(rp) :: M0g
        integer :: s

        !Set user hack functions
        !NOTE: Need silly double value for GNU
        eHack  => NULL()
        tsHack => NULL()
        Model%HackE => eHack
        Model%HackStep => tsHack

        !Get defaults from input deck

        !Density for magnetosphere/wind
        call inpXML%Set_Val(Rho0 ,"prob/Rho0",1.0_rp)
        call inpXML%Set_Val(P0   ,"prob/P0"  ,0.001_rp)
        
        !Set magnetosphere parameters
        call setMagsphere(Model,inpXML)
        P0 = P0/gP0 !Scale to magsphere units

        ! deallocate default BCs
        call WipeBCs(Model,Grid)

        !Set BCs (spherical, RPT)
        allocate(IonInnerBC_T       :: Grid%externalBCs(INI )%p)
        allocate(WindBC_T           :: Grid%externalBCs(OUTI)%p)
        allocate(lfmInBC_T          :: Grid%externalBCs(INJ )%p)
        allocate(lfmOutBC_T         :: Grid%externalBCs(OUTJ)%p)
        allocate(periodicInnerKBC_T :: Grid%externalBCs(INK )%p)
        allocate(periodicOuterKBC_T :: Grid%externalBCs(OUTK)%p)

        !Setup fields
        !Use cutoff dipole
        call genCutDipole(Model,Grid,State,inpXML)

        !Map IC to grid
        Wxyz => GasIC
        if (Model%doMultiF) then
            !Set fluid 1
            call GasIC2State(Model,Grid,State,Wxyz,1)
            !Set fluid 2 to plasmasphere
            Wxyz => PSphereIC
            call GasIC2State(Model,Grid,State,Wxyz,2)
            !Null remaining fluids
            do s=3,Model%nSpc
                State%Gas(:,:,:,DEN      ,s) = dFloor
                State%Gas(:,:,:,MOMX:MOMZ,s) = 0.0
                State%Gas(:,:,:,ENERGY   ,s) = pFloor
            enddo
            !Accumulate
            call State2Bulk(Model,Grid,State)
        else
            call GasIC2State(Model,Grid,State,Wxyz)
        endif

        !Set DT bounds
        Grid%isDT = Grid%is
        Grid%ieDT = Grid%ie
        Grid%jsDT = Grid%js
        Grid%jeDT = Grid%je
        Grid%ksDT = Grid%ks
        Grid%keDT = Grid%ke

        !Correction to E (from solar wind or ionosphere)        
        if (Grid%hasLowerBC(IDIR) .or. Grid%hasUpperBC(IDIR)) then
           !Set user hack functions
           !NOTE: Need silly double value for GNU
           eHack  => EFix
           Model%HackE => eHack
           Model%HackPredictor => PredFix
        end if

        !Setup perstep function for everybody
        tsHack => PerStep
        Model%HackStep => tsHack
        
        !Local functions
        !NOTE: Don't put BCs here as they won't be visible after the initialization call

        contains
            subroutine GasIC(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P

                real(rp) :: r,lambda,cL
                real(rp) :: phi,L,Lx,Ly
                real(rp) :: M

                r = sqrt(x**2.0+y**2.0+z**2.0)
                lambda = asin(z/r)
                cL = cos(lambda)
                L = r/(cL*cL)

                phi = atan2(y,x)
                Lx = L*cos(phi)
                Ly = L*sin(phi)

                
                if ( (L <= Lc) .and. doNewIC ) then
                    !P = max( Psk(Lx,Ly),P0 )
                    !P = max( pwolf(L),P0 )/gP0
                    M = RampDown(L,Lc-dLc,Lc)
                    D = max(M*DInner,Rho0)
                    P = P0
                else
                    P = P0
                    D = Rho0
                endif

                Vx = 0.0
                Vy = 0.0
                Vz = 0.0
            end subroutine GasIC

            subroutine PSphereIC(x,y,z,D,Vx,Vy,Vz,P)
                real(rp), intent(in) :: x,y,z
                real(rp), intent(out) :: D,Vx,Vy,Vz,P
                real(rp) :: r,lambda,cL,L

                r = sqrt(x**2.0+y**2.0+z**2.0)
                lambda = asin(z/r)
                cL = cos(lambda)
                L = r/(cL*cL)
                D = max(dFloor,psphD(L))
                P = P0
                Vx = 0.0
                Vy = 0.0
                Vz = 0.0

            end subroutine PSphereIC

    end subroutine initUser

    subroutine postBCInitUser(Model,Grid,State)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State

        SELECT type(pWind=>Grid%externalBCs(OUTI)%p)
            TYPE IS (WindBC_T)
                if (associated(pWind%getWind)) then
                    write(*,*) 'Using solar wind BC from file ...'
                else
                    write(*,*) 'No solar wind file provided/found ...'
                    stop
                endif
            CLASS DEFAULT
                write(*,*) 'Could not find Wind BC in IC'
                stop
        END SELECT

    end subroutine postBCInitUser

    subroutine PerStep(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        type(State_T), intent(inout) :: State

        integer :: i,j,k
        real(rp) :: dF

        !Call ingestion function
        if (Model%doSource) then
            call MagsphereIngest(Model,Gr,State)
        endif

        !Call cooling function/s
        call ChillOut(Model,Gr,State)

        !Do some nudging at the outermost cells to hit solar wind
        SELECT type(pWind=>Gr%externalBCs(OUTI)%p)
            TYPE IS (WindBC_T)
                if (Gr%hasUpperBC(IDIR)) then
                   call NudgeSW(pWind,Model,Gr,State)
                endif
        END SELECT

    end subroutine PerStep

    !Fixes electric field before application
    subroutine EFix(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        type(State_T), intent(inout) :: State

        !call ChkEFieldLFM(Model,Gr,State)
        
        !Fix inner shells
        SELECT type(iiBC=>Gr%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                if (Gr%hasLowerBC(IDIR)) then
                    call IonEFix(Model,Gr,State,iiBC%inEijk)
                endif
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in remix IC'
                stop
        END SELECT

        !Fix outer shells
        SELECT type(pWind=>Gr%externalBCs(OUTI)%p)
            TYPE IS (WindBC_T)
                if (Gr%hasUpperBC(IDIR)) then
                   call WindEFix(pWind,Model,Gr,State)
                end if
            CLASS DEFAULT
                write(*,*) 'Could not find Wind BC in remix IC'
                stop
        END SELECT
        
        !call FixEFieldLFM(Model,Gr,State)             
    end subroutine EFix

    !Fixes cell-centered fields in the predictor

    subroutine PredFix(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        type(State_T), intent(inout) :: State

        !Fix inner shells
        SELECT type(iiBC=>Gr%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                if (Gr%hasLowerBC(IDIR)) then
                    call IonPredFix(Model,Gr,State)
                endif
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in remix IC'
                stop
        END SELECT

        !Fix outer shells
        SELECT type(pWind=>Gr%externalBCs(OUTI)%p)
            TYPE IS (WindBC_T)
                if (Gr%hasUpperBC(IDIR)) then
                   call WindPredFix(pWind,Model,Gr,State)
                end if
            CLASS DEFAULT
                write(*,*) 'Could not find Wind BC in remix IC'
                stop
        END SELECT

    end subroutine PredFix

    !Ensure no flux through degenerate faces
    subroutine IonFlux(Model,Gr,gFlx,mFlx)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(inout) :: gFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(inout), optional :: mFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM,1:NDIM)

        integer :: n,s,j,k
        real(rp) :: igFlx(NVAR), imFlx(NDIM)

        !This is inner-most I tile
        if ( Gr%hasLowerBC(IDIR) .and. (.not. Model%doMultiF) ) then

            if (Model%doRing) then
                if (Model%Ring%doE) then
                    igFlx = sum(gFlx(Gr%is,Gr%je,Gr%ks:Gr%ke,1:NVAR,IDIR,BLK),dim=1)/Model%Ring%Np
                    imFlx = sum(mFlx(Gr%is,Gr%je,Gr%ks:Gr%ke,1:NDIM,IDIR    ),dim=1)/Model%Ring%Np
                    do k=Gr%ks,Gr%ke
                        gFlx(Gr%is,Gr%je,k,:,IDIR,BLK) = igFlx
                        mFlx(Gr%is,Gr%je,k,:,IDIR    ) = imFlx
                    enddo
                endif !doE

                if (Model%Ring%doS) then
                    igFlx = sum(gFlx(Gr%is,Gr%js,Gr%ks:Gr%ke,1:NVAR,IDIR,BLK),dim=1)/Model%Ring%Np
                    imFlx = sum(mFlx(Gr%is,Gr%js,Gr%ks:Gr%ke,1:NDIM,IDIR    ),dim=1)/Model%Ring%Np
                    do k=Gr%ks,Gr%ke
                        gFlx(Gr%is,Gr%js,k,:,IDIR,BLK) = igFlx
                        mFlx(Gr%is,Gr%js,k,:,IDIR    ) = imFlx
                    enddo
                endif !doE

            endif !doRing

            !Now loop over inner sphere (only need active since we're only touching I fluxes)
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(j,k)
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je
                    !Trap for outward mass flux
                    if (gFlx(Gr%is,j,k,DEN,IDIR,BLK) > 0) then
                        gFlx(Gr%is,j,k,DEN   ,IDIR,BLK) = min( 0.0,gFlx(Gr%is,j,k,DEN   ,IDIR,BLK) )
                        gFlx(Gr%is,j,k,ENERGY,IDIR,BLK) = min( 0.0,gFlx(Gr%is,j,k,ENERGY,IDIR,BLK) )
                    endif
                    
                enddo
            enddo !K loop

        endif !Inner i-tile and not MF

    end subroutine IonFlux

    !Initialization for Ion Inner BC
    subroutine InitIonInner(bc,Model,Grid,State,xmlInp)
        class(IonInnerBC_T), intent(inout) :: bc
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(in) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        integer :: PsiShells
        procedure(HackE_T), pointer :: eHack

        !Are we on the inner (REMIX) boundary
        if (Grid%hasLowerBC(IDIR)) then
            call xmlInp%Set_Val(PsiShells,"/remix/grid/PsiShells",5)

            !Create holders for coupling electric field
            allocate(bc%inExyz(1:PsiShells,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM))
            allocate(bc%inEijk(1:PsiShells+1,Grid%jsg:Grid%jeg+1,Grid%ksg:Grid%keg+1,1:NDIM))
            bc%inExyz = 0.0
            bc%inEijk = 0.0
            eHack  => EFix
            Model%HackE => eHack
            Model%HackFlux => IonFlux
        endif
    end subroutine InitIonInner


    !Inner-I BC for ionosphere
    subroutine IonInner(bc,Model,Grid,State)
        class(IonInnerBC_T), intent(inout) :: bc
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        real(rp) :: Rin,llBC,dA
        real(rp), dimension(NDIM) :: Bd,Exyz,Veb,rHat
        integer :: ig,ip,idip,j,k,jp,kp,n,np,d

        !Get inner radius and low-latitude
        Rin = norm2(Grid%xyz(Grid%is,Grid%js,Grid%ks,:))
        llBC = 90.0 - rad2deg*asin(sqrt(Rion/Rin)) !co-lat -> lat


        !$OMP PARALLEL DO default(shared) &
        !$OMP private(ig,ip,idip,j,k,jp,kp,n,np,d) &
        !$OMP private(Bd,Exyz,Veb,rHat,dA)
        do k=Grid%ksg,Grid%keg+1
            do j=Grid%jsg,Grid%jeg+1

                !Loop inward over ghosts
                do n=1,Model%Ng
                    ig = Grid%is-n
                    ip = Grid%is+n-1

                    !Map n=[1,ng] -> inExyz
                    !ASSUMING PsiSt=-3 if you're nudging, so n=[1,ng]->[4,...,1]
                    np = Model%nG-n+1 !Mapping into 1,4 of inExyz

                    !Do cell-centered stuff
                    if (isCellCenterG(Model,Grid,ig,j,k)) then
                        !Map to active ip,jp,kp (i=Grid%is => ip=Grid%is)
                        call lfmIJK(Model,Grid,Grid%is-1,j,k,idip,jp,kp)

                        !Get dipole value
                        Bd = VecDipole(Grid%xyzcc(idip,jp,kp,:))
                        !Get remix field
                        Exyz = bc%inExyz(np,jp,kp,:)

                        !ExB velocity
                        Veb = cross(Exyz,Bd)/dot_product(Bd,Bd)
                        !Remove radial component of velocity
                        rHat = normVec(Grid%xyzcc(idip,jp,kp,:))
                        Veb = Vec2Para(Veb,rHat)

                        !Now do spherical wall BC
                        call SphereWall(Model,State%Gas(ig,j,k,:,:),State%Gas(ip,jp,kp,:,:),Veb)

                    endif !Cell-centered

                !Now do face fluxes
                    !Loop over face directions
                    do d=IDIR,KDIR
                        if ( isLowLat(Grid%xfc(ig,j,k,:,d),llBC) ) then
                            !State%magFlux(ig,j,k,d) = 0.0
                            call lfmIJK(Model,Grid,ig,j,k,ip,jp,kp)
                            dA = Grid%face(ig,j,k,d)/Grid%face(Grid%is,jp,kp,d)
                            State%magFlux(ig,j,k,d) = dA*State%magFlux(Grid%is,jp,kp,d)
                        else
                            call lfmIJK(Model,Grid,ig,j,k,ip,jp,kp)
                            dA = Grid%face(ig,j,k,d)/Grid%face(Grid%is,jp,kp,d)
                            State%magFlux(ig,j,k,d) = dA*State%magFlux(Grid%is,jp,kp,d)
                        endif
                    enddo
                enddo !n loop (ig)
            enddo !j loop
        enddo !k loop

        contains 
            function isLowLat(xyz,llBC)
                real(rp), dimension(NDIM), intent(in) :: xyz
                real(rp), intent(in) :: llBC
                logical :: isLowLat

                real(rp) :: invlat

                invlat = rad2deg*InvLatitude(xyz)
                isLowLat = (invlat<=llBC)
            end function isLowLat

    end subroutine IonInner

    !Correct predictor Bxyz
    subroutine IonPredFix(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: n,ip,ig,ix,jp,kp,j,k

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,ip,ig,ix,j,k)        
        do k=Grid%ksg,Grid%keg
            do j=Grid%js,Grid%je
                do n=1,Model%Ng
                    ip = Grid%is
                    ig = Grid%is-n

                    call lfmIJK(Model,Grid,ig,j,k,ix,jp,kp)
                    State%Bxyz(ig,j,k,:) = State%Bxyz(ip,jp,kp,:)

                enddo !n loop
            enddo !j loop
        enddo !k loop

    end subroutine IonPredFix

    ! !Inner-I BC for ionosphere
    ! subroutine IonInner(bc,Model,Grid,State)
    !     class(IonInnerBC_T), intent(inout) :: bc
    !     type(Model_T), intent(in) :: Model
    !     type(Grid_T), intent(in) :: Grid
    !     type(State_T), intent(inout) :: State

    !     integer :: i,j,k,ip,jp,kp,ig,n,np
    !     logical :: isLL
    !     real(rp) :: rc,xc,yc,zc,Vr,invlat
    !     real(rp) :: xcg,ycg,zcg,xcd,ycd,zcd
    !     real(rp) :: Rin,llBC !Shared
    !     real(rp), dimension(NDIM) :: Exyz,Veb,dB,Bd,rHatG,rHatP,Vxyz,Vmir
    !     real(rp), dimension(NVAR) :: pW,pCon,gW,gCon

    !     !Get inner radius and low-latitude
    !     Rin = norm2(Grid%xyz(Grid%is,Grid%js,Grid%ks,:))
    !     llBC = 90.0 - rad2deg*asin(sqrt(Rion/Rin)) !co-lat -> lat

    !     !write(*,*) 'Rin / llbc = ',Rin,llBC

    !     !i-boundaries (IN)
    !     !$OMP PARALLEL DO default(shared) &
    !     !$OMP private(i,j,k,ip,jp,kp,ig,n,np,isLL) &
    !     !$OMP private(rc,xc,yc,zc,Vr,invlat) &
    !     !$OMP private(xcg,ycg,zcg,xcd,ycd,zcd) &
    !     !$OMP private(Exyz,Veb,Bd,dB,rHatG,rHatP,Vxyz,Vmir) &
    !     !$OMP private(pW,pCon,gW,gCon)
    !     do k=Grid%ksg,Grid%keg
    !         do j=Grid%jsg,Grid%jeg
    !             !Map to active ip,jp,kp (i=Grid%is => ip=Grid%is)
    !             call lfmIJK(Model,Grid,Grid%is,j,k,ip,jp,kp)

    !             !Loop inward over ghosts
    !             do n=1,Model%Ng
    !                 ig = Grid%is-n
    !                 ip = Grid%is+n-1
    !                 !Map n=[1,ng] -> inExyz
    !                 !ASSUMING PsiSt=-3 if you're nudging, so n=[1,ng]->[4,...,1]
    !                 np = Model%nG-n+1 !Mapping into 1,4 of inExyz

    !             !-------
    !             !Get geometry for this ghost and matching physical

    !                 !NOTE: Using j/k instead of jp/kp to deal with double-corner sign flip
    !                 call cellCenter(Grid,ig,j ,k ,xcg,ycg,zcg)
    !                 rHatG = normVec([xcg,ycg,zcg])

    !                 call cellCenter(Grid,ip,jp,kp,xc,yc,zc)
    !                 rHatP = normVec([xc,yc,zc])

    !                 invlat = rad2deg*InvLatitude([xcg,ycg,zcg]) !Convert to degrees

    !                 if (invlat<=llBC) then
    !                     isLL = .true.
    !                 else
    !                     isLL = .false.
    !                 endif

    !             !-------
    !             !Get E, dipole/perturbation and ExB velocity/Mirror velocity

    !                 !Get velocity from i-reflected active cell
    !                 Vmir = State%Gas(ip,jp,kp,MOMX:MOMZ,BLK)/max(State%Gas(ip,jp,kp,DEN,BLK),dFloor)
    !                 Exyz = bc%inExyz(np,jp,kp,:)

    !                 !Choose which dipole ExB speed to use, true ghost value is much faster
    !                 xcd = Grid%xyzcc(Grid%is-1,jp,kp,XDIR)
    !                 ycd = Grid%xyzcc(Grid%is-1,jp,kp,YDIR)
    !                 zcd = Grid%xyzcc(Grid%is-1,jp,kp,ZDIR)
    !                 call Dipole(xcd,ycd,zcd,Bd(XDIR),Bd(YDIR),Bd(ZDIR))
                    
    !                 dB = State%Bxyz(ip,jp,kp,:)
    !                 !ExB velocity
    !                 Veb = cross(Exyz,Bd)/dot_product(Bd,Bd)

    !                 !Use ExB (w/o radial) and mirror
    !                 Vxyz = Veb - rHatP*dot_product(rHatP,Veb) !- rHatP*dot_product(rHatP,Vmir)
    !             !-------
    !             !Set ghost hydro quantities
    !                 !Let density float
    !                 call SphereWall(Model,State%Gas(ig,j,k,:,:),State%Gas(ip,jp,kp,:,:),Vxyz)
    !                 !Now do polar outflow if testing
    !                 if (Model%doMultiF .and. (invlat>=70) .and. (Model%nSpc>2)) then
    !                     gW(DEN) = 100.0
    !                     gW(VELX:VELZ) = 0.2*rHatP + Veb - rHatP*dot_product(rHatP,Veb)
    !                     gW(PRESSURE) = 1.0e-3
    !                     call CellP2C(Model,gW,gCon)
    !                     State%Gas(ig,j,k,:,3) = gCon
    !                     !Reset bulk
    !                     call MultiF2Bulk(Model,State%Gas(ig,j,k,:,:))
    !                 endif

    !             !-------
    !             !Now handle magnetic quantities (perturbation field)

    !                 !Using rHatP to hold coefficients for singularity geometry
    !                 rHatP = 1.0
    !                 if ( (Model%Ring%doS) .and. (j < Grid%js) ) then
    !                     rHatP(JDIR:KDIR) = -1.0
    !                 endif
    !                 if ( (Model%Ring%doE) .and. (j >= Grid%je+1) ) then
    !                     rHatP(JDIR:KDIR) = -1.0
    !                 endif

    !                 State%magFlux(ig,j,k,IDIR) = rHatP(IDIR)*Grid%face(ig,j,k,IDIR)*dot_product(dB,Grid%Tf(ig,j,k,NORMX:NORMZ,IDIR))
    !                 State%magFlux(ig,j,k,JDIR) = rHatP(JDIR)*Grid%face(ig,j,k,JDIR)*dot_product(dB,Grid%Tf(ig,j,k,NORMX:NORMZ,JDIR))
    !                 State%magFlux(ig,j,k,KDIR) = rHatP(KDIR)*Grid%face(ig,j,k,KDIR)*dot_product(dB,Grid%Tf(ig,j,k,NORMX:NORMZ,KDIR))
    !                 if (j == Grid%jeg) then
    !                     State%magFlux(ig,j+1,k,IDIR:KDIR) = 0.0
    !                     State%magFlux(ig,j+1,k,JDIR)      = Grid%face(ig,j+1,k,JDIR)*dot_product(dB,Grid%Tf(ig,j+1,k,NORMX:NORMZ,JDIR))
    !                 endif

    !             enddo !n
    !         enddo
    !     enddo

    !     !NOTE: Currently just assuming we have full K
    !     if ( Grid%hasLowerBC(KDIR) .and. Grid%hasUpperBC(KDIR) ) then
    !         State%magFlux(Grid%is-Model%Ng:Grid%is-1,:,Grid%ke+1,IDIR:KDIR) = 0.0
    !         State%magFlux(Grid%is-Model%Ng:Grid%is-1,:,Grid%ke+1,KDIR) = State%magFlux(Grid%is-Model%Ng:Grid%is-1,:,Grid%ks,KDIR)
    !     else
    !         write(*,*) 'IonInner not implemented for K decomposition'
    !         stop
    !     endif                  

    ! end subroutine IonInner
    
end module uservoltic
