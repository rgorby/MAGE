!Toy gamera/remix magnetosphere

module uservoltic
    use gamtypes
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

    !Various global would go here
    real(rp) :: Rho0,P0
    real(rp) :: RhoW0,PrW0,VxW,BzW

    !Running values for BCs
    real(rp) :: T0  = 60.0
    real(rp) :: dCS = 0.0

    !doHeavy = Add plasmasphere
    logical :: doHeavy = .false.

    !doCool = Apply cooling function
    logical :: doCool       = .true.
    
    logical :: newMix = .false.


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
        call inpXML%Set_Val(RhoW0,"prob/RhoW",5.0_rp)
        call inpXML%Set_Val(P0   ,"prob/P0"  ,0.001_rp)
        call inpXML%Set_Val(PrW0 ,"prob/PrW" ,0.48_rp)


        !Use plasmasphere model for initial density
        call inpXML%Set_Val(doHeavy,"prob/doHeavy",doHeavy)

        !Set magnetosphere parameters
        call setMagsphere(Model,inpXML)


        !Solar wind
        call inpXML%Set_Val(VxW,"prob/Vx",4.0_rp)
        VxW = abs(VxW) !Assume positive value

        
        call inpXML%Set_Val(BzW,"prob/BzW",0.0_rp)
        
        ! deallocate default BCs
        deallocate(Grid%ExternalBCs(INI )%p)
        deallocate(Grid%ExternalBCs(OUTI)%p)
        deallocate(Grid%ExternalBCs(INJ )%p)
        deallocate(Grid%ExternalBCs(OUTJ)%p)
        deallocate(Grid%ExternalBCs(INK )%p)
        deallocate(Grid%ExternalBCs(OUTK)%p)

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

        !Set MG bounds
        Grid%isMG = Grid%is
        Grid%ieMG = Grid%ie
        Grid%jsMG = Grid%js
        Grid%jeMG = Grid%je
        Grid%ksMG = Grid%ks
        Grid%keMG = Grid%ke

        !Correction to E (from solar wind or ionosphere)        
        if ( (Model%Ri == Model%NumRi) .or. (Model%Ri == 1) ) then
           !Set user hack functions
           !NOTE: Need silly double value for GNU
           eHack  => EFix
           Model%HackE => eHack
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

                real(rp) :: r,lambda,cL,L

                r = sqrt(x**2.0+y**2.0+z**2.0)
                lambda = asin(z/r)
                cL = cos(lambda)
                L = r/(cL*cL)
                if (doHeavy) then
                    D = max(Rho0,psphD(L))
                else
                    D = Rho0
                endif

                P = P0
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
                    write(*,*) 'Using solar wind BC from subroutine ...'
                    pWind%getWind => SolarWindTS
                endif
            CLASS DEFAULT
                write(*,*) 'Could not find Wind BC in remix IC'
                stop
        END SELECT

    end subroutine postBCInitUser

    subroutine PerStep(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        type(State_T), intent(inout) :: State

        integer :: i

        !Call ingestion function
        if (Model%doSource) then
            call MagsphereIngest(Model,Gr,State)
        endif

        !Call cooling function/s
        if (doCool) call ChillOut(Model,Gr,State)
        
    end subroutine PerStep

    !Fixes electric field before application
    subroutine EFix(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        type(State_T), intent(inout) :: State

        integer :: i,j,k,kp
        real(rp) :: MaxEjp,MaxEjm,Ei,Ej,Ek

        !Fix inner shells
        SELECT type(iiBC=>Gr%externalBCs(INI)%p)
            TYPE IS (IonInnerBC_T)
                if (Model%Ri == 1) then
                    call IonEFix(Model,Gr,State,iiBC%inEijk)
                endif
            CLASS DEFAULT
                write(*,*) 'Could not find Ion Inner BC in remix IC'
                stop
        END SELECT

        !Fix outer shells
        SELECT type(pWind=>Gr%externalBCs(OUTI)%p)
            TYPE IS (WindBC_T)
                if (Model%Ri == Model%NumRi) then
                   call WindEFix(pWind,Model,Gr,State)
                end if
            CLASS DEFAULT
                write(*,*) 'Could not find Wind BC in remix IC'
                stop
        END SELECT

    end subroutine EFix

    !Ensure no flux through degenerate faces
    subroutine IonFlux(Model,Gr,gFlx,mFlx)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(inout) :: gFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NVAR,1:NDIM,BLK:Model%nSpc)
        real(rp), intent(inout), optional :: mFlx(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM,1:NDIM)

        integer :: n,s,j,k
        real(rp) :: igFlx(NVAR,2), imFlx(NDIM,2)
        real(rp) :: Rin,llBC,invlat

        !Ignore everything else
        return

        !This is inner-most I tile
        if ( (Model%Ri == 1) .and. (.not. Model%doMultiF) ) then
            !Get inner radius and low-latitude
            Rin = norm2(Gr%xyz(Gr%is,Gr%js,Gr%ks,:))
            llBC = 90.0 - rad2deg*asin(sqrt(Rion/Rin)) !co-lat -> lat

            !Now loop over inner sphere (only need active since we're only touching I fluxes)
            !$OMP PARALLEL DO default(shared) &
            !$OMP private(j,k,invlat)
            do k=Gr%ks,Gr%ke
                do j=Gr%js,Gr%je

                    !Only inward (negative) mass flux
                    gFlx(Gr%is,j,k,DEN,IDIR,BLK) = min( 0.0,gFlx(Gr%is,j,k,DEN,IDIR,BLK) )

                enddo
            enddo !K loop

        endif !Inner i-tile and not MF

    end subroutine IonFlux

    !Put BCs here for global access
    !Solar wind values
    subroutine SolarWindTS(windBC,Model,t,Rho,Pr,V,B)
        class(WindBC_T), intent(inout) :: windBC
        type(Model_T), intent(in) :: Model
        real(rp), intent(in) :: t
        real(rp), intent(out) :: Rho,Pr
        real(rp), dimension(NDIM), intent(out) :: V, B

        integer :: imfNS
        real(rp) :: vScl

        if (t <= T0) then
            imfNS = 0.0
        else if (t <= 3*T0) then
            imfNS = -1
        else if (t <= 6*T0) then
            imfNS =  1
        else
            imfNS = -1
        endif

        Rho = RhoW0
        Pr = PrW0

        V = 0
        B = 0
        vScl = 1.0
        B(ZDIR) = imfNS*BzW
        V(XDIR) = -vScl*VxW


    end subroutine SolarWindTS

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
        if (Model%Ri == 1) then
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

        integer :: i,j,k,ip,jp,kp,ig,n,np
        logical :: isLL
        real(rp) :: rc,xc,yc,zc,Vr,invlat
        real(rp) :: Rin,llBC !Shared
        real(rp), dimension(NDIM) :: Exyz,Veb,dB,Bd,rHatG,rHatP,Vxyz,Vmir
        real(rp), dimension(NVAR) :: pW,pCon,gW,gCon

        !Get inner radius and low-latitude
        Rin = norm2(Grid%xyz(Grid%is,Grid%js,Grid%ks,:))
        llBC = 90.0 - rad2deg*asin(sqrt(Rion/Rin)) !co-lat -> lat

        !write(*,*) 'Rin / llbc = ',Rin,llBC

        !i-boundaries (IN)
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,ip,jp,kp,ig,n,np,isLL) &
        !$OMP private(rc,xc,yc,zc,Vr,invlat) &
        !$OMP private(Exyz,Veb,Bd,dB,rHatG,rHatP,Vxyz,Vmir) &
        !$OMP private(pW,pCon,gW,gCon)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                !Map to active ip,jp,kp (i=Grid%is => ip=Grid%is)
                call lfmIJK(Model,Grid,Grid%is,j,k,ip,jp,kp)

                !Loop inward over ghosts
                do n=1,Model%Ng
                    ig = Grid%is-n
                    ip = Grid%is+n-1
                    !Map n=[1,ng] -> inExyz
                    !ASSUMING PsiSt=-3 if you're nudging, so n=[1,ng]->[4,...,1]
                    np = Model%nG-n+1 !Mapping into 1,4 of inExyz

                !-------
                !Get geometry for this ghost and matching physical

                    !call cellCenter(Grid,ig,jp,kp,xc,yc,zc)
                    !NOTE: Using j/k instead of jp/kp to deal with double-corner sign flip
                    call cellCenter(Grid,ig,j ,k ,xc,yc,zc)
                    rHatG = normVec([xc,yc,zc])

                    call cellCenter(Grid,ip,jp,kp,xc,yc,zc)
                    rHatP = normVec([xc,yc,zc])

                    invlat = rad2deg*InvLatitude([xc,yc,zc]) !Convert to degrees

                    if (invlat<=llBC) then
                        isLL = .true.
                    else
                        isLL = .false.
                    endif

                !-------
                !Get E, dipole/perturbation and ExB velocity/Mirror velocity

                    !Get velocity from i-reflected active cell
                    Vmir = State%Gas(ip,jp,kp,MOMX:MOMZ,BLK)/max(State%Gas(ip,jp,kp,DEN,BLK),dFloor)
                    Exyz = bc%inExyz(np,jp,kp,:)
                    call Dipole(xc,yc,zc,Bd(XDIR),Bd(YDIR),Bd(ZDIR))
                    dB = State%Bxyz(ip,jp,kp,:)
                    !Using ExB everywhere
                    Veb = cross(Exyz,Bd)/dot_product(Bd,Bd)
                    !Use ExB (w/o radial) and mirror
                    Vxyz = Veb - rHatP*dot_product(rHatP,Veb) - rHatP*dot_product(rHatP,Vmir)

                    
                !-------
                !Set ghost hydro quantities
                    !Let density float
                    call SphereWall(Model,State%Gas(ig,j,k,:,:),State%Gas(ip,jp,kp,:,:),Vxyz)
                    !Now do polar outflow if testing
                    if (Model%doMultiF .and. (invlat>=70) .and. (Model%nSpc>2)) then
                        gW(DEN) = 100.0
                        gW(VELX:VELZ) = 0.2*rHatP + Veb - rHatP*dot_product(rHatP,Veb)
                        gW(PRESSURE) = 1.0e-3
                        call CellP2C(Model,gW,gCon)
                        State%Gas(ig,j,k,:,3) = gCon
                        !Reset bulk
                        call MultiF2Bulk(Model,State%Gas(ig,j,k,:,:))
                    endif

                !-------
                !Now handle magnetic quantities
                    if (isLL) then
                        !In low-lat enforce full dipole
                        State%Bxyz   (ig,j,k,:) = 0.0
                        State%magFlux(ig,j,k,:) = 0.0

                    else
                        !Mirror fluxes to minimize gradient
                        State%Bxyz(ig,j,k,:) = dB
                        State%magFlux(ig,j,k,IDIR) = State%magFlux(ip,jp,kp,IDIR)
                        State%magFlux(ig,j,k,JDIR) = State%magFlux(ip,jp,kp,JDIR)
                        State%magFlux(ig,j,k,KDIR) = State%magFlux(ip,jp,kp,KDIR)
                    endif

                enddo !n
            enddo
        enddo

    end subroutine IonInner
    
end module uservoltic
