!Various utilities for magnetospheric runs

module msphutils
    use kdefs
    use gamtypes
    use volttypes
    use gamutils
    use math, magRampDown => PenticRampDown
    use gridutils
    use output
    use background
    use multifluid
    use earthhelper

    implicit none

!Set planet scalings here
!gx0 [m]
!gv0 [m/s]
!gD0 = 1.67e-21 ! 1 AMU/cc [kg/m3]
!gP0 [nPa]
!gB0 [nT]
!gT0 [s]
!gG0=0 ! Grav acceleration [m/s2]

!Example: Earth scaling
!x0 = 6.38e6 m [1 RE]
!rho0 = 1.67e-21 kg/m^3 [1 particle/cc]
!v0 = 10e5 m/s [100 km/s]
!B0 = 4.58 nT
!P0 = 1.67e-2 [nPa]
!t0 = 63.8 second
!M0 = -0.31*1.0e+5/B0

!Private module variables
    real(rp), private :: tScl !Needed for time output, TODO: Fix this
    real(rp), private :: Rion ! Planetary ionosphere radius
    !Chill out parameters
    logical , private :: doLFMChill = .true. !Do LFM-style chilling
    logical , private :: doGAMChill = .true. !Do GAM-style chilling
    real(rp), private :: RhoCO = 1.0e-3 ! Number density
    real(rp), private :: CsCO  = 1.0e-2  ! Cs chillout, m/s
    real(rp), private :: cLim = 1.5 ! Cool when sound speed is above cLim*Ca
    
    
    !Dipole cut values
    !real(rp), private :: rCut=4.5, lCut=3.5 !LFM values
    real(rp), private :: rCut=16.0,lCut=8.0
    real(rp), private :: dcScl=4.0 !Tailward dipole cut = dcScl x rCut, dawn/dusk = dcScl x rCut/2
    real(rp), private :: xSun,xTail,yMax
    real(rp), private :: x0,Aq,Bq,sInner

    !Grav/rotation values
    real(rp), private :: GM0  = 0.0 !Gravitational force coefficient
    real(rp), private :: Psi0 = 0.0 ! corotation potential coef
    real(rp), private :: M0   = 0.0 !Magnetic moment

    !Ingestion switch
    logical, private :: doIngest = .true.

    contains

    !Set magnetosphere parameters
    !Planet ID optional
    subroutine setMagsphere(Model,xmlInp,pStrO)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: xmlInp
        character(len=*), intent(in), optional :: pStrO

        character(len=strLen) :: pID !Planet ID string
        real(rp) :: M0g,rScl
        real(rp) :: gx0,gv0,gD0,gP0,gB0,gT0,gG0
        logical :: doCorot
        
        !Set some defaults
        gD0 = 1.67e-21 ! 1 AMU/cc [kg/m3]
        gG0 = 0.0 ! Grav acceleration [m/s2]

        if (present(pStrO)) then
            pID = pStrO
        else
            call xmlInp%Set_Val(pID,"prob/planet","Earth")
        endif

        Rion = 0.0 !Default
        select case (trim(toUpper(pID)))

        case("Earth","earth","EARTH")
            gx0 = REarth  ! [m]
            gv0 = 100.e3  ! [m/s]
            gG0 = 9.807   ! [m/s2]
            call xmlInp%Set_Val(M0g,"prob/M0",EarthM0g) !Mag moment [gauss]
            !Using corotation potential for Earth
            Psi0 = EarthPsi0 !kV
            Rion = RionE*1.e6/gx0 ! Radius of ionosphere in code units (RionE defined in kdefs in 1000km)
            Model%doGrav = .true.
        case("Saturn","saturn","SATURN")
            gx0 = RSaturnXE*REarth  ! [m]
            gv0 = 100.0e+3    ! [m/s]
            gG0 = 10.43       ! [m/s2]
            call xmlInp%Set_Val(M0g,"prob/M0",SaturnM0g) !Mag moment [gauss]
            Psi0 = 137.83*92 !kV
            Rion = 1.01 !Assuming
            Model%doGrav = .true.
        case("JUPITER")
            gx0 = RJupiterXE*REarth  ! [m]
            gv0 = 100.0e+3    ! [m/s]
            gG0 = 24.79       ! [m/s2]
            call xmlInp%Set_Val(M0g,"prob/M0",JupiterM0g) !Mag moment [gauss]
            Psi0 = -2.5*1702.9*92.0 !kV
            Rion = 1.01 !Assuming
            Model%doGrav = .true.
        case("Mercury","mercury","MERCURY")
            gx0 = RMercuryXE*REarth  ! [m]
            gv0 = 100.0e+3    ! [m/s]
            gG0 = 0.0       ! [m/s2]
            call xmlInp%Set_Val(M0g,"prob/M0",MercuryM0g) !Mag moment [gauss]
            Psi0 = 0.0 !kV
            Rion = 1.05 !Assuming
            Model%doGrav = .false.
        case("Neptune","NEPTUNE")
            gx0 = RNeptuneXE*REarth
            gv0 = 100.0e+3    ! [m/s]
            gG0 = 11.15       ! [m/s2]
            call xmlInp%Set_Val(M0g,"prob/M0",NeptuneM0g) !Mag moment [gauss]
            Psi0 = 10.024*92.0 !kV
            Rion = 1.01 !Assuming
            Model%doGrav = .false.
        case("Other","other","OTHER") ! Defaults to Earth values
            call xmlInp%Set_Val(gx0,"prob/x0",REarth)            ! [m]
            call xmlInp%Set_Val(gv0,"prob/v0",100.e3)            ! [m/s]
            call xmlInp%Set_Val(gG0,"prob/G0",9.807)             ! [m/s2] 
            call xmlInp%Set_Val(M0g,"prob/M0",EarthM0g)          ! [gauss]
            call xmlInp%Set_Val(Psi0,"prob/Psi0",EarthPsi0)      ! [kV]
            call xmlInp%Set_Val(Rion,"prob/Rion",RionE*1.e6/gx0) 
            call xmlInp%Set_Val(Model%doGrav,"prob/doGrav",.true.)
        end select

        !Set some chilling parameters
        !If we're using source term set default as chilling off
        if (Model%doSource) then
            doGAMChill = .false.
            doLFMChill = .false.
        endif
        !Now read XML parameters w/ default options
        call xmlInp%Set_Val(doGAMChill,"chill/doGAMChill",doGAMChill)
        call xmlInp%Set_Val(doLFMChill,"chill/doLFMChill",doLFMChill)
        if (doLFMChill) then
            call xmlInp%Set_Val(RhoCO,"chill/RhoCO",RhoCO)
        endif
        if (doGAMChill) then
            call xmlInp%Set_Val(cLim,"chill/cLim",cLim)
        endif

        !Whether to ignore ingestion (if set)
        if (Model%doSource) then
            call xmlInp%Set_Val(doIngest,"source/doIngest",.true.)
        endif
        
        call xmlInp%Set_Val(doCorot,"prob/doCorot",.true.)
        if (.not. doCorot) then
            !Zero out corotation potential
            Psi0 = 0.0
        endif

        gT0 = gx0/gv0 !Set time scaling
        gB0 = sqrt(Mu0*gD0)*gv0*1.0e+9 !T->nT
        gP0 = gD0*gv0*gv0*1.0e+9 !P->nPa
        M0  = -M0g*1.0e+5/gB0 !Magnetic moment
        GM0 = gG0*gx0/(gv0*gv0)

        Model%isMagsphere = .true.
        Model%MagM0 = M0
        
        !Add gravity if required
        if (Model%doGrav) then
            !Force spherical gravity (zap non-radial components)
            Model%doSphGrav = .true.
            Model%Phi => PhiGrav
        endif

        !Change console output pointer
        tScl = gT0
        timeString => magsphereTime
        
        if (Model%isLoud) then
            write(*,*) '---------------'
            write(*,*) 'Magnetosphere normalization'
            write(*,*) 'T0   [s]    = ', gT0
            write(*,*) 'x0   [m]    = ', gx0
            write(*,*) 'v0   [m/s]  = ', gv0
            write(*,*) 'P0   [nPa]  = ', gP0
            write(*,*) 'B0   [nT]   = ', gB0
            write(*,*) 'g    [m/s2] = ', gG0
            write(*,*) 'psi0 [kV]   = ', Psi0
            write(*,*) '---------------'
        endif
        
        !Save scaling to gUnits_T structure in Model
        Model%Units%gT0 = gT0
        Model%Units%gx0 = gx0
        Model%Units%gv0 = gv0
        Model%Units%gD0 = gD0
        Model%Units%gP0 = gP0
        Model%Units%gB0 = gB0
        Model%Units%gG0 = gG0

        !Add normalization/labels to output slicing
        Model%gamOut%tScl = gT0
        Model%gamOut%dScl = 1.0
        Model%gamOut%vScl = gv0*1.0e-3 !km/s
        Model%gamOut%pScl = gP0
        Model%gamOut%bScl = gB0

        Model%gamOut%uID = trim(toUpper(pID))
        Model%gamOut%tID = 's'
        Model%gamOut%dID = '#/cc'
        Model%gamOut%vID = 'km/s'
        Model%gamOut%pID = 'nPa'
        Model%gamOut%bID = 'nT'

        !Reinterpret pressure floor as nPa
        pFloor = pFloor/gP0
    end subroutine

    subroutine magsphereTime(T,tStr)
        real(rp), intent(in) :: T
        character(len=strLen), intent(out) :: tStr

        if (abs(T*tScl)>60.0) then
            write(tStr,'(f9.3,a)' ) T*tScl/60.0, ' [min]'
        else
            write(tStr,'(es9.2,a)') T*tScl     , ' [sec]'
        endif

    end subroutine magsphereTime
    
    !Just return module private magnetic moment
    function MagMoment() result(M)
        real(rp) :: M
        M = M0
    end function MagMoment

    !Just return module private ionospheric radius
    function RadIonosphere() result(R)
        real(rp) :: R
        R = Rion
    end function RadIonosphere

    !Fix inner shell electric fields (ijk) to remix values
    subroutine IonEFix(Model,Gr,State,inEijk)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State
        real(rp), intent(in) :: inEijk(PsiSh+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,1:NDIM)

        integer :: is0,j,k

        if (.not. Gr%hasLowerBC(IDIR)) return


        !Multiply by edge length to turn field into EMF
        !Check if Eijk has first shell
        if ( (PsiSh+PsiSt) >= Gr%is ) then
            !Find index in inEijk of first active shell
            is0 = 1-PsiSt+1

            !$OMP PARALLEL DO default(shared)
            do k=Gr%ksg,Gr%keg
                do j=Gr%jsg,Gr%jeg
                    State%Efld(Gr%is,j,k,JDIR) = Gr%edge(Gr%is,j,k,JDIR)*inEijk(is0,j,k,JDIR)
                    State%Efld(Gr%is,j,k,KDIR) = Gr%edge(Gr%is,j,k,KDIR)*inEijk(is0,j,k,KDIR)
                enddo
            enddo !k loop
        endif

    end subroutine IonEFix

    !Apply spherical wall boundary over all species
    !Take ghost conserved quantities and conjugate point conserved quantities
    subroutine SphereWall(Model,gU,pU,V,Dopt)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: gU(NVAR,BLK:Model%nSpc)
        real(rp), intent(in)    :: pU(NVAR,BLK:Model%nSpc)
        real(rp), intent(in) :: V(NDIM)
        real(rp), intent(in), optional :: Dopt

        integer :: s,s0,sE
        real(rp), dimension(NVAR) :: pW,pCon,gW,gCon
        real(rp), dimension(0:Model%nSpc) :: RhoMin
        real(rp) :: D
        logical :: doWall

        RhoMin(BLK) = 0.0
        if (Model%doMultiF) then
            !Don't do bulk
            s0 = 1
            sE = Model%nSpc
            RhoMin(1:Model%nSpc) = Spcs(:)%dVac
        else
            !Only do bulk
            s0 = BLK
            sE = BLK
        endif

        do s=s0,sE
            !Test if active cell for this species is good
            D = pU(DEN,s)
            doWall = ( D >= RhoMin(s) )
            !Get conserved quantities of active, map to prim
            pCon = pU(:,s)
            call CellC2P(Model,pCon,pW)

            !Now, if Dopt is set use that for density if BLK or 1st species
            if ( present(Dopt) .and. (s0<=1) ) then
                pW(DEN) = Dopt
            endif

            if (doWall) then
                !Set prim values of ghost
                gW(DEN)       = pW(DEN)
                gW(PRESSURE)  = pW(PRESSURE)
                gW(VELX:VELZ) = V
            else 
                !For now, do the same thing even if conjugate cell is under threshold
                !Probably necessary to get velocity correct
                gW(DEN)       = pW(DEN)
                gW(PRESSURE)  = pW(PRESSURE)
                gW(VELX:VELZ) = V
            endif
            
            !Now go back to con and set value
            call CellP2C(Model,gW,gCon)
            gU(:,s) = gCon
        enddo

        if (Model%doMultiF) then
            call MultiF2Bulk(Model,gU)
        endif

    end subroutine SphereWall

    !Cooling function to deal with anomalous heating
    subroutine ChillOut(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k
        integer :: s,s0,sE
        real(rp), dimension(NVAR) :: pW, pCon
        real(rp), dimension(0:Model%nSpc) :: RhoMin
        real(rp) :: D,P,CsC,Pc,Leq,tau
        logical :: doChill

        !Check if there's no chill and get out
        if ( (.not. doLFMChill) .and. (.not. doGAMChill) ) then
            return
        endif

        RhoMin(BLK) = 0.0
        if (Model%doMultiF) then
            !Don't do bulk
            s0 = 1
            sE = Model%nSpc
            RhoMin(1:Model%nSpc) = Spcs(:)%dVac
        else
            !Only do bulk
            s0 = BLK
            sE = BLK
        endif

        !$OMP PARALLEL DO default(shared) collapse(3) &
        !$OMP private(s,i,j,k,pW,pCon,D,P,CsC,Pc,Leq,tau,doChill)
        do s=s0,sE
            do k=Grid%ksg,Grid%keg
                do j=Grid%jsg,Grid%jeg
                    do i=Grid%isg,Grid%ieg

                        !Check species
                        pCon = State%Gas(i,j,k,:,s)
                        D = pCon(DEN)

                    !Check for low densities
                        if (s>BLK) then
                            doChill = ( D <= RhoCO ) .and. ( D >= RhoMin(s) )
                        else
                            doChill = ( D <= RhoCO )
                        endif

                        !If density is low keep things chill by setting sound speed
                        if (doChill .and. doLFMChill) then
                            pCon = State%Gas(i,j,k,:,s)
                            call CellC2P(Model,pCon,pW)
                            !Set pressure to ensure Cs = CsCO
                            CsC = CsCO/Model%Units%gv0 !Cs in code units
                            P = pW(DEN)*CsC*CsC/Model%gamma
                            pW(PRESSURE) = max(P,pFloor)
                            call CellP2C(Model,pW,pCon)
                            State%Gas(i,j,k,:,s) = pCon
                        endif

                    !Check for too fast sound speed

                        !Get sound speed
                        if (s>BLK) then
                            if ( D >= RhoMin(s) ) then
                                call CellPress2Cs(Model,pCon,CsC)
                            else
                                CsC = 0.0 !Ignore vacuum
                            endif
                        else
                            !Bulk fluid
                            call CellPress2Cs(Model,pCon,CsC)
                        endif

                        !If sound speed is faster than "light", chill the fuck out
                        doChill = Model%doBoris .and. (CsC>cLim*Model%Ca) .and. (doGAMChill)
                        if (doChill .and. Model%doSource) then
                            !If this is a pressure ingestion region, then let the pressure go wild
                            if (Grid%Gas0(i,j,k,IMPR ,BLK)>TINY) doChill = .false.
                        endif !doChill and doSource check
                        if (doChill) then
                            call CellC2P(Model,pCon,pW)
                            P = pW(PRESSURE) !Cell pressure

                            !Find target pressure w/ sound speed = cLim*Ca
                            Pc = pW(DEN)*(cLim*Model%Ca)**2.0/Model%gamma
                            !Calculate cooling rate, L/CsC ~ lazy bounce timescale
                            Leq = DipoleL(Grid%xyzcc(i,j,k,:))
                            !Convert to m/(m/s)/(code time)
                            tau = (Leq*Model%Units%gx0)/(CsC*Model%Units%gv0*Model%Units%gT0)

                            !Cool pressure to target w/ timescale tau
                            pW(PRESSURE) = P - (Model%dt/tau)*(P-Pc)

                            !Go back to conserved vars and save
                            call CellP2C(Model,pW,pCon)
                            State%Gas(i,j,k,:,s) = pCon
                        endif

                    enddo !i loop
                enddo
            enddo
        enddo !Species loop

                    
        if (Model%doMultiF) then
            call State2Bulk(Model,Grid,State)
        endif
                    
    end subroutine ChillOut

    !Create a cut dipole
    subroutine genCutDipole(Model,Grid,State,xmlInp)
        type(Model_T), intent(inout) :: Model
        type(Grid_T), intent(inout) :: Grid
        type(State_T), intent(inout) :: State
        type(XML_Input_T), intent(in) :: xmlInp

        real(rp), dimension(:,:,:,:), allocatable :: bFlux0, bFluxT
        procedure(VectorField_T), pointer :: Axyz

        integer :: i,j,k
        real(rp), dimension(8,NDIM) :: xyzC
        real(rp) :: rMax,rMin,MagP

        !Get values for initial field cutoffs

        ! !LFM values
        ! call xmlInp%Set_Val(xSun  ,"prob/xMax",20.0_rp  )
        ! call xmlInp%Set_Val(yMax  ,"prob/yMax",75.0_rp  )
        ! call xmlInp%Set_Val(xTail ,"prob/xMin",-185.0_rp)
        
        !Gamera values
        call xmlInp%Set_Val(xSun  ,"prob/xMax",20.0_rp  )
        call xmlInp%Set_Val(yMax  ,"prob/yMax",75.0_rp  )
        call xmlInp%Set_Val(xTail ,"prob/xMin",-225.0_rp)

        call xmlInp%Set_Val(sInner,"prob/sIn" ,0.96_rp  )

        !Get cut dipole values
        call xmlInp%Set_Val(rCut ,"background/rCut" ,rCut)
        call xmlInp%Set_Val(lCut ,"background/lCut" ,lCut)
        call xmlInp%Set_Val(dcScl,"background/dcScl",dcScl)

        !Calculate some derived quantities/alloc arrays
        Aq = 0.5*(xSun-xTail)
        x0 = Aq - xSun
        Bq = yMax

        call allocGridVec(Model,Grid,bFlux0,.true.,NDIM)
        call allocGridVec(Model,Grid,bFluxT,.true.,NDIM)
        bFlux0 = 0.0
        bFluxT = 0.0

        !Calculate flux from B0_xyz components (cut-off dipole)
        Model%doBackground = .true.

        Model%B0 => cutDipole
        Axyz     => cutDipole

        call AddB0(Model,Grid,Model%B0)

        !Be careful and forcibly zero out cut dipole forces near Earth
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyzC,rMin,rMax)
        do k=Grid%ks, Grid%ke
            do j=Grid%js, Grid%je
                do i=Grid%is, Grid%ie
                    !Get 8 cell corners
                    call cellCoords(Model,Grid,i,j,k,xyzC)
                    rMin = minval(norm2(xyzC,dim=2))
                    rMax = maxval(norm2(xyzC,dim=2))
                    !Force hard zero outside of dipole cut region (can be non-zero due to quadrature error)
                    if (rMin + TINY < rCut       ) Grid%dpB0(i,j,k,:) = 0.0
                    !if (rMax - TINY > rCut + lCut) Grid%dpB0(i,j,k,:) = 0.0

                enddo
            enddo
        enddo

        call VectorField2Flux(Model,Grid,State,Axyz)
        bFlux0(:,:,:,:) = State%magFlux(:,:,:,:) !bFlux0 = B0

        !Now use vector potential (elliptical cut-off) to set total initial field
        Axyz => VP_Init
        call VectorPot2Flux(Model,Grid,State,Axyz)
        bFluxT(:,:,:,:) = State%magFlux(:,:,:,:)

        !Subtract off background face-fluxes to get dB
        State%magFlux = bFluxT - bFlux0

    end subroutine genCutDipole

    
    !-----------------------------
    !Vector potential/vector fields for dipole/cut-dipoles
    subroutine VP_Init(x,y,z,Ax,Ay,Az)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Ax,Ay,Az

        real(rp) :: r,M,s
        r = sqrt(x**2.0+y**2.0+z**2.0)
        call VP_Dipole(x,y,z,Ax,Ay,Az)

        s = sqrt( ((x+x0)/Aq)**2.0 + (y**2.0+z**2.0)/Bq**2.0 )
        M = magRampDown(s,sInner,1.0-sInner)
        !Dipole moment already included in VP_Dipole
        Ax = M*Ax
        Ay = M*Ay
        Az = M*Az                
    end subroutine VP_Init

    subroutine VP_Dipole(x,y,z,Ax,Ay,Az)
        
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Ax,Ay,Az

        real(rp), dimension(NDIM) :: A,m,r,rhat
        m = [0,0,1]
        r = [x,y,z]
        rhat = r/norm2(r)

        A = M0*cross(m,rhat)/(dot_product(r,r))
        Ax = A(XDIR)
        Ay = A(YDIR)
        Az = A(ZDIR)
    end subroutine VP_Dipole

    subroutine Dipole(x,y,z,Ax,Ay,Az)
        
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Ax,Ay,Az

        real(rp), dimension(NDIM) :: m,r,A
        real(rp) :: rad

        r = [x,y,z]
        rad = norm2(r)
        m = [0.0_rp,0.0_rp,M0]

        A = 3*dot_product(m,r)*r/rad**5.0 - m/rad**3.0

        Ax = A(XDIR)
        Ay = A(YDIR)
        Az = A(ZDIR)

    end subroutine Dipole

    !Silly vector wrapper to make a dipole function
    function VecDipole(xyz) result(Bd)
        real(rp), intent(in) :: xyz(NDIM)
        real(rp), dimension(NDIM) :: Bd
        call Dipole(xyz(XDIR),xyz(YDIR),xyz(ZDIR),Bd(XDIR),Bd(YDIR),Bd(ZDIR))

    end function VecDipole

    subroutine cutDipole(x,y,z,Ax,Ay,Az)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Ax,Ay,Az
   
        real(rp) :: r,M,phid,rScl,yp
        r = sqrt(x**2.0 + y**2.0 + z**2.0)

        call Dipole(x,y,z,Ax,Ay,Az)
        !Get mollifier term
        yp = sqrt(y**2.0 + z**2.0)
        phid = atan2(yp,x)*180.0/PI
        rScl = 1.0 + RampUp(phid,0.0_rp,180.0_rp)*(dcScl-1.0) !Between 1,dcScl
        M = magRampDown(r,rCut*rScl,lCut*rScl)
        !M = magRampDown(r,rCut,lCut)

        Ax = M*Ax
        Ay = M*Ay
        Az = M*Az

    end subroutine cutDipole

    !Ingest density/pressure information from Grid%Gas0
    !Treat Gas0 as target value
    subroutine MagsphereIngest(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State

        integer :: i,j,k
        real(rp), dimension(NVAR) :: pW, pCon

        real(rp) :: Tau,dRho,dP,Pmhd,Prcm
        real(rp), dimension(NDIM) :: Mxyz,Vxyz,B
        real(rp), dimension(NDIM,NDIM) :: Lam,Laminv

        logical  :: doIngestIJK,doInD,doInP

        if (Model%doMultiF) then
            write(*,*) 'Source ingestion not implemented for multifluid, you should do that'
            stop
        endif

        if (Model%t<=0) return

        if (.not. doIngest) return

       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP private(i,j,k,doInD,doInP,doIngestIJK,pCon,pW) &
       !$OMP private(Tau,dRho,dP,Pmhd,Prcm,Mxyz,Vxyz,B) &
       !$OMP private(Lam,Laminv)

        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie
                    doInD = (Gr%Gas0(i,j,k,IMDEN,BLK)>TINY)
                    doInP = (Gr%Gas0(i,j,k,IMPR ,BLK)>TINY)
                    doIngestIJK = doInD .or. doInP

                    if (.not. doIngestIJK) cycle

                    pCon = State%Gas(i,j,k,:,BLK)
                    call CellC2P(Model,pCon,pW)
                    Pmhd = pW(PRESSURE)

                    !Calculate semi-rel momentum (identical to momentum if no boris correction)
                    B = State%Bxyz(i,j,k,:)
                    call Mom2Rel(Model,pW(DEN),B,Lam)
                    Mxyz = matmul(Lam,pCon(MOMX:MOMZ)) !semi-relativistic momentum

                    !Get timescale, taking directly from Gas0
                    Tau = Gr%Gas0(i,j,k,IMTSCL,BLK)
                                        
                    if (Tau<Model%dt) Tau = Model%dt !Unlikely to happen

                    if (doInD) then
                        dRho = Gr%Gas0(i,j,k,IMDEN,BLK) - pW(DEN)
                        pW(DEN) = pW(DEN) + (Model%dt/Tau)*dRho
                    endif

                    if (doInP) then
                        Prcm = Gr%Gas0(i,j,k,IMPR,BLK)
                        !Assume already wolf-limited or not
                        dP = Prcm - Pmhd
                        pW(PRESSURE) = pW(PRESSURE) + (Model%dt/Tau)*dP
                    endif

                    !Get new velocity, start w/ updated inverse matrix
                    call Rel2Mom(Model,pW(DEN),B,Laminv)
                    Vxyz = matmul(Laminv,Mxyz)/max(pW(DEN),dFloor)
                    pW(VELX:VELZ) = Vxyz

                    !Now put back
                    call CellP2C(Model,pW,pCon)
                    State%Gas(i,j,k,:,BLK) = pCon
                enddo
            enddo
        enddo
                    
    end subroutine MagsphereIngest

    !Set gPsi (corotation potential)
    subroutine CorotationPot(Model,Grid,gPsi)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp)  , intent(inout) :: gPsi(1:PsiSh+1,Grid%js:Grid%je+1,Grid%ks:Grid%ke+1)

        integer :: i,iG,j,k

        real(rp), dimension(NDIM) :: xyz
        real(rp) :: L

        if (abs(Psi0)<=TINY) then
            return
        endif

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,iG,j,k,xyz,L)        
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=1,PsiSh+1
                    !i = Grid%is+PsiSt-1
                    iG = i+PsiSt-1 !Global i index for Grid access
                    xyz = Grid%xyz(iG,j,k,:)
                    L = DipoleL(xyz)
                    gPsi(i,j,k) = gPsi(i,j,k) - Psi0/L
                enddo
            enddo
        enddo !K loop

    end subroutine CorotationPot

    subroutine PhiGrav(x,y,z,pot)
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: pot

        real(rp) :: rad
        rad = sqrt(x**2.0 + y**2.0 + z**2.0)
        pot = -GM0/rad
    end subroutine PhiGrav
      
end module msphutils
