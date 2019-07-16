!Various utilities for magnetospheric runs

module msphutils
    use kdefs
    use types
    use gamutils
    use math, magRampDown => CubicRampDown
    use output
    use gioH5

    use gridutils
    use xml_input
    use strings
    use background
    use multifluid
    use ingestpress

    implicit none
!Earth scaling
!x0 = 6.38e6 m [1 RE]
!rho0 = 1.67e-21 kg/m^3 [1 particle/cc]
!v0 = 10e5 m/s [100 km/s]
!B0 = 4.58 nT
!P0 = 1.67e-2 [nPa]
!t0 = 63.8 second
!M0 = -0.31*1.0e+5/B0

    real(rp) :: Rion ! planetary ionosphere radius
    real(rp) :: gx0 ! [m]
    real(rp) :: gv0 ! [m/s]
    real(rp) :: gD0 = 1.67e-21 ! 1 AMU/cc [kg/m3]
    real(rp) :: gP0 ! [nPa]
    real(rp) :: gB0 ! [nT]
    real(rp) :: gT0 ! [s]
    real(rp) :: gG0=0 ! Grav acceleration [m/s2]
    real(rp) :: M0 !Default magnetic moment [gauss]
    real(rp) :: GM0 !Gravitational force coefficient
    real(rp) :: Psi0=0. ! corotation potential coef
    
    !real(rp) :: rCut=4.5, lCut=3.5 !LFM values
    real(rp) :: rCut=4.0, lCut=5.0
    real(rp) :: xSun,xTail,yMax
    real(rp) :: x0,Aq,Bq,sInner

    integer :: JpSh = 1 !Number of cc current shells for CMI
    integer :: JpSt = 2 !First shell (i) to calculate current
    logical :: doRingFAC = .true. !Do some ring-processing on currents before sending to remix

    !integer :: PsiSh = 1 !Number of *SHELLS* getting nodes at, ie PsiSh+1 = # i nodes
    !integer :: PsiSt = 1 !Starting shell of potential

    !Get 5 shells, all ghosts plus first active
    integer :: PsiSh = 5 !Number of *SHELLS* getting nodes at, ie PsiSh+1 = # i nodes
    integer :: PsiSt = -3 !Starting shell of potential

    !Chill out parameters
    real(rp) :: RhoCO = 1.0e-3 ! Number density
    real(rp) :: CsCO  = 1.0e-2  ! Cs chillout, m/s

    contains

    !Set magnetosphere parameters
    !Planet ID optional
    subroutine setMagsphere(Model,xmlInp,pStrO)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: xmlInp
        character(len=strLen), intent(in), optional :: pStrO

        character(len=strLen) :: pID !Planet ID string
        real(rp) :: M0g,rScl

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
            Psi0 = 0*92.0 !kV
            Rion = RionE*1.e6/gx0 ! Radius of ionosphere in code units (RionE defined in kdefs in 1000km)
            Model%doGrav = .false.
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

        end select

        gT0 = gx0/gv0 !Set time scaling
        gB0 = sqrt(Mu0*gD0)*gv0*1.0e+9 !T->nT
        gP0 = gD0*gv0*gv0*1.0e+9 !P->nPa
        M0 = -M0g*1.0e+5/gB0 !Magnetic moment
        GM0 = gG0*gx0/(gv0*gv0)

        !Add gravity if required
        if (Model%doGrav) then
            Model%Phi => PhiGrav
        endif

        !Change console output pointer
        timeString => magsphereTime

        write(*,*) '---------------'
        write(*,*) 'Magnetosphere normalization'
        write(*,*) 'T0 [s]    = ', gT0
        write(*,*) 'x0 [m]    = ', gx0
        write(*,*) 'v0 [m/s]  = ', gv0
        write(*,*) 'P0 [nPa]  = ', gP0
        write(*,*) 'B0 [nT]   = ', gB0
        write(*,*) 'g  [m/s2] = ', gG0
        write(*,*) '---------------'

        !Add normalization/labels to output slicing
        gamOut%tScl = gT0
        gamOut%dScl = 1.0
        gamOut%vScl = gv0*1.0e-3 !km/s
        gamOut%pScl = gP0
        gamOut%bScl = gB0

        gamOut%uID = trim(toUpper(pID))
        gamOut%tID = 's'
        gamOut%dID = '#/cc'
        gamOut%vID = 'km/s'
        gamOut%pID = 'nPa'
        gamOut%bID = 'nT'

    end subroutine

    subroutine magsphereTime(T,tStr)
        real(rp), intent(in) :: T
        character(len=strLen), intent(out) :: tStr

        if (T*gT0>60.0) then
            write(tStr,'(f9.3,a)' ) T*gT0/60.0, ' [min]'
        else
            write(tStr,'(es9.2,a)') T*gT0     , ' [sec]'
        endif

    end subroutine magsphereTime
    
    !Convert electic potential from ionosphere to E fields for inner BCs
    subroutine Ion2MHD(Model,Grid,gPsi,inEijk,inExyz,pSclO)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), intent(inout) :: gPsi  (1:PsiSh+1,Grid%js:Grid%je+1,Grid%ks:Grid%ke+1)
        real(rp), intent(inout) :: inEijk(1:PsiSh+1,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), intent(inout) :: inExyz(1:PsiSh,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), intent(in), optional :: pSclO
        integer :: i,j,k,iG
        integer :: NumP
        real(rp) :: ijkDet, Mijk(NDIM,NDIM)
        real(rp), dimension(NDIM) :: iCC,jCC,kCC,ccEijk
        real(rp) :: ionScl

        !Set Eijk fields
        inEijk = 0.0
        inExyz = 0.0

        !Set scaling if present
        ionScl = 1.0
        if (present(pSclO)) then
            ionScl = pSclO
        endif

        NumP = Grid%ke-Grid%ks+1

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,iG) 
        do i=1,PsiSh+1
            do k=Grid%ks,Grid%ke+1
                do j=Grid%js,Grid%je+1

                    iG = i+PsiSt-1 !Global i index for Grid access

                    if (i <= PsiSh) then
                        inEijk(i,j,k,IDIR) = -( gPsi(i+1,j,k) - gPsi(i,j,k) )/Grid%edge(iG,j,k,IDIR)
                    endif
                    if (j <= Grid%je) then
                        inEijk(i,j,k,JDIR) = -( gPsi(i,j+1,k) - gPsi(i,j,k) )/Grid%edge(iG,j,k,JDIR)
                    endif
                    if (k <= Grid%ke) then
                        inEijk(i,j,k,KDIR) = -( gPsi(i,j,k+1) - gPsi(i,j,k) )/Grid%edge(iG,j,k,KDIR)
                    endif
                    !Scale field
                    inEijk(i,j,k,:) = inEijk(i,j,k,:)/ionScl

                enddo
            enddo
        enddo

        !$OMP PARALLEL DO default(shared)
        do i=1,PsiSh+1
            !Correct for pole
            if (Model%Ring%doS) then
                inEijk(i,Grid%js,:,KDIR) = 0.0
                inEijk(i,Grid%js,:,IDIR) = sum(inEijk(i,Grid%js,Grid%ks:Grid%ke,IDIR))/NumP
            endif
            if (Model%Ring%doE) then
                inEijk(i,Grid%je+1,:,JDIR) = 0.0
                inEijk(i,Grid%je+1,:,KDIR) = 0.0
                inEijk(i,Grid%je+1,:,IDIR) = sum(inEijk(i,Grid%je+1,Grid%ks:Grid%ke,IDIR))/NumP
            endif
            !Enforce periodicity constraints (needed for differencing in next step)
            inEijk(i,:,Grid%ke+1,IDIR:KDIR) = inEijk(i,:,Grid%ks,IDIR:KDIR)

        enddo

        !Now turn nc-IJK FIELDS into cc-XYZ fields
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,k,iG,ccEijk,iCC,jCC,kCC,Mijk,ijkDet)
        do i=1,PsiSh
            do k=Grid%ks,Grid%ke
                do j=Grid%js,Grid%je
                    iG = i+PsiSt-1 !Global i index into grid arrays
                    !Get cell-centered XYZ fields
                  
                    ccEijk(IDIR) = 0.25*(inEijk(i,j,k,IDIR)+inEijk(i,j,k+1,IDIR)+inEijk(i,j+1,k,IDIR)+inEijk(i,j+1,k+1,IDIR))
                    ccEijk(JDIR) = 0.25*(inEijk(i,j,k,JDIR)+inEijk(i,j,k+1,JDIR)+inEijk(i+1,j,k,JDIR)+inEijk(i+1,j,k+1,JDIR))
                    ccEijk(KDIR) = 0.25*(inEijk(i,j,k,KDIR)+inEijk(i+1,j,k,KDIR)+inEijk(i,j+1,k,KDIR)+inEijk(i+1,j+1,k,KDIR))

                    !Get cell-centered ijk vectors @ cell-center
                    iCC = ijkVec(Model,Grid,iG,j,k,IDIR)
                    jCC = ijkVec(Model,Grid,iG,j,k,JDIR)
                    kCC = ijkVec(Model,Grid,iG,j,k,KDIR)

                    !Get inverse matrix
                    call ijkMatrix(Model,Grid,iCC,jCC,kCC,Mijk)

                    !Determinant of ijk matrix
                    ijkDet = iCC(XDIR)*Mijk(1,1) + jCC(XDIR)*Mijk(1,2) + kCC(XDIR)*Mijk(1,3)

                    inExyz(i,j,k,XDIR) = (Mijk(1,1)*ccEijk(IDIR) + Mijk(1,2)*ccEijk(JDIR) + Mijk(1,3)*ccEijk(KDIR))/ijkDet
                    inExyz(i,j,k,YDIR) = (Mijk(2,1)*ccEijk(IDIR) + Mijk(2,2)*ccEijk(JDIR) + Mijk(2,3)*ccEijk(KDIR))/ijkDet
                    inExyz(i,j,k,ZDIR) = (Mijk(3,1)*ccEijk(IDIR) + Mijk(3,2)*ccEijk(JDIR) + Mijk(3,3)*ccEijk(KDIR))/ijkDet

                enddo
            enddo
            !Recalculate inner-most ring (degenerate metric)
            if (Model%doRing) then
                if (Model%Ring%doS) then
                    call FixRAVec_S(Model,Grid,inExyz(i,Grid%js:Grid%js+1,Grid%ks:Grid%ke,1:NDIM))
                endif
                if (Model%Ring%doE) then
                    call FixRAVec_E(Model,Grid,inExyz(i,Grid%je-1:Grid%je,Grid%ks:Grid%ke,1:NDIM))
                endif
            endif

        enddo !Shell loop
        
    end subroutine Ion2MHD

    !Fix inner shell electric fields (ijk) to remix values
    subroutine IonEFix(Model,Gr,State,inEijk)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(inout) :: State
        real(rp), intent(inout) :: inEijk(PsiSh+1,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM)

        integer :: is0,j,k

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

    !Hack to remove anomalous heating
    subroutine ChillOut(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k
        integer :: s,s0,sE
        real(rp), dimension(NVAR) :: pW, pCon
        real(rp) :: D,PrCO,CsC
        real(rp), dimension(0:Model%nSpc) :: RhoMin
        logical :: doChill

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
        !$OMP private(s,i,j,k,pW,pCon,D,doChill)
        do s=s0,sE
            do k=Grid%ksg,Grid%keg
                do j=Grid%jsg,Grid%jeg
                    do i=Grid%isg,Grid%ieg

                        D = State%Gas(i,j,k,DEN,s)
                        if (s>BLK) then
                            !Non-bulk species
                            doChill = ( D <= RhoCO ) .and. ( D >= RhoMin(s) )
                        else
                            doChill = ( D <= RhoCO )
                        endif

                        if (doChill) then
                            !Chill the fuck out
                            pCon = State%Gas(i,j,k,:,s)
                            call CellC2P(Model,pCon,pW)
                            !Set pressure to ensure Cs = CsCO
                            CsC = CsCO/gv0 !Cs in code units
                            PrCO = pW(DEN)*CsC*CsC/Model%gamma

                            pW(PRESSURE) = max(PrCO,pFloor)
                            call CellP2C(Model,pW,pCon)
                            State%Gas(i,j,k,:,s) = pCon

                        endif

                    enddo
                enddo
            enddo
        enddo

        if (Model%doMultiF) then
            call State2Bulk(Model,Grid,State)
        endif
        
    end subroutine ChillOut

    !Different hack for cooling
    !If Cs>Ca convert internal energy to mass using E=mc2, w/ c=Boris speed

    subroutine SuperChill(Model,Grid,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        type(State_T), intent(inout) :: State

        integer :: i,j,k
        real(rp), dimension(NVAR) :: pW, pCon
        real(rp) :: Cs0,alpha,Rho0,P0,e0,dRho,dE
        alpha = Model%gamma*(Model%gamma-1)

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,pW,pCon,Cs0,Rho0,P0,e0,dRho,dE)
        do k=Grid%ksg,Grid%keg
            do j=Grid%jsg,Grid%jeg
                do i=Grid%isg,Grid%ieg
                    !Check sound speed
                    pCon = State%Gas(i,j,k,:,BLK)
                    call CellC2P(Model,pCon,pW)
                    Rho0 = pW(DEN)
                    P0 = pW(PRESSURE)
                    Cs0 = sqrt(Model%gamma*P0/Rho0)
                    if (Cs0 > Model%Ca) then
                        
                        e0 = P0/(Model%gamma-1)
                        dRho = (alpha*e0/(1.0+alpha))*( (Cs0**2.0 - Model%Ca**2.0)/(Cs0*Model%Ca)**2.0)
                        dE = dRho*Model%Ca**2.0

                        !Set new primitive variables (modify velocity to conserve momentum)
                        pW(DEN) = Rho0+dRho
                        pW(VELX:VELZ) = pW(VELX:VELZ)*(Rho0/(Rho0+dRho))
                        pW(PRESSURE) = (Model%gamma-1)*(e0-dE)
                        !Convert back to conserved and store
                        call CellP2C(Model,pW,pCon)
                        State%Gas(i,j,k,:,BLK) = pCon
                        
                    endif

                enddo
            enddo
        enddo
    end subroutine SuperChill

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
        real(rp) :: rI(8),rMax,rMin,MagP

        !Get values for initial field cutoffs

        !LFM values
        !call xmlInp%Set_Val(xSun  ,"prob/xMax",20.0_rp  )
        call xmlInp%Set_Val(xSun  ,"prob/xMax",25.0_rp  )
        call xmlInp%Set_Val(xTail ,"prob/xMin",-185.0_rp)
        call xmlInp%Set_Val(yMax  ,"prob/yMax",80.0_rp  )
        call xmlInp%Set_Val(sInner,"prob/sIn" ,0.96_rp  )

        !Get cut dipole values
        call xmlInp%Set_Val(rCut,"prob/rCut",rCut)
        call xmlInp%Set_Val(lCut,"prob/lCut",lCut)
        
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
    !Calculate currents on inner-most active radial shells
    !NOTE: Only using perturbation field as we are assuming B0 is curl-free in inner region
    !Using IonG inner shells
    !See gridutils/bFld2Jxyz for full commented calculation
    !TODO: Better wrap this in routines to avoid replicated code
    subroutine GetShellJ(Model,Grid,Bxyz,Jxyz)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp), intent(in)  :: Bxyz(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp), intent(inout) :: Jxyz(1:JpSh,Grid%js:Grid%je,Grid%ks:Grid%ke,1:NDIM)

        integer :: ip0,ip1,ip2,ip3
        integer :: jp0,jp1,jp2,jp3
        integer :: kp0,kp1,kp2,kp3
        integer :: j,k,n,d,iG
        real(rp), dimension(NDIM) :: dl,bEdge !Edge vectors
        real(rp) :: bInt(1:JpSh+2,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        real(rp) :: JdS (1:JpSh+1,Grid%jsg:Grid%jeg,Grid%ksg:Grid%keg,1:NDIM)
        integer :: iR,NumP,jm,jp
        real(rp) :: Jbar

        Jxyz = 0.0

        !Start by integrating B.dl along cell edges
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(iG,n,j,k,d,dl,bEdge) &
        !$OMP private(ip0,ip1,ip2,ip3,jp0,jp1,jp2,jp3,kp0,kp1,kp2,kp3)
        do d=1,NDIM
            do k=Grid%ksg,Grid%keg
                do j=Grid%jsg,Grid%jeg
                    do n=1,JpSh+2
                        iG = JpSt+n-1

                        !Set edge vector
                        !Use ghost->active mapping to get ijk's for each 4-point average
                        select case(d)
                        case(IDIR)
                            dl = Grid%xyz(iG+1,j,k,:) - Grid%xyz(iG,j,k,:)
                            call ijk2Active(Model,Grid,iG  ,j  ,k  ,ip0,jp0,kp0)
                            call ijk2Active(Model,Grid,iG  ,j-1,k  ,ip1,jp1,kp1)
                            call ijk2Active(Model,Grid,iG  ,j  ,k-1,ip2,jp2,kp2)
                            call ijk2Active(Model,Grid,iG  ,j-1,k-1,ip3,jp3,kp3)
                        case(JDIR)
                            dl = Grid%xyz(iG,j+1,k,:) - Grid%xyz(iG,j,k,:)
                            call ijk2Active(Model,Grid,iG  ,j  ,k  ,ip0,jp0,kp0)
                            call ijk2Active(Model,Grid,iG-1,j  ,k  ,ip1,jp1,kp1)
                            call ijk2Active(Model,Grid,iG  ,j  ,k-1,ip2,jp2,kp2)
                            call ijk2Active(Model,Grid,iG-1,j  ,k-1,ip3,jp3,kp3)

                        case(KDIR)
                            dl = Grid%xyz(iG,j,k+1,:) - Grid%xyz(iG,j,k,:)
                            call ijk2Active(Model,Grid,iG  ,j  ,k  ,ip0,jp0,kp0)
                            call ijk2Active(Model,Grid,iG-1,j  ,k  ,ip1,jp1,kp1)
                            call ijk2Active(Model,Grid,iG  ,j-1,k  ,ip2,jp2,kp2)
                            call ijk2Active(Model,Grid,iG-1,j-1,k  ,ip3,jp3,kp3)
                        end select
                        
                        bEdge = 0.25*( Bxyz(ip0,jp0,kp0,:) + Bxyz(ip1,jp1,kp1,:) &
                                     +Bxyz(ip2,jp2,kp2,:) + Bxyz(ip3,jp3,kp3,:) )
                        bInt(n,j,k,d) = dot_product(bEdge,dl)

                    enddo
                enddo
            enddo
        enddo !IJK loop

        JdS = 0.0
        !Turn edge integrals into current flux through cell faces
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,j,k)
        do k=Grid%ksg,Grid%keg-1
            do j=Grid%jsg,Grid%jeg-1
                do n=1,JpSh+1
                    JdS(n,j,k,IDIR) = bInt(n,j,k,JDIR) + bInt(n,j+1,k,KDIR) - bInt(n,j,k+1,JDIR) - bInt(n,j,k,KDIR)
                    JdS(n,j,k,JDIR) = bInt(n,j,k,KDIR) + bInt(n,j,k+1,IDIR) - bInt(n+1,j,k,KDIR) - bInt(n,j,k,IDIR)
                    JdS(n,j,k,KDIR) = bInt(n,j,k,IDIR) + bInt(n+1,j,k,JDIR) - bInt(n,j+1,k,IDIR) - bInt(n,j,k,JDIR)
                enddo
            enddo
        enddo

        Jxyz = 0.0
        !Now turn surface fluxes into cell-centered XYZ using flux2field
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(n,j,k,iG)
        do k=Grid%ks,Grid%ke
            do j=Grid%js,Grid%je
                do n=1,JpSh
                    iG = JpSt+n-1
                    Jxyz(n,j,k,XDIR:ZDIR) = ( JdS(n+1,j  ,k  ,IDIR)*Grid%xfc(iG+1,j  ,k  ,:,IDIR) &
                                             +JdS(n  ,j+1,k  ,JDIR)*Grid%xfc(iG  ,j+1,k  ,:,JDIR) &
                                             +JdS(n  ,j  ,k+1,KDIR)*Grid%xfc(iG  ,j  ,k+1,:,KDIR) &
                                             -JdS(n  ,j  ,k  ,IDIR)*Grid%xfc(iG  ,j  ,k  ,:,IDIR) &
                                             -JdS(n  ,j  ,k  ,JDIR)*Grid%xfc(iG  ,j  ,k  ,:,JDIR) &
                                             -JdS(n  ,j  ,k  ,KDIR)*Grid%xfc(iG  ,j  ,k  ,:,KDIR) )/Grid%volume(iG,j,k)

                enddo
            enddo
        enddo

        !Do some smoothing on currents
        if (doRingFAC .and. Model%doRing) then
            
            NumP = Model%Ring%Np
            do n=1,JpSh
                if (Model%Ring%doS) then
                    call FixRAVec_S(Model,Grid,Jxyz(n,Grid%js:Grid%js+1,Grid%ks:Grid%ke,1:NDIM))
                endif

                if (Model%Ring%doE) then
                    call FixRAVec_E(Model,Grid,Jxyz(n,Grid%je-1:Grid%je,Grid%ks:Grid%ke,1:NDIM))
                endif

            enddo !JpSh loop

        endif !doRingFAC

    end subroutine GetShellJ
    
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

    subroutine cutDipole(x,y,z,Ax,Ay,Az)
        
        real(rp), intent(in) :: x,y,z
        real(rp), intent(out) :: Ax,Ay,Az
   
        real(rp) :: r,M
        r = sqrt(x**2.0 + y**2.0 + z**2.0)

        call Dipole(x,y,z,Ax,Ay,Az)
        M = magRampDown(r,rCut,lCut)

        Ax = M*Ax
        Ay = M*Ay
        Az = M*Az

    end subroutine cutDipole


    !Takes i,j,k cell index and returns active cell ip,jp,kp of mirror
    !Map in i,k,j order
    subroutine ijk2Active(Model,Grid,i,j,k,ip,jp,kp)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        integer, intent(in) :: i,j,k
        integer, intent(out) :: ip,jp,kp

        integer :: Np,Np2

        Np  = Grid%Nkp
        Np2 = Grid%Nkp/2

        !Start w/ i index, do mirror back into active
        if (i < Grid%is) then
            ip = Grid%is + (Grid%is-i) - 1
        elseif (i > Grid%ie) then
            ip = Grid%ie - (i-Grid%ie) + 1
        else
            ip = i
        endif

        !Next do k, map via periodicity
        if (k < Grid%ks) then
            kp = Grid%ke - (Grid%ks-k) + 1
        elseif (k > Grid%ke) then
            kp = Grid%ks + (k-Grid%ke) - 1
        else
            kp = k
        endif

        !Finally do j
        if (j < Grid%js) then
            jp = Grid%js + (Grid%js-j) - 1
            kp = k+Np2
            if (kp>Np) kp = kp-Np
        elseif (j > Grid%je) then
            jp = Grid%je - (j-Grid%je) + 1
            kp = k+Np2
            if (kp>Np) kp = kp-Np
        else
            jp = j
        endif

    end subroutine ijk2Active

    subroutine Heat(Model,Grid,State)
      type(Model_T), intent(in) :: Model
      type(Grid_T), intent(in) :: Grid
      type(State_T), intent(inout) :: State

      integer :: i,j,k
      real(rp), dimension(NVAR) :: pW, pCon

!      if (all(State%eqPres.eq.0.0))  then
!         ! TODO pack this and the next loop together (Maybe?)
!         call getEqPress(Model,Grid,State)
!      endif
      
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,k,pW,pCon)
      do k=Grid%ks,Grid%ke
         do j=Grid%js,Grid%je
            do i=Grid%is,Grid%ie
               pCon = State%Gas(i,j,k,:,BLK)
               call CellC2P(Model,pCon,pW)
               ! only heat -- don't cool
               pW(PRESSURE) = pW(PRESSURE)+Model%hRate/Model%hTau*Model%dt*max(0.,(State%eqPres(i,j,k)-pW(PRESSURE)))

               ! always on if heating
               if (Model%doPsphere) then
!                  pW(DEN) = pW(DEN)+ State%eqDen(i,j,k)
                  pW(DEN) = pW(DEN)+Model%hRate/Model%hTau*Model%dt*max(0.,(State%eqDen(i,j,k)-pW(DEN)))
               end if

               call CellP2C(Model,pW,pCon)
               State%Gas(i,j,k,:,BLK) = pCon
            enddo
         enddo
      enddo    
    end subroutine Heat

    !Set gPsi (corotation potential)
    subroutine CorotationPot(Model,Grid,gPsi)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Grid
        real(rp)  , intent(inout) :: gPsi(1:PsiSh+1,Grid%js:Grid%je+1,Grid%ks:Grid%ke+1)

        integer :: n,i,iG,j,k

        real(rp), dimension(NDIM) :: xyz
        real(rp) :: r, lambda

        if (abs(Psi0)<=TINY) then
            return
        endif
        
        do k=Grid%ks,Grid%ke+1
            do j=Grid%js,Grid%je+1
                do i=1,PsiSh+1
                    !i = Grid%is+PsiSt-1
                    iG = i+PsiSt-1 !Global i index for Grid access

                    xyz = Grid%xyz(iG,j,k,:)
                    r = norm2(xyz)
                    lambda = acos(xyz(ZDIR)/r)

                    gPsi(i,j,k) = gPsi(i,j,k)+Psi0*cos(lambda)*cos(lambda)/r                    
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

    !Calculate invariant latitude (RADIANS) for x,y,z vector (in Rx)
    function InvLatitude(r) result(invlat)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: invlat

        real(rp) :: z,rad,lat,Leq

        z = r(ZDIR)
        rad = norm2(r)

        lat = abs( asin(z/rad) )
        Leq = rad/( cos(lat)*cos(lat) )
        invlat = abs(acos(sqrt(1.0/Leq)))

    end function InvLatitude
  
end module msphutils
