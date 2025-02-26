!Routines to handle source term ingestion in magnetosphere runs

module msphingest
    
    use kdefs
    use gamtypes
    use imaghelper
    use earthhelper
    use planethelper
    use volttypes
    use gamutils
    use geopack
    use multifluid
    use gridutils

    implicit none

    !Ingestion switch
    logical , private :: doIngest = .true. !Whether to ignore ingestion value, ie 1-way coupling to RCM
    logical , private :: doAppetizer = .false. !Whether to ingest plasmasheet values during spinup
    real(rp), private :: dtAppetizer = 300 ![s], fiducial timescale for plasmasheet ingestion

    
    !Parameters for appetizer
    type Gas0App_T
        real(rp) :: dz = 5.0 !Wedge around equator to use
        real(rp) :: tScl !How to scale ingestion timescale
    end type Gas0App_T

    type(Gas0App_T), private :: Gas0App

    contains

    !Set ingestion parameters
    subroutine setIngestion(Model,xmlInp,pID)
        type(Model_T), intent(inout) :: Model
        type(XML_Input_T), intent(in) :: xmlInp
        character(len=*) , intent(in) :: pID

        real(rp) :: t0
        character(len=strLen) :: wID !Wind ID string

        if (.not. Model%doSource) return

        !Whether to ignore ingestion (if set)
        call xmlInp%Set_Val(doIngest,"source/doIngest",.true.)
        call xmlInp%Set_Val(doAppetizer,"/Kaiju/voltron/imag/doInit",.false.)

        if (doAppetizer) then
            call xmlInp%Set_Val(wID,"wind/tsfile","NONE")
            t0 = TINY !Just setting value as T=+0
            call xmlInp%Set_Val(Gas0App%dz ,"source0/dz" ,Gas0App%dz)
            call xmlInp%Set_Val(dtAppetizer,"source0/dt0",dtAppetizer)
            call SetTM03(Model,wID,t0)

        endif

        !------
        contains
            !TODO: Properly handle GSM rotation
            subroutine SetTM03(Model,wID,t0)
                type(Model_T), intent(in) :: Model
                character(len=*), intent(in) :: wID
                real(rp), intent(in) :: t0

                real(rp) :: D0,P0,Tau0,xyz(NDIM)
                logical  :: isIn

                !Setup TM03
                call InitTM03(wID,t0)

                !Setup ingestion timescale
                xyz = [-10.0_rp-TINY,0.0_rp,0.0_rp]
                call Appetizer_TM03(Model,xyz,D0,P0,Tau0,isIn)
                Gas0App%tScl = dtAppetizer/(Tau0*Model%Units%gT0)

            end subroutine SetTM03

    end subroutine setIngestion


!-----
    !Ingest density/pressure information from Grid%Gas0
    !Treat Gas0 as target value
    subroutine MagsphereIngest(Model,Gr,State)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr
        type(State_T), intent(inout) :: State

        integer :: i,j,k
        real(rp), dimension(NVAR) :: pCon,pConRC,pConPS,pW
        real(rp), dimension(NDIM) :: E,B,Veb
        real(rp) :: D0,P0,Tau,rcD,rcP,psD,psP
        logical  :: doIngestIJK,doIngestRC,doIngestCOLD

    !Do traps
        if (.not. doIngest) return
        if ( (Model%t<=0) .and. (.not. doAppetizer) ) return !You'll spoil your appetite

        if ( Model%doMultiF .and. (Model%nSpc<COLDFLUID) ) then
            write(*,*) "Not enough fluids to hold cold"
            stop
        endif

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,E,B,Veb,Tau)     &
        !$OMP private(D0,P0,rcD,rcP,psD,psP) &
        !$OMP private(doIngestIJK,doIngestRC,doIngestCOLD) &
        !$OMP private(pCon,pConRC,pConPS,pW)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie

                !Do relevant things for both multiF and singleF
                    !Get local field
                    B = State%Bxyz(i,j,k,:)
                    if (Model%doBackground) then
                        B = B + Gr%B0(i,j,k,:)
                    endif

                    !Get ingestion timescale
                    Tau = Gr%Gas0(i,j,k,IM_TSCL)
                    if (Tau<Model%dt) Tau = Model%dt !Unlikely to happen

                    !Get ingestion frame
                    if (Model%t < 0) then
                        !Use current fluid frame
                        pCon = State%Gas(i,j,k,:,BLK)
                        call CellC2P(Model,pCon,pW)
                        Veb = pW(VELX:VELZ)
                    else
                        !Use ExB
                        E = Gr%Gas0(i,j,k,IONEX:IONEZ) !Ionospheric E field
                        Veb = cross(E,B)/dot_product(B,B)
                    endif

                !Single fluid
                    if (.not. Model%doMultiF) then
                        D0 = Gr%Gas0(i,j,k,IM_D_RING) + Gr%Gas0(i,j,k,IM_D_COLD)
                        P0 = Gr%Gas0(i,j,k,IM_P_RING) + Gr%Gas0(i,j,k,IM_P_COLD)

                        pCon = State%Gas(i,j,k,:,BLK)
                        call Ingest2Con(Model,pCon,doIngestIJK,D0,P0,Veb,Tau)

                        if (doIngestIJK) then
                            State%Gas(i,j,k,:,BLK) = pCon
                        endif
                !Multi-F
                    else
                        !Try RC ingestion
                        rcD = Gr%Gas0(i,j,k,IM_D_RING)
                        rcP = Gr%Gas0(i,j,k,IM_P_RING)
                        pConRC = State%Gas(i,j,k,:,RCFLUID)
                        call Ingest2Con(Model,pConRC,doIngestRC,rcD,rcP,Veb,Tau)

                        !Try plasmasphere ingestion
                        psD = Gr%Gas0(i,j,k,IM_D_COLD)
                        psP = Gr%Gas0(i,j,k,IM_P_COLD)
                        pConPS = State%Gas(i,j,k,:,COLDFLUID)
                        call Ingest2Con(Model,pConPS,doIngestCOLD,psD,psP,Veb,Tau)

                        !Finish up
                        if (doIngestRC) then
                            State%Gas(i,j,k,:,RCFLUID) = pConRC
                        endif

                        if (doIngestCOLD) then
                            State%Gas(i,j,k,:,COLDFLUID) = pConPS
                        endif

                        if (doIngestRC .or. doIngestCOLD) then
                            !Rebuild BLK
                            call MultiF2Bulk(Model,State%Gas(i,j,k,:,:),B)
                        endif

                    endif

                enddo
            enddo
        enddo !k


    !---
        contains
            !Do some shenanigans to keep hot sound speed low where it doesn't matter
            subroutine ClampHot(Model,psD,psP,rcD,rcP)
                type(Model_T), intent(in)    :: Model
                real(rp)     , intent(in)    :: psD,psP,rcD
                real(rp)     , intent(inout) :: rcP

                real(rp) :: CsT,CsRC,rcP0,xFac

                if (rcP > psP) return !There's pressure here don't bother
                !If still here then this region can't even compete with stupid plasmasphere heat
                !Get that nonsense outta here

                rcP0 = rcP
                xFac = rcP/psP ! < 1
                if (Model%doBoris) then
                    CsT = xFac*Model%Ca
                else
                    CsT = HUGE !Not sure what to do for not boris
                endif

                CsRC = sqrt(Model%gamma*rcP/rcD) !RC fluid sound speed
                if (CsT < CsRC) then
                    !Reset rcP to hit target sound speed
                    rcP = max( rcD*(CsT**2.0)/Model%gamma,100.0*pFloor )
                    if (rcP > rcP0) rcP=rcP0
                endif
                
            end subroutine ClampHot

    end subroutine MagsphereIngest

    !Take conserved state (pCon) and ingest in place
    subroutine Ingest2Con(Model,pCon,doIngest,D0,P0,V0,Tau)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: pCon(NVAR)
        real(rp), intent(in)    :: D0,P0,Tau
        real(rp), intent(in)    :: V0(NDIM)
        logical , intent(out)   :: doIngest

        real(rp) :: pW(NVAR)
        real(rp) :: Dmhd,Pmhd,dRho,dP
        real(rp), dimension(NDIM) :: Mxyz,Vmhd,dMom

        !Decide if there's anything to eat
        doIngest = (D0>dFloor) .and. (P0>pFloor)
        if (.not. doIngest) return

        !Now we feast
        call CellC2P(Model,pCon,pW)

        Dmhd = pW(DEN)
        Pmhd = pW(PRESSURE)
        Mxyz = pCon(MOMX:MOMZ) !Momentum
        Vmhd = pW(VELX:VELZ)

        dRho = (D0-Dmhd)*(Model%dt/Tau)
        dP   = (P0-Pmhd)*(Model%dt/Tau)
        !Calculate change in momentum
        if (dRho > 0) then
            !If gaining mass do it in the ionospheric convection frame
            dMom = dRho*V0 !Change in momentum
        else
            !Lose mass in MHD frame
            dMom = dRho*Vmhd
        endif

        !Now calculate new primitive values
        pW(DEN) = pW(DEN) + dRho
        pW(VELX:VELZ) = ( Mxyz + dMom )/(pW(DEN)) !New M / new D
        pW(PRESSURE) = pW(PRESSURE) + dP

        !Now put back
        call CellP2C(Model,pW,pCon)
        
    end subroutine Ingest2Con

    subroutine CellPlasmaBC(Model,gU,pU,Gas0,Veb,rHat,B)
        type(Model_T), intent(in) :: Model
        real(rp), intent(inout) :: gU(NVAR,BLK:Model%nSpc)
        real(rp), intent(in)    :: pU(NVAR,BLK:Model%nSpc)
        real(rp), intent(in)    :: Gas0(NVARVOLTSRC)
        real(rp), dimension(NDIM), intent(in) :: Veb,rHat,B
        logical :: isIMCold,isIMRing
        real(rp), dimension(NVAR) :: pW,gW,pCon
        real(rp) :: Vr


        if ( (Model%t>0) .and. (Model%doSource) .and. (doIngest) ) then
            !Test for imag information
            if (Model%doMultiF) then
                !Use vacuum thresholds
                isIMRing = (Gas0(IM_D_RING) >= Spcs(RCFLUID)  %dVac)
                isIMCold = (Gas0(IM_D_COLD) >= Spcs(COLDFLUID)%dVac)
            else
                !Not MF, just use raw floors
                isIMRing = (Gas0(IM_D_RING) >= dFloor) .and. (Gas0(IM_P_RING) >= pFloor)
                isIMCold = (Gas0(IM_D_COLD) >= dFloor) .and. (Gas0(IM_P_COLD) >= pFloor)
            endif

        else
            !No useful info
            isIMRing = .false.
            isIMCold = .false.
        endif

        !General strategy: When we have informed values in ghost (ie, raiju plasmasphere) use full Veb
        !If we have good active cell info, use diode on Vr (allow inflow but not inflow)
        !If we got nothing, use Veb w/ Vr=0.0

        Vr = dot_product(Veb,rHat)

        if (Model%doMultiF) then

        !Do cold fluid
            gW(:) = 0.0 !Ghost prim quantities
            call CellC2P(Model,pU(:,COLDFLUID),pW) !Active cell prim quantities

            if (isIMCold) then
                !Have good source info in ghost
                gW(DEN)       = Gas0(IM_D_COLD)
                gW(PRESSURE)  = Gas0(IM_P_COLD)
                gW(VELX:VELZ) = Veb
            elseif (pW(DEN) > Spcs(COLDFLUID)%dVac) then
                !Have good active cell info
                gW(DEN)       = pW(DEN)     
                gW(PRESSURE)  = pW(PRESSURE)
                gW(VELX:VELZ) = Vec2Perp(Veb,rHat) + min(Vr,0.0)*rHat
            else
                !Don't got nothing
                gW(DEN)       = dFloor
                gW(PRESSURE)  = pFloor
                gW(VELX:VELZ) = Vec2Perp(Veb,rHat)
            endif
            !Have ghost prim values, convert to con and store
            call CellP2C(Model,gW,gU(:,COLDFLUID))

        !Do 'other' fluid
            gW(:) = 0.0 !Ghost prim quantities
            call CellC2P(Model,pU(:,RCFLUID),pW) !Active cell prim quantities

            if (isIMRing) then
                !Have good source info in ghost
                gW(DEN)       = Gas0(IM_D_RING)
                gW(PRESSURE)  = Gas0(IM_P_RING)
                gW(VELX:VELZ) = Veb
            elseif (pW(DEN) > Spcs(RCFLUID)%dVac) then
                !Have good active cell info
                gW(DEN)       = pW(DEN)     
                gW(PRESSURE)  = pW(PRESSURE)
                gW(VELX:VELZ) = Vec2Perp(Veb,rHat) + min(Vr,0.0)*rHat
            else
                !Got nothin
                gW(DEN)       = dFloor
                gW(PRESSURE)  = pFloor
                gW(VELX:VELZ) = Vec2Perp(Veb,rHat)
            endif
            !Have ghost prim values, convert to con and store
            call CellP2C(Model,gW,gU(:,RCFLUID))

        !Now set BLK
            call MultiF2Bulk(Model,gU,B)
        else
            !Single fluid
            gW(:) = 0.0
            call CellC2P(Model,pU(:,BLK),pW) !Active cell prim quantities
            if (isIMRing .or. isIMCold) then
                gW(DEN)      = Gas0(IM_D_COLD) + Gas0(IM_D_RING)
                gW(PRESSURE) = Gas0(IM_P_COLD) + Gas0(IM_P_RING)
                gW(VELX:VELZ) = Veb
            else
                !No data, just use zero grad based on active mhd info
                gW(DEN)      = pW(DEN)     
                gW(PRESSURE) = pW(PRESSURE)
                gW(VELX:VELZ) = Vec2Perp(Veb,rHat) + min(Vr,0.0)*rHat
            endif

            call CellP2C(Model,gW,gU(:,BLK))

        endif !MF v single

    end subroutine CellPlasmaBC

    !Loads Gas0 w/ t<0 ingestion values
    !Note: tSW is seconds elapsed (not gamera time units)
    subroutine LoadSpinupGas0(Model,Gr,tSW)
        type(Model_T), intent(in)    :: Model
        type(Grid_T) , intent(inout) :: Gr
        real(rp)     , intent(in)    :: tSW

        integer :: i,j,k
        real(rp) :: D0,P0,Tau
        logical  :: doIngestIJK,doInD,doInP
        real(rp), dimension(NDIM) :: xyzSM,xyzGSM

        if ( (Model%t<=0) .and. (.not. doAppetizer) ) return !You'll spoil your appetite

        call UpdateTM03(tSW)
        !Start by setting geopack for transformation
        call mjdRecalc( TM03_MJD() )

       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP private(i,j,k,doInD,doInP,doIngestIJK)  &
       !$OMP private(D0,P0,Tau,xyzSM,xyzGSM)
        do k=Gr%ksg,Gr%keg
            do j=Gr%jsg,Gr%jeg
                do i=Gr%isg,Gr%ieg
                    !Get GSM coordinates
                    xyzSM = Gr%xyzcc(i,j,k,XDIR:ZDIR)
                    call SM2GSW(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGSM(XDIR),xyzGSM(YDIR),xyzGSM(ZDIR))
                    
                    call Appetizer_TM03(Model,xyzGSM,D0,P0,Tau,doIngestIJK)

                    doInD = D0>dFloor
                    doinP = P0>pFloor
                    doIngestIJK = doInD .and. doInP .and. doIngestIJK
                    Gr%Gas0(i,j,k,:) = 0.0

                    if (doIngestIJK) then
                        Gr%Gas0(i,j,k,IM_D_RING)  = D0
                        Gr%Gas0(i,j,k,IM_P_RING)  = P0
                        Gr%Gas0(i,j,k,IM_TSCL  )  = Tau*Gas0App%tScl
                        !Leave D/P_COLD 0
                    endif

                enddo
            enddo
        enddo

    end subroutine LoadSpinupGas0
!-----

    !TM03 as an appetizer, provide D/P [code units] for a given XYZ and ingestion timescale [code]
    !isIn is whether the value is edible
    subroutine Appetizer_TM03(Model,xyzIN,D0,P0,Tau0,isIn)
        type(Model_T), intent(in) :: Model
        real(rp), intent(in)    :: xyzIN(NDIM)
        real(rp), intent(out)   :: D0,P0,Tau0
        logical , intent(inout) :: isIn

        real(rp) :: R,rho
        real(rp) :: D,P,Tau,Tev,Cs
        real(rp) :: xyz(NDIM)

        !Initialize
        D0 = 0.0
        P0 = 0.0
        Tau0 = 0.0
        isIn = .false.

        if (.not. inShueMP(xyzIn)) return

        R = norm2(xyzIN)
        rho = norm2([xyzIN(XDIR),xyzIN(YDIR)])

        if (xyzIN(XDIR)>0) then
            !Map back to Shue MP at terminator
            call ShueMP2Terminator(xyzIN,xyz)
        else
            xyz = xyzIN
        endif

        isIn = inShueMP(xyz) .and. (rho>10) .and. (xyz(XDIR)<=0) &
               .and. (xyz(XDIR)>-50) .and. (abs(xyz(ZDIR))<Gas0App%dz)

        call EvalTM03(xyz,D,P,isIn)

        if (.not. isIn) return

        !Get timescale [s], sonic bounce
        !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
        Tev = (1.0e+3)*DP2kT(D,P) !Temp in eV
        Cs = 9.79*sqrt((5.0/3)*TeV)

        Tau = (DipoleL(xyz)*Model%Units%gx0*1.0e-3)/Cs

        !Now have D,P,Tau in physical units. Convert back to code
        D0   = D !Magnetosphere is already in #/cc
        P0   = P/Model%Units%gP0
        Tau0 = Tau/Model%Units%gT0

    end subroutine Appetizer_TM03

end module msphingest
