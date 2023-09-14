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

	implicit none

    !Ingestion switch
    logical , private :: doIngest = .true. !Whether to ignore ingestion value, ie 1-way coupling to RCM
    logical , private :: doAppetizer = .false. !Whether to ingest plasmasheet values during spinup
    real(rp), private :: dtAppetizer = 300 ![s], fiducial timescale for plasmasheet ingestion

    
    !Parameters for appetizer
    type Gas0App_T
        logical  :: doInit = .true.
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
        real(rp), dimension(NVAR) :: pW, pCon

        real(rp) :: Tau,dRho,dP,Pmhd,Prcm
        real(rp), dimension(NDIM) :: Mxyz,Vxyz
        real(rp) :: D0,P0
        logical  :: doIngestIJK,doInD,doInP

    !Do traps
        if (.not. doIngest) return
        if ( (Model%t<=0) .and. (.not. doAppetizer) ) return !You'll spoil your appetite

    	!if ( (Model%t<0) .and. (Model%ts>1) .and. doAppetizer .and. TM03%doInit) then
        !For now redoing every 100 timesteps to make sure voltron doesn't overwrite (bit silly)
        !TODO: Remove this after overhauling source ingestion
        if ( (Model%t<0) .and. (modulo(Model%ts,100) == 0) .and. doAppetizer ) then
    		!Load TM03 into Gas0 for ingestion during spinup
    		call LoadSpinupGas0(Model,Gr)
    	endif

        if (Model%doMultiF) then
            write(*,*) 'Source ingestion not implemented for multifluid, you should do that'
            stop
        endif

    !Now real ingestion work
       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP private(i,j,k,doInD,doInP,doIngestIJK,pCon,pW) &
       !$OMP private(Tau,dRho,dP,Pmhd,Prcm,Mxyz,Vxyz,D0,P0)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie

            		D0 = Gr%Gas0(i,j,k,IMDEN,BLK)
            		P0 = Gr%Gas0(i,j,k,IMPR ,BLK)
            		Tau = Gr%Gas0(i,j,k,IMTSCL,BLK)
            		doInD = D0>dFloor
            		doinP = P0>pFloor
            		doIngestIJK = doInD .or. doInP

                    if (.not. doIngestIJK) cycle

                    pCon = State%Gas(i,j,k,:,BLK)
                    call CellC2P(Model,pCon,pW)
                    Pmhd = pW(PRESSURE)
                    Mxyz = pCon(MOMX:MOMZ) !Classical momentum
                    Vxyz = pW(VELX:VELZ) !Velocity pre-ingestion

                    if (Tau<Model%dt) Tau = Model%dt !Unlikely to happen

                    if (doInD) then
                        dRho = D0 - pW(DEN)
                        pW(DEN) = pW(DEN) + (Model%dt/Tau)*dRho
                        if (dRho>0) then
                            !If we're gaining mass, conserve momentum
                            !Don't increase speed if mass decreases
                            Vxyz = Mxyz/pW(DEN)
                        endif
                    endif

                    if (doInP) then
                        Prcm = P0
                        !Assume already wolf-limited or not
                        dP = Prcm - Pmhd
                        pW(PRESSURE) = pW(PRESSURE) + (Model%dt/Tau)*dP
                    endif

                    !Preserve velocity during ingestion, ie ingesting in the moving frame
                    pW(VELX:VELZ) = Vxyz

                    !Now put back
                    call CellP2C(Model,pW,pCon)
                    State%Gas(i,j,k,:,BLK) = pCon
                enddo
            enddo
        enddo
                    
    end subroutine MagsphereIngest

    !Loads Gas0 w/ t<0 ingestion values
    subroutine LoadSpinupGas0(Model,Gr)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(inout) :: Gr

        integer :: i,j,k
        real(rp) :: D0,P0,Tau
        logical  :: doIngestIJK,doInD,doInP
        real(rp), dimension(NDIM) :: xyzSM,xyzGSM

        !Start by setting geopack for transformation
	    call mjdRecalc( TM03_MJD() )

       !$OMP PARALLEL DO default(shared) collapse(2) &
       !$OMP private(i,j,k,doInD,doInP,doIngestIJK)  &
       !$OMP private(D0,P0,Tau,xyzSM,xyzGSM)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is,Gr%ie
                	!Get GSM coordinates
                	xyzSM = Gr%xyzcc(i,j,k,XDIR:ZDIR)
                	call SM2GSW(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGSM(XDIR),xyzGSM(YDIR),xyzGSM(ZDIR))
            		
            		call Appetizer_TM03(Model,xyzGSM,D0,P0,Tau,doIngestIJK)

            		doInD = D0>dFloor
            		doinP = P0>pFloor
            		doIngestIJK = doInD .and. doInP .and. doIngestIJK
            		if (doIngestIJK) then
            			Gr%Gas0(i,j,k,IMDEN,BLK)  = D0
            			Gr%Gas0(i,j,k,IMPR ,BLK)  = P0
            			Gr%Gas0(i,j,k,IMTSCL,BLK) = Tau*Gas0App%tScl
            		else
            			Gr%Gas0(i,j,k,IMDEN,BLK)  = 0.0
            			Gr%Gas0(i,j,k,IMPR ,BLK)  = 0.0
            			Gr%Gas0(i,j,k,IMTSCL,BLK) = 0.0
            		endif

               	enddo
            enddo
        enddo

        !We're done now
        Gas0App%doInit = .false.
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
			!Map back to Shue MP at dusk
			call ShueMP2Dusk(xyzIN,xyz)
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
