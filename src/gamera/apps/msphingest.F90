!Routines to handle source term ingestion in magnetosphere runs

module msphingest
	
	use kdefs
    use gamtypes
    use earthhelper
    use volttypes
    use gamutils
    use geopack

	implicit none

    !Ingestion switch
    logical, private :: doIngest = .true. !Whether to ignore ingestion value, ie 1-way coupling to RCM
    logical, private :: doAppetizer = .false. !Whether to ingest plasmasheet values during spinup
    logical, private :: dtAppetizer = 300 ![s], fiducial timescale for plasmasheet ingestion
    !Coefficients for TM03
    real(rp), dimension(11), parameter, private :: AN = [-0.15930, 0.60805,0.50555,0.07959,0.27462,0.03611, &
    											         -0.03419,-0.79347,1.16222,0.47559,0.71166]

    real(rp), dimension(14), parameter, private :: AP = [ 0.05701,0.52401,0.09075,0.52720,0.07814,-4.42226,-1.53325, &
    											         -1.21666,2.53964,0.31988,0.75433,1.04862,-0.07390,1.01540]

    

    !Parameters for appetizer, ie solar wind values for TM03
    type TM03_T
    	!Vectors should be in GSM
    	!#/cc,km/s,nPa,nT
        real(rp) :: dSW,vSW,PdynSW,BxSW,BySW,BzSW,MJD0
        logical  :: doInit = .true.
        real(rp) :: dz = 5.0 !Wedge around equator to use
        real(rp) :: tScl !How to scale ingestion timescale
    end type TM03_T

	type(TM03_T), private :: TM03

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
        call xmlInp%Set_Val(doAppetizer,"/voltron/imag/doInit",.false.)

        if (doAppetizer) then
        	call xmlInp%Set_Val(wID,"wind/tsfile","NONE")
        	t0 = TINY !Just setting value as T=+0
        	call xmlInp%Set_Val(TM03%dz,"TM03/dz",TM03%dz)
        	call xmlInp%Set_Val(dtAppetizer,"TM03/dt0",dtAppetizer)
        	call SetTM03(Model,wID,t0)

        endif

        !------
        contains
        	!TODO: Properly handle GSM rotation
        	subroutine SetTM03(Model,wID,t0)
        		type(Model_T), intent(in) :: Model
	            character(len=*), intent(in) :: wID
	            real(rp), intent(in) :: t0

	            real(rp) :: D,Vx,Vy,Vz,Bx,By,Bz
	            real(rp) :: D0,P0,Tau0,xyz(NDIM)
	            logical  :: isIn

	            !Get values from SW file
	            D = GetSWVal("D",wID,t0) !#/cc
	            Vx = GetSWVal("Vx",wID,t0)*1.0e-3 !km/s
	            Vy = GetSWVal("Vy",wID,t0)*1.0e-3 !km/s
	            Vz = GetSWVal("Vz",wID,t0)*1.0e-3 !km/s
	            Bx = GetSWVal("Bx",wID,t0) !nT
	            By = GetSWVal("By",wID,t0) !nT
	            Bz = GetSWVal("Bz",wID,t0) !nT
	            TM03%MJD0 = GetSWVal("MJD",wID,t0)
	            !Start storing, want GSM
	            TM03%dSW  = D
	            TM03%vSW  = norm2([Vx,Vy,Vz])
	            TM03%PdynSW = PV2PDyn(TM03%dSW,TM03%vSW)
	            !Do B field transform

	            call mjdRecalc(TM03%MJD0)
	            call SM2GSW(Bx,By,Bz,TM03%BxSW,TM03%BySW,TM03%BzSW)

	            !Setup ingestion timescale
	            xyz = [-10.0_rp-TINY,0.0_rp,0.0_rp]
	            call Appetizer_TM03(Model,xyz,D0,P0,Tau0,isIn)
	            TM03%tScl = dtAppetizer/(Tau0*Model%Units%gT0)

        	end subroutine SetTM03

        	!Very inefficiently get values from solar wind file
	        function GetSWVal(vID,fID,t0) result(qSW)
	            character(len=*), intent(in) :: vID,fID
	            real(rp), intent(in) :: t0
	            real(rp) :: qSW

	            type(TimeSeries_T) :: tsQ
	            tsQ%wID = trim(fID)
	            call tsQ%initTS(trim(vID),doLoudO=.false.)
	            qSW = tsQ%evalAt(t0)
	        end function GetSWVal

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

    	if ( (Model%t<0) .and. (Model%ts>1) .and. doAppetizer .and. TM03%doInit) then
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
                                        
                    if (Tau<Model%dt) Tau = Model%dt !Unlikely to happen

                    if (doInD) then
                        dRho = D0 - pW(DEN)
                        pW(DEN) = pW(DEN) + (Model%dt/Tau)*dRho
                    endif

                    if (doInP) then
                        Prcm = P0
                        !Assume already wolf-limited or not
                        dP = Prcm - Pmhd
                        pW(PRESSURE) = pW(PRESSURE) + (Model%dt/Tau)*dP
                    endif

                    Vxyz = Mxyz/max(pW(DEN),dFloor) !Conserve classical momentum
                    !Don't allow mass ingestion to speed things up
                    if ( norm2(Vxyz) <= norm2(pW(VELX:VELZ)) ) then
                        !Conserve momentum if it doesn't increase speed
                        pW(VELX:VELZ) = Vxyz
                    endif

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
	    call mjdRecalc(TM03%MJD0)

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
            			Gr%Gas0(i,j,k,IMTSCL,BLK) = Tau*TM03%tScl
            		else
            			Gr%Gas0(i,j,k,IMDEN,BLK)  = 0.0
            			Gr%Gas0(i,j,k,IMPR ,BLK)  = 0.0
            			Gr%Gas0(i,j,k,IMTSCL,BLK) = 0.0
            		endif

               	enddo
            enddo
        enddo

        !We're done now
        TM03%doInit = .false.
    end subroutine LoadSpinupGas0
!-----
	!TM03 plasma sheet model, adapted from Matlab code by Meg Noah adapted from TM03

	!TODO: Properly handle SM vs GSM
	!Vectors should properly be in GSM
	function inShueMP(xyz)
		real(rp), intent(in) :: xyz(NDIM)
		logical :: inShueMP

		real(rp) :: R,theta,CosT,R0,alpha,Rm

		R = norm2(xyz)
		theta = acos(xyz(XDIR)/R)
		CosT = cos(theta)
		
		R0 = (10.22+1.29*tanh(0.184*(TM03%BzSW+8.14)))*(TM03%PdynSW)**(-1.0/6.6)
		alpha = (0.58-0.007*TM03%BzSW)*(1.0+0.024*log(TM03%PdynSW))
		RM = R0*(2.0/max(1.0+CosT,TINY))**(alpha)
		if (R<RM) then
			inShueMP = .true.
		else
			inShueMP = .false.
		endif

	end function inShueMP

	!Map xyz coordinates to dusk along Shue MP contours
	subroutine ShueMP2Dusk(xyz,xyzD)
		real(rp), intent(in)  :: xyz (NDIM)
		real(rp), intent(out) :: xyzD(NDIM)

		real(rp) :: R,theta,CosT,R0,alpha,rScl
		xyzD = 0.0

		R = norm2(xyz)
		theta = acos(xyz(XDIR)/R)
		CosT = cos(theta)
		alpha = (0.58-0.007*TM03%BzSW)*(1.0+0.024*log(TM03%PdynSW))

		!Solve for R0
		rScl = (2.0/max(1.0+CosT,TINY))**(alpha)
		R0 = R/rScl

		xyzD = [0.0_rp,R0*(2.0**alpha),xyz(ZDIR)]

	end subroutine ShueMP2Dusk

	!TM03 as an appetizer, provide D/P [code units] for a given XYZ and ingestion timescale [code]
	!isIn is whether the value is edible
	subroutine Appetizer_TM03(Model,xyzIN,D0,P0,Tau0,isIn)
		type(Model_T), intent(in) :: Model
		real(rp), intent(in)    :: xyzIN(NDIM)
		real(rp), intent(out)   :: D0,P0,Tau0
		logical , intent(inout) :: isIn

		real(rp) :: R,clockang,Fst,Pst,rho,rho10,Nst,BNst,BSst,phi
		real(rp) :: D,P,Tau,Tev,Cs
		real(rp) :: xyz(NDIM)

		!Initialize
		D0 = 0.0
		P0 = 0.0
		Tau0 = 0.0

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
			   .and. (xyz(XDIR)>-50) .and. (abs(xyz(ZDIR))<TM03%dz)

		if (.not. isIn) return

		!If still here, get to work
		clockang = katan2(TM03%BySW,TM03%BzSW)
		Fst = norm2([TM03%BySW,TM03%BzSW])*sin(clockang/2)/5.0
		Pst = TM03%PdynSW/3.0
		
		rho10 = rho/10.0
		Nst = TM03%dSW/10.0
		BNst = TM03%BzSW/5.0
		if (BNst<0) BNst = 0.0
		BSst = (-TM03%BzSW/5.0)*(TM03%vSW/500.0)
		if (BSst<0) BSst = 0.0

		phi = katan2(xyz(YDIR),-xyz(XDIR))

		!Get TS D [#/cc], P [nPa]
		D = (AN(1) + AN(2)*(Nst**AN(10)) + AN(3)*BNst + AN(4)*BSst)*(RHO10**AN(8)) + &
		    (AN(5)*(Nst**AN(11)) + AN(6)*BNst + AN(7)*BSst)*(RHO10**AN(9))*(sin(phi)**2.0)


		P =  AP(1)*(RHO10**AP(6)) + &
		     AP(2)*(Pst**AP(11))*(RHO10**AP(7)) + &
		     AP(3)*(Fst**AP(12))*(RHO10**AP(8)) + &
		    (AP(4)*(Pst**AP(13))*exp(-AP( 9)*RHO10) + &
		     AP(5)*(Fst**AP(14))*exp(-AP(10)*RHO10))*(sin(phi)**2.0)

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
