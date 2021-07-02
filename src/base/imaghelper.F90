!Various useful empirical inner magnetosphere models
module imaghelper
    use kdefs
    use kronos
    use earthhelper
    use math
    use geopack

    implicit none

    !Stuff to calculate TM03 plasma sheet model
    !TM03 plasma sheet model, adapted from Matlab code by Meg Noah adapted from TM03

    type TM03_T
    	!Vectors should be in GSM
    	!#/cc,km/s,nPa,nT
        real(rp) :: dSW,vSW,PdynSW,BxSW,BySW,BzSW,MJD0
        logical :: isInit = .false. !Has it been initialized
        type(TimeSeries_T) :: tsD,tsVx,tsVy,tsVz,tsBx,tsBy,tsBz,tsMJD
    end type TM03_T

	type(TM03_T), private :: TM03
    !Coefficients for TM03
    real(rp), dimension(11), parameter, private :: AN = [-0.15930, 0.60805,0.50555,0.07959,0.27462,0.03611, &
    											         -0.03419,-0.79347,1.16222,0.47559,0.71166]

    real(rp), dimension(14), parameter, private :: AP = [ 0.05701,0.52401,0.09075,0.52720,0.07814,-4.42226,-1.53325, &
    											         -1.21666,2.53964,0.31988,0.75433,1.04862,-0.07390,1.01540]
  
    contains

    !Initialize TM03 from solar wind file (wID = filename)
    subroutine InitTM03(wID,t0)
    	character(len=*), intent(in) :: wID
    	real(rp), intent(in) :: t0

    	if (TM03%isInit) return !Already initialized

    	write(*,*) "Initializing TM03 model ..."

    	!Get time series for TM03 object
    	call SetupTS("D" ,wID,TM03%tsD )
    	call SetupTS("Vx",wID,TM03%tsVx)
    	call SetupTS("Vy",wID,TM03%tsVy)
    	call SetupTS("Vz",wID,TM03%tsVz)
    	call SetupTS("Bx",wID,TM03%tsBx)
    	call SetupTS("By",wID,TM03%tsBy)
    	call SetupTS("Bz",wID,TM03%tsBz)

    	call SetupTS("MJD",wID,TM03%tsMJD)

    	TM03%isInit = .true.
    	call UpdateTM03(t0)

    	contains
    		subroutine SetupTS(vID,fID,tsQ)
    			character(len=*)  , intent(in)    :: vID,fID
    			type(TimeSeries_T), intent(inout) :: tsQ

    			tsQ%wID = trim(fID)
    			call tsQ%initTS(trim(vID),doLoudO=.false.)
    		end subroutine SetupTS
    end subroutine InitTM03

    !Just return current MJD of TM model
    function TM03_MJD()
    	real(rp) :: TM03_MJD

    	TM03_MJD = TM03%MJD0

    end function TM03_MJD

    !Update time for TM03
    subroutine UpdateTM03(t0)
    	real(rp), intent(in) :: t0
    	real(rp) :: D,Vx,Vy,Vz,Bx,By,Bz

    	if (.not. TM03%isInit) then
    		write(*,*) "Attempting to update uninitialized TM03 model, bailing ..."
    		stop
    	endif
    	D  = TM03%tsD %evalAt(t0)        !#/cc
    	Vx = TM03%tsVx%evalAt(t0)*1.0e-3 !km/s
    	Vy = TM03%tsVy%evalAt(t0)*1.0e-3 !km/s
    	Vz = TM03%tsVz%evalAt(t0)*1.0e-3 !km/s
    	Bx = TM03%tsBx%evalAt(t0)        !nT
    	By = TM03%tsBy%evalAt(t0)        !nT
    	Bz = TM03%tsBz%evalAt(t0)        !nT

    	TM03%MJD0 = TM03%tsMJD%evalAt(t0)
        !Start storing, want GSM
        TM03%dSW  = D
        TM03%vSW  = norm2([Vx,Vy,Vz])
        TM03%PdynSW = PV2PDyn(TM03%dSW,TM03%vSW)

        !Do B field transform
        call mjdRecalc(TM03%MJD0)
        call SM2GSW(Bx,By,Bz,TM03%BxSW,TM03%BySW,TM03%BzSW)

    end subroutine UpdateTM03

	!Vectors should properly be in GSM
	!sScl is safety factor, >1 means less likely to be called out
	function inShueMP(xyz,sSclO)
		real(rp), intent(in) :: xyz(NDIM)
		logical :: inShueMP
		real(rp), intent(in), optional :: sSclO
		real(rp) :: R,R0,RM,sScl
		
		if (present(sSclO)) then
			sScl = sSclO
		else
			sScl = 1.0
		endif
		R = norm2(xyz)
		call ShueHelper(xyz,R0,RM)
		if (R<RM*sScl) then
			inShueMP = .true.
		else
			inShueMP = .false.
		endif

	end function inShueMP

	function inShueMP_SM(xyzSM,sSclO)
		real(rp), intent(in) :: xyzSM(NDIM)
		logical :: inShueMP_SM
		real(rp), intent(in), optional :: sSclO

		real(rp), dimension(NDIM) :: xyzGSM
		call SM2GSW(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGSM(XDIR),xyzGSM(YDIR),xyzGSM(ZDIR))
		if (present(sSclO)) then
			inShueMP_SM = inShueMP(xyzGSM,sSclO)
		else
			inShueMP_SM = inShueMP(xyzGSM)
		endif
		
	end function inSHueMP_SM

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

	subroutine ShueHelper(xyz,R0,RM)
		real(rp), intent(in ) :: xyz (NDIM)
		real(rp), intent(out) :: R0,RM
		real(rp) :: R,theta,CosT,alpha

		R = norm2(xyz)
		theta = acos(xyz(XDIR)/R)
		CosT = cos(theta)
		
		R0 = (10.22+1.29*tanh(0.184*(TM03%BzSW+8.14)))*(TM03%PdynSW)**(-1.0/6.6)
		alpha = (0.58-0.007*TM03%BzSW)*(1.0+0.024*log(TM03%PdynSW))
		RM = R0*(2.0/max(1.0+CosT,TINY))**(alpha)

	end subroutine ShueHelper

	!Return density [#/cc] and pressure [nPa] from TM03
	!Assuming XYZ are in GSM
	subroutine EvalTM03(xyzIN,D,P,isIn)
		real(rp), intent(in)    :: xyzIN(NDIM)
		real(rp), intent(out)   :: D,P
		logical , intent(out)   :: isIn

		real(rp) :: R,clockang,Fst,Pst,rho,rho10,Nst,BNst,BSst,phi
		real(rp) :: xyz(NDIM)

		D = 0.0
		P = 0.0
		isIn = .false.

		R = norm2(xyzIN)
		rho = norm2([xyzIN(XDIR),xyzIN(YDIR)])

		if (xyzIN(XDIR)>0) then
			!Map back to Shue MP at dusk
			call ShueMP2Dusk(xyzIN,xyz)
		else
			xyz = xyzIN
		endif

		isIn = inShueMP(xyz) .and. (rho>10) .and. (xyz(XDIR)<=0) &
			   .and. (xyz(XDIR)>-50)

		if (.not. isIn) return

		!----
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

	end subroutine EvalTM03

	!Evaluate for SM coordinates incoming
	subroutine EvalTM03_SM(xyzSM,D,P,isIn)
		real(rp), intent(in)    :: xyzSM(NDIM)
		real(rp), intent(out)   :: D,P
		logical , intent(out)   :: isIn

		real(rp), dimension(NDIM) :: xyzGSM

		call SM2GSW(xyzSM(XDIR),xyzSM(YDIR),xyzSM(ZDIR),xyzGSM(XDIR),xyzGSM(YDIR),xyzGSM(ZDIR))

		call EvalTM03(xyzGSM,D,P,isIn)

	end subroutine EvalTM03_SM
	
 end module imaghelper
