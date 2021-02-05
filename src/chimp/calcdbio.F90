module calcdbio

	use chmpdefs
	use chmpunits
	use chmpfields
	use ebtypes
	use ioH5
	use xml_input
	use files
	use clocks
    use calcdbutils

	implicit none

    character(len=strLen), private :: dbOutF
    integer, parameter, private :: MAXDBVS = 20
    integer, private :: i0 !Shell to start at
    real(rp), private :: rMax = 25.0 !Radius of magnetospheric ball to integrate over [Re]

    contains

    subroutine initDBio(Model,ebState,inpXML,NumP)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(XML_Input_T), intent(in) :: inpXML
        integer, intent(inout) :: NumP

        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        integer :: i,j,k,NLat,NLon,Nz
        real(rp) :: z,R,lat,phi
        write(dbOutF,'(2a)') trim(adjustl(Model%RunID)),'.deltab.h5'

        associate( ebGr=>ebState%ebGr )
        call CheckAndKill(dbOutF)

        call inpXML%Set_Val(NLat,'Grid/NLat',45) !Number of latitudinal cells
        call inpXML%Set_Val(NLon,'Grid/NLon',90) !Number of longitudinal cells
        call inpXML%Set_Val(Nz  ,'Grid/Nz'  , 2) !Number of longitudinal cells
        call inpXML%Set_Val(i0,'CalcDB/i0',1)
        call inpXML%Set_Val(rMax,'CalcDB/rMax',rMax)

        NumP = NLat*NLon*Nz

        !Store sizes
        gGr%NLat = NLat
        gGr%NLon = NLon
        gGr%Nz   = Nz

        !Calculate corner/center
        allocate(gGr%xyzI(NLat+1,NLon+1,Nz+1,NDIM))
        allocate(gGr%xyzC(NLat  ,NLon  ,Nz  ,NDIM))

        !Do corners
        do k=1,Nz+1
        	z = -0.5*dzGG + (k-1)*dzGG !km above ground
        	z = (1.0e+5)*z !cm above ground
        	R = 1.0 + z/Re_cgs !Assuming Earth here

        	do j=1,NLon+1
        		phi = (j-1)*360.0/NLon 
        		do i=1,NLat+1
        			lat = -90.0 + (i-1)*180.0/NLat
        			gGr%xyzI(i,j,k,XDIR) = R*cos(lat*PI/180.0)*cos(phi*PI/180.0)
        			gGr%xyzI(i,j,k,YDIR) = R*cos(lat*PI/180.0)*sin(phi*PI/180.0)
        			gGr%xyzI(i,j,k,ZDIR) = R*sin(lat*PI/180.0)

        		enddo !i
        	enddo !j
        enddo

        !Do centers
        do k=1,Nz
        	do j=1,NLon
        		do i=1,NLat
        			gGr%xyzC(i,j,k,:) = 0.125*(  gGr%xyzI(i,j,k,:)     + gGr%xyzI(i+1,j,k,:) &
                                                 + gGr%xyzI(i,j+1,k,:)   + gGr%xyzI(i,j,k+1,:) &
                                                 + gGr%xyzI(i+1,j+1,k,:) + gGr%xyzI(i+1,j,k+1,:) &
                                                 + gGr%xyzI(i,j+1,k+1,:) + gGr%xyzI(i+1,j+1,k+1,:) )
        		enddo
        	enddo
        enddo

        end associate

        !Write grid
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",gGr%xyzI(:,:,:,XDIR))
        call AddOutVar(IOVars,"Y",gGr%xyzI(:,:,:,YDIR))
        call AddOutVar(IOVars,"Z",gGr%xyzI(:,:,:,ZDIR))

        call AddOutVar(IOVars,"Xcc",gGr%xyzC(:,:,:,XDIR))
        call AddOutVar(IOVars,"Ycc",gGr%xyzC(:,:,:,YDIR))
        call AddOutVar(IOVars,"Zcc",gGr%xyzC(:,:,:,ZDIR))

        call WriteVars(IOVars,.true.,dbOutF)
        call ClearIO(IOVars)

    end subroutine initDBio

    subroutine initRM(Model,ebState,rmState)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(rmState_T), intent(inout) :: rmState
        
        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        integer :: Ni,Nj
        character(len=strLen) :: rmF !Remix file
        write(rmF,'(2a)') trim(adjustl(Model%RunID)),'.mix.h5'
        write(*,*) 'Initializing w/ ', trim(rmF)

        call ClearIO(IOVars)
        call AddInVar(IOVars,"X")
        call AddInVar(IOVars,"Y")
        call ReadVars(IOVars,.true.,rmF)

        Ni = IOVars(1)%dims(1) - 1
        Nj = IOVars(2)%dims(2) - 1

        write(*,*) 'Shape = ', IOVars(1)%dims(:)

        !Initialize 4 hemispheres (N/S i1/i2)
        call initHemi(rmState%rmN1,Ni,Nj)
        call initHemi(rmState%rmN2,Ni,Nj)
        call initHemi(rmState%rmS1,Ni,Nj)
        call initHemi(rmState%rmS2,Ni,Nj)

        !Now do main rmState
        rmState%Np  = Ni !Phi cells
        rmState%Nth = Nj !Theta cells

        allocate(rmState%nFac (Ni,Nj))
        allocate(rmState%nPot (Ni,Nj))
        allocate(rmState%nSigP(Ni,Nj))
        allocate(rmState%nSigH(Ni,Nj))
        allocate(rmState%sFac (Ni,Nj))
        allocate(rmState%sPot (Ni,Nj))
        allocate(rmState%sSigP(Ni,Nj))
        allocate(rmState%sSigH(Ni,Nj))

        call hemi2rm(rmState,0.0_rp,0.0_rp) !Zero out main state arrays

        contains

            subroutine initHemi(rmHemi,Ni,Nj)
                type(rmHemi_T), intent(inout) :: rmHemi
                integer, intent(in) :: Ni,Nj

                rmHemi%nStp = -1 !Not yet set
                rmHemi%time = 0.0
                rmHemi%Np  = Ni !Phi cells
                rmHemi%Nth = Nj !Theta cells

                allocate(rmHemi%xFac (Ni,Nj))
                allocate(rmHemi%xPot (Ni,Nj))
                allocate(rmHemi%xSigP(Ni,Nj))
                allocate(rmHemi%xSigH(Ni,Nj))

                rmHemi%xFac  = 0.0
                rmHemi%xPot  = 0.0
                rmHemi%xSigP = 0.0
                rmHemi%xSigH = 0.0

            end subroutine initHemi
    end subroutine initRM

    !Update remix data
    subroutine updateRemix(Model,ebState,t,rmState)
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)   :: ebState
        real(rp), intent(in) :: t
        type(rmState_T), intent(inout)   :: rmState

        character(len=strLen) :: rmF !Remix file    
        integer :: i1,i2
        real(rp) :: w1,w2

        !TODO: Remove redundant code here
        write(rmF,'(2a)') trim(adjustl(Model%RunID)),'.mix.h5'
        call findSlc(ebState%ebTab,t,i1,i2)

        !Read 4 hemispheres
        call readHemi(rmState%rmN1,rmF,ebState%ebTab,i1,NORTH)
        call readHemi(rmState%rmS1,rmF,ebState%ebTab,i1,SOUTH)

        call readHemi(rmState%rmN2,rmF,ebState%ebTab,i2,NORTH)
        call readHemi(rmState%rmS2,rmF,ebState%ebTab,i2,SOUTH)

        !Now fill in remix main state for this t
        call GetTWgts(Model,ebState,t,w1,w2)
        call hemi2rm(rmState,w1,w2)

        !TODO: Fill in fac and ionospheric grid data
        call facGridUpdate(Model,ebState,rmState)

        contains

        !Read hemisphere data
        !NOTE: nStp here refers to array index in ebTab, not necessarily Step#X
        subroutine readHemi(rmHemi,rmF,ebTab,nStp,nsID)
            type(rmHemi_T), intent(inout) :: rmHemi
            character(len=strLen), intent(in) :: rmF
            type(ebTab_T), intent(in) :: ebTab
            integer, intent(in) :: nStp,nsID

            character(len=strLen) :: hID,gStr
            type(IOVAR_T), dimension(MAXDBVS) :: IOVars

            !Check to see if we need to read
            if (rmHemi%nStp == nStp) then
                !We've already read this data
                return
            endif
            !Otherwise get the data
            !Which hemisphere?
            if (nsID == NORTH) then
                hID = "NORTH"
            else
                hID = "SOUTH"
            endif
            gStr = trim(ebState%ebTab%gStrs(nStp))

            rmHemi%time = ebTab%times(nStp)
            rmHemi%nStp = nStp

            call ClearIO(IOVars)
            call AddInVar(IOVars,"Field-aligned current" // hID)
            call AddInVar(IOVars, "Pedersen conductance" // hID)
            call AddInVar(IOVars,     "Hall conductance" // hID)
            call AddInVar(IOVars,            "Potential" // hID)
            call ReadVars(IOVars,.true.,rmF,gStr)

            call IOArray2DFill(IOVars,"Field-aligned current" // hID,rmHemi%xFac )
            call IOArray2DFill(IOVars, "Pedersen conductance" // hID,rmHemi%xSigP)
            call IOArray2DFill(IOVars,     "Hall conductance" // hID,rmHemi%xSigH)
            call IOArray2DFill(IOVars,            "Potential" // hID,rmHemi%xPot )

        end subroutine readHemi

    end subroutine updateRemix

    !Calculate and output delta-B data
    subroutine writeDB(Model,ebState,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        real(rp), dimension(:,:,:,:), allocatable :: dbMAG_xyz,dbMAG_rtp
        real(rp) :: mjd

        mjd = MJDAt(ebState%ebTab,Model%t)
        allocate(dbMAG_xyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Magnetospheric delta-B
        allocate(dbMAG_rtp(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Magnetospheric delta-B in spherical vectors

        call Tic("CalcMagDB")
        call CalcMagDB(Model,ebState,dbMAG_xyz)
        call Toc("CalcMagDB")
        
        
        call xyz2rtp(dbMAG_xyz,dbMAG_rtp)

        call ClearIO(IOVars)
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"MJD",mjd)
        call AddOutVar(IOVars,"dBx_M" ,dbMAG_xyz(:,:,:,XDIR),"nT")
        call AddOutVar(IOVars,"dBy_M" ,dbMAG_xyz(:,:,:,YDIR),"nT")
        call AddOutVar(IOVars,"dBz_M" ,dbMAG_xyz(:,:,:,ZDIR),"nT")
        !Write out spherical vectors (XDIR:ZDIR = RDIR,TDIR,PDIR)
        call AddOutVar(IOVars,"dBr_M" ,dbMAG_rtp(:,:,:,XDIR),"nT")
        call AddOutVar(IOVars,"dBt_M" ,dbMAG_rtp(:,:,:,YDIR),"nT")
        call AddOutVar(IOVars,"dBp_M" ,dbMAG_rtp(:,:,:,ZDIR),"nT")

        call WriteVars(IOVars,.true.,dbOutF,gStr)
        call ClearIO(IOVars)
        
        
    end subroutine writeDB

    !Convert XYZ to RTP vectors
    subroutine xyz2rtp(dbXYZ,dbRTP)
        real(rp), dimension(:,:,:,:), intent(in ) :: dbXYZ
        real(rp), dimension(:,:,:,:), intent(out) :: dbRTP

        integer :: i,j,k
        real(rp) :: rad,theta,phi,dBx,dBy,dBz
        dbRTP = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,rad,theta,phi,dBx,dBy,dBz)
        do k=1,gGr%Nz
            do j=1,gGr%NLon
                do i=1,gGr%NLat
                    rad  = norm2(gGr%xyzC(i,j,k,:))
                    theta = acos(gGr%xyzC(i,j,k,ZDIR)/rad)
                    phi  = atan2(gGr%xyzC(i,j,k,YDIR),gGr%xyzC(i,j,k,XDIR))
                    dBx = dbXYZ(i,j,k,XDIR)
                    dBy = dbXYZ(i,j,k,YDIR)
                    dBz = dbXYZ(i,j,k,ZDIR)

                    !Radial,Theta,Phi system
                    dbRTP(i,j,k,1) =  dbX*sin(theta)*cos(phi) + dBy*sin(theta)*sin(phi) + dBz*cos(theta)
                    dbRTP(i,j,k,2) =  dbX*cos(theta)*cos(phi) + dBy*cos(theta)*sin(phi) - dBz*sin(theta)
                    dbRTP(i,j,k,3) = -dbX           *sin(phi) + dBy           *cos(phi) 

                enddo
            enddo
        enddo

    end subroutine xyz2rtp
    
    !Get time weights for time t (assuming proper bracketing)
    subroutine GetTWgts(Model,ebState,t,w1,w2)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), intent(in)  :: t
        real(rp), intent(out) :: w1,w2

        real(rp) :: dt
        !Start by getting time weight
        if (ebState%doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = ebState%eb2%time-ebState%eb1%time
            w1 = (ebState%eb2%time-t)/dt
            w2 = (t-ebState%eb1%time)/dt
        endif

    end subroutine GetTWgts

    !Get ground delta-B at grid points
    subroutine CalcMagDB(Model,ebState,dbMAG)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), dimension(:,:,:,:), intent(inout) :: dbMAG

        real(rp), dimension(:,:,:,:), allocatable :: Jxyz
        real(rp) :: w1,w2,dV,B0,r3
        real(rp), dimension(NDIM) :: x0,xCC,ddB
        integer :: iG,jG,kG,iM,jM,kM
        
        B0 = bScale()
        call GetTWgts(Model,ebState,Model%t,w1,w2)

    	allocate(Jxyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM))
        !$OMP PARALLEL WORKSHARE
    	Jxyz = w1*ebState%eb1%Jxyz + w2*ebState%eb2%Jxyz
        !$OMP END PARALLEL WORKSHARE
        
    	dbMAG = 0.0
    	associate( ebGr=>ebState%ebGr )
    	!Big-ass loop, loop over grid cells we want dB at and then loop over MHD grid contribution
    	!iG,jG,kG = ground grid
    	!iM,jM,kM = MHD grid

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(iG,jG,kG,iM,jM,kM) &
        !$OMP private(dV,r3,x0,xCC,ddB)
    	do kG=1,gGr%Nz
    		do jG=1,gGr%NLon
    			do iG=1,gGr%NLat
    				x0 = gGr%xyzC(iG,jG,kG,:) !Cell center of ground grid

    				do kM=ebGr%ks,ebGr%ke
    					do jM=ebGr%js,ebGr%je
    						do iM=ebGr%is+i0-1,ebGr%ie
    							xCC = ebGr%xyzcc(iM,jM,kM,:) !MHD grid cell center
    							if ( norm2(xCC) > rMax ) cycle
                                
                                r3 = norm2(x0-xCC)**3.0
    							dV = ebGr%dV(iM,jM,kM)
    							!Get differential contribution
                                !Avoid array temporary
    							!ddB = -cross(Jxyz(iM,jM,kM,:),xCC)/r3
                                ddB(XDIR) = -( Jxyz(iM,jM,kM,YDIR)*xCC(ZDIR) - Jxyz(iM,jM,kM,ZDIR)*xCC(YDIR) )/r3
                                ddB(YDIR) = -( Jxyz(iM,jM,kM,ZDIR)*xCC(XDIR) - Jxyz(iM,jM,kM,XDIR)*xCC(ZDIR) )/r3
                                ddB(ZDIR) = -( Jxyz(iM,jM,kM,XDIR)*xCC(YDIR) - Jxyz(iM,jM,kM,YDIR)*xCC(XDIR) )/r3

                                !Pulling out overall scaling
    							!dbMAG(iG,jG,kG,:) = dbMAG(iG,jG,kG,:) + B0*dV*ddB/(4.0*PI)
                                dbMAG(iG,jG,kG,:) = dbMAG(iG,jG,kG,:) + dV*ddB

    						enddo !iM
    					enddo !jM
    				enddo !kM

    			enddo !iG
    		enddo !jG
    	enddo !kG

        !Moving overall scaling factor to secondary loop to only apply once
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(iG,jG,kG)
        do kG=1,gGr%Nz
            do jG=1,gGr%NLon
                do iG=1,gGr%NLat
                    dbMAG(iG,jG,kG,:) = B0*dbMAG(iG,jG,kG,:)/(4.0*PI)
                enddo !iG
            enddo !jG
        enddo !kG

    	end associate

    end subroutine CalcMagDB

    !Lazy function to return scaling, should be replaced
    function bScale() result(B0)
    	real(rp) :: B0
    	real(rp) :: mu0,Mp,x0,u0,t0,d0,p0
    	mu0 = 4*PI*1e-7
    	Mp = 1.67e-27 ![kg]
		x0 = 1*6.38e6 ![m]   - RE
		u0 = 100e3    ![m/s] - 100 km/s
		t0 = x0/u0 ![s]   -
		d0 = Mp*1e6 ! [kg/m^3] - 1 particle/cc
		p0 = d0*u0*u0 ![N/m^2]
		B0 = sqrt(mu0*d0*u0*u0)*1e9 ! [nT]
    end function bScale

end module calcdbio