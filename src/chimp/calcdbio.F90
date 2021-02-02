module calcdbio

	use chmpdefs
	use chmpunits
	use chmpfields
	use ebtypes
	use ioH5
	use xml_input
	use files
	use clocks
	implicit none

    !Remix holders
    type rmState_T
        real(rp) :: time !CHIMP units
        real(rp), dimension(:,:,:), allocatable :: XY
        integer :: Np,Nth !Remix cap sizes
        integer :: i1=-1,i2=-1 !Bracketing step numbers
        !Data arrays are size Np,Nth (lon,lat)
        real(rp), dimension(:,:), allocatable :: nFac,nSigP,nSigH,nPot
        real(rp), dimension(:,:), allocatable :: sFac,sSigP,sSigH,sPot
    
    end type rmState_T

    character(len=strLen), private :: dbOutF
    integer, parameter, private :: MAXDBVS = 20
    !Parameters for output DB grid, 3D thin shell (lat/lon/height)
    integer, private :: NLat,NLon,Nz
    real(rp), private :: dz = 60.0 !Default height spacing [km]
    real(rp), dimension(:,:,:,:), allocatable, private :: xyzI,xyzC !Corner/Center points
    integer, private :: i0 !Shell to start at
    real(rp), private :: rMax = 25.0 !Radius of magnetospheric ball to integrate over [Re]

    contains

    subroutine initDBio(Model,ebState,inpXML,NumP)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(XML_Input_T), intent(in) :: inpXML
        integer, intent(inout) :: NumP

        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        integer :: i,j,k
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

        !Calculate corner/center
        allocate(xyzI(NLat+1,NLon+1,Nz+1,NDIM))
        allocate(xyzC(NLat  ,NLon  ,Nz  ,NDIM))


        !Do corners
        do k=1,Nz+1
        	z = -0.5*dz + (k-1)*dz !km above ground
        	z = (1.0e+5)*z !cm above ground
        	R = 1.0 + z/Re_cgs !Assuming Earth here

        	do j=1,NLon+1
        		phi = (j-1)*360.0/NLon 
        		do i=1,NLat+1
        			lat = -90.0 + (i-1)*180.0/NLat
        			xyzI(i,j,k,XDIR) = R*cos(lat*PI/180.0)*cos(phi*PI/180.0)
        			xyzI(i,j,k,YDIR) = R*cos(lat*PI/180.0)*sin(phi*PI/180.0)
        			xyzI(i,j,k,ZDIR) = R*sin(lat*PI/180.0)

        		enddo !i
        	enddo !j
        enddo

        !Do centers
        do k=1,Nz
        	do j=1,NLon
        		do i=1,NLat
        			xyzC(i,j,k,:) = 0.125*(  xyzI(i,j,k,:)     + xyzI(i+1,j,k,:) &
                                           + xyzI(i,j+1,k,:)   + xyzI(i,j,k+1,:) &
                                           + xyzI(i+1,j+1,k,:) + xyzI(i+1,j,k+1,:) &
                                           + xyzI(i,j+1,k+1,:) + xyzI(i+1,j+1,k+1,:) )
        		enddo
        	enddo
        enddo

        end associate

        !Write grid
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",xyzI(:,:,:,XDIR))
        call AddOutVar(IOVars,"Y",xyzI(:,:,:,YDIR))
        call AddOutVar(IOVars,"Z",xyzI(:,:,:,ZDIR))

        call AddOutVar(IOVars,"Xcc",xyzC(:,:,:,XDIR))
        call AddOutVar(IOVars,"Ycc",xyzC(:,:,:,YDIR))
        call AddOutVar(IOVars,"Zcc",xyzC(:,:,:,ZDIR))

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

    end subroutine initRM

    !Update remix data
    subroutine updateRemix(Model,ebState,t,rmState)
        type(chmpModel_T), intent(in)    :: Model
        type(ebState_T)  , intent(in)   :: ebState
        real(rp), intent(in) :: t
        type(rmState_T), intent(inout)   :: rmState

        character(len=strLen) :: rmF !Remix file
        character(len=strLen) :: gStr
        integer :: i1,i2
        type(IOVAR_T), dimension(MAXDBVS) :: IOVars

        !TODO: Remove redundant code here
        write(rmF,'(2a)') trim(adjustl(Model%RunID)),'.mix.h5'

        call findSlc(ebState%ebTab,t,i1,i2)

        !Right now just lazily reading one slice and storing
        !TODO: Fix this to properly interpolate
        gStr = trim(ebState%ebTab%gStrs(i1))

        call ClearIO(IOVars)
        !Northern hemisphere
        call AddInVar(IOVars,"Field-aligned current NORTH")
        call AddInVar(IOVars, "Pedersen conductance NORTH")
        call AddInVar(IOVars,     "Hall conductance NORTH")
        call AddInVar(IOVars,            "Potential NORTH")
        call AddInVar(IOVars,"Field-aligned current SOUTH")
        call AddInVar(IOVars, "Pedersen conductance SOUTH")
        call AddInVar(IOVars,     "Hall conductance SOUTH")
        call AddInVar(IOVars,            "Potential SOUTH")

        call ReadVars(IOVars,.true.,rmF,gStr)

        write(*,*) 'Done reading ...'
        !Pull data into arrays
        call IOArray2DFill(IOVars,"Field-aligned current NORTH",rmState%nFac )
        call IOArray2DFill(IOVars, "Pedersen conductance NORTH",rmState%nSigP)
        call IOArray2DFill(IOVars,     "Hall conductance NORTH",rmState%nSigH)
        call IOArray2DFill(IOVars,            "Potential NORTH",rmState%nPot )

        call IOArray2DFill(IOVars,"Field-aligned current SOUTH",rmState%sFac )
        call IOArray2DFill(IOVars, "Pedersen conductance SOUTH",rmState%sSigP)
        call IOArray2DFill(IOVars,     "Hall conductance SOUTH",rmState%sSigH)
        call IOArray2DFill(IOVars,            "Potential SOUTH",rmState%sPot )

 
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
        allocate(dbMAG_xyz(NLat,NLon,Nz,NDIM)) !Magnetospheric delta-B
        allocate(dbMAG_rtp(NLat,NLon,Nz,NDIM)) !Magnetospheric delta-B in spherical vectors

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
        do k=1,Nz
            do j=1,NLon
                do i=1,NLat
                    rad = norm2(xyzC(i,j,k,:))
                    theta = acos(xyzC(i,j,k,ZDIR)/rad)
                    phi = atan2(xyzC(i,j,k,YDIR),xyzC(i,j,k,XDIR))
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
    
    !Get ground delta-B at grid points
    subroutine CalcMagDB(Model,ebState,dbMAG)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        real(rp), dimension(:,:,:,:), intent(inout) :: dbMAG

        real(rp), dimension(:,:,:,:), allocatable :: Jxyz
        real(rp) :: w1,w2,dt
        real(rp) :: dV,B0,r3
        real(rp), dimension(NDIM) :: x0,xCC,ddB
        integer :: iG,jG,kG,iM,jM,kM
        B0 = bScale()

        !Start by getting time weight
        if (ebState%doStatic) then
            w1 = 1.0
            w2 = 0.0
        else
            dt = ebState%eb2%time-ebState%eb1%time
            w1 = (ebState%eb2%time-Model%t)/dt
            w2 = (Model%t-ebState%eb1%time)/dt
        endif

    	allocate(Jxyz(NLat,NLon,Nz,NDIM))
    	Jxyz = w1*ebState%eb1%Jxyz + w2*ebState%eb2%Jxyz

    	dbMAG = 0.0
    	associate( ebGr=>ebState%ebGr )
    	!Big-ass loop, loop over grid cells we want dB at and then loop over MHD grid contribution
    	!iG,jG,kG = ground grid
    	!iM,jM,kM = MHD grid

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(iG,jG,kG,iM,jM,kM) &
        !$OMP private(dV,r3,x0,xCC,ddB)
    	do kG=1,Nz
    		do jG=1,NLon
    			do iG=1,NLat
    				x0 = xyzC(iG,jG,kG,:) !Cell center of ground grid

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
        do kG=1,Nz
            do jG=1,NLon
                do iG=1,NLat
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