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
    use calcdbcore
    use ebinterp
    
	implicit none

    character(len=strLen), private :: dbOutF
    integer, parameter, private :: MAXDBVS = 20
    logical, private :: doParInT = .false. !// in time
    integer, private :: NumB = 0

    contains

    subroutine initDBio(Model,ebState,gGr,inpXML,NumP)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)    :: ebState
        type(grGrid_T), intent(inout) :: gGr
        type(XML_Input_T), intent(in) :: inpXML
        integer, intent(inout) :: NumP

        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        character(len=strLen) :: cID

        integer :: i,j,k,NLat,NLon,Nz,dOut
        real(rp) :: z,R,lat,phi
        real(rp) :: dtB,T0

        associate( ebGr=>ebState%ebGr )
        !Equate dtout/dt since the difference doesn't matter here
        Model%dtOut = Model%dt

    !Read info from XML
        call inpXML%Set_Val(NLat,'Grid/NLat',45) !Number of latitudinal cells
        call inpXML%Set_Val(NLon,'Grid/NLon',90) !Number of longitudinal cells
        call inpXML%Set_Val(Nz  ,'Grid/Nz'  , 2) !Number of longitudinal cells
        call inpXML%Set_Val(gGr%rMax,'CalcDB/rMax',gGr%rMax)
        call inpXML%Set_Val(gGr%doGEO,'Grid/doGEO',gGr%doGEO) !Whether to do GEO on ground
        if (gGr%doGEO) then
            cID = "GEO"
        else
            cID = "SM"
        endif

    !Possible // in time
        call inpXML%Set_Val(NumB,'parintime/NumB',NumB)
        if (NumB > 1) then
            doParInT = .true.
            if ( (Model%Nblk>NumB) .or. (Model%Nblk<1) ) then
                write(*,*) "This block outside of acceptable bounds"
                write(*,*) "Block = ",Model%Nblk
                write(*,*) "Bounds = ",1,NumB
                write(*,*) "Bailing ..."
                stop
            endif
            !Reset time bounds
            T0 = Model%T0
            dtB = (Model%tFin-Model%T0)/NumB
            write(*,*) '------'
            write(*,*) 'Resetting T0/TFin = ',Model%T0*oTScl,Model%tFin*oTScl
            write(*,*) 'Using block ', Model%Nblk
            Model%T0 = (Model%Nblk-1)*dtB + T0
            Model%tFin = Model%T0 + dtB
            write(*,*) 'To        T0/TFin = ',Model%T0*oTScl,Model%tFin*oTScl
            if (Model%Nblk < NumB) then
                !Cut off a bit from TFin to avoid overlap w/ start of next
                Model%tFin = Model%tFin-0.01*dtB
            endif
            !Get step# offset
            !NOTE: Assuming here nice divisibility
            dOut = nint(dtB/Model%dtOut)
            Model%nOut = (Model%Nblk-1)*dOut + 0
            write(*,*) 'Offsetting Step# by ', Model%nOut

            write(dbOutF,'(a,a,I0.4,a)') trim(adjustl(Model%RunID)),'.',Model%Nblk,'.deltab.h5'
            write(*,*) '------'

        else
            doParInT = .false.
            NumB = 0
            write(dbOutF,'(2a)') trim(adjustl(Model%RunID)),'.deltab.h5'
        endif
        write(*,*) "Writing output to ", trim(dbOutF)
    
    !Setup output file
        call CheckAndKill(dbOutF)

    !Create ground grid
        NumP = NLat*NLon*Nz

        !Store sizes
        gGr%NLat = NLat
        gGr%NLon = NLon
        gGr%Nz   = Nz

        !Calculate corner/center
        allocate(gGr%GxyzI (NLat+1,NLon+1,Nz+1,NDIM))
        allocate(gGr%GxyzC (NLat  ,NLon  ,Nz  ,NDIM))
        allocate(gGr%SMxyzC(NLat  ,NLon  ,Nz  ,NDIM))

        !Do corners
        do k=1,Nz+1
        	z = -0.5*dzGG + (k-1)*dzGG !km above ground
        	z = (1.0e+5)*z !cm above ground
        	R = 1.0 + z/Re_cgs !Assuming Earth here

        	do j=1,NLon+1
        		phi = (j-1)*360.0/NLon 
        		do i=1,NLat+1
        			lat = -90.0 + (i-1)*180.0/NLat
        			gGr%GxyzI(i,j,k,XDIR) = R*cos(lat*PI/180.0)*cos(phi*PI/180.0)
        			gGr%GxyzI(i,j,k,YDIR) = R*cos(lat*PI/180.0)*sin(phi*PI/180.0)
        			gGr%GxyzI(i,j,k,ZDIR) = R*sin(lat*PI/180.0)

        		enddo !i
        	enddo !j
        enddo

        !Do centers
        do k=1,Nz
        	do j=1,NLon
        		do i=1,NLat
        			gGr%GxyzC(i,j,k,:) = 0.125*( gGr%GxyzI(i  ,j  ,k  ,:) + gGr%GxyzI(i+1,j  ,k  ,:) &
                                               + gGr%GxyzI(i  ,j+1,k  ,:) + gGr%GxyzI(i  ,j  ,k+1,:) &
                                               + gGr%GxyzI(i+1,j+1,k  ,:) + gGr%GxyzI(i+1,j  ,k+1,:) &
                                               + gGr%GxyzI(i  ,j+1,k+1,:) + gGr%GxyzI(i+1,j+1,k+1,:) )
                    gGr%SMxyzC(i,j,k,:) = gGr%GxyzC(i,j,k,:) !Just set SM=G for now
        		enddo
        	enddo
        enddo

        end associate

        !Allocate db holders (XYZ)
        allocate(gGr%dbMAG_xyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Magnetospheric delta-B
        allocate(gGr%dbION_xyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Ionospheric delta-B
        allocate(gGr%dbFAC_xyz(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !FAC delta-B
        !Spherical holders
        allocate(gGr%dbMAG_rtp(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Magnetospheric delta-B
        allocate(gGr%dbION_rtp(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !Ionospheric delta-B
        allocate(gGr%dbFAC_rtp(gGr%NLat,gGr%NLon,gGr%Nz,NDIM)) !FAC delta-B

        gGr%dbMAG_xyz = 0.0
        gGr%dbION_xyz = 0.0
        gGr%dbFAC_xyz = 0.0
        gGr%dbMAG_rtp = 0.0
        gGr%dbION_rtp = 0.0
        gGr%dbFAC_rtp = 0.0

        !Write grid
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",gGr%GxyzI(:,:,:,XDIR))
        call AddOutVar(IOVars,"Y",gGr%GxyzI(:,:,:,YDIR))
        call AddOutVar(IOVars,"Z",gGr%GxyzI(:,:,:,ZDIR))

        call AddOutVar(IOVars,"Xcc",gGr%GxyzC(:,:,:,XDIR))
        call AddOutVar(IOVars,"Ycc",gGr%GxyzC(:,:,:,YDIR))
        call AddOutVar(IOVars,"Zcc",gGr%GxyzC(:,:,:,ZDIR))

        call AddOutVar(IOVars,"CoordinatesID",cID)

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
        write(rmF,'(2a)') trim(adjustl(ebState%ebTab%bStr)),'.mix.h5'

        write(*,*) 'Initializing w/ ', trim(rmF)

        call ClearIO(IOVars)
        call AddInVar(IOVars,"X")
        call AddInVar(IOVars,"Y")
        call ReadVars(IOVars,.true.,rmF)

        Ni = IOVars(1)%dims(1) - 1
        Nj = IOVars(2)%dims(2) - 1

        !Initialize 4 hemispheres (N/S i1/i2)
        call initHemi(rmState%rmN1,Ni,Nj)
        call initHemi(rmState%rmN2,Ni,Nj)
        call initHemi(rmState%rmS1,Ni,Nj)
        call initHemi(rmState%rmS2,Ni,Nj)

        !Now do main rmState
        rmState%Np  = Ni !Phi cells
        rmState%Nth = Nj !Theta cells

        allocate(rmState%XY(Ni+1,Nj+1,XDIR:YDIR))

        allocate(rmState%nFac (Ni,Nj))
        allocate(rmState%nPot (Ni,Nj))
        allocate(rmState%nSigP(Ni,Nj))
        allocate(rmState%nSigH(Ni,Nj))
        allocate(rmState%sFac (Ni,Nj))
        allocate(rmState%sPot (Ni,Nj))
        allocate(rmState%sSigP(Ni,Nj))
        allocate(rmState%sSigH(Ni,Nj))

        call IOArray2DFill(IOVars,"X",rmState%XY(:,:,XDIR))
        call IOArray2DFill(IOVars,"Y",rmState%XY(:,:,YDIR))


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
        write(rmF,'(2a)') trim(adjustl(ebState%ebTab%bStr)),'.mix.h5'
        
        call findSlc(ebState%ebTab,t,i1,i2)

        !Read 4 hemispheres
        call readHemi(rmState%rmN1,rmF,ebState%ebTab,i1,NORTH)
        call readHemi(rmState%rmS1,rmF,ebState%ebTab,i1,SOUTH)

        call readHemi(rmState%rmN2,rmF,ebState%ebTab,i2,NORTH)
        call readHemi(rmState%rmS2,rmF,ebState%ebTab,i2,SOUTH)

        !Now fill in remix main state for this t
        call GetTWgts(Model,ebState,t,w1,w2)
        call hemi2rm(rmState,w1,w2)

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
            call AddInVar(IOVars,"Field-aligned current " // hID)
            call AddInVar(IOVars, "Pedersen conductance " // hID)
            call AddInVar(IOVars,     "Hall conductance " // hID)
            call AddInVar(IOVars,            "Potential " // hID)
            call ReadVars(IOVars,.true.,rmF,gStr)

            !Pull from arrays
            call IOArray2DFill(IOVars,"Field-aligned current " // hID,rmHemi%xFac )
            call IOArray2DFill(IOVars, "Pedersen conductance " // hID,rmHemi%xSigP)
            call IOArray2DFill(IOVars,     "Hall conductance " // hID,rmHemi%xSigH)
            call IOArray2DFill(IOVars,            "Potential " // hID,rmHemi%xPot )

        end subroutine readHemi

    end subroutine updateRemix

    !Calculate and output delta-B data
    subroutine writeDB(Model,ebState,gGr,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(grGrid_T), intent(inout) :: gGr
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXDBVS) :: IOVars
        real(rp) :: mjd
        real(rp), dimension(:,:,:,:), allocatable :: dbRTP

        write(*,*) 'Writing ', trim(gStr)

        mjd = MJDAt(ebState%ebTab,Model%t)
        
        !Do conversion to spherical
        call xyz2rtp(gGr,gGr%dbMAG_xyz,gGr%dbMAG_rtp)
        call xyz2rtp(gGr,gGr%dbION_xyz,gGr%dbION_rtp)
        call xyz2rtp(gGr,gGr%dbFAC_xyz,gGr%dbFAC_rtp)

        !Get total perturbation
        allocate(dbRTP(gGr%NLat,gGr%NLon,gGr%Nz,NDIM))
        !$OMP PARALLEL WORKSHARE
        dbRTP = gGr%dbMAG_rtp + gGr%dbION_rtp + gGr%dbFAC_rtp
        !$OMP END PARALLEL WORKSHARE

        call ClearIO(IOVars)
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"MJD",mjd)

    !Write out spherical vectors (XDIR:ZDIR = RDIR,TDIR,PDIR)
        if (.not. Model%doSlim) then
            call AddOutVar(IOVars,"dBrM" ,gGr%dbMAG_rtp(:,:,:,XDIR),"nT")
            call AddOutVar(IOVars,"dBtM" ,gGr%dbMAG_rtp(:,:,:,YDIR),"nT")
            call AddOutVar(IOVars,"dBpM" ,gGr%dbMAG_rtp(:,:,:,ZDIR),"nT")

            call AddOutVar(IOVars,"dBrI" ,gGr%dbION_rtp(:,:,:,XDIR),"nT")
            call AddOutVar(IOVars,"dBtI" ,gGr%dbION_rtp(:,:,:,YDIR),"nT")
            call AddOutVar(IOVars,"dBpI" ,gGr%dbION_rtp(:,:,:,ZDIR),"nT")

            call AddOutVar(IOVars,"dBrF" ,gGr%dbFAC_rtp(:,:,:,XDIR),"nT")
            call AddOutVar(IOVars,"dBtF" ,gGr%dbFAC_rtp(:,:,:,YDIR),"nT")
            call AddOutVar(IOVars,"dBpF" ,gGr%dbFAC_rtp(:,:,:,ZDIR),"nT")
        endif
        
        call AddOutVar(IOVars,"dBr" ,dbRTP(:,:,:,XDIR),"nT")
        call AddOutVar(IOVars,"dBt" ,dbRTP(:,:,:,YDIR),"nT")
        call AddOutVar(IOVars,"dBp" ,dbRTP(:,:,:,ZDIR),"nT")

        call WriteVars(IOVars,.true.,dbOutF,gStr)
        call ClearIO(IOVars)
        
    end subroutine writeDB

    !Convert XYZ to RTP vectors
    subroutine xyz2rtp(gGr,dbXYZ,dbRTP)
        type(grGrid_T), intent(in) :: gGr
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
                    rad  = norm2(gGr%GxyzC(i,j,k,:))
                    theta = acos(gGr%GxyzC(i,j,k,ZDIR)/rad)
                    phi  = atan2(gGr%GxyzC(i,j,k,YDIR),gGr%GxyzC(i,j,k,XDIR))
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
    

end module calcdbio