module psdio
    use chmpdefs
    use chmpunits
    use ebtypes
    use psdtypes
    use psdutils
    use particleio
    use lineio
    use files
    use pdfuns

    implicit none
    character(len=strLen) :: ps4OutF,psOutF,pseqOutF,wgtOutF

    integer, parameter :: MAXPSVS = 20
    integer, parameter :: PSQNUM = 1 !Min TP per cell for quality

    logical :: doInitPSIO = .false. !Has IO been initialized
    logical :: doLogKCyl = .true. !Log(K) for K-cylinders
    logical :: doPSLines = .false.
    logical, parameter :: doEBKEQ = .true.

    contains

    subroutine fOutPSD(Model,psGr,psPop)
        type(chmpModel_T), intent(inout) :: Model
        type(PSEq_T), intent(in) :: psGr
        type(psdPop_T), intent(in) :: psPop
        character(len=strLen) :: gStr

        integer :: Nfl

        !write(*,*) 'Writing PSD output at t = ', oTScl*Model%t
        !Do init if needed
        if (.not. doInitPSIO) then
            call InitPSIO(Model,psGr)
            if (doPSLines) call initFLio(Model)
            doInitPSIO = .true.
        endif

        !Write current slice
        Nfl = size(psGr%bLns)
        !write(*,*) 'Writing N lines, N = ',Nfl
        write(gStr,'(A,I0)') "Step#", Model%nOut
        call writePSD(Model,psGr,psPop,gStr)
        if (doPSLines) then
            write(*,*) 'PSLines Not yet implemented!'
            stop
            !call writeLines(Model,psGr%bLns,Nfl)
        endif
        
        !Setup for next output
        Model%tOut = Model%tOut + Model%dtOut
        Model%nOut = Model%nOut + 1

    end subroutine fOutPSD

    subroutine writePSD(Model,psGr,psPop,gStr)
        type(chmpModel_T), intent(inout) :: Model
        type(PSEq_T), intent(in) :: psGr
        type(psdPop_T), intent(in) :: psPop
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXPSVS) :: IOVars
        !3D vars (r,p,K)
        real(rp), dimension(:,:,:), allocatable :: cP,jP,fP,dG
        !2D vars (r,p)
        real(rp), dimension(:,:), allocatable :: isC,dvB,cP2,nQ
        real(rp), dimension(:,:,:), allocatable :: Pxyz

        integer :: ir,ip,ik,ia,idx(NVARPS)
        real(rp) :: ds,da,pMag

        associate(Nr=>psGr%Nr,Np=>psGr%Np,Nk=>psGr%Nk,Na=>psGr%Na)

    !Do 4D variables if desired
        if (Model%doFat) then
            call ClearIO(IOVars)
            call AddOutVar(IOVars,"Ntp",1.0_rp*psPop%nPSD)
            call AddOutVar(IOVars,"fPSD",psPop%fPSD)
            call AddOutVar(IOVars,"dG",psGr%dG)
            call AddOutVar(IOVars,"time",oTScl*Model%t)
            call WriteVars(IOVars,.true.,ps4OutF,gStr)
        endif
    !Do 3D variables
        allocate(cP(Nr,Np,Nk))
        allocate(fP(Nr,Np,Nk))
        allocate(jP(Nr,Np,Nk))
        allocate(dG(Nr,Np,Nk))
        fP = 0.0
        jP = 0.0

        !Sum over alpha to get K-Cyls
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(ia,ir,ip,ik,idx,da,ds,pMag)
        do ik=1,Nk
            do ip=1,Np
                do ir=1,Nr
                    cP(ir,ip,ik) = 1.0*sum(psPop%nPSD(ir,ip,ik,:))
                    dG(ir,ip,ik) = sum(psGr%dG(ir,ip,ik,:))

                    !Integrate PSD over alpha,psi (i.e. intensity)
                    do ia=1,Na
                        idx = [ir,ip,ik,ia]
                        da = dGamma(Model,psGr,idx,PSALPHA)
                        ds = dGamma(Model,psGr,idx,PSPSI)
                        jP(ir,ip,ik) = jP(ir,ip,ik) + ds*da*psPop%fPSD(ir,ip,ik,ia)
                    enddo
                    pMag = psMagP(Model,psGr,idx) !ia of idx doesn't matter
                    !Correct for /steradian
                    jP(ir,ip,ik) = jP(ir,ip,ik)/(4*PI)
                    fP(ir,ip,ik) = jP(ir,ip,ik)/(pMag*pMag)
                enddo
            enddo
        enddo

        

        call ClearIO(IOVars)
        call AddOutVar(IOVars,"Ntp",cP)
        call AddOutVar(IOVars,"fPSD",fP)
        call AddOutVar(IOVars,"jPSD",jP)
        call AddOutVar(IOVars,"dG",dG)
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call WriteVars(IOVars,.true.,psOutF,gStr)

    !Do 2D variables
        allocate(dvB(Nr,Np))
        allocate(isC(Nr,Np))
        allocate(cP2(Nr,Np))
        allocate(nQ (Nr,Np)) !Quality factor

        !$OMP PARALLEL DO default(shared) collapse(2)
        do ip=1,Np
            do ir=1,Nr
                dvB(ir,ip) = maxval(psGr%dvB(ir,ip,:))
                if (psGr%isClosed(ir,ip)) then
                    isC(ir,ip) = 1.0
                else
                    isC(ir,ip) = 0.0
                endif
                cP2(ir,ip) = sum(psPop%nPSD(ir,ip,:,:))
                nQ(ir,ip) = count( psPop%nPSD(ir,ip,:,:) > PSQNUM )/(1.0*Na*Nk)

            enddo
        enddo
        !Get PSD pressure contributions
        allocate(Pxyz(Nr,Np,NDIM))
        Pxyz = 0.0

        call CalcP(Model,psGr,psPop,Pxyz)
        
        !For now, just not deal with flow speed stuff
        Pxyz(:,:,XDIR) = Pxyz(:,:,YDIR)

        !Prep and write values
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"Ntp",cP2)
        call AddOutVar(IOVars,"isC",isC)
        call AddOutVar(IOVars,"dvB",dvB)
        call AddOutVar(IOVars,"nQ" ,nQ )
        
        call AddOutVar(IOVars,"D" ,      psGr%Qrp(:,:,DEN))
        call AddOutVar(IOVars,"Vx",oVScl*psGr%Qrp(:,:,VELX))
        call AddOutVar(IOVars,"Vy",oVScl*psGr%Qrp(:,:,VELY))
        call AddOutVar(IOVars,"Vz",oVScl*psGr%Qrp(:,:,VELZ))
        call AddOutVar(IOVars,"Pg",      psGr%Qrp(:,:,PRESSURE))
        call AddOutVar(IOVars,"Pk",sum(Pxyz,dim=3)/3.0)
        call AddOutVar(IOVars,"Pxy",0.5*(Pxyz(:,:,XDIR)+Pxyz(:,:,YDIR)))
        call AddOutVar(IOVars,"Pz",Pxyz(:,:,ZDIR))
        call AddOutVar(IOVars,"kT",psGr%kTeq)
        call AddOutVar(IOVars,"Vr",oVScl*psGr%Vreq)
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call WriteVars(IOVars,.true.,pseqOutF,gStr)

        !Write some text to screen
        !Find index to demarcate geosynchronous orbit
        ir = maxloc(psGr%rI,dim=1,mask=6.5>=psGr%rI)

        write(*,*) '----------------------------'
        write(*,*) 'Writing PSD output at t = ', oTScl*Model%t
        write(*,*) 'Total Domain'
        write(*,*) '   Normalization = ', sum(psGr%dG*psPop%fPSD)
        write(*,*) '   TPs           = ', sum(psPop%nPSD)
        write(*,*) 'Inside Geo'
        write(*,*) '   Normalization = ', sum(psGr%dG(1:ir,:,:,:)*psPop%fPSD(1:ir,:,:,:))
        write(*,*) '   TPs           = ', sum(psPop%nPSD(1:ir,:,:,:))
        write(*,*) '----------------------------'
        end associate
    end subroutine writePSD

    subroutine initPSIO(Model,psGr)
        type(chmpModel_T), intent(in) :: Model
        type(PSEq_T), intent(in) :: psGr

        type(IOVAR_T), dimension(MAXPSVS) :: IOVars
        integer :: ir,ip,ik
        real(rp), dimension(:,:,:), allocatable :: X3,Y3,K3
        real(rp), dimension(:,:), allocatable :: X2,Y2

        associate(Nr=>psGr%Nr,Np=>psGr%Np,Nk=>psGr%Nk,Na=>psGr%Na)
        
        !3D grid data (R,P,K)
        write(psOutF,'(2a)') trim(adjustl(Model%RunID)),'.ps.h5'
        call CheckAndKill(psOutF)

        allocate(X3(Nr+1,Np+1,Nk+1))
        allocate(Y3(Nr+1,Np+1,Nk+1))
        allocate(K3(Nr+1,Np+1,Nk+1))

        !Create grid
        do ik=1,Nk+1
            do ip=1,Np+1
                do ir=1,Nr+1
                    X3(ir,ip,ik) = psGr%rI(ir)*cos(psGr%pI(ip))
                    Y3(ir,ip,ik) = psGr%rI(ir)*sin(psGr%pI(ip))
                    K3(ir,ip,ik) = log10(psGr%kI(ik))
                enddo
            enddo
        enddo

        !Write grid
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",X3)
        call AddOutVar(IOVars,"Y",Y3)
        call AddOutVar(IOVars,"Z",K3)
        
        call WriteVars(IOVars,.true.,psOutF)

        !Do 2D part
        write(pseqOutF,'(2a)') trim(adjustl(Model%RunID)),'.pseq.h5'
        call CheckAndKill(pseqOutF)

        allocate(X2(Nr+1,Np+1))
        allocate(Y2(Nr+1,Np+1))
        X2 = X3(:,:,1)
        Y2 = Y3(:,:,1)
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",X2)
        call AddOutVar(IOVars,"Y",Y2)
        call WriteVars(IOVars,.true.,pseqOutF)

        !Do 4D part if necessary
        if (Model%doFat) then
            write(ps4OutF,'(2a)') trim(adjustl(Model%RunID)),'.ps4.h5'
            call CheckAndKill(ps4OutF)

            call ClearIO(IOVars)
            call AddOutVar(IOVars,"X",X3)
            call AddOutVar(IOVars,"Y",Y3)
            call AddOutVar(IOVars,"Z",K3)
            call AddOutVar(IOVars,"K",psGr%kI)
            call AddOutVar(IOVars,"A",psGr%aI)
            call WriteVars(IOVars,.true.,ps4OutF)
        endif

        end associate

    end subroutine initPSIO

    !Update population to time t
    subroutine updatePop(Model,psPop,ebSt,t)
        type(chmpModel_T), intent(in) :: Model
        type(psdPop_T), intent(inout) :: psPop
        type(ebState_T), intent(in)  :: ebSt
        real(rp), intent(in) :: t

        character(len=strLen) :: fIn,gStr
        integer :: n,nStp,np,i0,i1
        type(IOVAR_T), dimension(MAXTPVS) :: IOVars
        real(rp) :: tpScl

        !Calculate step to read
        nStp = nint( (t-psPop%T0)/psPop%dtStp )
        if (nStp<0) nStp = 0

        write(*,*) 'Reading H5p step ', nStp

        write(gStr,'(A,I0)') "Step#", nStp

        !TP timing data is always in seconds
        !Use specific scaling (may disagree) with MHD data scaling
        tpScl = 1.0*vc_cgs/L0 !in2s = 1

        i0 = 1
        !Loop over files, and read this step
        do n=psPop%ns,psPop%ne
            write(fIn,'(a,a,I0.6,a)') trim(adjustl(psPop%popid)),'.',n,'.h5part'
            call AddInVar(IOVars,"xeq")
            call AddInVar(IOVars,"yeq")

            if (doEBKEQ) then
                !Read eb kinetic energy (kin energy in ExB frame)
                call AddInVar(IOVars,"ebKeq")
            else
                !Read equatorial kinetic energy
                call AddInVar(IOVars,"Keq")
            endif

            call AddInVar(IOVars,"Aeq")
            call AddInVar(IOVars,"Teq")
            call AddInVar(IOVars,"isIn")
            call ReadVars(IOVars,.true.,fIn,gStr)


            np = IOVars(1)%dims(1) !Number of particles in this file
            i1 = i0+np-1
            psPop%TPs(i0:i1,PSRAD)   = sqrt(IOVars(1)%data**2.0 + IOVars(2)%data**2.0)
            psPop%TPs(i0:i1,PSPHI)   = atan2(IOVars(2)%data,IOVars(1)%data)
            psPop%TPs(i0:i1,PSKINE)  = IOVars(3)%data
            psPop%TPs(i0:i1,PSALPHA) = (1/rad2deg)*IOVars(4)%data

            psPop%tx(i0:i1) = tpScl*IOVars(5)%data
            !Set boolean for TP alive
            psPop%isIn(i0:i1) = (IOVars(6)%data>0.5)
            
            !Update starting point
            i0 = i0+np

            !Probably don't need to clear every time
            call ClearIO(IOVars) 
        enddo

        !Now we have all the particles, could re-project if necessary

    end subroutine updatePop

    !Write out weights
    subroutine fOutWgts(Model,psGr,psPop)
        type(chmpModel_T), intent(in) :: Model
        type(PSEq_T), intent(in) :: psGr
        type(psdPop_T), intent(in) :: psPop
        
        type(IOVAR_T), dimension(MAXPSVS) :: IOVars

        write(wgtOutF,'(2a)') trim(adjustl(Model%RunID)),'.wgt.h5'
        call CheckAndKill(wgtOutF)

        call ClearIO(IOVars)
        call AddOutVar(IOVars,"wgt",psPop%wgt)
        
        call WriteVars(IOVars,.true.,wgtOutF)

    end subroutine fOutWgts
    
end module psdio