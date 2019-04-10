module sliceio

    use chmpdefs
    use chmpunits
    use ebtypes
    use streamline
    use ioH5
    use xml_input

    implicit none

    character(len=strLen) :: ebOutF
    integer, parameter :: MAXEBVS = 30
    !Parameters for output EB grid (2D equatorial)
    !Fixme: Clean up and generalize
    integer  :: Nx1 = 128, Nx2 = 256
    real(rp) :: xSun = 12.5,xTail=-20.0,yM=20.0 !Default Cartesian slice
    real(rp) :: dx0=0.05
    real(rp) :: z0=0.05
    real(rp), dimension(:,:), allocatable :: xxi,yyi,xxc,yyc
    real(rp), dimension(:,:,:), allocatable :: B02D
    logical :: doXY = .true. !Do XY or XZ slice
    logical :: doVelGC = .false.

    contains

    subroutine initEBio(Model,ebState,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(XML_Input_T), intent(in) :: inpXML

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars
        real(rp) :: dx1,dx2,Rin,Rout,x1,x2
        integer :: i,j,ijk(NDIM)
        integer :: iS,iE,n,Npow
        real(rp) :: xcc(NDIM)
        real(rp), dimension(:,:,:), allocatable :: ijkXY
        character(len=strLen) :: idStr
        
        write(ebOutF,'(2a)') trim(adjustl(Model%RunID)),'.eb.h5'

        associate( ebGr=>ebState%ebGr )

        call CheckAndKill(ebOutF)


        !Figure out grid slice type
        call inpXML%Set_Val(doXY,'slice/doXY',.true.)
        call inpXML%Set_Val(idStr,'slice/grType',"XY")
        call inpXML%Set_Val(Npow,'slice/Npow',0) !Number of times to double

        select case(trim(toUpper(idStr)))
            case("XY")
                call inpXML%Set_Val(xSun,'slice/xSun',xSun)
                call inpXML%Set_Val(xTail,'slice/xTail',xTail)
                call inpXML%Set_Val(yM,'slice/yM',yM)
                Nx1 = nint((xSun-xTail)/dx0)
                Nx2 = nint(2*yM/dx0)

                call inpXML%Set_Val(Nx1,'slice/Nx1',Nx1)
                call inpXML%Set_Val(Nx2,'slice/Nx2',Nx2)
            case("RP")
                Nx1 = ebState%ebGr%Nip
                Nx2 = 2*ebState%ebGr%Njp
                call inpXML%Set_Val(Nx1,'slice/Nx1',Nx1)
                call inpXML%Set_Val(Nx2,'slice/Nx2',Nx2)
                !Do polar grid, take Rin/Rout from sunward line
                Rin  = ebGr%xyz(ebGr%is+1,ebGr%js,ebGr%ks,XDIR)
                Rout = ebGr%xyz(ebGr%ie-1,ebGr%js,ebGr%ks,XDIR)

                call inpXML%Set_Val(Rin,'slice/Rin',Rin)
                call inpXML%Set_Val(Rout,'slice/Rout',Rout)
                write(*,*) 'Radial grid bounds = ', Rin,Rout

            case("LFM2D")
                !Create 2D LFM slice (w/ full 2pi)
                call inpXML%Set_Val(xSun,'slice/xSun',xSun)
                !Find first cell that's "in" according to set locator
                i = ebGr%is
                xcc = ebGr%xyz(i,ebGr%js,ebGr%ks,:)
                do while ( .not. inDomain(xcc,Model,ebGr) )
                    i = i+1
                    xcc = ebGr%xyz(i,ebGr%js,ebGr%ks,:)
                enddo
                iS = i
                iE = minloc(abs(ebGr%xyz(:,ebGr%js,ebGr%ks,XDIR)-xSun),dim=1)
                write(*,*) 'i-Start/End/Total = ',iS,iE,ebGr%Nip
                
                Nx1 = iE-iS+1
                Nx2 = 2*ebState%ebGr%Njp
        end select

        !Create 2D grid, wait to allocate rest until after embiggening
        allocate(xxi(Nx1+1,Nx2+1))
        allocate(yyi(Nx1+1,Nx2+1))

        !Create grid
        select case(trim(toUpper(idStr)))
        !---------
        case("XY")
            
            dx1 = (xSun-xTail)/Nx1
            dx2 = 2*yM/Nx2
            do j=1,Nx2+1
                do i=1,Nx1+1
                    x1 = xTail + (i-1)*dx1
                    x2 = -yM + (j-1)*dx2
                    xxi(i,j) = x1
                    yyi(i,j) = x2

                enddo
            enddo
        !---------
        case("RP")            
            dx1 = ( log10(Rout)-log10(Rin) )/Nx1
            dx2 = (2*PI-  0)/Nx2
            do j=1,Nx2+1
                do i=1,Nx1+1
                    !x1 = Rin + (i-1)*dx1
                    x1 = 10**( log10(Rin) + (i-1)*dx1 )
                    x2 = 0.0 + (j-1)*dx2
                    xxi(i,j) = x1*cos(x2)
                    yyi(i,j) = x1*sin(x2)
                enddo
            enddo
        !---------
        case("LFM2D")
            !Do upper half plane
            do j=1,(Nx2/2)+1
                do i=1,Nx1+1
                    xxi(i,j) = ebGr%xyz(i+iS-1,j,ebGr%ks,XDIR)
                    yyi(i,j) = ebGr%xyz(i+iS-1,j,ebGr%ks,YDIR)
                enddo
            enddo
            !Reflect for bottom
            do j=1,Nx2/2
                do i=1,Nx1+1
                    xxi(i,j+1+Nx2/2) =  ebGr%xyz(i+iS-1,Nx2/2 - j+1,ebGr%ks,XDIR)
                    yyi(i,j+1+Nx2/2) = -ebGr%xyz(i+iS-1,Nx2/2 - j+1,ebGr%ks,YDIR)
                enddo
            enddo
        end select

        !Embiggen grid (i.e. upscale)
        do n=1,Npow
            write(*,*) 'Embiggening grid ...'
            write(*,*) 'Before: ',Nx1,Nx2
            call Embiggen(xxi,yyi,Nx1,Nx2)
            write(*,*) 'After : ',Nx1,Nx2
        enddo

        !Allocate remaining holders
        allocate(xxc(Nx1,Nx2),yyc(Nx1,Nx2))
        !Hold mappings from cell center of 2D grid to original XYZ grid cells
        allocate(ijkXY(Nx1,Nx2,NDIM))

        !Calculate cell-centers/background field
        allocate(B02D(Nx1,Nx2,NDIM))
        do j=1,Nx2
            do i=1,Nx1
                xxc(i,j) = 0.25*( xxi(i,j) + xxi(i+1,j) + xxi(i,j+1) + xxi(i+1,j+1) )
                yyc(i,j) = 0.25*( yyi(i,j) + yyi(i+1,j) + yyi(i,j+1) + yyi(i+1,j+1) )
                if (doXY) then
                    xcc(XDIR) = xxc(i,j)
                    xcc(YDIR) = yyc(i,j)
                    xcc(ZDIR) = 0.0
                else
                    !Assuming XZ, x->x and y->z
                    xcc(XDIR) = xxc(i,j)
                    xcc(YDIR) = 0.0
                    xcc(ZDIR) = yyc(i,j)
                endif
                
                if (inDomain(xcc,Model,ebGr)) then
                    B02D(i,j,:) = oBScl*Model%B0(xcc)
                    !Get 2D->3D mapping info for debugging
                    call locate(xcc,ijk,Model,ebGr)
                    ijkXY(i,j,:) = 1.0*ijk + Map2ezp(xcc,ijk,Model,ebGr)
                else
                    B02D(i,j,:) = 0.0
                    ijkXY(i,j,:)= 0.0
                endif
            enddo
        enddo

        !Write grid
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",xxi)
        call AddOutVar(IOVars,"Y",yyi)
        call AddOutVar(IOVars,"Bx0",B02D(:,:,XDIR))
        call AddOutVar(IOVars,"By0",B02D(:,:,YDIR))
        call AddOutVar(IOVars,"Bz0",B02D(:,:,ZDIR))
        call AddOutVar(IOVars,"iG",ijkXY(:,:,IDIR))
        call AddOutVar(IOVars,"jG",ijkXY(:,:,JDIR))
        call AddOutVar(IOVars,"kG",ijkXY(:,:,KDIR))            

        call WriteVars(IOVars,.true.,ebOutF)
        call ClearIO(IOVars)

        end associate
    end subroutine initEBio

    !Do output slice for eb data
    !Scale output data using values set in chimpunits
    subroutine writeEB(Model,ebState,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars
        real(rp), dimension(:,:,:), allocatable :: dB2D,E2D,Q,J2D,MagEQ
        real(rp), dimension(:,:,:), allocatable :: Vcrv,Vgrd
        real(rp), dimension(:,:), allocatable :: OCb,dvB,bS,Vr,Lb

        integer :: i,j
        real(rp), dimension(NDIM) :: xp,xm,dB,Ep,Em,Bp,Bm
        real(rp) :: MagB,MagJ,oVGScl
        real(rp), dimension(NVARMHD) :: Qij
        type(gcFields_T) :: gcFields
        associate( ebGr=>ebState%ebGr,ebTab=>ebState%ebTab,eb1=>ebState%eb1,eb2=>ebState%eb2 )


        allocate( dB2D(Nx1,Nx2,NDIM))
        allocate(  E2D(Nx1,Nx2,NDIM))
        allocate(  J2D(Nx1,Nx2,NDIM))
        allocate(MagEQ(Nx1,Nx2,NDIM))
        allocate( Vcrv(Nx1,Nx2,NDIM))
        allocate( Vgrd(Nx1,Nx2,NDIM))
        allocate(OCb(Nx1,Nx2))
        allocate(dvB(Nx1,Nx2))
        allocate(bS (Nx1,Nx2))
        allocate(Vr (Nx1,Nx2))
        allocate(Lb (Nx1,Nx2))
        dB2D  = 0.0
        E2D   = 0.0
        J2D   = 0.0
        MagEQ = 0.0
        OCb = 0.0
        dvB = 0.0
        bS  = 0.0
        Vr  = 0.0
        Lb  = 0.0
        Vcrv = 0.0
        Vgrd = 0.0

        if (Model%doFat) then
            doVelGC = .true.
        endif

        if (Model%doMHD) then
            allocate(Q(Nx1,Nx2,NVARMHD))
            Q = 0.0
        endif
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,xp,xm,Bp,Bm,Ep,Em,dB,Qij,gcFields,MagB,MagJ)
        do j=1,Nx2
            do i=1,Nx1
                !Straddle slice plane
                !Get fields at x,y,z0
                if (doXY) then
                    xp = [xxc(i,j),yyc(i,j), z0]
                    xm = [xxc(i,j),yyc(i,j),-z0]
                else
                    xp = [xxc(i,j), z0,yyc(i,j)]
                    xm = [xxc(i,j),-z0,yyc(i,j)]
                endif

                call ebFields(xp,Model%t,Model,ebState,Ep,Bp)
                call ebFields(xm,Model%t,Model,ebState,Em,Bm)

                !Background already scaled to output units
                db = oBScl*0.5*(Bp+Bm)-B02D(i,j,:)
                dB2D(i,j,:) = db
                E2D (i,j,:) = oEScl*0.5*(Em+Ep)

                !Get MHD vars if requested
                if (Model%doMHD) then
                    !Qij = mhdInterp(xp,Model%t,Model,ebState)
                    !Q(i,j,:) = Qij
                    Qij = 0.5*(mhdInterp(xp,Model%t,Model,ebState) + mhdInterp(xm,Model%t,Model,ebState))
                    Q(i,j,:) = Qij
                    Vr(i,j) = (xxc(i,j)*Qij(VELX) + yyc(i,j)*Qij(VELY))/norm2([xxc(i,j),yyc(i,j)])
                endif
                !Do only point on slice plane
                if (doXY) then
                    xp = [xxc(i,j),yyc(i,j),0.0_rp]
                else
                    xp = [xxc(i,j),0.0_rp,yyc(i,j)]
                endif    

                !Get current
                call ebFields(xp,Model%t,Model,ebState,Ep,Bp,gcFields=gcFields)
                J2D(i,j,:) = Jac2Curl(gcFields%JacB)

                !Get B lengthscale
                MagB = norm2(db+B02D(i,j,:))
                MagJ = oBScl*sqrt(sum(gcFields%JacB**2.0))
                Lb(i,j) = MagB/max(MagJ,TINY)

                if (doVelGC) then
                    !Get drift terms
                    call VelGC(Bp,gcFields%JacB,Vgrd(i,j,:),Vcrv(i,j,:))
                endif
                if (Model%doTrc) then
                    !Get field line topology stuff
                    call SliceFL(Model,ebState,xp,Model%t,OCb(i,j),dvB(i,j),bS(i,j),MagEQ(i,j,:))
                endif
            enddo
        enddo

        call ClearIO(IOVars)

        !-------------------
        !Variables to always output
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"dBx",db2D(:,:,XDIR))
        call AddOutVar(IOVars,"dBy",db2D(:,:,YDIR))
        call AddOutVar(IOVars,"dBz",db2D(:,:,ZDIR))
        call AddOutVar(IOVars,"Lb" ,Lb  (:,:)     )

        if (Model%doTrc) then
            call AddOutVar(IOVars,"OCb",OCb)
            call AddOutVar(IOVars,"dvB",dvB)
            call AddOutVar(IOVars,"bS",bS)
            call AddOutVar(IOVars,"zBEQ",MagEQ(:,:,ZDIR))
        endif

        if (.not. Model%doSlim) then
            call AddOutVar(IOVars,"Ex" , E2D(:,:,XDIR))
            call AddOutVar(IOVars,"Ey" , E2D(:,:,YDIR))
            call AddOutVar(IOVars,"Ez" , E2D(:,:,ZDIR))

            call AddOutVar(IOVars,"Jx" , J2D(:,:,XDIR))
            call AddOutVar(IOVars,"Jy" , J2D(:,:,YDIR))
            call AddOutVar(IOVars,"Jz" , J2D(:,:,ZDIR))
        endif        

        if (Model%doMHD) then
            call AddOutVar(IOVars,"Vx" , oVScl*Q(:,:,VELX))
            call AddOutVar(IOVars,"Vy" , oVScl*Q(:,:,VELY))
            call AddOutVar(IOVars,"Vz" , oVScl*Q(:,:,VELZ))
            call AddOutVar(IOVars,"Vr" , oVScl*Vr)
            call AddOutVar(IOVars,"D"  ,       Q(:,:,DEN))
            call AddOutVar(IOVars,"P"  ,       Q(:,:,PRESSURE))
        endif

        !Convert gradient drift to km/s/keV
        oVGScl = 0.5*3.1*(1.0e+2)*4.6*(1.0e+3)/(6380.0)/oBScl

        if (doVelGC) then
            call AddOutVar(IOVars,"Vgx", oVGScl*Vgrd(:,:,XDIR))
            call AddOutVar(IOVars,"Vgy", oVGScl*Vgrd(:,:,YDIR))
            call AddOutVar(IOVars,"Vgz", oVGScl*Vgrd(:,:,ZDIR))
            call AddOutVar(IOVars,"Vcx", oVGScl*Vcrv(:,:,XDIR))
            call AddOutVar(IOVars,"Vcy", oVGScl*Vcrv(:,:,YDIR))
            call AddOutVar(IOVars,"Vcz", oVGScl*Vcrv(:,:,ZDIR))
        endif

        call WriteVars(IOVars,.true.,ebOutF,gStr)
        call ClearIO(IOVars)
        deallocate(dB2D,E2D)

        end associate
    end subroutine writeEB

    subroutine SliceFL(Model,ebState,x0,t,OCb,dvB,S,xEQ)
        real(rp), intent(in) :: x0(NDIM),t
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        
        real(rp), intent(out) :: dvB,S,OCb,xEQ(NDIM)

        type(fLine_T) :: bTrc
        real(rp) :: bEQ
        OCb = 0
        dvB = 0.0
        S = 0.0
        xEQ = 0.0

        if (.not. inDomain(x0,Model,ebState%ebGr)) return
        !Trace field line
        call genStream(Model,ebState,x0,t,bTrc)

        !Get diagnostics
        OCb = 1.0*FLTop(Model,ebState%ebGr,bTrc)
        if (OCb > 0) then
            dvB = FLVol    (Model,ebState%ebGr,bTrc)
            S   = FLEntropy(Model,ebState%ebGr,bTrc)
            call FLEq(Model,bTrc,xEQ,bEQ)
        endif
        !write(*,*) 'FL size = ', bTrc%Nm+bTrc%Np+1
    end subroutine SliceFL

    !Calculate drift velocities
    subroutine VelGC(B,JacB,Vgrd,Vcrv)
        real(rp), dimension(NDIM), intent(in) :: B
        real(rp), dimension(NDIM,NDIM), intent(in) :: JacB
        real(rp), dimension(NDIM), intent(out) :: Vgrd, Vcrv

        real(rp) :: MagB
        real(rp), dimension(NDIM) :: bhat,gMagB,cbhat,kappa,alpha

        real(rp), dimension(NDIM,NDIM) :: gbhat
        Vgrd = 0.0
        Vcrv = 0.0

        MagB = norm2(B)
        bhat = normVec(B)

        if (MagB < TINY) then
            return
        endif
        !gMagB = gradient(|B|), vector
        !      = bhat \cdot \grad \vec{B}
        gMagB = VdT(bhat,JacB)

        !Curl of bhat
        !Curl(bhat) = (1/MagB)*( Curl(B) + bhat \cross gMagB)
        cbhat = (1/MagB)*( Jac2Curl(JacB) + cross(bhat,gMagB) )

        !Get grad bhat
        gbhat = ( MagB*JacB - Dyad(gMagB,B) )/(MagB**2.0)

        kappa = cross(cbhat,bhat)
        !alpha = bhat cross bhat dot grad bhat
        alpha = cross(bhat,VdT(bhat,gbhat))
        Vgrd = cross(B,gMagB)/(MagB**3.0)
        !write(*,*) 'Vcrv = ', -2*cross(kappa,B)/(MagB**2.0), 2*alpha/MagB
        !Vcrv = -2*cross(kappa,B)/(MagB**2.0)
        Vcrv =  2*alpha/MagB

    end subroutine VelGC

    !Double grid from corners
    subroutine Embiggen(xxi,yyi,Nx1,Nx2)
        real(rp), allocatable, intent(inout) :: xxi(:,:),yyi(:,:)
        integer, intent(inout) :: Nx1,Nx2

        real(rp), dimension(:,:), allocatable :: xxO,yyO
        integer :: i,j,ip,jp,nNx1,nNx2

        write(*,*) shape(xxi)
        !Start by saving old corners
        allocate(xxO(Nx1+1,Nx2+1))
        allocate(yyO(Nx1+1,Nx2+1))

        xxO = xxi
        yyO = yyi

        !Now reallocate xxi/yyi to bigger size
        nNx1 = 2*Nx1
        nNx2 = 2*Nx2

        deallocate(xxi)
        deallocate(yyi)
        allocate(xxi(nNx1+1,nNx2+1))
        allocate(yyi(nNx1+1,nNx2+1))
        !Embed old points into new grid
        !$OMP PARALLEL DO default(shared) private(i,j,ip,jp)
        do j=1,Nx2+1
            do i=1,Nx1+1
                ip = 2*i-1
                jp = 2*j-1
                xxi(ip,jp) = xxO(i,j)
                yyi(ip,jp) = yyO(i,j)
            enddo
        enddo

        !Create I-midpoints
        !$OMP PARALLEL DO default(shared) private(i,j,ip,jp)
        do i=1,Nx1
            xxi(2*i,:) = 0.5*( xxi(2*i-1,:) + xxi(2*i+1,:) )
            yyi(2*i,:) = 0.5*( yyi(2*i-1,:) + yyi(2*i+1,:) )
        enddo
        !Create J-midpoints
        !$OMP PARALLEL DO default(shared) private(i,j,ip,jp)
        do j=1,Nx2
            xxi(:,2*j) = 0.5*( xxi(:,2*j-1) + xxi(:,2*j+1) )
            yyi(:,2*j) = 0.5*( yyi(:,2*j-1) + yyi(:,2*j+1) )
        enddo

        !$OMP PARALLEL DO default(shared) private(i,j,ip,jp)
        do j=1,Nx2
            do i=1,Nx1
                ip = 2*i
                jp = 2*j
                !Create new I-J midpoints
                xxi(ip,jp) = 0.25*( xxi(ip-1,jp-1) + xxi(ip-1,jp+1) + xxi(ip+1,jp-1) + xxi(ip+1,jp+1) )
                yyi(ip,jp) = 0.25*( yyi(ip-1,jp-1) + yyi(ip-1,jp+1) + yyi(ip+1,jp-1) + yyi(ip+1,jp+1) )
            enddo
        enddo

        !Finally, change sizing information
        Nx1 = nNx1
        Nx2 = nNx2
    end subroutine Embiggen
end module sliceio