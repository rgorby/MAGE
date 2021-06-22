module chopio

    use chmpdefs
    use chmpunits
    use ebtypes
    use parintime
    use ioH5
    use xml_input
    use ebinterp
    use files
    use streamline
    use ebtabutils
    
    implicit none

    character(len=strLen) :: eb3DOutF
    integer, parameter :: MAXEBVS = 30
    
    integer  :: Nx1 = 64, Nx2 = 64, Nx3 = 64

    real(rp) :: xyzMax = 10.0 !Default bound
    real(rp), dimension(:,:,:), allocatable :: xxi,yyi,zzi,xxc,yyc,zzc
    integer, private :: NumB = 0

    contains

    subroutine initEB3Dio(Model,ebState,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(XML_Input_T), intent(in) :: inpXML

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars
        character(len=strLen) :: idStr
        real(rp) :: x1Max,x1Min,x2Max,x2Min,x3Max,x3Min
        real(rp) :: x1,x2,x3
        real(rp) :: dx1,dx2,dx3
        integer :: i,j,k
        integer :: iS,iE
        real(rp) :: xcc(NDIM)

    !Check for time parallelism
        call InitParInTime(Model,inpXML,"eb3",eb3DOutF)
    !Setup output file
        associate( ebGr=>ebState%ebGr )
        
        call CheckAndKill(eb3DOutF)

        !Figure out information for chop region
        call inpXML%Set_Val(idStr,'chop/grType',"XYZ")
        !Get 3D bounds, use xyzMax as default max
        call inpXML%Set_Val(x1Max,'chop/x1Max',xyzMax)
        call inpXML%Set_Val(x2Max,'chop/x2Max',xyzMax)
        call inpXML%Set_Val(x3Max,'chop/x3Max',xyzMax)

        !Use negative of maxes as default mins
        call inpXML%Set_Val(x1Min,'chop/x1Min',-x1Max)
        call inpXML%Set_Val(x2Min,'chop/x2Min',-x2Max)
        call inpXML%Set_Val(x3Min,'chop/x3Min',-x3Max)

        call inpXML%Set_Val(Nx1,'chop/Nx1',Nx1)
        call inpXML%Set_Val(Nx2,'chop/Nx2',Nx2)
        call inpXML%Set_Val(Nx3,'chop/Nx3',Nx3)

        !for LFM grid, take all i-shells within x1Max, in the Sun direction
        select case(trim(toUpper(idStr)))
        !---------
        case("LFM")
            !Find first cell that's "in" according to set locator
            i = ebGr%is
            xcc = ebGr%xyz(i,ebGr%js,ebGr%ks,:)
            do while ( .not. inDomain(xcc,Model,ebGr) )
                i = i+1
                xcc = ebGr%xyz(i,ebGr%js,ebGr%ks,:)
            enddo
            iS = i
            iE = minloc(abs(ebGr%xyz(:,ebGr%js,ebGr%ks,XDIR)-x1Max),dim=1)
            write(*,*) 'i-Start/End/Total = ',iS,iE,ebGr%Nip
            
            Nx1 = iE-iS+1
            Nx2 = ebState%ebGr%Njp
            Nx3 = ebState%ebGr%Nkp
        
        end select

        !Allocate chop grids
        allocate(xxi(Nx1+1,Nx2+1,Nx3+1))
        allocate(yyi(Nx1+1,Nx2+1,Nx3+1))
        allocate(zzi(Nx1+1,Nx2+1,Nx3+1))
        
        allocate(xxc(Nx1  ,Nx2  ,Nx3  ))
        allocate(yyc(Nx1  ,Nx2  ,Nx3  ))
        allocate(zzc(Nx1  ,Nx2  ,Nx3  ))

        !Create corners
        select case(trim(toUpper(idStr)))
        !---------
        case("XYZ")
            dx1 = (x1Max-x1Min)/Nx1
            dx2 = (x2Max-x2Min)/Nx2
            dx3 = (x3Max-x3Min)/Nx3
            !$OMP PARALLEL DO
            do k=1,Nx3+1
                do j=1,Nx2+1
                    do i=1,Nx1+1
                        xxi(i,j,k) = x1Min + (i-1)*dx1
                        yyi(i,j,k) = x2Min + (j-1)*dx2
                        zzi(i,j,k) = x3Min + (k-1)*dx3
                    enddo
                enddo
            enddo
        !---------
        ! Spherical
        case("RTP")            
            dx1 = ( log10(x1Max)-log10(x1Min) )/Nx1
            dx2 = (x2Max-x2Min)*deg2rad/Nx2
            dx3 = (x3Max-x3Min)*deg2rad/Nx3
            !$OMP PARALLEL DO
            do k=1, Nx3+1
                do j=1,Nx2+1    
                    do i=1,Nx1+1
                        !x1 = Rin + (i-1)*dx1
                        x1 = 10**( log10(x1Min) + (i-1)*dx1 )
                        x2 = x2Min*deg2rad + (j-1)*dx2
                        x3 = x3Min*deg2rad + (k-1)*dx3
                        xxi(i,j,k) = x1*cos(x2)*sin(x3)
                        yyi(i,j,k) = x1*sin(x2)*sin(x3)
                        zzi(i,j,k) = x1*cos(x3)
                    enddo
                enddo
            enddo
        !---------
        ! only x1Max is used to set how far downtail you chop
        case("LFM")
            !$OMP PARALLEL DO
            do k=1,Nx3+1
                do j=1,Nx2+1
                    do i=1,Nx1+1
                        xxi(i,j,k) = ebGr%xyz(i+iS-1,j,k,XDIR)
                        yyi(i,j,k) = ebGr%xyz(i+iS-1,j,k,YDIR)
                        zzi(i,j,k) = ebGr%xyz(i+iS-1,j,k,ZDIR)
                    enddo
                enddo
            enddo

            case default
                write(*,*) 'Unknown chop type, exiting ...'
                stop

        end select
        !Calculate cell centers (just simple 8 point average)
        !$OMP PARALLEL DO
        do k=1,Nx3
            do j=1,Nx2
                do i=1,Nx1

                    xxc(i,j,k) = 0.125*( xxi(i,j,k  ) + xxi(i+1,j,k  ) + xxi(i,j+1,k  ) + xxi(i+1,j+1,k  ) + &
                                       & xxi(i,j,k+1) + xxi(i+1,j,k+1) + xxi(i,j+1,k+1) + xxi(i+1,j+1,k+1)    )
                    yyc(i,j,k) = 0.125*( yyi(i,j,k  ) + yyi(i+1,j,k  ) + yyi(i,j+1,k  ) + yyi(i+1,j+1,k  ) + &
                                       & yyi(i,j,k+1) + yyi(i+1,j,k+1) + yyi(i,j+1,k+1) + yyi(i+1,j+1,k+1)    )
                    zzc(i,j,k) = 0.125*( zzi(i,j,k  ) + zzi(i+1,j,k  ) + zzi(i,j+1,k  ) + zzi(i+1,j+1,k  ) + &
                                       & zzi(i,j,k+1) + zzi(i+1,j,k+1) + zzi(i,j+1,k+1) + zzi(i+1,j+1,k+1)    )


                enddo
            enddo
        enddo

        !Write grid
        call ClearIO(IOVars)
        call AddOutVar(IOVars,"X",xxi)
        call AddOutVar(IOVars,"Y",yyi)
        call AddOutVar(IOVars,"Z",zzi)

        call WriteVars(IOVars,.true.,eb3DOutF)

        call ClearIO(IOVars)  
        end associate     
    end subroutine initEB3Dio

    !Do output slice for eb data
    !Scale output data using values set in chimpunits
    subroutine writeEB3D(Model,ebState,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars

        real(rp), dimension(:,:,:,:), allocatable :: B,Q,E,J3
        real(rp), dimension(NVARMHD) :: Qijk
        real(rp), dimension(NDIM) :: Bijk,Eijk,xyz
        type(gcFields_T) :: gcFields
        real(rp), dimension(NDIM,NDIM) :: jB
        integer :: i,j,k

        real(rp) :: oJScl

        !Data for tracing
        type(ebTrc_T), dimension(:,:,:), allocatable :: ebTrcIJK

        allocate(B (Nx1,Nx2,Nx3,NDIM))
        allocate(E (Nx1,Nx2,Nx3,NDIM))
        allocate(J3(Nx1,Nx2,Nx3,NDIM))
        allocate(Q (Nx1,Nx2,Nx3,NVARMHD))

        if (Model%doTrc) then
            !Create an ebTrc for each point on slice
            allocate(ebTrcIJK(Nx1,Nx2,Nx3))
        endif

        B = 0.0
        E = 0.0
        J = 0.0
        Q = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyz,Bijk,Eijk,Qijk,gcFields,jB)
        do k=1,Nx3
            do j=1,Nx2
                do i=1,Nx1
                    xyz = [xxc(i,j,k),yyc(i,j,k),zzc(i,j,k)]

                    call ebFields(xyz,Model%t,Model,ebState,Eijk,Bijk,gcFields=gcFields)
                    B(i,j,k,:) = Bijk
                    E(i,j,k,:) = Eijk

                    jB = gcFields%JacB

                    J3(i,j,k,:) = Jac2Curl(jB)

                    if (Model%doMHD) then
                        Qijk = mhdInterp(xyz,Model%t,Model,ebState)
                        Q(i,j,k,:) = Qijk
                    endif

                    if (Model%doTrc) then
                        !Get field line topology stuff
                        call SliceFL(Model,ebState,xyz,Model%t,ebTrcIJK(i,j,k))
                    endif

                enddo
            enddo
        enddo

        !Calculate output J (current) scaling
        !Mu0 [Tm/A] 
        ! Curl(B) * oBScl/L0 => nT/cm x (1.0e-9)/(1.0e-2) => T/m
        !oJScl => A/m^2
        oJScl = ( (1.0e-7)*oBScl/L0 )/Mu0
        
        !Now write out
        call ClearIO(IOVars)

        !-------------------
        !Variables to always output
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"MJD",MJDAt(ebState%ebTab,Model%t))
        call AddOutVar(IOVars,"Bx"  ,oBScl*B(:,:,:,XDIR))
        call AddOutVar(IOVars,"By"  ,oBScl*B(:,:,:,YDIR))
        call AddOutVar(IOVars,"Bz"  ,oBScl*B(:,:,:,ZDIR))

        call AddOutVar(IOVars,"Ex"  ,oEScl*E(:,:,:,XDIR))
        call AddOutVar(IOVars,"Ey"  ,oEScl*E(:,:,:,YDIR))
        call AddOutVar(IOVars,"Ez"  ,oEScl*E(:,:,:,ZDIR))

        call AddOutVar(IOVars,"Jx"  ,oJScl*J3(:,:,:,XDIR))
        call AddOutVar(IOVars,"Jy"  ,oJScl*J3(:,:,:,YDIR))
        call AddOutVar(IOVars,"Jz"  ,oJScl*J3(:,:,:,ZDIR))

        if (Model%doMHD) then
            call AddOutVar(IOVars,"Vx"  ,oVScl*Q(:,:,:,VELX    ))
            call AddOutVar(IOVars,"Vy"  ,oVScl*Q(:,:,:,VELY    ))
            call AddOutVar(IOVars,"Vz"  ,oVScl*Q(:,:,:,VELZ    ))
            call AddOutVar(IOVars,"D"   ,      Q(:,:,:,DEN     ))
            call AddOutVar(IOVars,"P"   ,      Q(:,:,:,PRESSURE))
        endif

        if (Model%doTrc) then
            !Field line tracing metrics
            call AddOutVar(IOVars,"OCb" ,ebTrcIJK(:,:,:)%OCb )

            !Equator and end-points
            call AddOutVar(IOVars,"xBEQ",ebTrcIJK(:,:,:)%MagEQ(XDIR))
            call AddOutVar(IOVars,"yBEQ",ebTrcIJK(:,:,:)%MagEQ(YDIR))
            call AddOutVar(IOVars,"zBEQ",ebTrcIJK(:,:,:)%MagEQ(ZDIR))

            call AddOutVar(IOVars,"xP",ebTrcIJK(:,:,:)%xEPp(XDIR))
            call AddOutVar(IOVars,"yP",ebTrcIJK(:,:,:)%xEPp(YDIR))
            call AddOutVar(IOVars,"zP",ebTrcIJK(:,:,:)%xEPp(ZDIR))

            call AddOutVar(IOVars,"xM",ebTrcIJK(:,:,:)%xEPm(XDIR))
            call AddOutVar(IOVars,"yM",ebTrcIJK(:,:,:)%xEPm(YDIR))
            call AddOutVar(IOVars,"zM",ebTrcIJK(:,:,:)%xEPm(ZDIR))

        endif

        call WriteVars(IOVars,.true.,eb3DOutF,gStr)
        call ClearIO(IOVars)
        
    end subroutine writeEB3D
end module chopio
