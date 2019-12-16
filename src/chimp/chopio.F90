module chopio

    use chmpdefs
    use chmpunits
    use ebtypes
    use ioH5
    use xml_input
    use ebinterp
    use files
    
    implicit none

    character(len=strLen) :: eb3DOutF
    integer, parameter :: MAXEBVS = 30
    
    integer  :: Nx1 = 64, Nx2 = 64, Nx3 = 64

    real(rp) :: xyzMax = 10.0 !Default bound
    real(rp), dimension(:,:,:), allocatable :: xxi,yyi,zzi,xxc,yyc,zzc

    contains

    subroutine initEB3Dio(Model,ebState,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(ebState_T), intent(in)   :: ebState
        type(XML_Input_T), intent(in) :: inpXML

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars
        character(len=strLen) :: idStr
        real(rp) :: x1Max,x1Min,x2Max,x2Min,x3Max,x3Min
        real(rp) :: dx1,dx2,dx3
        integer :: i,j,k

        write(eb3DOutF,'(2a)') trim(adjustl(Model%RunID)),'.eb.h5'

        
        call CheckAndKill(eb3DOutF)

        !Figure out information for chop region
        call inpXML%Set_Val(idStr,'chop/grType',"XYZ")
        select case(trim(toUpper(idStr)))
            case("XYZ")
                !Get 3D bounds, use xyzMax as default max
                call inpXML%Set_Val(x1Max,'chop/xMax',xyzMax)
                call inpXML%Set_Val(x2Max,'chop/yMax',xyzMax)
                call inpXML%Set_Val(x3Max,'chop/zMax',xyzMax)

                !Use negative of maxes as default mins
                call inpXML%Set_Val(x1Min,'chop/xMin',-x1Max)
                call inpXML%Set_Val(x2Min,'chop/yMin',-x2Max)
                call inpXML%Set_Val(x3Min,'chop/zMin',-x3Max)

                call inpXML%Set_Val(Nx1,'chop/Nx1',Nx1)
                call inpXML%Set_Val(Nx2,'chop/Nx2',Nx2)
                call inpXML%Set_Val(Nx3,'chop/Nx3',Nx3)
            case default
                write(*,*) 'Unknown chop type, exiting ...'
                stop
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
    end subroutine initEB3Dio

    !Do output slice for eb data
    !Scale output data using values set in chimpunits
    subroutine writeEB3D(Model,ebState,gStr)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T), intent(in)   :: ebState
        character(len=strLen), intent(in) :: gStr

        type(IOVAR_T), dimension(MAXEBVS) :: IOVars

        real(rp), dimension(:,:,:,:), allocatable :: B,Q,E
        real(rp), dimension(NVARMHD) :: Qijk
        real(rp), dimension(NDIM) :: Bijk,Eijk,xyz
        integer :: i,j,k

        allocate(B(Nx1,Nx2,Nx3,NDIM))
        allocate(E(Nx1,Nx2,Nx3,NDIM))
        allocate(Q(Nx1,Nx2,Nx3,NVARMHD))

        B = 0.0
        E = 0.0
        Q = 0.0

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyz,Bijk,Eijk,Qijk)
        do k=1,Nx3
            do j=1,Nx2
                do i=1,Nx1
                    xyz = [xxc(i,j,k),yyc(i,j,k),zzc(i,j,k)]

                    call ebFields(xyz,Model%t,Model,ebState,Eijk,Bijk)
                    B(i,j,k,:) = Bijk
                    E(i,j,k,:) = Eijk

                    if (Model%doMHD) then
                        Qijk = mhdInterp(xyz,Model%t,Model,ebState)
                        Q(i,j,k,:) = Qijk
                    endif

                enddo
            enddo
        enddo


        !Now write out
        call ClearIO(IOVars)

        !-------------------
        !Variables to always output
        call AddOutVar(IOVars,"time",oTScl*Model%t)
        call AddOutVar(IOVars,"Bx"  ,oBScl*B(:,:,:,XDIR))
        call AddOutVar(IOVars,"By"  ,oBScl*B(:,:,:,YDIR))
        call AddOutVar(IOVars,"Bz"  ,oBScl*B(:,:,:,ZDIR))

        call AddOutVar(IOVars,"Ex"  ,oEScl*E(:,:,:,XDIR))
        call AddOutVar(IOVars,"Ey"  ,oEScl*E(:,:,:,YDIR))
        call AddOutVar(IOVars,"Ez"  ,oEScl*E(:,:,:,ZDIR))
        if (Model%doMHD) then
            call AddOutVar(IOVars,"Vx"  ,oVScl*Q(:,:,:,VELX    ))
            call AddOutVar(IOVars,"Vy"  ,oVScl*Q(:,:,:,VELY    ))
            call AddOutVar(IOVars,"Vz"  ,oVScl*Q(:,:,:,VELZ    ))
            call AddOutVar(IOVars,"D"   ,      Q(:,:,:,DEN     ))
            call AddOutVar(IOVars,"P"   ,      Q(:,:,:,PRESSURE))
        endif

        call WriteVars(IOVars,.true.,eb3DOutF,gStr)
        call ClearIO(IOVars)
        
    end subroutine writeEB3D
end module chopio