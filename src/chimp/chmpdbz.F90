!Various possible background field options
module chmpdbz
    use chmpdefs
    use chmpunits
    use ebtypes
    use math
    use xml_input
    implicit none

    contains

    subroutine setBackground(Model,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(inout) :: inpXML


        !Model%B0 => NullB0
        Model%B0 => DipoleB0
        Model%JacB0 => DipoleJ

    end subroutine setBackground
    !----------------------------
    !Dipole options
    function DipoleB0(r) result(Bxyz)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: Bxyz(NDIM)
        real(rp), dimension(NDIM) :: m
        real(rp) :: rad
        
        rad = norm2(r)
        if (rad<=TINY) then
            Bxyz = 0.0
        else
            m = [0.0_rp,0.0_rp,MagM0]
            Bxyz = 3*dot_product(m,r)*r/rad**5.0 - m/rad**3.0
        endif

    end function DipoleB0

    !Jacobian for dipole
    !TODO: Simplify this for hard-coded z-dipole
    function DipoleJ(r) result(Jb0)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: Jb0(NDIM,NDIM)

        real(rp) :: rad,rm5
        real(rp), dimension(NDIM) :: zhat,B1,B2,B3,M1
        real(rp), dimension(NDIM,NDIM) :: Tx,Ty,Tz

        rad = norm2(r)
        if (rad<=TINY) then
            Jb0 = 0.0
            return
        endif

        associate(x=>r(XDIR),y=>r(YDIR),z=>r(ZDIR))

        Tx(XDIR,:) = [2*x,y,z]
        Tx(YDIR,:) = [y,0.0_rp,0.0_rp]
        Tx(ZDIR,:) = [z,0.0_rp,0.0_rp]

        Ty(XDIR,:) = [0.0_rp,x,0.0_rp]
        Ty(YDIR,:) = [x,2*y,z]
        Ty(ZDIR,:) = [0.0_rp,z,0.0_rp]

        Tz(XDIR,:) = [0.0_rp,0.0_rp,x]
        Tz(YDIR,:) = [0.0_rp,0.0_rp,y]
        Tz(ZDIR,:) = [x,y,2*z]

        zhat = [0.0,0.0,1.0]

        rm5 = rad**(-5.0)
        B1 = 3*MagM0*r*dot_product(r,zhat)*rm5
        B2 = -zhat*MagM0/(rad**3.0)
        B3 = (-5*B1 -3*B2)/(rad**2.0)
        M1 = 3*MagM0*zhat*rm5

        Jb0(:,XDIR) = B3*r(XDIR) + matmul(Tx,M1)
        Jb0(:,YDIR) = B3*r(YDIR) + matmul(Ty,M1)
        Jb0(:,ZDIR) = B3*r(ZDIR) + matmul(Tz,M1)
        
        end associate
    end function DipoleJ


    !----------------------------
    !Null options, use by default
    function NullB0(r) result(Bxyz)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: Bxyz(NDIM)
        
        Bxyz = 0        
    end function NullB0
    function NullJacB0(r) result(Jb0)
        real(rp), intent(in) :: r(NDIM)
        real(rp) :: Jb0(NDIM,NDIM)

        Jb0 = 0.0
    end function NullJacB0

end module chmpdbz