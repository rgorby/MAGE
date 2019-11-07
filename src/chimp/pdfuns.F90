!PSD initial condition routines
!Returns PDF in (keV*s)^(-3)
module pdfuns
    use chmpdefs
    use chmpunits
    use psdtypes
    use xml_input
    use strings

    implicit none

    procedure(PDFun_T), pointer :: fPSD0 => NULL()
    real(rp) :: kappa=3.0

    contains

    !Set f0 information
    subroutine SetPSD0(Model,inpXML)
        type(chmpModel_T), intent(inout) :: Model
        type(XML_Input_T), intent(in)    :: inpXML
        
        character(len=strLen) :: f0Str

        call inpXML%Set_Val(f0Str,"population/f0","Max")

        select case (trim(toUpper(f0Str)))
        !-------
        case("MAXWELLIAN","MAX")
            fPSD0 => fMaxwellian
        case("KAPPA","KAP")
            fPSD0 => fKappa
            call inpXML%Set_Val(kappa,"population/k0",3.0_rp)
        case default
            write(*,*) '<Unknown f0, using Maxwellian>'
            fPSD0 => fMaxwellian
        end select

    end subroutine SetPSD0

    function fMaxwellian(Model,L,phi,K,alpha,n0,kT0) result(fD)
        type(chmpModel_T), intent(in) :: Model
        real(rp), intent(in) :: L,phi,K,alpha
        real(rp), intent(in), optional :: n0,kT0
        real(rp) :: fD

        real(rp) :: e0,fScl
        e0 = Model%m0*mec2*1.0e+3 !Particle rest mass energy [keV]
        fScl = (2*PI)**(-1.5)
        fScl = fScl*(vc_cgs**3.0)*(e0**(-1.5))
        fD = fScl*(n0)*(kT0**(-1.5))*exp(-K/kT0)

    end function fMaxwellian

    function fKappa(Model,L,phi,K,alpha,n0,kT0) result(fD)
        type(chmpModel_T), intent(in) :: Model
        real(rp), intent(in) :: L,phi,K,alpha
        real(rp), intent(in), optional :: n0,kT0
        real(rp) :: fD

        real(rp) :: e0,K0,gScl,kScl,Ak,fPow
        e0 = Model%m0*mec2*1.0e+3 !Particle rest mass energy [keV]
        K0 = kT0*(kappa-1.5)/kappa
        gScl = gamma(kappa)/gamma(kappa-0.5)
        kScl = sqrt(kappa)*(2.0*e0*PI*K0)**1.5
        Ak = n0*(vc_cgs**3.0)*gScl/kScl
        fPow = (1.0 + (K/(kappa*K0)))**(-kappa-1.0)
        fD = Ak*fPow
        
    end function fKappa
    
end module pdfuns
