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
    real(rp), private :: kappa=3.0

    !See src/voltron/sstimag.F90 for example on reading HDF5
    type PSDIN_T
        integer :: Nr,Np,Nk,Na
        real(rp), dimension(:,:), allocatable :: X,Y,xxc,yyc
        real(rp), dimension(:,:,:,:), allocatable :: Fpsd
        character(len=strLen) :: filename
    end type PSDIN_T
    
    type(PSDIN_T), private :: PSDInput
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
        case("RBPSD")
            fPSD0 => fRBPSD
        case("HDF5IN")
            fPSD0 => fHDF5
            !Read filename of hdf5 from XML
            call inpXML%Set_Val(PSDInput%filename,"population/f0data","psd.h5")
            !Initialize data object
            !call InitPSDIn()

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

    !Function for radiation belt trapped population
    function fRBPSD(Model,L,phi,K,alpha,n0,kT0) result(fD)
        type(chmpModel_T), intent(in) :: Model
        real(rp), intent(in) :: L,phi,K,alpha
        real(rp), intent(in), optional :: n0,kT0
        real(rp) :: fD

        real(rp) :: AScl,B,a,n
        real(rp) :: L0,Lm,Flk,lScl,nXP

        AScl = 4.54238496204e+26 !4.5e5 ![keV^-3 cm^-1 s^-3 sr^-1]
        B = -3.5
        a = 0.5
        n = 2.0
        
        if (K >= 1000.0) then
            L0 = 3.5
            Lm = 4.5
        else
            L0 = 4.75-1.25*(K*1.0e-3)
            Lm = 6.37-1.88*(K*1.0e-3)
        endif

        if (L < L0) then
            Flk = 0.0
        else
            lScl = (L-L0)/(Lm-L0)
            nXP = (2.0/n)*(1.0-lScl**n)
            Flk  = lScl**2.0*exp(nXP)
        endif

        fD = AScl*(K**B)*(sin(alpha)**a)*Flk

    end function fRBPSD

    !Function for taking PSD from a file
    function fHDF5(Model,L,phi,K,alpha,n0,kT0) result(fD)
        type(chmpModel_T), intent(in) :: Model
        real(rp), intent(in) :: L,phi,K,alpha
        real(rp), intent(in), optional :: n0,kT0
        real(rp) :: fD

        fD = 0.0

    end function fHDF5
end module pdfuns
