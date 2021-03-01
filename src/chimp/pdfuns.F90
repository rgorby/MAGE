!PSD initial condition routines
!Returns PDF in (keV*s)^(-3)
module pdfuns
    use chmpdefs
    use chmpunits
    use psdtypes
    use xml_input
    use strings
    use ioh5

    implicit none

    procedure(PDFun_T), pointer :: fPSD0 => NULL()
    real(rp), private :: kappa=3.0

    !See src/voltron/sstimag.F90 for example on reading HDF5
    type PSDIN_T
        integer :: Nl,Np,Nk,Na
        real(rp), dimension(:), allocatable :: Li,Ki,Ai,Pi
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
            !Initialize data object
            call InitPSDIn(PSDInput,inpXML)

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
    !Currently for reading in PSD(K,L) in units of (keV*s)^(-3)
    !FIXME: need to update to include alpha and phi the file
    function fHDF5(Model,L,phi,K,alpha,n0,kT0) result(fD)
        type(chmpModel_T), intent(in) :: Model
        real(rp), intent(in) :: L,phi,K,alpha
        real(rp), intent(in), optional :: n0,kT0
        integer :: Nl,Np,Nk,Na
        real(rp), dimension(:), allocatable :: Li,Ki,Ai,Pi ! values at the interface
        real(rp), dimension(:), allocatable :: Lc,Kc,Ac,Pc ! cell centered values
        real(rp), dimension(:,:,:,:), allocatable :: psdH5
        integer :: iL,iK,iP,iA
        integer :: n
        logical :: inL,inP,inK,inA,inPS
        real(rp) :: fD,sinN,paScl

        Li   = PSDInput%Li
        Ki   = PSDInput%Ki
        psdH5 = PSDInput%Fpsd

        !Check if within bounds of psd file
        inL = ( L >= minval(Li) ) .and. ( L <= maxval(Li) )
        inK = ( K >= minval(Ki) ) .and. ( K <= maxval(Ki) )
        inPS = inL .and. inK !.and. inP .and. inA
        
        if (.not. inPS) then
            fD = 0 ! not in domain of file, set PSD to zero
            return
        endif

        ! Calculating cell center values
        allocate(Kc(PSDInput%Nk))
        do n=1,PSDInput%Nk
            Kc(n) = 0.5*(Ki(n)+Ki(n+1))
        enddo

        allocate(Lc(PSDInput%Nl))
        do n=1,PSDInput%Nl
            Lc(n) = 0.5*(Li(n)+Li(n+1))
        enddo

        !Now we know we're in this PS grid
        iL = maxloc(Lc,dim=1,mask=Lc<=L)
        iK = maxloc(Kc,dim=1,mask=Kc<=K)

        ! scaling PSD for equatorial pitch angles to others with sin^n(alpha) fit
        sinN = 0.8
        paScl = sin(alpha)**sinN

        fD = psdH5(iL,iK,1,1)*paScl

    end function fHDF5

    subroutine InitPSDIn(PSDInput,inpXML)
        class(PSDIN_T), intent(inout) :: PSDInput
        type(XML_Input_T), intent(in) :: inpXML

        integer, parameter :: NIOVAR = 6
        type(IOVAR_T), dimension(NIOVAR) :: IOVars
        character(len=strLen) :: psdFile
        integer :: dims(2) ! update when add higher dimensions
        integer :: Nl,Np,Nk,Na,Ndim
        real(rp), dimension(:,:), allocatable :: fLK

        !Read filename of hdf5 from XML
        call inpXML%Set_Val(psdFile,"population/f0data","psd.h5")  
        call CheckFileOrDie(psdFile,"Error opening PSD initial condition file")
        PSDInput%filename = psdFile

        !Read file
        call ClearIO(IOVars)
        call AddInVar(IOVars,"fPSD")
        call AddInVar(IOVars,"Ki")
        call AddInVar(IOVars,"Li")
        call ReadVars(IOVars,.false.,psdFile) !Don't use io precision

        ! FIXME: only compatible with PSD(K,L) need to extend to include alpha & phi
        Ndim = IOVars(1)%Nr
        if ( Ndim /= 2) then
            write(*,*) 'Currently only support input files in the form PSD(L,K) and not higher dimensions'
            stop
        endif
        
        dims = IOVars(1)%dims(1:Ndim)
        Nk = IOVars(2)%N-1
        Nl = IOVars(3)%N-1
        Na = 1
        Np = 1

        if ( Nl /=  dims(1) .or. Nk /= dims(2)) then
            write(*,*) 'Dimensions of PSD file are not compatible with each other'
            stop
        endif

        allocate(PSDInput%Fpsd(Nl,Nk,Na,Np))

        PSDInput%Nk = Nk
        PSDInput%Nl = Nl
        PSDInput%Na = Na
        PSDInput%Np = Np

        PSDInput%Fpsd(:,:,1,1) = reshape(IOVars(1)%data,dims)
        PSDInput%Ki = IOVars(2)%data
        PSDInput%Li = IOVars(3)%data 

    end subroutine InitPSDIn
end module pdfuns
