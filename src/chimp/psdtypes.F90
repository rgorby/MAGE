!Main data types for PSD calculations

module psdtypes
    use chmpdefs
    use streamline
    
    implicit none

    real(rp), parameter :: dGTINY = 1.0e-16
    
    enum, bind(C)
        enumerator :: PSRAD=1,PSPHI,PSKINE,PSALPHA,PSPSI
    endenum

    integer, parameter :: NVARPS = 4 !Number of PS variables

    !Particle distribution function type
    !Background density and temperature (n0,kT0) in [#/cc] and [keV]
    !Returns PDF in (keV*s)^(-3)

    abstract interface
        function PDFun_T(Model,n0,kT0,K,alpha) result(fD)
            import :: rp,chmpModel_T
            type(chmpModel_T), intent(in) :: Model
            real(rp), intent(in) :: n0,kT0,K,alpha
            real(rp) :: fD
        end function PDFun_T
    end interface


    !Equatorial phase space
    !4D: R,phi,K,alpha
    !Note: Volume element is still dx3dp3, using flux tube volume
    type PSEq_T
        real(rp) :: time=0 !Time
        integer :: Nr,Np,Nk,Na !Cells in each dimension
        real(rp) :: dimBds(4,2) !Min/Max of each dimension
        !Cell centers & interfaces in each dimension
        real(rp), dimension(:), allocatable :: rC,pC,kC,aC
        real(rp), dimension(:), allocatable :: rI,pI,kI,aI
        real(rp), dimension(:), allocatable :: rD,pD,kD,aD

        logical :: doShape !Use shape functions when depositing weight

        !Flux tube volume, dVb = |Nr,Np,Na|
        real(rp), allocatable :: dVb(:,:,:)
        logical, allocatable :: isClosed(:,:)
        real(rp) :: rCrit = 4.0 !Radius for field closure

        !Full volume element, dG = |Nr,Np,Nk,Na|
        !dGx ~ L0^3
        !dGp ~ keV^3/c^3
        !dG ~ (keV*s)^3
        real(rp), allocatable :: dG(:,:,:,:) 

        !Holder for field lines traced from equator
        type(fLine_T), allocatable :: bLns(:,:)

        !Equatorial MHD values (primitive)
        !Den (#/cc), Vxyz (cm/s), Pressure (nPa)
        real(rp), dimension(:,:,:), allocatable :: Qrp
        !Equatorial kT (keV)
        real(rp), dimension(:,:), allocatable :: Vreq,kTeq      
    end type PSEq_T

    !PSD test particle population
    type psdPop_T
        real(rp) :: kTScl=1.0 !Parameter to scale temperature (ie MHD->electron)
        
        !File info
        !Assuming form popid.(I0.6).h5part
        !With I0.6=ns,ne w/ 0 padding
        character(len=strLen) :: popid
        integer :: ns,ne !Start/stop blocks

        !TP info
        integer :: NumTP !Total number of test particles
        real(rp) :: T0 !First time info in data file
        real(rp) :: dtStp !Time between Steps in h5p input files
        !TP values, size (NumTP,NVARPS)
        real(rp), dimension(:,:), allocatable :: TPs
        logical, dimension(:), allocatable :: isWgt,isIn
        real(rp), dimension(:), allocatable :: wgt,tx

        real(rp) :: dTau = 0.0 !Injection timescale (code time)
        real(rp) :: dShell = 0.0 !Thickness of injection shell (L0)

        !Phase space density for this population (on a given PS-grid)
        !fPSD = |Nr,Np,Nk,Na|
        real(rp), dimension(:,:,:,:), allocatable :: fPSD
        integer , dimension(:,:,:,:), allocatable :: nPSD
    end type psdPop_T

end module psdtypes
