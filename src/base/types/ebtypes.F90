!Main CHIMP field/grid data types

module ebtypes
    use chmpdefs
    
    implicit none

    !Primitive variables
    integer, parameter :: NVARMHD = 5

    !Data necessary to update fields, ie time->field data file mapping
    !File: bStr + gStr (base+group)
    type ebTab_T
        integer :: N !Number of time slices
        character(len=strLen) :: bStr
        character(len=strLen), dimension(:), allocatable :: gStrs
        real(rp), dimension(:), allocatable :: times
        real(rp), dimension(:), allocatable :: MJDs

        !Information for decomposed data
        logical :: isMPI = .false.,hasMJD=.false.
        integer :: Ri,Rj,Rk
        integer :: dNi,dNj,dNk
    end type ebTab_T

    !Holds single slice of field data, attached to grid data in ebState
    type ebField_T
        real(rp), dimension(:,:,:,:), allocatable :: dB,E !Fields
        real(rp), dimension(:,:,:,:), allocatable :: W !Primitive MHD variables
        real(rp) :: time !Time in code units for this slice
    end type ebField_T

    !Grid data for eb fields
    type ebGrid_T
        integer :: Ni,Nj,Nk !Grid dimensions (w/ ghosts)
        integer :: Nip,Njp,Nkp !Active cells (w/o ghosts)
        integer :: is,ie,js,je,ks,ke !Grid indices (w/o ghosts)
        integer :: isg,ieg,jsg,jeg,ksg,keg !Grid indices (w ghosts)
        integer :: GrID=CARTGRID !Grid ID (see chmpdefs)
        !4D cell corners, cell centers
        !(isg:ieg,jsg:jeg,ksg:keg,1:NDIM)
        real(rp), dimension(:,:,:,:), allocatable :: xyz,xyzcc,B0cc
        !Cell-centered Jacobians.  Txi: (xyz),ijk / Tix: (ijk),xyz
        !Ie, Txi(XDIR,JDIR) = dx/dj
        !(isg:ieg,jsg:jeg,ksg:keg,1:NDIM,1:NDIM)
        real(rp), dimension(:,:,:,:,:), allocatable :: Txi,Tix
    end type ebGrid_T

    !Main field data structure
    type ebState_T
        type(ebGrid_T)  :: ebGr
        type(ebField_T) :: eb1,eb2 !Two slices of field data
        type(ebTab_T)   :: ebTab
        logical :: doStatic = .false.
    end type ebState_T
    
    !Type used for tracing field lines
    type ebTrcApp_T
        type(chmpModel_T) :: ebModel
        type(ebState_T)   :: ebState
    end type ebTrcApp_T

    !Data holder for doing field line tracing at point
    type ebTrc_T
        real(rp) :: OCb !Topology
        real(rp) :: dvB,bD,bP,bS !Flux-tube volume, averaged density/pressure, integrated entropy
        real(rp), dimension(NDIM) :: MagEQ, xEPm,xEPp !xyz of equator/ -/+ field endpoints
        real(rp) :: bMin !Minimum B (@ equator)
    end type ebTrc_T

    !Fields and derivatives necessary for GC update
    type gcFields_T
        real(rp), dimension(NDIM) :: DotE,DotB
        real(rp), dimension(NDIM,NDIM) :: JacE,JacB

    end type gcFields_T

    integer, parameter :: MaxFL = 5000 !Reduced for multi-threading speed
    integer, parameter :: NumVFL = NVARMHD !Number of field line variables (other than |B|)

    !Streamline variable
    type lnVar_T
        character(len=strLen) :: idStr !Variable name
        real(rp), allocatable, dimension(:) :: V !Same spacing as xyz in main streamline structure (-Nm:Np)
        real(rp) :: V0 !Value @ "equator"
    end type lnVar_T

    !Individual streamline
    type fLine_T
        integer :: Nm,Np
        real(rp), dimension(NDIM) :: x0 !Seed point
        real(rp), allocatable, dimension(:,:) :: xyz

        !xyz(-Nm:Np,:), w/ xyz(0,:) = seed point
        type(lnVar_T), dimension(0:NumVFL) :: lnVars

        !Localization data, ie ijk of each node of field line (not set yet)
        integer, allocatable, dimension(:,:) :: ijk

        !Whether this is a degenerate line (seed point not in domain)
        logical :: isGood = .false.
    end type fLine_T
        
    contains


    !Allocates ebField_T given ebGr
    subroutine allocEB(Model,ebGr,ebF)
        type(chmpModel_T), intent(inout) :: Model
        type(ebGrid_T), intent(in)  :: ebGr
        type(ebField_T), intent(inout) :: ebF

        ebF%time = 0.0
        if (.not. allocated(ebF%dB)) then
            allocate(ebF%dB(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        endif

        if (.not. allocated(ebF%E)) then
            allocate(ebF%E(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NDIM))
        endif

        ebF%dB = 0.0
        ebF%E  = 0.0
        if (Model%doMHD .and. (.not. allocated(ebF%W))) then
            allocate(ebF%W(ebGr%isg:ebGr%ieg,ebGr%jsg:ebGr%jeg,ebGr%ksg:ebGr%keg,NVARMHD))
            ebF%W = 0.0
        endif
        
    end subroutine allocEB
    
end module ebtypes
