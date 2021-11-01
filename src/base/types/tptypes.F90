!Main CHIMP particle data types

module tptypes
    use chmpdefs

    implicit none
!Constant sizes
    integer, parameter :: NVARTP = 7 !Number of TP dynamical variables

!Numbered accessors
    !Full orbit (FO) dynamical variables
    enum, bind(C) 
        enumerator :: XPOS=1, YPOS, ZPOS, PXFO,PYFO,PZFO, PSITP
    endenum

    !Guiding center (GC) dynamical variables
    !Using (FO) terms for spatial position of GC
    enum, bind(C) 
        enumerator :: P11GC=ZPOS+1,MUGC,GAMGC
    endenum 

    !Equatorial variables
    enum, bind(C)
        !(EQX,EQY,EQTIME) = position & time of last EQ-X (equatorial crossing)
        !EQKEV,EQALP = Energy [kev] & pitch angle [radians] @ EQ-X
        !EQKEB = Energy [kev] in ExB frame @ EQ-X
        enumerator :: EQX=1,EQY,EQTIME,EQKEV,EQALP,EQKEB
    endenum
    integer, parameter :: NVAREQ = 6

!Individual particle data
    !Goal is to keep this as slim as possible
    !q is the evolved variables for each particle
    !Full orbit (FO): x,y,z,px,py,pz
    !Guiding center (GC): x,y,z (of GC in ExB frame), p11,gamma,psi
    type prt_T
        real(rp) :: Q(NVARTP),Qeq(NVAREQ)
        logical :: isGC=.false., isIn=.false.,isInit=.false.
        real(rp) :: T0p=0,Tfp !This particles birthday/death
        integer :: id
        integer :: Ngc=0,Nfo=0 !Number of steps
        integer :: ijk0(NDIM) !Last known location in eb grid
        real(rp) :: ddt !Size of next substep
        real(rp) :: alpha=0.0
        real(rp) :: dAwpi=0,dKwpi=0 !change in pitch angle and energy due to wpi
        real(rp) :: xj=0.0,yj=0.0 !normalized resonant wave frequency and wave # if wpi present
        integer :: OCb=0 !Topology, [0,1,2] = Open/Clopen/Closed
        real(rp) :: Nwpi=0 !to keep track of number of wpi undergone
    end type prt_T

!Particle collection data
    !Assuming all particles in one tpState are same species
    !Extra derived quantities/metrics for particles stored here
    type tpState_T
        integer :: Np !Number of particles
        integer :: NpT !Number of currently active particles
        type(prt_T), dimension(:), allocatable :: TPs
    end type tpState_T

end module tptypes