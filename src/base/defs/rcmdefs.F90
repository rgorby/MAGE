! Constants for RCM

module rcmdefs

    ! preventing ip and rp from polluting rcm due to use association
    use kdefs, ONLY: kip => ip, krp => rp

    implicit none

    private kip,krp

    INTEGER, parameter :: RCMIOVARS = 50

    INTEGER, parameter :: ICONWRITERESTART = 31337
    INTEGER, parameter :: ICONWRITEOUTPUT  = ICONWRITERESTART + 1
    INTEGER, parameter :: ICONRESTART      = ICONWRITERESTART - 1
    INTEGER, parameter :: RCMCOLD      = 0
    INTEGER, parameter :: RCMELECTRON  = 1
    INTEGER, parameter :: RCMPROTON    = 2
    INTEGER, parameter :: RCMOXYGEN    = 3
    INTEGER, parameter :: RCMVIBRANIUM = 4678371489657242888
    INTEGER, parameter :: RCMNUMFLAV = 2 !Number of RCM flavors
    INTEGER, PARAMETER :: isize = RCMSIZEI !RCM grid size in colatitude
    INTEGER, PARAMETER :: jsize = RCMSIZEJ !RCM grid size in longitude
    INTEGER, PARAMETER :: ksize = RCMSIZEK !RCM grid size in lambda
    INTEGER, PARAMETER :: jwrap = RCMWRAPJ !RCM wrapper cells in j
    INTEGER, PARAMETER :: nptmax = 50000
    INTEGER, PARAMETER :: iesize = 2 !Number of species flavors
    INTEGER, PARAMETER :: ncoeff = 5
    INTEGER, PARAMETER :: kcsize = ksize !No idea why
    LOGICAL, PARAMETER :: asci_flag = .FALSE.
    LOGICAL, PARAMETER :: isGAMRCM = .TRUE. !Whether running coupled to Gamera
    LOGICAL, PARAMETER :: doQuietRCM = .TRUE.
    LOGICAL, PARAMETER :: doDiskWrite = .FALSE.
    integer(kip), parameter :: RCMTOPCLOSED=-1,RCMTOPOPEN=+1,RCMTOPNULL=0
    REAL(krp) :: DenPP0 = 0.0 !Defining plasmasphere density cutoff, [#/cc]
    REAL(krp), PARAMETER :: machine_tiny = 1.0e-32
    REAL(krp), PARAMETER :: tiote_RCM = 4.0

    enum, bind(C)
      enumerator :: ELOSS_FDG=1,ELOSS_SS,ELOSS_C05,ELOSS_C19,ELOSS_DW !Choice of electron loss model
    end enum

    REAL(krp), PARAMETER :: bMin_C_DEF  = 1.0 ![nT], default min allowable field strength
    !For min-B value see Ohtani+Motoba 17
    REAL(krp), PARAMETER :: wImag_C_DEF = 0.10 !Default min allowable RCM "weight" (see rcmimag)

    !Dumb clawpack hard-coded values
    REAL(krp), PARAMETER :: CLAW_MAXCFL = 0.95
    REAL(krp), PARAMETER :: CLAW_REGCFL = 0.80
    
    !Standard config (CTU + Superbee)
    INTEGER  , PARAMETER :: ICLAW_TRANSORDER = 2 !2 is standard
    INTEGER  , PARAMETER :: ICLAW_LIMITER   = +2 !Superbee
    
    !Toy config (Dim split + MC)
    !INTEGER  , PARAMETER :: ICLAW_TRANSORDER = -1 !2 is standard, maybe try -1?
    !INTEGER  , PARAMETER :: ICLAW_LIMITER = +4 !MC Limiter

    !Enumerators for rcm boundary shape
    enum, bind(C)
        enumerator :: RCMDOMELLIPSE=1,RCMDOMCONTOUR,RCMDOMNONE
    endenum     
end module rcmdefs
