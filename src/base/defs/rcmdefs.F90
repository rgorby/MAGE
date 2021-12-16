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
    INTEGER, parameter :: RCMELECTRON = 1
    INTEGER, parameter :: RCMPROTON   = 2
    INTEGER, parameter :: RCMOXYGEN   = 3
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
    REAL(krp) :: DenPP0 = 5.0 !Defining plasmasphere density cutoff, [#/cc]
    REAL(krp) :: PSPHKT = 1.0e-3 !Characteristic temperature for plasmasphere [keV]
    REAL(krp), PARAMETER :: machine_tiny = 1.0e-32

    enum, bind(C)
      enumerator :: ELOSS_FDG=1,ELOSS_SS,ELOSS_C05,ELOSS_C19,ELOSS_DW !Choice of electron loss model
    end enum

end module rcmdefs
