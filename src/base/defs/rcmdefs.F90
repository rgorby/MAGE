! Constants for RCM

module rcmdefs
    use kdefs, ONLY: ip,rp
    implicit none

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
    integer(ip), parameter :: RCMTOPCLOSED=-1,RCMTOPOPEN=+1,RCMTOPNULL=0
    real(rp) :: DenPP0 = 1.0e-6 !Defining plasmasphere density cutoff, [#/cc]
end module rcmdefs