    INTEGER (iprec), PARAMETER :: &
       isize = 200, &
       jsize = 101, &
       jwrap =   3, &
       ksize  = 090,&
       nptmax = 50000, &
       iesize =   2, &
       ncoeff =   5, &
       kcsize = ksize
    LOGICAL, PARAMETER :: asci_flag = .FALSE.
    LOGICAL, PARAMETER :: isGAMRCM = .TRUE. !Whether running coupled to Gamera
    LOGICAL, PARAMETER :: doQuietRCM = .TRUE.
    !If you came here looking for doClaw95, use the xml option
    !rcm/clawpack/doKaiClaw="T"