    INTEGER (iprec), PARAMETER :: &
       isize = 200, &
       jsize = 101, &
       jwrap =   3, &
       ksize  = 090,  kcsize = 090, &
       nptmax = 50000, &
       iesize =   2, &
       ncoeff =   5
    LOGICAL, PARAMETER :: asci_flag = .FALSE.
    LOGICAL, PARAMETER :: isGAMRCM = .TRUE. !Whether running coupled to Gamera
    LOGICAL, PARAMETER :: doQuietRCM = .TRUE.
    LOGICAL, PARAMETER :: doClaw95 = .TRUE.
