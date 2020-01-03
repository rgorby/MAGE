    INTEGER (iprec), PARAMETER :: &
       isize = 200, &
       jsize = 101, &
       jwrap =   3, &
        ksize  = 090,  kcsize = 090, &
       nptmax = 50000, &
       iesize =   2, &
       ncoeff =   5
    LOGICAL :: asci_flag = .FALSE.
    LOGICAL :: isGAMRCM = .TRUE. !Whether running coupled to Gamera
