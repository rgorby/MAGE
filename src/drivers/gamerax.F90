!Driver for Gamera (uncoupled)

program gamerax
    use gamapp
    use drivertypes
    use usergamic

    implicit none

    type(AppDriver_T) :: appDriver
    type(GameraApp_T), allocatable :: gApp

    ! create instance of gamera app, and perform any needed configuration
    allocate(gApp)
    gApp%gOptions%userInitFunc => initUser

    ! use the single-app helper function to move the gamera app into the driver structure
    ! note that after this point, gamApp will no longer be allocated
    allocate(appDriver%appPointers(1))
    call move_alloc(gApp, appDriver%appPointers(1)%p)

    ! run the gamera app
    call appDriver%RunApps()

end program gamerax

