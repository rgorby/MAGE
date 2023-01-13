

module sifgrids
    use ioh5
    use xml_input
    use shellgrid

    use sifdefs
    use siftypes
    use sifSpeciesHelper

    implicit none

    contains

!------
! Spatial grid stuff
!------


    subroutine sifGenUniSphGrid(Grid, iXML)
        type(sifGrid_T)  , intent(inout) :: Grid
        type(XML_Input_T), intent(in)    :: iXML

        real(rp), dimension(:), allocatable :: Theta
        real(rp), dimension(:), allocatable :: Phi
        real(rp) :: dTheta, dPhi, thetaL, thetaU
        integer :: Nt,Np,Ng,Ngn,Ngs,Nge,Ngw
        integer :: i


        call iXML%Set_Val(thetaL , "grid/ThetaL", 50.0)
            !! Lower lat boundary [deg]
        call iXML%Set_Val(thetaU , "grid/ThetaU", 80.0)
            !! Upper lat boundary [deg]
        call iXML%Set_Val(Nt, "grid/Nt", 30 )  ! 1 deg resolution
        call iXML%Set_Val(Np, "grid/Np", 360)  ! 1 deg resolution
        call iXML%Set_Val(Ng, "grid/Ng", 4  )  ! Number of ghosts, in every direction for now

        ! Turn degrees into radians
        thetaL = thetaL*deg2rad
        thetaU = thetaU*deg2rad

        ! Probably need to change
        Ngn = Ng
        Ngs = Ng
        Nge = Ng
        Ngw = Ng

        ! Allocate arrays
        allocate(Theta(Nt))
        allocate(Phi(Np))

        ! Create uniform grids
        dTheta = (thetaU-thetaL)/(Nt-1)
        dPhi = 2*PI/(Np-1)

        do i=1,Nt
            Theta(i) = thetaL + (i-1)*dTheta
        enddo

        do i=1,Np
            Phi(i) = (i-1)*dPhi
        enddo

        !write(*,*)"Orig theta:",Theta
        !write(*,*)"Orig phi:",Phi

        call GenShellGrid(Grid%shGrid,Theta,Phi,Ngn,Ngs,Nge,Ngw)

    end subroutine sifGenUniSphGrid

    subroutine sifGenGridFromShGrid(Grid, shGrid)
        type(sifGrid_T)  , intent(inout) :: Grid
        type(ShellGrid_T), intent(inout) :: shGrid

        !!TODO
    end subroutine sifGenGridFromShGrid

    subroutine finalizeLLGrid(Grid)
        !! Use a fully-created shell grid to allocate and populate the rest of the grid parameters
        type(sifGrid_T), intent(inout) :: Grid


        associate(shGr=>Grid%shGrid)
            ! First allocate remaining arrays
            !allocate(Grid%llfc(shGr%isg:shGr%ieg, shGr%isg:shGr%ieg, 2, 2))  ! Face-centered
            !allocate(Grid%llfc(shGr%isg:shGr%ieg, shGr%isg:shGr%ieg, 2   ))  ! Cell-centered
            allocate(Grid%iBnd(shGr%is :shGr%ie ))  ! i/lat boundary for valid domain

            Grid%iBnd = 0

        end associate

    end subroutine finalizeLLGrid


!------
! Lambda grid stuff
!------

!! Maybe should leave just spatial grid stuff in sifGrids and move lambda stuff to a lambdaUtils

    

    subroutine initLambdaGrid(Model, Grid, configfname)
        type(sifModel_T)  , intent(inout) :: Model
        type(sifGrid_T), intent(inout) :: Grid
        character(len=strLen), intent(in) :: configfname

        ! First read in speies information from config file
        call populateSpeciesFromConfig(Model, Grid, configfname)

        ! Now prepare our alamc grid, tell each Species what its range is in alamc dimension
        !! TODO
    end subroutine initLambdaGrid

end module sifgrids