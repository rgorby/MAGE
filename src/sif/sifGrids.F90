

module sifgrids
    use ioh5
    use xml_input
    use sifdefs
    use siftypes
    use shellgrid


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

    subroutine populateSpeciesFromConfig(Grid, configfname, nSpc)
        type(sifGrid_T), intent(inout) :: Grid
        character(len=strLen), intent(in) :: configfname
        integer, intent(in) :: nSpc

        integer(HID_T) :: h5fId
        integer :: herr, i, flav
        character(len=strLen) :: gStr
        integer(HID_T) :: gId
        logical :: gExist, isEnd
        type(IOVAR_T), dimension(2) :: IOVars ! Just grabbing lami and we're done
        type(IOVAR_T) :: IOV


        call CheckFileOrDie(configfname,"Unable to open file")
        call h5open_f(herr) !Setup Fortran interface
        ! Open file
        call h5fopen_f(trim(configfname), H5F_ACC_RDONLY_F, h5fId, herr)
        
        ! Make sure Species group is there
        gStr = "Species"
        call h5lexists_f(h5fId, gStr, gExist, herr)
        if (.not. gExist) then
            write(*,*) "This config file not structured for SIF, use genSIF.py. Good day."
            stop
        endif
        
        ! If still here, start populating Grid%spc
        allocate(Grid%spc(nSpc))

        associate(spc=>Grid%spc)
        
        ! TODO: Think through edge cases that will cause errors

        isEnd = .false.
        do i=1,nSpc
            write(gStr, '(A,I0)') "Species/",i-1  ! Species indexing in config starts at 0
            call h5lexists_f(h5fId,gStr,gExist,herr)
            if (.not. gExist) then
                write(*,*) "ERROR in sifGrids.F90:populateSpeciesFromConfig"
                write(*,'(A,I0,A,I0)') "  Expected ",nSpc," species but only found ",i-1
            else
                write(*,*)"Found spc group ",trim(gStr)
                
                call h5gopen_f(h5fId,trim(gStr),gId,herr)
                ! TODO: Figure out how to get Name string from attrs
                spc(i)%N = readIntHDF(gId, "N")
                spc(i)%flav = readIntHDF(gId, "flav")
                spc(i)%fudge = readRealHDF(gId, "fudge")
                
                write(*,*)spc(i)%N
                write(*,*)spc(i)%flav
                write(*,*)spc(i)%fudge

                ! Now get our alami values
                allocate(spc(i)%alami(spc(i)%N+1))
                
                ! Do this manually, can't use IOVARS/ReadVars cause we already have the group open
                IOV%toRead = .true.
                IOV%idStr = "alami"
                IOV%vType = IONULL
                IOV%scale = 1.0
                call ReadHDFVar(IOV, gId)
                spc(i)%alami = IOV%data
                write(*,*)spc(i)%alami

                call h5gclose_f(gId,herr)
                
            endif
        enddo

        ! TODO: Put in a check for if there are more species in the file than what we read
        !       Warn user if so, list the ones we are using
        !       Maybe stop if certain override XML parameter isn't set
        ! TODO: Once we establish 0=psph, 1=ele, 2=proton, etc, use this info to ignore psphere in config if doPsphere=F

        end associate

        call h5fclose_f(h5fId,herr)
        call h5close_f(herr) !Close intereface

    end subroutine populateSpeciesFromConfig

    subroutine initLambdaGrid(Grid, configfname, nSpc)
        type(sifGrid_T)  , intent(inout) :: Grid
        character(len=strLen), intent(in) :: configfname
        integer, intent(in) :: nSpc

        ! First read in speies information from config file
        call populateSpeciesFromConfig(Grid, configfname, nSpc)

        ! Now prepare our alamc grid, tell each Species what its range is in alamc dimension

    end subroutine initLambdaGrid

end module sifgrids