

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

    subroutine initLambdaGrid(Grid, configfname, nSpc)
        type(sifGrid_T)  , intent(inout) :: Grid
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
        i = 0
        do while (.not. isEnd)
            write(gStr, '(A,I0)') "Species/",i
            call h5lexists_f(h5fId,gStr,gExist,herr)
            if (.not. gExist) then
                isEnd = .true.
            else
                write(*,*)"Found spc group ",trim(gStr)
                
                call h5gopen_f(h5fId,trim(gStr),gId,herr)
                ! TODO: Figure out how to get Name string from attrs
                spc(i+1)%N = readIntHDF(gId, "N")
                spc(i+1)%flav = readIntHDF(gId, "flav")
                spc(i+1)%fudge = readRealHDF(gId, "fudge")
                
                !write(*,*)spc(i+1)%N
                !write(*,*)spc(i+1)%flav
                !write(*,*)spc(i+1)%fudge

                ! Now get our alami values
                allocate(spc(i+1)%alami(spc(i+1)%N+1))
                
                ! Do this manually, can't use IOVARS/ReadVars cause we already have the group open
                ! TODO: Abstract this away a little bit
                IOV%toRead = .true.
                IOV%idStr = "alami"
                IOV%vType = IONULL
                IOV%scale = 1.0
                call ReadHDFVar(IOV, gId)
                spc(i+1)%alami = IOV%data
                write(*,*)spc(i+1)%alami

                call h5gclose_f(gId,herr)

                !call ClearIO(IOVars) !Reset IO chain
                !call AddInVar(IOVars,"alami")
                !write(*,*)"a"
                !call ReadVars(IOVars,.false.,configfname, gStr)
                !write(*,*)"b"
                !write(*,*)IOVars(1)%idStr
                !write(*,*)IOVars(1)%data
                !call IOArray1DFill(IOVars, "alami", spc(i+1)%alami)
                !write(*,*)"c"
!
                !write(*,*)spc(i+1)%alami

                
            endif
            i = i+1
        enddo


        
        end associate

        ! Read lambda interfaces from file
        !call ClearIO(IOVars) !Reset IO chain
        !call AddInVar(IOVars,"alamc")
        !call ReadVars(IOVars,.false.,configfname)

        !! TODO: sifconfig.h5 should have lami(Nk+1), flav(Nk), etc...
        !! We get Nk as IOVars(lami)%N-1, then allocate and initialize accordingly
        !Nk = IOVars(FindIO(IOVars, "alamc"))%N

        !allocate(Grid%lamc(Nk))

        !call IOArray1DFill(IOVars, "alamc", Grid%lamc)

        !write(*,*)Grid%lamc

        call h5fclose_f(h5fId,herr)
        call h5close_f(herr) !Close intereface

    end subroutine initLambdaGrid

end module sifgrids