! Initialize kaimag
module sifstarter
    use sifdefs
    use siftypes
    use sifgrids
    use xml_input
    use planethelper

    implicit none






    contains

!------
! Main Initialization Routines
!------

    ! Sets up Model, but Grid and State must be set up separately
    ! Its up to a higher being to determine how we get our grid
    ! After we have a grid, we can initialize our first state
    subroutine sifInitModel(Model, Grid, planet, iXML)
        type(sifModel_T), intent(inout) :: Model
        type(sifGrid_T) , intent(inout) :: Grid
        type(planet_T), intent(in)  :: planet
        type(XML_Input_T), intent(in) :: iXML
 
        character(len=strLen) :: xmlStr
        type(XML_Input_T) :: xmlInp
        write(*,*) "sifInitModel is starting"



        !! NOT SET HERE:
        ! nG, nB, t0, tFin, dt, fixedTimestep

        ! Set some settings
        call iXML%Set_Val(Model%isMPI, "mpi/isMPI",.false.)
        if (Model%isMPI) then
            write(*,*) "MPI not implemented for SIF yet, dying."
            stop
        endif

        call iXML%Set_Val(Model%isRestart, "Kaiju/Gamera/restart/doRes",.false.)
        call iXML%Set_Val(Model%isLoud, "debug/isLoud",.false.)
        call iXML%Set_Val(Model%writeGhosts, "debug/writeGhosts",.false.)

        ! Plasmasphere settings
        call iXML%Set_Val(Model%doPlasmasphere, "plasmasphere/doPsphere",.false.)

        ! Lambda channel settings
        call iXML%Set_Val(Model%doDynamicLambdaRanges, "lambdas/dynamicRanges",.false.)


        ! Set planet params
        !! This should only be kept for as long as planet_T doesn't contain pointers
        !! In this current case, there should be a full copy to our own planet params
        Model%planet = planet

        ! Generate grid
        call sifInitGrid(Grid, iXML)

    end subroutine sifInitModel

    

    subroutine sifInitGrid(Grid, iXML)
        !type(sifModel_T), intent(inout) :: Model
        type(sifGrid_T), intent(inout) :: Grid


        type(XML_Input_T) :: iXML


        ! Set grid params
        ! Grid settings
        call iXML%Set_Val(Grid%gType, "grid/gType",G_UNISPH)
        call iXML%Set_Val(Grid%Ng , "grid/Ng", 4  )  ! Ghost cells
        call iXML%Set_Val(Grid%Nip, "grid/Np", 90 )  ! Phi
        !! Ignoring lambda grid for now, add in later
        call iXML%Set_Val(Grid%Njp, "grid/Nt", 360)  ! Theta
        call iXML%Set_Val(Grid%Nkp, "grid/Nl", 140)  ! Lambda
        call iXML%Set_Val(Grid%latBndL, "grid/LatL", 50.0) ! Lat boundary, lower
        call iXML%Set_Val(Grid%latBndU, "grid/LatU", 80.0) ! Lat boundary, upper

        ! Finalize grid definition
        Grid%Ni = Grid%Nip + 2*Grid%Ng
        Grid%Nj = Grid%Njp + 2*Grid%Ng

        Grid%is = 1; Grid%ie = Grid%Nip
        Grid%js = 1; Grid%je = Grid%Njp

        Grid%isg = Grid%is-Grid%Ng
        Grid%ieg = Grid%ie+Grid%Ng
        Grid%jsg = Grid%js-Grid%Ng
        Grid%jeg = Grid%je+Grid%Ng

        call allocGrid(Grid)

        ! Now fill in out coordinates depending on specification
        select case(Grid%gType)
            case(G_UNISPH)
                ! pass
            case(G_VOLTRON)
                ! pass
            case DEFAULT
                write(*,*) "SIF Received invalid grid definition: ",Grid%gType
                write(*,*) " Dying."
                stop
        end select

    end subroutine sifInitGrid





!------
! Defaults
!------

    ! This will make a simple Maxwellian distribution in a dipole field
    subroutine sifInitStateDefault(State, D, T, r)
        type(sifState_T), intent(inout) :: State
        real(rp), intent(in) :: D, T, r  ! density, temperature, r value
        !! TODO
    end subroutine sifInitStateDefault


end module sifstarter