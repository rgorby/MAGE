! Initialize kaimag
module sifstarter
    use sifdefs
    use siftypes
    use sifgrids
    use shellgrid
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
    subroutine sifInitModel(Model, planet, iXML)
        type(sifModel_T) , intent(inout) :: Model
        type(planet_T)   , intent(in)    :: planet
        type(XML_Input_T), intent(in)    :: iXML
 
        write(*,*) "sifInitModel is starting"



        !! NOT SET HERE:
        ! nG, nB, t0, tFin, dt, fixedTimestep

        ! Set some settings
        call iXML%Set_Val(Model%configFName, "sim/config","sifconfig.h5")

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
        ! Determine number of species. First set default, then read from xml to overwrite if present
        if (Model%doPlasmasphere) then
            Model%nSpc = 3
        else
            Model%nSpc = 2
        endif
        call iXML%Set_Val(Model%nSpc, "sim/nSpc",Model%nSpc)


        ! Lambda channel settings
        call iXML%Set_Val(Model%doDynamicLambdaRanges, "lambdas/dynamicRanges",.false.)


        ! Set planet params
        !! This should only be kept for as long as planet_T doesn't contain pointers
        !! In this current case, there should be a full copy to our own planet params
        Model%planet = planet

        ! Set up timing
        !!TODO
        ! If we are running stand-alone, look for timing info inside SIF XML block
        ! If voltron is present, look for timing information there

        

    end subroutine sifInitModel

    

    subroutine sifInitGrid(Model, Grid, iXML, shGridO)
        type(sifModel_T)  , intent(inout) :: Model
        type(sifGrid_T)  , intent(inout) :: Grid
        type(XML_Input_T), intent(in)   :: iXML
        type(ShellGrid_T), optional, intent(in)    :: shGridO

        character(len=strLen) :: tmpStr
        type(ShellGrid_T) :: shGrid

        ! Set grid params
        ! Grid settings
        call iXML%Set_Val(tmpStr, "grid/gType","UNISPH")


        ! Fill out Grid object depending on chosen method
        select case(tmpStr)
            case("UNISPH")
                Grid%gType = G_UNISPH
                ! Generate our own grid from scratch
                call sifGenUniSphGrid(Grid, iXML)
            case("SHGRID")
                Grid%gType = G_SHGRID
                ! Then we should be receiving a predefined ShellGrid that Voltron has set up
                if(present(shGridO)) then
                    shGrid = shGridO
                    call sifGenGridFromShGrid(Grid, shGrid)
                else
                    write(*,*) "SIF expecting a ShellGrid_T but didn't receive one. Dying."
                endif
            case DEFAULT
                write(*,*) "SIF Received invalid grid definition: ",Grid%gType
                write(*,*) " Dying."
                stop
        end select

        ! Finalize the spatial part of the grid
        call finalizeLLGrid(Grid)

        ! Now handle lambda grid
        !!!TODO
        call initLambdaGrid(Grid, Model%configFName, Model%nSpc)

    end subroutine sifInitGrid





!------
! Defaults
!------

    ! This will make a simple Maxwellian distribution in a dipole field
    subroutine sifInitStateDefault(Model, Grid, State, D, T, r)
        type(sifModel_T), intent(inout) :: Model
        type(sifGrid_T) , intent(inout) :: Grid
        type(sifState_T), intent(inout) :: State
        real(rp), intent(in) :: D, T, r  ! density, temperature, r value
        !! TODO
    end subroutine sifInitStateDefault


end module sifstarter