! Initialize kaimag
module sifStarter
	use sifdefs
	use siftypes
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
	subroutine sifInitModel(Model, planet, iXML, optFilename)
		type(sifModel_T), intent(inout) :: Model
		type(planet_T), intent(in)  :: planet
		type(XML_Input_T), intent(out), optional :: iXML
        character(len=*) , intent(in) , optional :: optFilename
 
        character(len=strLen) :: xmlStr
        type(XML_Input_T) :: xmlInp
        write(*,*) "kmInitModel is starting"

        !Create input XML object
        if (present(optFilename)) then
            xmlStr = optFilename
        else
            call getIDeckStr(xmlStr)  ! Returns string based on what's listed in command-line arguments
        endif
        
        xmlInp = New_XML_Input(trim(xmlStr),"Kaiju/SIF",.true.)
        ! In Chimp/starter.F90/goApe: Why send it back if optional arg present


        !! NOT SET HERE:
        ! nG, nB, t0, tFin, dt, fixedTimestep

        ! Set some settings
        call xmlInp%Set_Val(Model%isRestart, "/Kaiju/Gamera/restart/doRes",.false.)
        call xmlInp%Set_Val(Model%isLoud, "debug/isLoud",.false.)
        call xmlInp%Set_Val(Model%writeGhosts, "debug/writeGhosts",.false.)

        ! Plasmasphere settings
        call xmlInp%Set_Val(Model%doPlasmasphere, "plasmasphere/doPsphere",.false.)

        ! Lambda channel settings
        call xmlInp%Set_Val(Model%doDynamicLambdaRanges, "lambdas/dynamicRanges",.false.)


        ! Set planet params
        !! This should only be kept for as long as planet_T doesn't contain pointers
        !! In this current case, there should be a full copy to our own planet params
        Model%planet = planet


	end subroutine

	subroutine sifGenLLGrid(Grid, latBndL, latBndU, Ni, Nj, Ng)
		type(sifGrid_T), intent(inout) :: Grid
		real(rp), intent(in) :: latBndL, latBndU
		integer, intent(in) :: Ni,Nj,Ng
		!! TODO
	end subroutine





!------
! Defaults
!------

	! This will make a simple Maxwellian distribution in a dipole field
	subroutine sifInitStateDefault(State, D, T, r)
		type(sifState_T), intent(inout) :: State
		real(rp), intent(in) :: D, T, r  ! density, temperature, r value
		!! TODO
	end subroutine


end module sifStarter