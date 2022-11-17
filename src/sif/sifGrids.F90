

module sifgrids
	use sifdefs
	use siftypes


	implicit none

	contains



	subroutine allocGrid(Grid)
		type(sifGrid_T), intent(inout) :: Grid

		allocate(Grid%latlon(Grid%isg:Grid%ieg+1,Grid%jsg:Grid%jeg+1,2))
        Grid%latlon = 0.0

        allocate(Grid%latloncc(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,2))
        Grid%latloncc = 0.0

        allocate(Grid%llfc(Grid%isg:Grid%ieg,Grid%jsg:Grid%jeg,2,2))
        Grid%llfc = 0.0

        ! TODO: Lambda channel allocation

	end subroutine allocGrid

	subroutine sifGenUniSphGrid(Grid)
		type(sifGrid_T), intent(inout) :: Grid


		

		
		! Allocate arrays

	end subroutine


end module sifgrids