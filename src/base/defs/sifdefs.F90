! kaimag definitions/constants

module sifdefs
	use kdefs

	implicit none

	!------
	! Enumerators
	!------

	! Grid settings
	enum, bind(C)
		enumerator :: G_UNISPH, G_SHGRID
	endenum

	! Topology
	enum, bind(C)
		enumerator :: SIFOPEN, SIFCLOSED
	endenum

	! Active/buffer/inactive cells
	enum, bind(C)
		enumerator :: SIFACTIVE, SIFBUFFER, SIFINACTIVE
	endenum

	!------
	! Defaults
	!------


end module sifdefs