! kaimag definitions/constants

module sifdefs
	use kdefs

	implicit none

	! Topology
	enum, bind(C)
		enumerator :: SIFOPEN, SIFCLOSED
	endenum

	! Active/buffer/inactive cells
	enum, bind(C)
		enumerator :: SIFACTIVE, SIFBUFFER, SIFINACTIVE
	endenum

end module sifdefs