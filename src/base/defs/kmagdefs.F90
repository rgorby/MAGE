! kaimag definitions/constants

module kmagdefs
	use kdefs

	implicit none

	! Topology
	enum, bind(C)
		enumerator :: KMOPEN, KMCLOSED
	endenum

	! Active/buffer/inactive cells
	enum, bind(C)
		enumerator :: KMACTIVE, KMBUFFER, KMINACTIVE
	endenum

end module kmagdefs