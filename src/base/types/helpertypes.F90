

module helpertypes
	use kdefs
	implicit none
    
    ! Types to hold planet-specific information
	type planet_T
        character(len=strLen) :: name
        real(rp) :: rp_m  ! Planet radius [m]
        real(rp) :: ri_m  ! Ionosphere radius [m]
        real(rp) :: grav  ! Gravity [m/s2]
        real(rp) :: magMoment  ! [Gauss]
        real(rp) :: psiCorot  ! [kV] Should maybe instead store rotation rate, and let corotation potential be decided elsewhere
        logical :: doGrav  
    end type planet_T


end module helpertypes