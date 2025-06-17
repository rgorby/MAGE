

module helpertypes
	use kdefs
	implicit none

    ! Type to hold info for dynamic load balancing
    type loadBal_T
        integer :: nL ! number of loads to balance
        integer :: nO ! number of outputs to split across loads
        real(rp), dimension(:), allocatable :: instantTimes
        real(rp), dimension(:), allocatable :: smoothLoads
        real(rp) :: hAlpha = 0.8 ! rate at which smoothed loads change with new data
        integer, dimension(:), allocatable :: balStartInd ! load balanced starting indices
    end type loadBal_T
    
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