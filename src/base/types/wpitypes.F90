!Main data types for WPI calculations

module wpitypes
    use chmpdefs
    use tptypes
    !use ebtypes

    implicit none

    !wave specific information
    type wave_T
        character(len=strLen) :: mode
        real(rp) :: s      ! notation from Summers 2005, s holds wave mode (L = -1, R = 1)
        real(rp) :: lam    ! lam holds particle info(-1 = e-, me/mp = p+) 
        real(rp) :: eMin=0 !minimum energy needed for particle to resonant
        integer :: id
    end type wave_T

    !wave model information
    !currently assumes gaussian wave distribution
    !will need to be updated to hold these variables as a function of location/allocatable arrays 
    type wModel_T
        character(len=strLen) :: model
        logical :: isWave=.false. !boolean if wave is present at particle location
        real(rp) :: xm,Dx !This is mean frequency and half width of gaussian wave band normalized by gyro frequency
        real(rp) :: B1 !wave amplitude

    end type wModel_T

end module wpitypes

