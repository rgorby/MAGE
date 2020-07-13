!Main data types for WPI calculations

module wpitypes
    use chmpdefs
    use tptypes

    implicit none

    integer, parameter :: NROOTS = 4 !Number of possible resonant roots

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

    !Group Velocity function type
    !Returns Vg for a given resonant wave in normalized units
    abstract interface
        function VgFun_T(wave,astar,x,y) result(Vg)
            import :: rp
            type(wave_T), intent(in) :: wave
            real(rp), intent(in) :: astar,x,y
            real(rp) :: Vg
        end function VgFun_T
    end interface

    !Resonant root function type
    !Returns resonant root in unitless variables (divided by gyrofrequency)
    abstract interface
        subroutine ResFun_T(Model,wave,prt,astar,x,y) 
            import :: rp,prt_T,chmpModel_T
            type(chmpModel_T), intent(in) :: Model
            type(prt_T), intent(in) :: prt
            type(wave_T), intent(in) :: wave
            real(rp), intent(in) :: astar
            real(rp), intent(inout) :: x,y
        end subroutine ResFun_T
    end interface

    !wave frequency function type
    !Returns wave spectral density for a given wave number
    abstract interface
        function wsFun_T(wModel,x) result(Ws)
            import :: rp
            type(wModel_T), intent(in) :: wave
            real(rp), intent(in) :: x
            real(rp) :: Ws
        end function wsFun_T
    end interface

end module wpitypes

