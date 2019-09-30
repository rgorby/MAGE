!Various simple routines for managing time/dates
module dates
    use kdefs

    implicit none

    real(rp), parameter, private :: mjdScl = 1.0/(60.0*60.0*24.0)

    contains

    !Return MJD given t as seconds elapsed from mjd0
    function T2MJD(t,mjd0) result(mjd)
        real(rp), intent(in) :: t,mjd0
        real(rp) :: mjd

        mjd = mjd0 + t*mjdScl
    end function T2MJD


end module dates
