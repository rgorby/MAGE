module timeHelpers

    use kdefs

    contains

    subroutine timeStrFmt(T,tStr, tSclO)
        !! Simple function to format time in either units of seconds or minutes
        real(rp), intent(in) :: T
            !! Time value. Assumed in seconds, or can be in code units if tSclO is provided
        character(len=strLen), intent(out) :: tStr
            !! Returned formatted time string
        real(rp), optional, intent(in) :: tSclO
            !! (Default 1.0) time to scale T by to get units of seconds

        real(rp) :: tScl

        if (present(tSclO)) then
            tScl = tSclO
        else
            tScl = 1.0_rp
        endif

        if (abs(T*tScl)>60.0) then
            write(tStr,'(f9.3,a)' ) T*tScl/60.0, ' [min]'
        else
            write(tStr,'(es9.2,a)') T*tScl     , ' [sec]'
        endif

    end subroutine timeStrFmt


end module timeHelpers