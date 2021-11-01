module ebtabutils
    use chmpdefs
    use chmpunits
    use ebtypes

    implicit none

    contains

    !Finds bounding slices from ebTab file
    !NOTE: findloc isn't supported by most gfortran versions, so this is a lazy workaround
    !When gfortran gets its act together, just comment out last few lines using findloc
    subroutine findSlc(ebTab,t,i1,i2)
        type(ebTab_T), intent(in) :: ebTab
        real(rp), intent(in) :: t
        integer, intent(out) :: i1,i2

        integer :: n

        !Old code
        ! i1 = findloc(ebTab%times .le. t,.true.,dim=1,back=.true.)
        ! i2 = findloc(ebTab%times .gt. t,.true.,dim=1)
        ! i1 = max(1,i1)
        ! i2 = min(ebTab%N,i2)
    
        !Work-around code        
        do n=1,ebTab%N
            if (ebTab%times(n)>t) exit
        enddo
        i1 = n-1
        i2 = i1+1
        i1 = max(1,i1)
        
        if (i2 == i1) i2=i1+1 !Possible if none of the tab slices are in range
        i2 = min(ebTab%N,i2)
        
        !write(*,*) 'i1 / i2 = ', i1,i2
        !write(*,*) 'T(i1) / T / T(i2) = ', oTScl*ebTab%times(i1),oTScl*t,oTScl*ebTab%times(i2)

    end subroutine findSlc

    !Convert time (chimp units) to MJD
    function MJDAt(ebTab,t)
        type(ebTab_T), intent(in) :: ebTab
        real(rp), intent(in) :: t
        real(rp) :: MJDAt

        real(rp) :: dt
        integer :: i1,i2
        if (.not. ebTab%hasMJD) then
            MJDAt = 0.0
            return
        endif
        call findSlc(ebTab,t,i1,i2)

        if (t>=ebTab%times(i1)) then
            dt = oTScl*(t-ebTab%times(i1)) !Seconds
            MJDAt = ebTab%MJDs(i1) + dt/(60.0*60.0*24.0)
            
        else
            MJDAt = ebTab%MJDs(i1)
        endif
    end function MJDAt
     
end module ebtabutils