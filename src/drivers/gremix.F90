!Driver for testing gremix

program gremixx
    use shellGrid

    implicit none

    real(rp) :: LowLatBoundary, HighLatBoundary
    real(rp), dimension(100) :: t,p
    integer :: i, Nt=100
    type(ShellGrid_T) :: shGr

    LowLatBoundary = 90.*deg2rad
    HighLatBoundary = 0.

    do i=1,Nt
        t(i) = (LowLatBoundary-HighLatBoundary)/(Nt-1)*(i-1) + HighLatBoundary
    end do

    do i=1,Nt
       p(i) = 2*pi/Nt*(i-1)+pi
    end do

    call GenShellGrid(shGr,t,p)
    
end program gremixx
