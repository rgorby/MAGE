!Driver for testing gremix

program gremixx
    use shellGrid
    use shellInterp

    implicit none

    real(rp) :: LowLatBoundary, HighLatBoundary,QQ
    real(rp), allocatable, dimension(:) :: t,p
    real(rp), allocatable, dimension(:,:) :: Q,Qinterp
    integer :: i,j, Nt=180, Np=360
    type(ShellGrid_T) :: shGr

    LowLatBoundary = PI !180.*deg2rad
    HighLatBoundary = 0.

    allocate(t(Nt))
    allocate(p(Np))
    allocate(Q(Nt-1,Np-1))
    allocate(Qinterp(Nt-3,Np-3))

    do i=1,Nt
        t(i) = (LowLatBoundary-HighLatBoundary)/(Nt-1)*(i-1) + HighLatBoundary
    end do

    do j=1,Np
       p(j) = 2*pi/(Np-1)*(j-1)
    end do

    do i=1,Nt-1
       do j=1,Np-1
          Q(i,j) = sin(t(i))*cos(p(j))
       end do
    end do

    call GenShellGrid(shGr,t,p)
    do i=2,Nt-2
       do j=2,Np-2
          call InterpShell(shGr,Q,0.5*(t(i)+t(i+1)),0.5*(p(j)+p(j+1)),Qinterp(i-1,j-1))
          write(*,*) 0.5*(t(i)+t(i+1))*rad2deg,0.5*(p(j)+p(j+1))*rad2deg,Qinterp(i-1,j-1)
       end do
    end do
    
end program gremixx
