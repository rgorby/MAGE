!Test driver for compilation debugging

program kaitoyx
    use kdefs
    use clocks
    use xml_input
    use ioH5
    use files

    implicit none

    integer, parameter :: N = 10000
    integer, parameter :: MAXIOVAR = 5

    integer :: i
    real(rp), dimension(:), allocatable :: X
    character(len=strLen) :: H5Out
    type(IOVAR_T), dimension(MAXIOVAR) :: IOVars

    H5Out = "kaitoy.h5"

    !Setup timers
    call initClocks()
    call Tic("OMEGA")

    !Create some toy data w/ OMP
    call Tic("GenData")
    allocate(X(N))

    !$OMP PARALLEL DO
    do i=1,N
        X(i) = REarth*i*sin(i*i*1.0)
    enddo
    call Toc("GenData")

    !Output data
    call Tic("IOData")
    call CheckAndKill(H5Out)

    call ClearIO(IOVars)
    call AddOutVar(IOVars,"X",X)
    call WriteVars(IOVars,.true.,H5Out)
    call Toc("IOData")

    call Toc("OMEGA")
    call printClocks()

end program kaitoyx