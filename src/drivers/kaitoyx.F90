!Test driver for compilation debugging

program kaitoyx
    use kdefs
    use clocks
    use xml_input
    use ioH5
    use files

    implicit none

    type ToyGr_T
        integer :: is,ie,js,je,ks,ke
    end type ToyGr_T


    integer, parameter :: N = 10000
    integer, parameter :: M = 500
    integer, parameter :: MAXIOVAR = 5

    integer :: i,j,k
    real(rp), dimension(:), allocatable :: X
    real(rp), dimension(:,:,:), allocatable :: Z
    real(rp) :: ijkQ
    character(len=strLen) :: H5Out
    type(IOVAR_T), dimension(MAXIOVAR) :: IOVars
    type(ToyGr_T) :: ToyGr

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

    !Do test with associate/OMP
    write(*,*) 'Doing 3D test ...'
    call Tic("GenData3D")
    allocate(Z(M,M,M))
    ToyGr%is = 1
    ToyGr%js = 1
    ToyGr%ks = 1
    ToyGr%ie = M
    ToyGr%je = M
    ToyGr%ke = M

    associate(Gr=>toyGr)
    !$OMP PARALLEL DO default(shared) collapse(2) &
    !$OMP private(i,j,k,ijkQ)
    do k=Gr%ks,Gr%ke
        do j=Gr%js,Gr%je
            do i=Gr%is,Gr%ie
                ijkQ = i*j*k*1.0
                Z(i,j,k) = REarth*sin(ijkQ)
            enddo
        enddo
    enddo
    end associate
    call Toc("GenData3D")
    write(*,*) 'Finished 3D test ...'

    !Output data
    call Tic("IOData")
    call CheckAndKill(H5Out)

    call ClearIO(IOVars)
    call AddOutVar(IOVars,"X",X)
    call AddOutVar(IOVars,"Z",Z)
    call WriteVars(IOVars,.true.,H5Out)
    call Toc("IOData")

    call Toc("OMEGA")
    call printClocks()

end program kaitoyx