!> Reads REMIX data from output files, performs interpolation in time, etc.
module remixReader

    use kdefs
    use ioH5
    use XML_Input
    use ioH5
    use shellGrid

    implicit none

    integer :: MAXIOVARS = 20


    type rmState_T
        ! Model stuff
        logical :: doPed, doHall
            !! TODO: used by calcdb. Maybe they shouldn't live here. We don't do anything with them inside this module

        ! Grid stuff
        real(rp), dimension(:,:,:), allocatable :: XY
            !! X/Y coordinates in 2D. Convenient since we need it for lots of 2D arrays
        type(shellGrid_T) :: shGr

        ! State stuff
        real(rp) :: time
            !! Current sim time according to whoever's using rmState
        integer :: i1=-1,i2=-1
            !! Input file step numbers bracketing time

        type(ShellGridVar_T), dimension(:) ::
    end type rmState_T


    contains


    subroutine initRM(ftag,inpXML,rmState)
        character(len=strLen), intent(in) :: ftag
            !! Filename tag in this format: <ftag>.mix.h5
        type(XML_Input_T), intent(in) :: inpXML
        type(rmState_T), intent(inout) :: rmState

        type(IOVAR_T), dimension(MAXIOVARS) :: IOVars
        real(rp), dimension(:), allocatable :: th1D, ph1D
            !! 1D arrays of theta and phi coordinates, derived from XY locations
        real(rp), dimension(:,:,:), allocatable :: tmpXY
            !! Remix data is coming in at (Np,Nt), we use this to read and then convert to (Nt,Np)
        real(rp), dimension(:,:), allocatable :: tmpXcn, tmpYcn
            !! Temporary X/Y in real mix corner format
        integer :: Nt,Np
            !! # of non-ghost cell centers
        character(len=strLen) :: rmF 
            !! Remix file
        integer :: j
        
        write(rmF,'(2a)') trim(adjustl(ftag)),'.mix.h5'

        write(*,*) 'Initializing RMState w/ ', trim(rmF)

        ! Start with grid
        call ClearIO(IOVars)
        call AddInVar(IOVars,"X")
        call AddInVar(IOVars,"Y")
        call ReadVars(IOVars,.true.,rmF)


        !! Grid conversion
        !! X/Y in file are grid corners around the data
        !! To get original X/Y coords, we need to take the average and retrieve the cell-centered coordinates
        !! This means that the variable data is now stored at cell corners
        !! But, the outputted variable data only goes from phi [0,2Pi-dphi]
        !! So in order to comply with shellGrid expectations, we need to copy first phi column to end of array as well
        Np = IOVars(FindIO(IOVars, "X"))%dims(1)
        Nt = IOVars(FindIO(IOVars, "X"))%dims(2)
        write(*,*)Nt,Np
        allocate(rmState%XY(Nt-1,Np,XDIR:YDIR))
        allocate(tmpXY(Np,Nt,XDIR:YDIR))
        call IOArray2DFill(IOVars,"X",tmpXY(:,:,XDIR))
        call IOArray2DFill(IOVars,"Y",tmpXY(:,:,YDIR))

        ! One day we will get grid in this format directly from the mix.h5 file
        ! But today is not that day
        call genInGrid(tmpXY(:,:,XDIR), tmpXY(:,:,YDIR), tmpXcn, tmpYcn)
        rmState%XY(:,:Np-1,XDIR) = transpose(tmpXcn)
        rmState%XY(:,:Np-1,YDIR) = transpose(tmpYcn)
        ! Wrap in j
        rmState%XY(:,Np,:) = rmState%XY(:,1,:)

        ! XY are in units of ionospheric radii (r = 1 Ri)
        ! So the theta and phi we calculate are also for the ionospheric grid
        allocate(th1D(Nt-1))
        allocate(ph1D(Np))
        th1D = asin(rmState%XY(:,1,XDIR))
        ph1D = atan2(rmState%XY(Nt-1,:,YDIR), rmState%XY(Nt-1,:,XDIR))

        ! Clean up phi for shellGrid generation
        do j=Np/2,Np
            if (ph1D(j) < 0) then
                ph1D(j) = ph1D(j) + 2.0_rp*PI
            endif
        enddo
        if (abs(ph1D(1)) < TINY) then
            ph1D(1) = 0.0
        endif
        !write(*,*) "th1D"
        !write(*,*) th1D
        !write(*,*) "ph1D"
        !write(*,*) ph1D*180.0_rp/PI
        call GenShellGrid(rmState%shGr, th1D, ph1D)

        ! Hooray we have a shellGrid now



    end subroutine initRM



!------
! Temporary Helpers
!------

    subroutine genInGrid(xc,yc,x,y)
        !! Mix h5 grid to "real" mix grid
        !! NOTE: This has been modified from the version in mixio
        real(rp), dimension(:,:),intent(in) :: xc,yc ! with corners 1/2-cell shifted from original    
        real(rp), dimension(:,:),allocatable,intent(out) :: x,y
    
        integer, dimension(2) :: dims
        integer :: Np, Nt
    
        dims = shape(xc)
        Np = dims(1)-1; Nt = dims(2)
    
        if (.not.allocated(x)) allocate(x(Np,Nt-1))
        if (.not.allocated(y)) allocate(y(Np,Nt-1))
    
        !x(:,2:Nt) = 0.25*(xc(2:Np+1,2:Nt)+xc(1:Np,2:Nt)+xc(2:Np+1,1:Nt-1)+xc(1:Np,1:Nt-1))
        !y(:,2:Nt) = 0.25*(yc(2:Np+1,2:Nt)+yc(1:Np,2:Nt)+yc(2:Np+1,1:Nt-1)+yc(1:Np,1:Nt-1))
        x = 0.25*(xc(2:Np+1,2:Nt)+xc(1:Np,2:Nt)+xc(2:Np+1,1:Nt-1)+xc(1:Np,1:Nt-1))
        y = 0.25*(yc(2:Np+1,2:Nt)+yc(1:Np,2:Nt)+yc(2:Np+1,1:Nt-1)+yc(1:Np,1:Nt-1))
    
        ! fix pole
        !x(:,1) = 0
        !y(:,1) = 0
      end subroutine genInGrid

end module remixReader