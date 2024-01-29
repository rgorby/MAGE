!> Reads REMIX data from output files, performs interpolation in time, etc.
module remixReader

    use kdefs
    use ioH5
    use XML_Input
    use ioH5
    use shellGrid

    implicit none

    integer :: MAXIOVARS = 20

    integer :: nShellVars = 8
    enum, bind(C)
        enumerator :: RMR_NFAC=1,RMR_NSIGP,RMR_NSIGH,RMR_NPOT, \
                        RMR_SFAC,RMR_SSIGP,RMR_SSIGH,RMR_SPOT
    end enum


    type rmState_T
        ! Model stuff
        logical :: doPed, doHall
            !! TODO: used by calcdb. Maybe they shouldn't live here. We don't do anything with them inside this module

        ! Grid stuff
        real(rp), dimension(:,:,:), allocatable :: XY
            !! X/Y coordinates in 2D. Convenient since we need it for lots of 2D arrays

        ! shGrid is a mix (haha) of Grid and State info. idk I just work here
        type(shellGrid_T) :: shGrid
            !! Also holds current "State" of remix variables within shellVars

        ! State stuff
        real(rp) :: time
            !! Current sim time according to whoever's using rmState
        integer :: i1=-1,i2=-1
            !! Input file step numbers bracketing time
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

        Np = IOVars(FindIO(IOVars, "X"))%dims(1) - 1
        Nt = IOVars(FindIO(IOVars, "X"))%dims(2) - 1
        write(*,*)Nt,Np
        allocate(rmState%XY(Np+1,Nt+1,XDIR:YDIR))
        call IOArray2DFill(IOVars,"X",rmState%XY(:,:,XDIR))
        call IOArray2DFill(IOVars,"Y",rmState%XY(:,:,YDIR))
        
        ! XY are in units of ionospheric radii (r = 1 Ri)
        ! So the theta and phi we calculate are also for the ionospheric grid
        allocate(th1D(Nt+1))
        allocate(ph1D(Np+1))
        th1D = asin(rmState%XY(1,:,XDIR))
        ph1D = atan2(rmState%XY(:,Nt,YDIR), rmState%XY(:,Nt,XDIR))
        do j=Np/2,Np+1
            if (ph1D(j) < 0) then
                ph1D(j) = ph1D(j) + 2.0_rp*PI
            endif
        enddo
        write(*,*) "th1D"
        write(*,*) th1D
        write(*,*) "ph1D"
        write(*,*) ph1D
        call GenShellGrid(rmState%shGrid, th1D, ph1D, nShellVarsO=nShellVars)

    end subroutine initRM

end module remixReader