!> Various data structures and routines to define grids on spherical shells
module shellGrid
    use kdefs
    use math
    
    ! note, this will conflict with mixdefs
    ! but NORTH, SOUTH are defined the same way
    enum, bind(C)
        enumerator :: NORTH=1,SOUTH,EAST,WEST
    end enum
    
    !> Data type for holding 2D spherical shell grid
    type ShellGrid_T
        integer :: nShellVars = 1 
            !! Number of shell grid functions
        type(shellGridFunction_T), dimension(:), allocatable :: shellVars
            !! Array of (Nt, Np) cell-centered variables
        
        integer :: Np,Nt 
            !! Number of lon/lat cells (phi/theta)
        ! xxc = centers (Nx)
        ! xx  = corners (Nx+1)
        ! th (theta) runs from north pole toward south
        ! Assuming lat \in -pi/2,pi/2 and lon \in [0,2pi]
        ! note, th/ph are colatitude/longitude; also, defining lat/latc for latitude
        real(rp), dimension(:), allocatable :: th, ph, thc, phc, lat, latc  
            !! [radians]
        logical :: doSP = .false., doNP = .false. 
            !! Whether active grid contains south/north pole, no ghosts in this case
        logical :: ghostSP = .false., ghostNP = .false. 
            !! Whether the last ghost corner is the south/north pole

        real(rp) :: minTheta, maxTheta, minPhi, maxPhi
            !! Theta and phi bounds of grid excluding ghost cells
        real(rp) :: minGTheta, maxGTheta, minGPhi, maxGPhi
            !! Theta and phi bounds of grid including ghost cells

        integer :: is=1,ie,js=1,je
            !! Local indices, active region
        integer :: isg,ieg,jsg,jeg
            !! Local indices, including ghosts

        integer :: Ngn=0, Ngs=0, Nge=1, Ngw=1
            !! Ghosts for north, south, east, and west boundaries. East=>larger phi. West => smaller phi.
            !! default to 0

        logical :: isPeriodic   
            !! Whether the low/high phi boundary is periodic
        logical :: isPhiUniform 
            !! Define this to speed up search for interpolation
    end type ShellGrid_T

    type ShellGridFunction_T
        real(rp), dimension(:,:), allocatable :: cellData
! corner data currently unused. If ever uncommented, need to allocate/init in initShellGridFunctions below
!        real(rp), dimension(:,:), allocatable :: cornerData 
        logical, dimension(4) :: bcsApplied 
            !! Flag indicating whether BCs were applied (ghosts filled) for [n,s,e,w] boundaries
    end type ShellGridFunction_T

    contains

    !> Create a shell grid data structure
    !> Takes Theta and Phi 1D arrays (uniform or not)
    !> Decides if isPeriodic and isPhiUniform based on the Phi array passed in
    subroutine GenShellGrid(shGr,Theta,Phi,nGhosts,nShellVarsO)
        type(ShellGrid_T), intent(inout) :: shGr
        real(rp), dimension(:), intent(in) :: Theta, Phi
        integer, optional, dimension(4), intent(in) :: nGhosts 
            !! How many ghosts on each side (n,s,e,w)
        integer, optional, intent(in) :: nShellVarsO
            !! Number of shell vars we should allocate space for

        integer :: i,j
        real(rp) :: delta
        real(rp), dimension(:), allocatable :: dphi

        ! Parse optional parameters
        if (present(nGhosts)) then
            shGr%Ngn = nGhosts(NORTH) ! otherwise, always 0 as set in ShellGrid type
            shGr%Ngs = nGhosts(SOUTH) ! otherwise, always 0 as set in ShellGrid type
            shGr%Nge = nGhosts(EAST) ! otherwise, always 1 as set in ShellGrid type
            shGr%Ngw = nGhosts(WEST) ! otherwise, always 1 as set in ShellGrid type
        end if

        if (present(nShellVarsO)) then
            shGr%nShellVars = nShellVarsO
        endif

        ! do some checks first
        if (.not.(isAscending(Theta))) then 
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Theta array must be ascending. Quitting..."
            stop
        end if

        if (.not.(isAscending(Phi))) then 
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Phi array must be ascending. Quitting..."
            stop
        end if

        if (any(Theta<0.).or.(any(Theta>PI))) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Theta array should be in the range [0,PI]. Quitting..."
            stop
        end if

        if (any(Phi<0.).or.(any(Phi>2*PI))) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Phi array should be in the range [0,2*PI]. Quitting..."
            stop
        end if

        ! decide if the grid is periodic
        if ( ( minval(Phi) <= TINY ).and.( maxval(Phi) >= 2*PI - TINY) ) then
            shGr%isPeriodic = .True.
        else
            shGr%isPeriodic = .False.
        end if
        
        if ( (shGr%isPeriodic).and.(shGr%Nge/=shGr%Ngw) ) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Periodic grid must have the same number of ghosts on low/high phi boundaries. Quitting..."
            stop
        end if

        if ( (.not.(shGr%isPeriodic)).and.(shGr%Nge/=shGr%Ngw) ) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Grids with Nge/=Ngw have not been implemented. Quitting..."
            stop
        end if
        
        ! Decide if this has north/south pole
        shGr%doNP = .false.
        shGr%doSP = .false.
        if ( maxval(Theta) >= (PI - TINY) ) then
            shGr%doSP = .true.
        endif

        if ( minval(Theta) <= (TINY) ) then
            shGr%doNP = .true.
        endif

        if ( (shGr%doNP).and.(shGr%Ngn/=0) ) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Grid containing the north pole can't have Ngn/=0. Quitting..."
            stop
        end if

        if ( (shGr%doSP).and.(shGr%Ngs/=0) ) then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Grid containing the south pole can't have Ngs/=0. Quitting..."
            stop
        end if

        ! define various array indices
        shGr%Np = size(Phi)-1
        shGr%Nt = size(Theta)-1

        shGr%is = 1; shGr%ie = shGr%Nt
        shGr%js = 1; shGr%je = shGr%Np

        shGr%isg = shGr%is - shGr%Ngn
        shGr%ieg = shGr%ie + shGr%Ngs
        shGr%jsg = shGr%js - shGr%Ngw
        shGr%jeg = shGr%je + shGr%Nge

        ! decide if the grid is uniform in Phi
        ! helps to speed up interpolation search
        if (allocated(dphi)) deallocate(dphi)
        allocate(dphi(shGr%Np))   ! helper
        dphi = Phi(2:shGr%Np+1) - Phi(1:shGr%Np)
        if ( all( dphi - dphi(1) <= TINY ) ) then
            shGr%isPhiUniform = .True.
        else 
            shGr%isPhiUniform = .False.
        end if
        deallocate(dphi)

        associate(is=>shGr%is, ie=>shGr%ie, js=>shGr%js, je=>shGr%je, &
            isg=>shGr%isg, ieg=>shGr%ieg, jsg=>shGr%jsg, jeg=>shGr%jeg)

        ! Nuke arrays if already allocated
        if (allocated(shGr%th))   deallocate(shGr%th)
        if (allocated(shGr%ph))   deallocate(shGr%ph)
        if (allocated(shGr%thc))  deallocate(shGr%thc)
        if (allocated(shGr%phc))  deallocate(shGr%phc)

        ! Create new arrays
        allocate(shGr%th  (isg:ieg+1))
        allocate(shGr%thc (isg:ieg  ))
        allocate(shGr%ph  (jsg:jeg+1))
        allocate(shGr%phc (jsg:jeg  ))

        ! Set grid coordinates
        shGr%th(is:ie+1) = Theta  ! note the arrays are conformable because of the index definitions above
        shGr%ph(js:je+1) = Phi

        ! Define ghost coordinates

        ! Do north, unless it's the pole (no ghosts for poles)
        if ( (.not.shGr%doNP).and.(shGr%Ngn/=0) ) then
            ! linearly extrapolate
            do i=1,shGr%Ngn
                shGr%th(is-i) = 2*shGr%th(is) - shGr%th(is+i)
            end do

            ! check if we got into negative thetas
            ! then just fill the space with same size cells down to theta=0
            if ( shGr%th(is-shGr%Ngn) < 0 ) then
                delta = shGr%th(is)/shGr%Ngn

                do i=1,shGr%Ngn
                    shGr%th(is-i) = shGr%th(is) - delta*i
                end do

                shGr%ghostNP = .true.
            end if
        end if

        ! Do south, unless it's the pole (no ghosts for poles)
        if ( (.not.shGr%doSP).and.(shGr%Ngs/=0) ) then
            ! linearly extrapolate
            do i=1,shGr%Ngs
                shGr%th(ie+1+i) = 2*shGr%th(ie+1) - shGr%th(ie+1-i)
            end do

            ! check if we got into thetas > pi
            ! then just fill the space with same size cells up to theta=pi
            if ( shGr%th(ie+1+shGr%Ngs) > PI ) then
                delta = (PI - shGr%th(ie+1))/shGr%Ngs
                do i=1,shGr%Ngs
                    shGr%th(ie+1+i) = shGr%th(ie+1) + delta*i
                end do

                shGr%ghostSP = .true.
            end if
        end if

        ! if non-periodic grid is needed, can implement ghosts on E/W ends similarly to N/S above
        ! but we assume the grid is always periodic
        if ( (.not.shGr%isPeriodic) )  then
            write(*,*) "Inside shell grid generator (GenShellGrid)."
            write(*,*) "Non-periodic grids are not implemented. Quitting..."
            stop
        endif

        ! Note, we already made sure that for a periodic grid Nge=Ngw, phi(js)=0 and phi(je+1)=pi
        do j=1,shGr%Nge
            shGr%ph(js-j)   = shGr%ph(je+1-j) - 2*PI
            shGr%ph(je+1+j) = shGr%ph(js+j) + 2*PI
        end do

        ! Set centers
        shGr%thc = ( shGr%th(isg+1:ieg+1) + shGr%th(isg:ieg) )/2.
        shGr%phc = ( shGr%ph(jsg+1:jeg+1) + shGr%ph(jsg:jeg) )/2.

        ! Note, the bounds below only include the active cells
        shGr%minTheta = minval(shGr%th(is:ie+1))
        shGr%maxTheta = maxval(shGr%th(is:ie+1))
        shGr%minPhi   = minval(shGr%ph(js:je+1))
        shGr%maxPhi   = maxval(shGr%ph(js:je+1))

        ! this includes ghosts
        shGr%minGTheta = minval(shGr%th)
        shGr%maxGTheta = maxval(shGr%th)
        shGr%minGPhi   = minval(shGr%ph)
        shGr%maxGPhi   = maxval(shGr%ph)

        end associate

        ! finally, initialize the shell grid functions
        ! nuke everything 
        if (allocated(shGr%shellVars))   deallocate(shGr%shellVars)
        allocate(shGr%shellVars(shGr%nShellVars)) 
        do i=1,shGr%nShellVars
            call initShellGridFunctions(shGr%shellVars(i))
        end do

        contains
            subroutine initShellGridFunctions(shellVar)
                type(shellGridFunction_T), intent(out) :: shellVar
                
                if (.not.allocated(shellVar%cellData)) allocate(shellVar%cellData(shGr%Nt,shGr%Np))
                shellVar%cellData = 0.  ! initialize to 0
                
                ! unset all BC's
                shellVar%bcsApplied = .false.
                
            end subroutine initShellGridFunctions

    end subroutine GenShellGrid
end module shellGrid
