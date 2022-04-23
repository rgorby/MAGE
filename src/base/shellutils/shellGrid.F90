! Various data structures and routines to define grids on spherical shells
module shellGrid
    use kdefs
    use math

    ! Data type for holding 2D spherical shell grid
    type ShellGrid_T
        integer :: Np,Nt ! Number of lon/lat cells (phi/theta)
        ! xxc = centers (Nx)
        ! xx  = corners (Nx+1)
        ! th (theta) runs from north pole toward south
        ! Assuming lat \in -pi/2,pi/2 and lon \in [0,2pi]
        ! note, th/ph are colatitude/longitude; also, defining lat/latc for latitude
        real(rp), dimension(:), allocatable :: th, ph, thc, phc, lat, latc  ! Radians
        logical :: doSP = .false., doNP = .false. ! Whether grid contains south/north pole

        ! Local indices, active region
        integer :: is=1,ie,js=1,je
        ! Local indices, including ghosts
        integer :: isg,ieg,jsg,jeg

        ! Ghosts for north, south, east, and west boundaries. East=>larger phi. West => smaller phi.
        ! default to 0
        integer :: Ngn=0, Ngs=0, Nge=1, Ngw=1

        logical :: isPeriodic   ! Whether the low/high phi boundary is periodic
    end type ShellGrid_T

    contains

    ! Create a shell grid data structure
    ! Takes Theta and Phi 1D arrays (uniform or not)
    ! Default the number of ghosts to 0
    subroutine GenShellGrid(shGr,Theta,Phi,Ngn,Ngs,Nge,Ngw)
        type(ShellGrid_T), intent(inout) :: shGr
        real(rp), dimension(:), intent(in) :: Theta, Phi
        integer, optional, intent(in) :: Ngn, Ngs, Nge, Ngw  ! how many ghosts on each side

        integer :: i,j
        integer :: Np, Nt
        real(rp) :: delta

        ! Parse optional parameters
        if (present(Ngn)) shGr%Ngn = Ngn ! otherwise, always 0 as set in ShellGrid type
        if (present(Ngs)) shGr%Ngs = Ngs ! otherwise, always 0 as set in ShellGrid type
        if (present(Ngw)) shGr%Ngw = Ngw ! otherwise, always 1 as set in ShellGrid type
        if (present(Nge)) shGr%Nge = Nge ! otherwise, always 1 as set in ShellGrid type

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

        associate(is=>shGr%is, ie=>shGr%ie, js=>shGr%js, je=>shGr%je, &
             isg=>shGr%isg, ieg=>shGr%ieg, jsg=>shGr%jsg, jeg=>shGr%jeg)

        ! Nuke arrays if already allocated
        if (allocated(shGr%th))   deallocate(shGr%th)
        if (allocated(shGr%ph))   deallocate(shGr%ph)
        if (allocated(shGr%thc))  deallocate(shGr%thc)
        if (allocated(shGr%phc))  deallocate(shGr%phc)
        if (allocated(shGr%lat))  deallocate(shGr%lat)
        if (allocated(shGr%latc)) deallocate(shGr%latc)

        ! Create new arrays
        allocate(shGr%th  (isg:ieg+1))
        allocate(shGr%thc (isg:ieg  ))
        allocate(shGr%lat (isg:ieg+1))
        allocate(shGr%latc(isg:ieg  ))
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

        !Set centers
        shGr%thc = ( shGr%th(isg+1:ieg+1) + shGr%th(isg:ieg) )/2.
        shGr%phc = ( shGr%ph(jsg+1:jeg+1) + shGr%ph(jsg:jeg) )/2.

        ! Define latitudes just in case
        shGr%lat  = PI/2 - shGr%th
        shGr%latc = PI/2 - shGr%thc

        end associate

        ! TODO: define isPhiUniform
    end subroutine GenShellGrid
end module shellGrid
