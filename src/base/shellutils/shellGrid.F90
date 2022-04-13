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
        real(rp) :: minTheta, maxTheta, minPhi, maxPhi

        ! Local indices, active region
        integer :: is=1,ie,js=1,je
        ! Local indices, including ghosts
        integer :: isg,ieg,jsg,jeg

        ! Ghosts for north, south, east, and west boundaries. East=>larger phi. West => smaller phi.
        ! default to 0
        integer :: Ngn=0, Ngs=0, Nge=0, Ngw=0

        logical ::  isPeriodic=.True.   ! Whether the low/high phi boundary is periodic
    end type ShellGrid_T

    contains

    ! Create a shell grid data structure
    ! Takes Theta and Phi 1D arrays (uniform or not)
    ! Default the number of ghosts to 0
    subroutine GenShellGrid(shGr,Theta,Phi,Ngn,Ngs,Nge,Ngw,isPeriodic)
        type(ShellGrid_T), intent(inout) :: shGr
        real(rp), dimension(:), intent(in) :: Theta, Phi
        integer, optional, intent(in) :: Ngn, Ngs, Nge, Ngw  ! how many ghosts on each side
        logical, optional, intent(in) :: isPeriodic

        integer :: Np, Nt

        ! Parse optional parameters
        if (present(isPeriodic)) shGr%isPeriodic = isPeriodic   ! otherwise, always periodic as set in ShellGrid type
        if (present(Ngn)) shGr%Ngn = Ngn ! otherwise, always 0 as set in ShellGrid type
        if (present(Ngs)) shGr%Ngs = Ngs ! otherwise, always 0 as set in ShellGrid type
        if (present(Ngw)) shGr%Ngw = Ngw ! otherwise, always 0 as set in ShellGrid type
        if (present(Nge)) shGr%Nge = Nge ! otherwise, always 0 as set in ShellGrid type

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

        if (any(Theta<0.).or.(any(Theta>PI/2))) then
           write(*,*) "Inside shell grid generator (GenShellGrid)."
           write(*,*) "Theta array should be in the range [0,PI/2]. Quitting..."
           stop
        end if

        if (any(Phi<0.).or.(any(Phi>2*PI))) then
           write(*,*) "Inside shell grid generator (GenShellGrid)."
           write(*,*) "Phi array should be in the range [0,2*PI]. Quitting..."
           stop
        end if

        if ( (shGr%isPeriodic).and.(shGr%Nge/=shGr%Ngw) ) then
           write(*,*) "Inside shell grid generator (GenShellGrid)."
           write(*,*) "Periodic grid must have the same number of ghosts on low/high phi boundaries. Quitting..."
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

        ! Nuke arrays if already allocated
        if (allocated(shGr%th))   deallocate(shGr%th)
        if (allocated(shGr%ph))   deallocate(shGr%ph)
        if (allocated(shGr%thc))  deallocate(shGr%thc)
        if (allocated(shGr%phc))  deallocate(shGr%phc)
        if (allocated(shGr%lat))  deallocate(shGr%lat)
        if (allocated(shGr%latc)) deallocate(shGr%latc)

        ! Create new arrays
        allocate(shGr%th  (shGr%Nt+shGr%Ngn+shGr%Ngs+1))
        allocate(shGr%ph  (shGr%Np+shGr%Ngw+shGr%Nge+1))
        allocate(shGr%thc (shGr%Nt+shGr%Ngn+shGr%Ngs))
        allocate(shGr%phc (shGr%Np+shGr%Ngw+shGr%Nge))
        allocate(shGr%lat (shGr%Nt+shGr%Ngn+shGr%Ngs+1))
        allocate(shGr%latc(shGr%Nt+shGr%Ngn+shGr%Ngs))


        ! !Set edges
        ! shGr%LatI(:) = iLats
        ! dphi = 2*PI/shGr%NLon
        ! do n=1,shGr%NLon+1
        !     shGr%LonI(n) = (n-1)*dphi
        ! enddo

        ! !Set centers
        ! shGr%LatC = ( shGr%LatI(2:shGr%NLat+1) + shGr%LatI(1:shGr%NLat) )/2.0
        ! shGr%LonC = ( shGr%LonI(2:shGr%NLon+1) + shGr%LonI(1:shGr%NLon) )/2.0

        ! !Decide if this has north/south pole
        ! shGr%doNP = .false.
        ! shGr%doSP = .false.
        ! if ( maxval(shGr%LatI) >= (PI/2.0 - TINY) ) then
        !     shGr%doNP = .true.
        ! endif

        ! if ( minval(shGr%LatI) <= (-PI/2.0 + TINY) ) then
        !     shGr%doSP = .true.
        ! endif

        ! if (shGr%doSP .or. shGr%doNP) then
        !     !Die for now
        !     write(*,*) "This routine does not yet support handling north/south pole"
        !     stop
        ! endif
        
        ! shGr%minLat = minval(shGr%LatI)
        ! shGr%maxLat = maxval(shGr%LatI)

    end subroutine GenShellGrid
end module shellGrid
