module mixinterp
  use mixdefs
  use mixtypes
  use math

  implicit none
  
  contains
    subroutine mix_interpolant(G)
      ! calculates and stores interpolant for grid G
      type(mixGrid_T),intent(inout) :: G
      real(rp), dimension(4,4) :: A,Ainv
      real(rp), dimension(4) :: B
      real(rp) :: dp ! delta phi to add to the right edge of the cell -- treats the periodic boundary
      integer :: i,j,jp1

      ! note, looping over cells -- not nodes -- here
      do i=1,G%Nt-1  
         do j=1,G%Np  ! we go to +1 point in phi to cover the gap between Np and 1
            if (j.ne.G%Np) then
               jp1=j+1
               dp = 0._rp
            else
               jp1=1
               dp = 2*pi
            endif

            A = reshape(&
                 (/&
                 1.0_rp,G%p(j,i),       G%t(j,i),     G%p(j,i)        *G%t(j,i),&
                 1.0_rp,G%p(jp1,i)  +dp,G%t(jp1,i),  (G%p(jp1,i)  +dp)*G%t(jp1,i),&
                 1.0_rp,G%p(j,i+1),     G%t(j,i+1),   G%p(j,i+1)      *G%t(j,i+1),&
                 1.0_rp,G%p(jp1,i+1)+dp,G%t(jp1,i+1),(G%p(jp1,i+1)+dp)*G%t(jp1,i+1)&
                 /)&
                 ,shape(A))
            
            ! exclude the poles to avoid inverting the bad matrix there
            if (.not.(  &
                 ( (i.eq.1).and.(G%t(1,1).le.TINY) ) .or. &
                 ( (i.eq.G%Nt-1).and.((pi-G%t(1,G%Nt)).le.TINY) ) &
                 )) then
               call matinv4(A,Ainv)
               G%Interpolant(j,i,:,:) = Ainv   ! for cell (j,i)
            end if
         enddo
      enddo
    end subroutine mix_interpolant

    subroutine mix_set_map(G1,G2,Map)
      type(mixGrid_T), intent(in) :: G1,G2   ! Interpolate from G1 to G2
      type(Map_T), Intent(out):: Map
      integer :: i1,j1,i2,j2,j,i
      real(rp), dimension(4) :: x  ! (1,p,t,pt) for point to which we interpolate

      if (.not.allocated(Map%M)) allocate(Map%M(G2%Np,G2%Nt,4))
      if (.not.allocated(Map%I1)) allocate(Map%I1(G2%Np,G2%Nt))
      if (.not.allocated(Map%J1)) allocate(Map%J1(G2%Np,G2%Nt))

      ! loop over points onto which we're interpolating
      do i2=1,G2%Nt
         do j2=1,G2%Np
            call mix_search(G1,G2%p(j2,i2),G2%t(j2,i2),j1,i1)
            Map%I1(j2,i2) = i1
            Map%J1(j2,i2) = j1
            x = (/1._rp, G2%p(j2,i2),G2%t(j2,i2), G2%p(j2,i2)*G2%t(j2,i2) /)

            ! coming out of mix_search: 
            ! i1=0 if the destination point (j2,i2) is above of the poleward boundary of the source grid (G1) or
            ! i1=G1%Nt if the destination point (j2,i2) is below the equatorward boundary of the source grid (G1)
            ! Set the Map to 0.0 there for now and see if we need a more nuanced treatement later
            if ( (i1.eq.0).or.(i1.eq.G1%Nt) ) then
              Map%M(j2,i2,:) = 0.0
              cycle
            end if

            ! treat the poles
            if ((i1.eq.1).and.(G1%t(1,1)).lt.TINY) then
               ! just average over the three vertices of the triangle
               ! could make it more fancy to interpolate using a0+a1*t+a2*t*p interpolant, but why bother???
               Map%M(j2,i2,:) = 1._rp/3._rp*(/0.5_rp, 0.5_rp, 1._rp, 1._rp/)
            else if ((i1.eq.G1%Nt-1).and.(pi-G1%t(1,G1%Nt-1)).lt.TINY) then
               Map%M(j2,i2,:) = 1._rp/3._rp*(/1._rp, 1._rp, 0.5_rp, 0.5_rp/)               
            else 
               Map%M(j2,i2,:) = matmul(G1%Interpolant(j1,i1,:,:),x)
            end if
         end do
      end do

    end subroutine mix_set_map

    subroutine mix_search(G,p,t,j,i)
      ! for a given point (p,t) finds (j,i) of the grid
      type(mixGrid_T),intent(in) :: G
      real(rp) :: p, t
      integer, intent(out) :: i,j

      ! FIXME: treat when outside the boundary
      j = minloc(abs(G%p(:,1)-p),1)
      i = minloc(abs(G%t(1,:)-t),1)

      ! correct for when we're close but just below this makes a big
      ! difference -- smoothes out the interpolation and removes
      ! spikes
      if ( (G%p(j,1)-p).gt.0 ) j=j-1
      if ( (G%t(1,i)-t).gt.0 ) i=i-1

      ! treat periodic
      if ( ((p-G%p(G%Np,1)).gt.0).and.(p.lt.2*pi) ) j=G%Np
    end subroutine mix_search

    subroutine mix_map_grids(Map,F1,F2)
      type(Map_T), Intent(in):: Map
      real(rp), dimension(:,:), intent(in) :: F1
      real(rp), dimension(:,:), allocatable, Intent(out):: F2
      integer :: i1,j1,j1p1,i2,j2,j,i,G2Np,G2Nt
      real(rp), dimension(4) :: F  ! vector of function F1 values at vertices surrounding the point on G2
      integer, dimension(3) :: dims

      ! get grid size from Map so we don't have to pass the grid itself
      dims = shape(Map%M); G2Nt = dims(2); G2Np = dims(1)
      
      if (.not.allocated(F2)) allocate(F2(G2Np,G2Nt))
      ! loop over points onto which we're interpolating
      do i2=1,G2Nt
         do j2=1,G2Np
            i1 = Map%I1(j2,i2)
            j1 = Map%J1(j2,i2)
            
            if (j1.eq.size(F1,1)) then
               j1p1=1
            else 
               j1p1=j1+1
            end if
            
            !FIXME: Something here causing overrun on index 2
            if (i1 == size(F1,2)) then
              F = (/ F1(j1,i1), F1(j1p1,i1), F1(j1,i1), F1(j1p1,i1)/)
            else
              F = (/ F1(j1,i1), F1(j1p1,i1), F1(j1,i1+1), F1(j1p1,i1+1)/)
            endif
            F2(j2,i2) = dot_product(Map%M(j2,i2,:),F)
         end do
      end do
    end subroutine mix_map_grids

end module mixinterp
