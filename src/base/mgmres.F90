module mgmres

    use kdefs
    
    implicit none

    real(rp), parameter :: delta = 1.0e-3

    contains

subroutine ax_cr ( n, nz_num, ia, ja, a, x, w )

  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) k2
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) x(n)

  w(1:n) = 0.0D+00

  do i = 1, n
    k1 = ia(i)
    k2 = ia(i+1) - 1
    w(i) = w(i) + dot_product ( a(k1:k2), x(ja(k1:k2)) )
  end do

  return
end

subroutine diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) ua(n)

  ua(1:n) = -1

  do i = 1, n
    do k = ia(i), ia(i+1) - 1
      if ( ja(k) == i ) then
        ua(i) = k
      end if
    end do
  end do

  return
end
subroutine ilu_cr ( n, nz_num, ia, ja, a, ua, l )

  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) iw(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jrow
  integer ( kind = 4 ) jw
  integer ( kind = 4 ) k
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) tl
  integer ( kind = 4 ) ua(n)
!
!  Copy A.
!
  l(1:nz_num) = a(1:nz_num)

  do i = 1, n
!
!  IW points to the nonzero entries in row I.
!
    iw(1:n) = -1

    do k = ia(i), ia(i+1) - 1
      iw(ja(k)) = k
    end do

    do j = ia(i), ia(i+1) - 1
      jrow = ja(j)
      if ( i <= jrow ) then
        exit
      end if
      tl = l(j) * l(ua(jrow))
      l(j) = tl
      do jj = ua(jrow) + 1, ia(jrow+1) - 1
        jw = iw(ja(jj))
        if ( jw /= -1 ) then
          l(jw) = l(jw) - tl * l(jj)
        end if
      end do
    end do

    ua(i) = j

    if ( jrow /= i ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a)' ) '  JROW ~= I'
      write ( *, '(a,i8)' ) '  JROW = ', jrow
      write ( *, '(a,i8)' ) '  I    = ', i
      stop
    end if

    if ( l(j) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ILU_CR - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step I = ', i
      write ( *, '(a,i8,a)' ) '  L(', j, ') = 0.0'
      stop
    end if

    l(j) = 1.0D+00 / l(j)

  end do

  l(ua(1:n)) = 1.0D+00 / l(ua(1:n))

  return
end

! !---------------
!     !A*x (compressed row)
!     subroutine ax_cr(n,nz_num,ia,ja,A,x,W)
!         integer , intent(in)  :: n,nz_num
!         integer , intent(in)  :: ia(n+1),ja(nz_num)
!         real(rp), intent(in)  :: A(nz_num)
!         real(rp), intent(in)  :: x(n)
!         real(rp), intent(inout) :: W(n)

!         integer :: i,k1,k2

!         W(1:n) = 0.0
!         !TODO: Add OMP bindings here
!         do i=1,n
!             k1 = ia(i)
!             k2 = ia(i+1) - 1
!             W(i) = W(i) + dot_product( A(k1:k2), x(ja(k1:k2)) )
!         enddo
!     end subroutine ax_cr
! !---------------
!     subroutine diagonal_pointer_cr(n,nz_num,ia,ja,ua)
!         integer, intent(in) :: n,nz_num
!         integer, intent(in) :: ia(n+1),ja(nz_num)
!         integer, intent(out) :: ua(n)

!         integer :: i,k

!         ua(1:n) = -1
!         do i=1,n
!             do k=ia(i),ia(i+1)-1
!                 if ( ja(k) == i ) then
!                     ua(i) = k
!                 endif
!             enddo
!         enddo

!     end subroutine diagonal_pointer_cr


! !---------------
!     !Incomplete LU factorization
!     subroutine ilu_cr(n,nz_num,ia,ja,A,ua,L)
!         integer, intent(in) :: n,nz_num
!         integer, intent(in) :: ia(n+1),ja(nz_num)
!         integer, intent(inout) :: ua(n)
!         real(rp), intent(in) :: A(nz_num)
!         real(rp), intent(out) :: L(nz_num)

!         integer :: i,j,jj,jrow,jw,k
!         integer :: iw(n)
!         real(rp) :: tl

!         L(1:nz_num) = A(1:nz_num)
!         do i=1,n
!             iw(1:n) = -1
!             do k=ia(i),ia(i+1)-1
!                 iw(ja(k)) = k
!             enddo

!             do j=ia(i),ia(i+1)-1
!                 jrow = ja(j)
!                 if ( i<= jrow ) exit
!                 tl = L(j)*L(ua(jrow))
!                 L(j) = tl
!                 do jj=ua(jrow)+1,ia(jrow+1)-1
!                     jw = iw(ja(jj))
!                     if ( jw /= -1) then
!                         L(jw) = L(jw) - tl*L(jj)
!                     endif
!                 enddo
!             enddo

!             ua(i) = j
!             L(j) = 1.0/L(j)
!         enddo
!         L(ua(1:n)) = 1.0/L(ua(1:n))
!     end subroutine ilu_cr

subroutine lus_cr ( n, nz_num, ia, ja, l, ua, r, z )

  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  real ( kind = 8 ) l(nz_num)
  real ( kind = 8 ) r(n)
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) w(n)
  real ( kind = 8 ) z(n)
!
!  Copy R in.
!
  w(1:n) = r(1:n)
!
!  Solve L * w = w where L is unit lower triangular.
!
  do i = 2, n
    do j = ia(i), ua(i) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
  end do
!
!  Solve U * w = w, where U is upper triangular.
!
  do i = n, 1, -1
    do j = ua(i) + 1, ia(i+1) - 1
      w(i) = w(i) - l(j) * w(ja(j))
    end do
    w(i) = w(i) / l(ua(i))
  end do
!
!  Copy Z out.
!
  z(1:n) = w(1:n)

  return
end

! !---------------
!     subroutine lus_cr(n,nz_num,ia,ja,L,ua,R,Z)
!         integer, intent(in) :: n,nz_num
!         integer, intent(in) :: ia(n+1),ja(nz_num),ua(n)
!         real(rp), intent(inout) :: L(nz_num)
!         real(rp), intent(in) :: R(n)
!         real(rp), intent(out) :: Z(n)

!         integer :: i,j
!         real(rp), dimension(n) :: W

!         !Copy in
!         W(1:n) = R(1:n)

!         !Solve L*W=W, L unit lower triangular
!         !TODO: Add OMP bindings here
!         do i=2,n
!             do j=ia(i),ua(i)-1
!                 W(i) = W(i) - L(j)*W(ja(j))
!             enddo
!         enddo

!         !Solve U*W=W, U upper triangular
!         do i=n,1,-1
!             do j=ua(i)+1,ia(i+1)-1
!                 W(i) = W(i) - L(j)*W(ja(j))
!             enddo
!             W(i) = W(i)/L(ua(i))
!         enddo

!         !Copy out
!         Z(1:n) = W(1:n)

!     end subroutine lus_cr

! subroutine mult_givens ( c, s, k, g )

!   implicit none

!   integer ( kind = 4 ) k

!   real ( kind = 8 ) c
!   real ( kind = 8 ) g(1:k+1)
!   real ( kind = 8 ) g1
!   real ( kind = 8 ) g2
!   real ( kind = 8 ) s

!   g1 = c * g(k) - s * g(k+1)
!   g2 = s * g(k) + c * g(k+1)

!   g(k)   = g1
!   g(k+1) = g2

!   return
! end

!---------------
    !Givens rotation
    subroutine mult_givens(c,s,k,g)
        real(rp), intent(in) :: c,s !cos/sin of rotation
        integer, intent(in) :: k
        real(rp), intent(inout) :: g(1:k+1)

        real(rp) :: g1,g2
        g1 = c*g(k) - s*g(k+1)
        g2 = s*g(k) - c*g(k+1)

        g(k  ) = g1
        g(k+1) = g2
    end subroutine mult_givens

! !---------------
!     !Preconditioned restarted GMRES
!     subroutine pmgmres_ilu_cr(n,nz_num,ia,ja,A,x,rhs,itr_max,mr,tol_abs,tol_rel)
!         integer, intent(in) :: n,nz_num
!         integer, intent(inout) :: ia(n+1),ja(nz_num)
!         real(rp), intent(inout) :: A(nz_num)
!         real(rp), intent(inout) :: x(n)
!         real(rp), intent(in) :: rhs(n)
!         integer, intent(in) :: itr_max,mr
!         real(rp), intent(in) :: tol_abs,tol_rel

!         real(rp) :: av,htmp,mu,rho,rho_tol
!         real(rp) :: r(n)
!         integer :: ua(n)
!         real(rp), dimension(mr+1) :: c,g,s,y
!         real(rp) :: H(mr+1,mr)
!         real(rp) :: L(ia(n+1)+1)
!         real(rp) :: V(n,mr+1)
!         integer :: i,itr,itr_used,j,k,k_copy


!         itr_used = 0
!         call rearrange_cr(n,nz_num,ia,ja,A)
!         call diagonal_pointer_cr(n,nz_num,ia,ja,ua)
!         call ilu_cr(n,nz_num,ia,ja,A,ua,L)

!         do itr=1,itr_max
!             call ax_cr(n,nz_num,ia,ja,A,x,r)
!             r(1:n) = rhs(1:n) - r(1:n)
!             call lus_cr(n,nz_num,ia,ja,L,ua,r,r)
!             rho = sqrt(dot_product(r,r))

!             if (itr == 1) then
!                 rho_tol = rho*tol_rel
!             endif

!             V(1:n,1) = r(1:n)/rho
!             g(1) = rho
!             g(2:mr+1) = 0.0
!             H(1:mr+1,1:mr) = 0.0

!             do k=1,mr
!                 k_copy = k
!                 call ax_cr (n,nz_num,ia,ja,A,V(1:n,k),V(1:n,k+1) )
!                 call lus_cr(n,nz_num,ia,ja,L,ua,V(1:n,k+1),V(1:n,k+1) )
!                 av = sqrt( dot_product( V(1:n,k+1), V(1:n,k+1) ) )
!                 do j=1,k
!                     H(j,k) = dot_product ( V(1:n,k+1), V(1:n,j) )
!                     V(1:n,k+1) = V(1:n,k+1) - V(1:n,j) * H(j,k)
!                 enddo

!                 H(k+1,k) = sqrt ( dot_product ( V(1:n,k+1), V(1:n,k+1) ) )

!                 if ( ( av + delta * H(k+1,k)) == av ) then
!                     do j=1,k
!                         htmp = dot_product ( V(1:n,k+1), V(1:n,j) )
!                         H(j,k) = H(j,k) + htmp
!                         V(1:n,k+1) = V(1:n,k+1) - htmp * V(1:n,j)
!                     enddo
!                     H(k+1,k) = sqrt ( dot_product ( V(1:n,k+1), V(1:n,k+1) ) )
!                 endif

!                 if ( abs(H(k+1,k)) >= TINY ) then
!                     V(1:n,k+1) = V(1:n,k+1)/H(k+1,k)
!                 endif

!                 if ( 1<k ) then
!                     y(1:k+1) = H(1:k+1,k)
!                     do j=1,k-1
!                         call mult_givens(c(j),s(j),j,y)
!                     enddo
!                     H(1:k+1,k) = y(1:k+1)
!                 endif

!                 mu = sqrt( H(k,k)**2.0 + H(k+1,k)**2.0)
!                 c(k) =  H(k  ,k)/mu
!                 s(k) = -H(k+1,k)/mu
!                 H(k,k) = c(k)*H(k,k) - s(k)*H(k+1,k)
!                 H(k+1,k) = 0.0
!                 call mult_givens(c(k),s(k),k,g)
!                 rho = abs(g(k+1))

!                 itr_used = itr_used+1
!                 if ( rho <= rho_tol .and. rho <= tol_abs ) then
!                     !Exit if done
!                     exit
!                 endif
!             enddo !K-loop
!             k = k_copy - 1
!             y(k+1) = g(k+1)/H(k+1,k+1)

!             do i=k,1,-1
!                 y(i) = ( g(i) - dot_product ( H(i,i+1:k+1), y(i+1:k+1) ) ) / H(i,i)
!             enddo

!             do i=1,n
!                 x(i) = x(i) + dot_product ( V(i,1:k+1), y(1:k+1) )
!             enddo

!             if ( (rho<=rho_tol) .and. (rho<=tol_abs) ) then
!                 exit
!             endif


!         enddo !Iteration loop
!     end subroutine pmgmres_ilu_cr

subroutine pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, &
  tol_abs, tol_rel )

  implicit none

  integer ( kind = 4 ) mr
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  real ( kind = 8 ) av
  real ( kind = 8 ) c(mr+1)
  real ( kind = 8 ), parameter :: delta = 1.0D-03
  real ( kind = 8 ) g(mr+1)
  real ( kind = 8 ) h(mr+1,mr)
  real ( kind = 8 ) htmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) itr
  integer ( kind = 4 ) itr_max
  integer ( kind = 4 ) itr_used
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) k_copy
  real ( kind = 8 ) l(ia(n+1)+1)
  real ( kind = 8 ) mu
  real ( kind = 8 ) r(n)
  real ( kind = 8 ) rho
  real ( kind = 8 ) rho_tol
  real ( kind = 8 ) rhs(n)
  real ( kind = 8 ) s(mr+1)
  real ( kind = 8 ) tol_abs
  real ( kind = 8 ) tol_rel
  integer ( kind = 4 ) ua(n)
  real ( kind = 8 ) v(n,mr+1);
  logical, parameter :: verbose = .false.
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(mr+1)

  itr_used = 0

  call rearrange_cr ( n, nz_num, ia, ja, a )

  call diagonal_pointer_cr ( n, nz_num, ia, ja, ua )

  call ilu_cr ( n, nz_num, ia, ja, a, ua, l )

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR'
    write ( *, '(a,i4)' ) '  Number of unknowns = ', n
  end if

  do itr = 1, itr_max

    call ax_cr ( n, nz_num, ia, ja, a, x, r )

    r(1:n) = rhs(1:n) - r(1:n)

    call lus_cr ( n, nz_num, ia, ja, l, ua, r, r )

    rho = sqrt ( dot_product ( r, r ) )

    if ( verbose ) then
      write ( *, '(a,i4,a,g14.6)' ) '  ITR = ', itr, '  Residual = ', rho
    end if

    if ( itr == 1 ) then
      rho_tol = rho * tol_rel
    end if

    v(1:n,1) = r(1:n) / rho

    g(1) = rho
    g(2:mr+1) = 0.0D+00

    h(1:mr+1,1:mr) = 0.0D+00

    do k = 1, mr

      k_copy = k

      call ax_cr ( n, nz_num, ia, ja, a, v(1:n,k), v(1:n,k+1) ) 

      call lus_cr ( n, nz_num, ia, ja, l, ua, v(1:n,k+1), v(1:n,k+1) )

      av = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      do j = 1, k
        h(j,k) = dot_product ( v(1:n,k+1), v(1:n,j) )
        v(1:n,k+1) = v(1:n,k+1) - v(1:n,j) * h(j,k)
      end do

      h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )

      if ( ( av + delta * h(k+1,k)) == av ) then
        do j = 1, k
          htmp = dot_product ( v(1:n,k+1), v(1:n,j) )
          h(j,k) = h(j,k) + htmp
          v(1:n,k+1) = v(1:n,k+1) - htmp * v(1:n,j)
        end do
        h(k+1,k) = sqrt ( dot_product ( v(1:n,k+1), v(1:n,k+1) ) )
      end if

      if ( h(k+1,k) /= 0.0D+00 ) then
        v(1:n,k+1) = v(1:n,k+1) / h(k+1,k)
      end if

      if ( 1 < k ) then
        y(1:k+1) = h(1:k+1,k)
        do j = 1, k - 1
          call mult_givens ( c(j), s(j), j, y )
        end do
        h(1:k+1,k) = y(1:k+1)
      end if

      mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )

      c(k) = h(k,k) / mu
      s(k) = -h(k+1,k) / mu
      h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
      h(k+1,k) = 0.0D+00
      call mult_givens ( c(k), s(k), k, g )

      rho = abs ( g(k+1) )

      itr_used = itr_used + 1

      if ( verbose ) then
        write ( *, '(a,i4,a,g14.6)' ) '  K = ', k, '  Residual = ', rho
      end if

      if ( rho <= rho_tol .and. rho <= tol_abs ) then
        exit
      end if

    end do

    k = k_copy - 1

    y(k+1) = g(k+1) / h(k+1,k+1)

    do i = k, 1, -1
      y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
    end do

    do i = 1, n
      x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
    end do

    if ( rho <= rho_tol .and. rho <= tol_abs ) then
      exit
    end if

  end do

  if ( verbose ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMGMRES_ILU_CR:'
    write ( *, '(a,i6)' ) '  Iterations = ', itr_used
    write ( *, '(a,g14.6)' ) '  Final residual = ', rho
  end if

  return
end


subroutine rearrange_cr ( n, nz_num, ia, ja, a )

!*****************************************************************************80
!
!! REARRANGE_CR sorts a sparse compressed row matrix.
!
!  Discussion:
!
!    This routine guarantees that the entries in the CR matrix
!    are properly sorted.
!
!    After the sorting, the entries of the matrix are rearranged in such
!    a way that the entries of each column are listed in ascending order
!    of their column values.
!
!    The matrix A is assumed to be stored in compressed row format.  Only
!    the nonzero entries of A are stored.  The vector JA stores the
!    column index of the nonzero value.  The nonzero values are sorted
!    by row, and the compressed row vector IA then has the property that
!    the entries in A and JA that correspond to row I occur in indices
!    IA(I) through IA(I+1)-1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 July 2007
!
!  Author:
!
!    Original C version by Lili Ju.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
!    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
!    Charles Romine, Henk van der Vorst,
!    Templates for the Solution of Linear Systems:
!    Building Blocks for Iterative Methods,
!    SIAM, 1994.
!    ISBN: 0898714710,
!    LC: QA297.8.T45.
!
!    Tim Kelley,
!    Iterative Methods for Linear and Nonlinear Equations,
!    SIAM, 2004,
!    ISBN: 0898713528,
!    LC: QA297.8.K45.
!
!    Yousef Saad,
!    Iterative Methods for Sparse Linear Systems,
!    Second Edition,
!    SIAM, 2003,
!    ISBN: 0898715342,
!    LC: QA188.S17.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the system.
!
!    Input, integer ( kind = 4 ) NZ_NUM, the number of nonzeros.
!
!    Input, integer ( kind = 4 ) IA(N+1), the compressed row indices.
!
!    Input/output, integer ( kind = 4 ) JA(NZ_NUM), the column indices.
!    On output, these may have been rearranged by the sorting.
!
!    Input/output, real ( kind = 8 ) A(NZ_NUM), the matrix values.  On output,
!    the matrix values may have been moved somewhat because of the sorting.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nz_num

  real ( kind = 8 ) a(nz_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia(n+1)
  integer ( kind = 4 ) i4temp
  integer ( kind = 4 ) ja(nz_num)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) r8temp

  do i = 1, n

    do k = ia(i), ia(i+1) - 2
      do l = k + 1, ia(i+1) - 1

        if ( ja(l) < ja(k) ) then
          i4temp = ja(l)
          ja(l)  = ja(k)
          ja(k)  = i4temp

          r8temp = a(l)
          a(l)   = a(k)
          a(k)   = r8temp
        end if

      end do
    end do

  end do

  return
end

! !---------------
!     !Sort a compressed row matrix
!     subroutine rearrange_cr(n,nz_num,ia,ja,A)
!         integer, intent(in) :: n,nz_num
!         integer, intent(in) :: ia(n+1)
!         integer, intent(inout) :: ja(nz_num)
!         real(rp), intent(inout) :: A(nz_num)

!         integer :: i,itemp,k,l
!         real(rp) :: rtemp

!         do i=1,n
!             do k=ia(i),ia(i+1)-2
!                 do l=k+1,ia(i+1)-1
!                     if ( ja(l) < ja(k) ) then
!                         itemp = ja(l)
!                         ja(l) = ja(k)
!                         ja(k) = itemp

!                         rtemp = A(l)
!                         A(l) = A(k)
!                         A(k) = rtemp
!                     endif
!                 enddo
!             enddo
!         enddo

!     end subroutine rearrange_cr

end module mgmres
