module mgmres

    use kdefs
    implicit none

    contains


!---------------
    !A*x (compressed row)
    subroutine ax_cr(n,nz_num,ia,ja,A,x,W)
        integer , intent(in)  :: n,nz_num
        integer , intent(in)  :: ia(n+1),ja(nz_num)
        real(rp), intent(in)  :: A(nz_num)
        real(rp), intent(in)  :: x(n)
        real(rp), intent(inout) :: W(n)

        integer :: i,k1,k2

        W(1:n) = 0.0
        !TODO: Add OMP bindings here
        do i=1,n
            k1 = ia(i)
            k2 = ia(i+1) - 1
            W(i) = W(i) + dot_product( A(k1:k2), x(ja(k1:k2)) )
        enddo
    end subroutine ax_cr
!---------------
    subroutine diagonal_pointer_cr(n,nz_num,ia,ja,ua)
        integer, intent(in) :: n,nz_num
        integer, intent(in) :: ia(n+1),ja(nz_num)
        integer, intent(out) :: ua(n)

        integer :: i,k

        ua(1:n) = -1
        do i=1,n
            do k=ia(i),ia(i+1)-1
                if ( ja(k) == i ) then
                    ua(i) = k
                endif
            enddo
        enddo

    end subroutine diagonal_pointer_cr
    

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


!---------------
    !Sort a compressed row matrix
    subroutine rearrange_cr(n,nz_num,ia,ja,A)
        integer, intent(in) :: n,nz_num
        integer, intent(in) :: ia(n+1)
        integer, intent(inout) :: ja(nz_num)
        real(rp), intent(inout) :: A(nz_num)

        integer :: i,itemp,k,l
        real(rp) :: rtemp

        do i=1,n
            do k=ia(i),ia(i+1)-2
                do l=k+1,ia(i+1)-1
                    if ( ja(l) < ja(k) ) then
                        itemp = ja(l)
                        ja(l) = ja(k)
                        ja(k) = itemp

                        rtemp = A(l)
                        A(l) = A(k)
                        A(k) = rtemp
                    endif
                enddo
            enddo
        enddo

    end subroutine rearrange_cr

end module mgmres
