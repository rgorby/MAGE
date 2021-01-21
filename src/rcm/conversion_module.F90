MODULE conversion_module
  USE rcm_precision
  USE rcmdefs
  USE kdefs, ONLY : qp,PI
  IMPLICIT NONE
  REAL(rprec), ALLOCATABLE :: bndloc_old(:),almmin(:),almmax(:),almdel(:),&
      eta_midnight(:)
  REAL(rprec), ALLOCATABLE :: x0(:,:),y0(:,:),z0(:,:)
  REAL(rprec), ALLOCATABLE :: x0_sm(:,:),y0_sm(:,:),z0_sm(:,:)
  REAL(rprec), ALLOCATABLE :: te(:,:),ti(:,:),to(:,:),&
      den(:,:),press(:,:),&
      deno(:,:),presso(:,:),&
      beta_average(:,:),wImag(:,:)

  REAL(rprec), ALLOCATABLE :: eeta_new(:,:,:)
  INTEGER(iprec), ALLOCATABLE :: iopen(:,:),imin_j_old(:),inner_bndy(:)

  contains

    !Calculates difference of erfs - diff of exps, i.e. Eqn B5 from Pembroke+ 2012
    function erfexpdiff(A,x,y) result(z)

      real(rp), intent(in) :: A,x, y
      real(rp) :: z

      !QUAD precision holders
      real(qp) :: xq,yq,zq,differf,diffexp

      xq = x
      yq = y
      !Replacing erf(x)-erf(y) w/ erfc to avoid flooring to zero
      differf = erfc(yq)-erfc(xq)
      diffexp = xq*exp(-xq**2.0) - yq*exp(-yq**2.0)
      diffexp = 2.0*diffexp/sqrt(PI)

      zq = A*(differf-diffexp)
      z = zq

    end function erfexpdiff

END MODULE conversion_module
