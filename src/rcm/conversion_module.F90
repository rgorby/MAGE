MODULE conversion_module
  USE rcm_precision
  USE rcmdefs
  USE kdefs, ONLY : qp,QPI,PI
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

  !Quad prec. parameters for erf difference
  real(qp), parameter, private :: p  =  0.3275911  , &
                                  a1 =  0.254829592, &
                                  a2 = -0.284496736, &
                                  a3 =  1.421413741, &
                                  a4 = -1.453152027, &
                                  a5 =  1.061405429
  contains

    !Calculates difference of erfs - diff of exps, i.e. Eqn B5 from Pembroke+ 2012
    function erfexpdiff(A,xp,xm) result(eta)
      real(rprec), intent(in) :: A,xp,xm
      real(rprec) :: eta

      real(qp) :: qp,qm,tp,tm,ep,em,erfdiff,expdiff,etaq
      qp = xp
      qm = xm


      tp = 1.0/(1.0+p*qp)
      tm = 1.0/(1.0+p*qm)
      ep = exp(-(qp**2.0))
      em = exp(-(qm**2.0))

      !Difference of erf's using Abramowitz & Stegun, 7.1.26
      !erfdiff(qp,qm) = erf(qp)-erf(qm)
      !erf(x) ~ 1 - (a1.t + a2.t^2 + a3.t^3 + a4.t^5 + a5.t^5)*exp(-x^2) + eps(x)
      !t = 1/(1+px)
      !|eps(x)| <= 1.5e-7
      !Explicitly remove canceling 1's for difference
      erfdiff = - (a1*tp + a2*(tp**2.0) + a3*(tp**3.0) + a4*(tp**4.0) + a5*(tp**5.0))*ep &
                + (a1*tm + a2*(tm**2.0) + a3*(tm**3.0) + a4*(tm**4.0) + a5*(tm**5.0))*em

      expdiff = 2.0*(qp*ep - qm*em)/sqrt(QPI) !Using quad prec PI
      etaq = A*(erfdiff-expdiff)
      eta = etaq !Cast back down to rp
    end function erfexpdiff


    ! function olderfexpdiff(A,x,y) result(z)

    !   real(rp), intent(in) :: A,x, y
    !   real(rp) :: z

    !   !QUAD precision holders
    !   real(qp) :: xq,yq,zq,differf,diffexp

    !   xq = x
    !   yq = y
    !   !Replacing erf(x)-erf(y) w/ erfc to avoid flooring to zero
    !   differf = erfc(yq)-erfc(xq)
    !   diffexp = xq*exp(-xq**2.0) - yq*exp(-yq**2.0)
    !   diffexp = 2.0*diffexp/sqrt(PI)

    !   zq = A*(differf-diffexp)
    !   z = zq

    ! end function olderfexpdiff

END MODULE conversion_module
