MODULE conversion_module
  USE rcm_precision
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
END MODULE conversion_module
