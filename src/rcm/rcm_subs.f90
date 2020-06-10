!
    MODULE Rcm_mod_subs
    use rcm_precision
    IMPLICIT NONE
    SAVE
!
!
    INTEGER, PARAMETER :: LUN = 11, LUN_2 = 12, LUN_3 = 13
!    INTEGER, PARAMETER :: iprec = SELECTED_INT_KIND (9)
!    INTEGER, PARAMETER :: rprec = SELECTED_REAL_KIND (6,37)
! use gamera precision
!    INTEGER, PARAMETER :: iprec = ip
!    INTEGER, PARAMETER :: rprec = rp
!
!
    INCLUDE 'rcm_include.h'
!
!
!      Define a number of universal useful constants and parameters:
!      Part 1 is machine-dependent parameters and they should not be changed
!      under any circumstances.
!      Part 2 is physical constants; these may require editing, in which case
!      all the code must be recompiled.
!
    REAL (RPREC), PARAMETER ::     &
!                 Part 1: machine-specific and mathematical parameters
                               zero         = 0.0_rprec, &
                               one          = 1.0_rprec, &
                               two          = 2.0_rprec, &
                               three        = 3.0_rprec, &
                               four         = 4.0_rprec, &
                               five         = 5.0_rprec, &
                               six          = 6.0_rprec, &
                               eight        = 8.0_rprec, &
                               half         = 0.5_rprec, &
                               qtr          = 0.25_rprec,&
                               machine_eps1 = EPSILON (1.0_rprec), &
                               machine_eps2 = machine_eps1*10_rprec, &
                               machine_tiny = TINY (one),&
                               machine_huge = HUGE (one),&
                               pi_two       = two * pi, &
                               pi_by_two    = pi / two, &
                               rtd          = 180.0_rprec/pi, &
                               dtr          = pi/180.0_rprec, &
                               rth          = 12.0_rprec / pi,&
                               htr          = one / rth      ,&
!
!                 Part 2: physical constants          ! EDITING ALLOWED HERE
                               xmass (2)    = (/ 9.1E-31_rprec, &
                                                 1.67E-27_rprec /), &
                               besu         = 3.0584E+4_rprec, &
                               signbe       = one, &
                               romeca       = zero, &
                               charge_e     = 1.6E-19_rprec, &
                               sgn (ksize)  = one
                  INTEGER (iprec) :: ie_el = 1, ie_hd = 2 ! coding for e and proton
!
!
!   Potential solver GMRESM tolerance:
    REAL (rprec) :: tol_gmres
!
!
!
!
!   This is a definition of the label structure, for I/O:
    TYPE :: label_def
       INTEGER (iprec)   :: intg (20)
       REAL (rprec)      :: real (20)
       CHARACTER(LEN=80) :: char
    END TYPE label_def
    TYPE (label_def) :: label
!
!
!   Define an ellipse:
    TYPE :: ellipse_def
       REAL(rprec) :: aa, bb, xx, yy
    END TYPE ellipse_def
!
    TYPE (ellipse_def) :: boundary (2)
!
!
!   Grid info:
    REAL (rprec) :: dlam, dpsi, Ri, Re, &
                    alpha (isize, jsize), &
                    beta  (isize, jsize), &
                    colat (isize, jsize), &
                    aloct (isize, jsize), &
                    bir   (isize, jsize), &
                    sini  (isize, jsize), &
                    vcorot(isize, jsize), &
                    fac   (isize, jsize)
    INTEGER (iprec) :: i1, i2, iint, j1, j2, jint, imin, imin_j(jsize), ibnd_type
!
!
    LOGICAL ::  L_move_plasma_grid = .TRUE.
    LOGICAL ::  L_doKaiClaw        = .FALSE.
!
!
!   Plasma on grid:
    REAL (rprec) :: alamc (kcsize), etac (kcsize), fudgec (kcsize), &
                    eeta (isize,jsize,kcsize), eeta_cutoff, cmax, &
                    eeta_avg (isize,jsize,kcsize)
    INTEGER (iprec) :: ikflavc (kcsize), i_advect, i_eta_bc, i_birk
    LOGICAL :: L_dktime
    INTEGER (iprec), PARAMETER :: irdk=18, inrgdk=13, isodk=2, iondk=2
    REAL (rprec) :: dktime (irdk, inrgdk, isodk, iondk), sunspot_number
 
     logical :: kill_fudge
!
!
!   Magnetic field:
    REAL (rprec) :: xmin (isize,jsize), ymin (isize,jsize), zmin (isize,jsize), &
                    bmin (isize,jsize), vm (isize,jsize), &
                    rmin (isize,jsize), pmin(isize,jsize),&
                    x1 (isize,jsize), x2 (isize,jsize), &
                    y1 (isize,jsize), y2 (isize,jsize), &
                    b1 (isize,jsize), b2 (isize,jsize), &
                    vm1(isize,jsize), vm2(isize,jsize), &
                    bndloc (jsize)
    INTEGER (iprec), ALLOCATABLE :: ibtime (:)
    REAL    (rprec) :: fstoff, fclps, fdst, fmeb, ftilt
    INTEGER (iprec) :: itype_bf
         ! itype_bf = 1 -- read time sequence from input files and interpolate in time
         ! itype_bf = 2 -- read single B-field configuration from input files           
         ! itype_bf = 3 -- expect that an external program will assign B-field arrays
         !                 including pmin and rmin, so do nothing.
!
!
!   Ionospheric quantities:
    REAL (rprec) :: qtped (isize,jsize), pedpsi (isize,jsize), &
                    qtplam(isize,jsize), pedlam (isize,jsize), &
                    qthall(isize,jsize), hall   (isize,jsize), &
                    ss (jsize), &
                    pwe (isize,jsize), pwn (isize,jsize), &
                    hwe (isize,jsize), hwn (isize,jsize), &
                    sw  (jsize), &
                    eflux (isize,jsize,iesize), eavg (isize,jsize,iesize)
    INTEGER (iprec) :: icond, nsmthi, nsmthj, iwind
    LOGICAL :: ifloor, icorrect
!
!
!   Magnetospheric quantities:
    REAL (rprec) :: v (isize,jsize), vpar (isize,jsize), vbnd (jsize), &
                    birk (isize,jsize), pvgamma (isize,jsize,iesize), &
                    pressrcm (isize,jsize), &
                    v_avg (isize,jsize), birk_avg (isize,jsize), &
                    densrcm(isize,jsize),denspsph(isize,jsize)
    INTEGER (iprec) :: ipcp_type, ipot
!
!
!   Input PCP drop and its current value:
    INTEGER (iprec), ALLOCATABLE :: ivtime (:)
    REAL    (rprec), ALLOCATABLE :: vinput (:), vinput_phase(:)
    REAL    (rprec)              :: vdrop,      vdrop_phase
!
!
 
    INCLUDE 'rcmdir.h'

    Logical :: IsCoupledExternally = .false.  ! flag to determine if RCM is standalone or not

    ! Variables for internal RCM timing:
    INTEGER(iprec) :: timer_start(10) = 0, timer_stop(10) = 0, count_rate
    REAL (rprec) :: timer_values (10)=0.0_rprec

!
    INTERFACE Gntrp
       MODULE PROCEDURE Gntrp_2d_ang
    END INTERFACE
 
    INTERFACE Interp
      MODULE PROCEDURE Interp_1d, Interp_2d, Interp_2d_of3d
    END INTERFACE

    INTERFACE Bjmod
       MODULE PROCEDURE Bjmod_int, Bjmod_real
    END INTERFACE

    INTERFACE Circle
       MODULE PROCEDURE Circle_1d, Circle_2d
    END INTERFACE
!
!
!
    CONTAINS
!
!
!
      SUBROUTINE Comput (jtime, dt )
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime
      REAL (rprec),    INTENT (IN) :: dt
!
      INTEGER (iprec) :: j
      REAL (rprec)  ::  a(3), b(3), dx(3), dy(3), deqdt
!
!
      IF (IsCoupledExternally) then
         vdrop = (MAXVAL(vbnd) - MINVAL(vbnd))/1.0E+3
         vdrop_phase = 0.0
      ELSE
         vdrop = Get_vdrop    (ivtime,  vinput,  jtime)
         vdrop_phase = Get_vdrop_phase    (ivtime,  vinput_phase,  jtime)
      END IF
      IF (i_eta_bc == 1) THEN
!        DO NOTHING
      ELSE IF (i_eta_bc == 2) THEN
!        DO NOTHING
      ELSE 
         STOP 'ILLEGAL VALUE OF I_eta_bc'
      END IF
      IF (i_birk == 1) THEN
         CALL Get_jbirk
      ELSE IF (i_birk == 2) THEN
         stop 'do not use'
      ELSE IF (i_birk ==3) THEN
         CALL Get_jbirk2

      ELSE
          STOP 'ILLEGAL VALUE OF BIRK'
      END IF
      CALL Get_vparallel ()
!
!
      IF (ibnd_type == 4) THEN
!        DO NOTHING, VBND IS ALREADY SET
         DO j = 1, jsize
            v (1:imin_j(j)-1,j) = vbnd (j)
         END DO
      ELSE
         STOP 'COMPUT: ibnd_type not implemented'
      END IF
!
!
      RETURN 
      END SUBROUTINE Comput                         
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      FUNCTION Get_vdrop (ivtime, vinput, jtime)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime, ivtime(:)
      REAL (rprec), INTENT (IN) :: vinput (:)
      REAL (rprec) :: Get_vdrop
!                                                                       
!-------------------------------------------------------------
!     Subroutine to specify total cross-polar-cap potential drop
!     vdrop (in kV) at time jtime.  this is accomplished by
!     interpolating vinput in time.
!     rws     3/20/97
!     If jtime <= ivtime(1) then vdrop = vinput(1)
!     If jtime >  ivtime(nvmax) then vdrop = vinput(nvmax)
!     all other cases--interpolated.
!-------------------------------------------------------------
!
      INTEGER (iprec) :: nv, nvmax
      REAL (rprec)    :: f
!
      nvmax = SIZE (vinput)
      DO nv = 1, nvmax 
         IF (jtime <= ivtime (nv) ) THEN 
            IF (nv == 1) THEN 
               Get_vdrop = vinput (1)
               RETURN 
            ELSE 
               f = REAL(jtime-ivtime(nv-1),rprec) / &
                   REAL(ivtime(nv)-ivtime(nv-1), rprec)
               Get_vdrop = (one - f) * vinput(nv-1) + f * vinput(nv)
               RETURN 
            END IF 
         END IF 
      END DO 
      Get_vdrop = vinput (nvmax)
!
      RETURN 
      END FUNCTION Get_vdrop
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      FUNCTION Get_vdrop_phase (ivtime, vinput_phase, jtime)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime, ivtime(:)
      REAL (rprec), INTENT (IN) :: vinput_phase (:)
      REAL (rprec) :: Get_vdrop_phase
!                                                                       
!-------------------------------------------------------------
!     Subroutine to specify total cross-polar-cap potential drop
!     vdrop (in kV) at time jtime.  this is accomplished by
!     interpolating vinput in time.
!     rws     3/20/97
!     If jtime <= ivtime(1) then vdrop = vinput(1)
!     If jtime >  ivtime(nvmax) then vdrop = vinput(nvmax)
!     all other cases--interpolated.
!-------------------------------------------------------------
!
      INTEGER (iprec) :: nv, nvmax
      REAL (rprec)    :: f
!
      nvmax = SIZE (vinput_phase)
      DO nv = 1, nvmax 
         IF (jtime <= ivtime (nv) ) THEN 
            IF (nv == 1) THEN 
               Get_vdrop_phase = vinput_phase (1)
               RETURN 
            ELSE 
               f = REAL(jtime-ivtime(nv-1),rprec) / &
                   REAL(ivtime(nv)-ivtime(nv-1), rprec)
               Get_vdrop_phase = (one - f) * vinput_phase(nv-1) + f * vinput_phase(nv)
               RETURN 
            END IF 
         END IF 
      END DO 
      Get_vdrop_phase = vinput_phase (nvmax)
!
      RETURN 
      END FUNCTION Get_vdrop_phase
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_jbirk ( )
      IMPLICIT NONE
!__________________________________________________________________________
!                                                                       
!  Program written by: r.w. spiro        
!  last update:
!     04-05-88          
!     01-29-96 frt                 - added ain,min_j arr
!  Algorithm by: r.a. wolf                                              
!                                                                       
!  This subroutine computes birk(i,j) given inner edge
!  locations 
!  modified 04-05-88 to include effects of gradients in eta.            
!  see raw document re including eeta in computation of jbirk           
!    dated feb 6, 1988.                                                 
!  birk is current density (2 hemispheres) in units of
!  microamp/m**2    
!
!  birk(i,j) here is the field-aligned current density per
!  unit of ionospheric area, so that it already includes
!  the factor sin(I); this is J_parallel*sin(I) in the RHS
!  of the Vasyliunas equation.
!
!  Issues with non-integer boundary (Stanislav's notes):
!  for BIRK from inner edge segments, this is not an issue
!  (except that if a segment is entirely outside the bndry,
!  then we don't compute its contribution); of course, we 
!  have to care about this somewhere else where motion of 
!  test particles is computed. For BIRK from gradients of 
!  EETA,  
!  removed edges 3/19 frt
!                                                                       
!______________________________________________________________________________
!
      REAL (rprec), PARAMETER :: cf1 = one / pi_two, &
                                 cf2 =  - (three/four)*( (two / pi) - half)
!
      INTEGER (iprec) :: i, j, k, kc, klbeg, klend, kl, klnext, &
                 ibmin, ibmax, jbmin, jbmax, jb1, jb2,  &
                 ig, jj, jindex, ib1, ib2
      REAL (rprec) :: detadi(isize,jsize), detadj(isize,jsize), &
                      dvmdi(isize,jsize), dvmdj(isize,jsize), dbirk, &
                      vmkl, vmnext, sum, b1, b2, x, y, el, umax, umin, ss, &
                      z, dg1, dg2, dg3, qmin, qmax, qn, qx,  &
                      denom, a1, a2, bjm, range, bim, gkl (5000)
!                                                                       
!
      birk (:,:) = zero
!
!
!     Compute J_parallel due to continuous channel:
!                                                                       
      dvmdi  = Deriv_i (vm, imin_j)
      dvmdj  = Deriv_j (vm, imin_j, j1, j2, 1.0E+25_rprec)
      WHERE (ABS(dvmdj) > 1.0E+24)  ! to prevent artificial inflows on bndy
          dvmdi = 0.0
          dvmdj = 0.0
      END WHERE
!
      DO kc = 1, kcsize
!
         detadi = Deriv_i (eeta (:,:,kc), imin_j)
         detadj = Deriv_j (eeta (:,:,kc), imin_j, j1, j2, 1.0E+32_rprec)
         WHERE (ABS(detadj) > 1.0E+31)
           detadi = 0.0
           detadj = 0.0
         END WHERE
!
         DO  j = j1, j2 
            DO  i = imin_j(j)+1, i2
!              IF (i < imin_j(j) + 3) CYCLE
!              IF (i <= imin_j(j-1) .or. i<=imin_j(j-2).or.&
!                 i<=imin_j(Bjmod_int(j-3,jwrap,jsize)) ) CYCLE
!              IF (i <= imin_j(j+1) .or. i<=imin_j(Bjmod_int(j+2,jwrap,jsize))&
!                  .or. i<=imin_j(Bjmod_int(j+3,jwrap,jsize))) CYCLE
               dbirk  = charge_e * signbe * ABS(alamc(kc)) * &
                        (detadj(i,j) * dvmdi(i,j) - detadi(i,j)*dvmdj(i,j)) / &
                        (alpha(i,j)*beta(i,j)*dlam*dpsi*Ri**2)
               birk (i, j) = birk (i, j) + dbirk 
            END DO 
         END DO 
      END DO 
!     pause
!     print *,birk
!        DO  j = j1, j2 
!           DO  i = imin_j(j)+1, i2
!              birk (i, j) = MIN (birk(i,j),10.)
!              birk (i, j) = MAX (birk(i,j),-10.)
!        END DO 
!     END DO 
!                                                                       
      CALL Circle (birk)
!
      RETURN 
      END SUBROUTINE Get_jbirk
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_jbirk2 ( )
      IMPLICIT NONE
!__________________________________________________________________________
!                                                                       
!  Program written by: r.w. spiro        
!  last update:
!     04-05-88          
!     01-29-96 frt                 - added ain,min_j arr
!     11-14-02 frt uses a new version on the grid-based scheme
!  Algorithm by: r.a. wolf                                              
!                                                                       
!  This subroutine computes birk(i,j) given inner edge
!  locations 
!  modified 04-05-88 to include effects of gradients in eta.            
!  see raw document re including eeta in computation of jbirk           
!    dated feb 6, 1988.                                                 
!  birk is current density (2 hemispheres) in units of
!  microamp/m**2    
!
!  birk(i,j) here is the field-aligned current density per
!  unit of ionospheric area, so that it already includes
!  the factor sin(I); this is J_parallel*sin(I) in the RHS
!  of the Vasyliunas equation.
!
!  Issues with non-integer boundary (Stanislav's notes):
!  for BIRK from inner edge segments, this is not an issue
!  (except that if a segment is entirely outside the bndry,
!  then we don't compute its contribution); of course, we 
!  have to care about this somewhere else where motion of 
!  test particles is computed. For BIRK from gradients of 
!  EETA,  
! 
!  This version uses a new way to compute jbirk from the cts channel - frt
!                                                                       
!______________________________________________________________________________
!
      REAL (rprec), PARAMETER :: cf1 = one / pi_two, &
                                 cf2 =  - (three/four)*( (two / pi) - half)
!
      INTEGER (iprec) :: i, j, k, kc, klbeg, klend, kl, klnext, &
                 ibmin, ibmax, jbmin, jbmax, jb1, jb2,  &
                 ig, jj, jindex, ib1, ib2
      REAL (rprec) :: detadi(isize,jsize), detadj(isize,jsize), &
                      dvmdi(isize,jsize), dvmdj(isize,jsize), dbirk, &
                      vmkl, vmnext, sum, b1, b2, x, y, el, umax, umin, ss, &
                      z, dg1, dg2, dg3, qmin, qmax, qn, qx,  &
                      denom, a1, a2, bjm, range, bim, gkl (5000)
      REAL (rprec) :: eeta2(isize,jsize,ksize),vm2(isize,jsize)
!                                                                       
!
      birk (:,:) = zero
!
!
!     Compute J_parallel due to continuous channel:
!                                                                       
! define new temporary work arrays eeta2 and vm2
      eeta2 = eeta
      vm2   = vm
      DO j=1,jsize
       DO i=1,imin_j(j)-1
        DO k=1,kcsize
         eeta2(i,j,k) = eeta(imin_j(j),j,k)
        END DO
         vm2(i,j) = vm(imin_j(j),j)
       END DO
      END DO

!    1          2            3
! i-1,j+1------i,j+1------i+1,j+1
!    |          |            |  
!    |          |            |  
!    |          |            |  
!    4          |            5
! i-1,j--------i,j--------i+1,j
!    |          |            |  
!    |          |            |  
!    |          |            |  
!    6          7            8
! i-1,j-1------i,j-1------+1,j-1 
! the basic equation is
! dbirk  is proportional to 
! (eta_1*(vm_2-vm_1)+eta_2*(vm_3-vm_1)+eta_3*(vm_5-vm_2)+eta_4*(vm_1-vm_6)
!  eta_5*(vm_8-vm_3)+eta_6*(vm_4-vm_7)+eta_7*(vm_6-vm_8)+eta_8*(vm_7-vm_5))/8
      DO kc = 1, kcsize
!
         DO  j = j1, j2 
           DO  i = imin_j(j)+1, i2

            dbirk =(eeta2(i-1,j+1,kc)*(vm2(i  ,j+1)-vm2(i-1,j  )) + &
                    eeta2(i  ,j+1,kc)*(vm2(i+1,j+1)-vm2(i-1,j+1)) + &
                    eeta2(i+1,j+1,kc)*(vm2(i+1,j  )-vm2(i  ,j+1)) + &
                    eeta2(i-1,j  ,kc)*(vm2(i-1,j+1)-vm2(i-1,j-1)) + &
                    eeta2(i+1,j  ,kc)*(vm2(i+1,j-1)-vm2(i+1,j+1)) + &
                    eeta2(i-1,j-1,kc)*(vm2(i-1,j  )-vm2(i  ,j-1)) + &
                    eeta2(i  ,j-1,kc)*(vm2(i-1,j-1)-vm2(i+1,j-1)) + &
                    eeta2(i+1,j-1,kc)*(vm2(i  ,j-1)-vm2(i+1,j  ))) /8.

               dbirk  = charge_e * signbe * ABS(alamc(kc)) * dbirk/ &
                        (alpha(i,j)*beta(i,j)*dlam*dpsi*Ri**2)
               birk (i, j) = birk (i, j) + dbirk 
           END DO 
         END DO 
      END DO 
!                                                                       
!print*,'zeroing birk'
      CALL Circle (birk)
!
      RETURN 
      END SUBROUTINE Get_jbirk2
!
!
!*************************************************************************

!
!
      SUBROUTINE Get_vparallel ()
      IMPLICIT NONE
!______________________________________________________________________________
!  last update: 
!     05-05-87       by:rws                              
!     02-10-96          frt - added arrays ain,min_j     
!                                                                       
!  Birk is sum of current densities into both hemispheres.              
!  (micro amp/m**2).  Before activating parallel potential drop        
!  we need to check if birk is being used correctly in
!  this routine.   
!
!  Stanislav: VPAR is computed inside IE loop (for both
!             negative and positive particles), and will
!             be the one for the largest IE value. Which
!             is nonsense.
!  Stanislav: this subroutine needs grid-based formulation
!             of plasma (EETA). Before it was done by
!             computing EETA for electrons from the inner
!             edges of electrons, then it was changed to
!             use directly grid-based population. In the
!             latter case, array PVEC returned by this
!             routine is the electron pressure (without
!             the factor of 2/3) and is the same as what
!             routine PV returns as array PVGAM. If the
!             electrons are on the grid only, as in my case,
!             then we call PV in rcm main program to compute
!             the ion pressure, and we use PVEC from this
!             routine for the electron pressure. (04/20/99)
!  Stanislav, may 18,99: make all loops over electrons only,
!             by using iedim_local and setting it to 1.
!______________________________________________________________________________
!
      INTEGER (iprec) :: i, j, ie, iedim_local, kc
      REAL (rprec)    :: en, ekt, therm, sum1 (iesize), sum2 (iesize)
!                                                                       
!
!                                  
      iedim_local = 2
!
      vpar  (:,:)   = zero
      eavg  (:,:,:) = zero
      eflux (:,:,:) = zero
!
      loop_j: DO j = j1, j2
      loop_i: DO i = imin_j(j), isize
!
!           For each grid point, clear sum1 and sum2:
!
            sum1 (1:iedim_local) = zero
            sum2 (1:iedim_local) = zero
!
!
!           Now for each grid point, consider all species
!           present at that grid point, and compute sum1 and
!           sum2 for positive and negative particles separately:
!
            GRID_BASED: DO kc = 1, kcsize
             IF ( ABS(alamc(kc))*vm(i,j) > 500.0_rprec) THEN
               IF (alamc (kc) < zero) THEN
                  ie = 1 
               ELSE 
                  ie = 2 
!                 STOP 'BALGN4: ie is 2'
               END IF
               sum1(ie) = sum1(ie) + eeta(i,j,kc)*fudgec(kc)
               sum2(ie) = sum2(ie) + eeta(i,j,kc)*fudgec(kc)*ABS(alamc(kc))
             END IF
            END DO GRID_BASED 
!
!           For positive and negative particles separately,
!           compute precipitating number flux, average energy,
!           and parallel potential drop:
!
            DO ie = 1, iedim_local 
!                                                                       
               IF (sum1 (ie) > 10.*machine_tiny) THEN  ! zero  sbao 07/2019
!
!                compute thermal electron current, field-aligned
!                potential drop, electron energy flux,
!                and average electron energy at (i,j):          
!
                  en    = sum1 (ie) * vm (i, j)**1.5 / 6.38E+21
                  ekt   = (two/three) * sum2 (ie) * vm (i,j) / sum1 (ie)
                  therm = 0.02675 * en * SQRT(ekt*xmass(1)/xmass(ie))
!
                  IF (therm < 1.E-30) THEN 
                     therm      = zero
                     vpar (i,j) = zero
                  ELSE 
                     IF (- birk (i, j) / therm > one) THEN
                        vpar (i,j) = ekt * (- birk (i,j) / therm - one)
                     ELSE 
                        vpar (i,j) = one
                     END IF
                     vpar(i,j) = MIN (vpar (i, j), 10000.0_rprec)
                  END IF 
!
!    !!!!!!!      ALERT: VPAR(I,J) IS SET TO 0 !!!!!!!!!!!!!!!!!
!
                  vpar (i, j) = zero
!
!
                  eflux(i,j,ie) = 0.002 * therm * &
                                 ( ekt + vpar(i,j) + half*vpar(i,j)**2/ekt)
                  eavg(i,j,ie) = two*(ekt+vpar(i,j)+half*vpar(i,j)**2 /ekt) / &
                                 (one + vpar (i, j) / ekt)
                  ! sbao 6/19 detect Nan 
                  if (ISNAN(eflux(i,j,ie)))then
                       if (.not. doQuietRCM) write(*,*)'eflux,i,j,therm,ekt,vpar,sum1,sum2,vm',eflux(i,j,ie),i,j,therm,ekt,vpar(i,j),sum1(ie),sum2(ie),vm(i,j)
                       eflux(i,j,ie) = 0.0
                       eavg(i,j,ie) = 0.0
                  end if

               ELSE 
!                                                                       
!                 Case fudge=0: we want eflux=0 and eavg=0 for no precipitation.
!
                  eflux (i, j, ie) = zero
                  eavg  (i, j, ie) = zero
!
               END IF 
               ! corrections to eavg at eflux(i,j) == 0.0     sbao 07/2019 
               ! == does not work well with real number, use lt threshold instead. ldong 04/2020
               IF (eflux(i,j,ie) .lt. 0.01) eavg(i,j,ie) = 0.0
               IF (eavg(i,j,ie) .lt. 0.01) eflux(i,j,ie) = 0.0

!                                                                       
            END DO
!
      END DO loop_i
      END DO loop_j 
!                                                                       
!
!
      CALL Circle (vpar)
      CALL Circle (eflux (:, :, ie_el))
      CALL Circle (eavg  (:, :, ie_el))
      CALL Circle (eflux (:, :, ie_hd))
      CALL Circle (eavg  (:, :, ie_hd))
!
      RETURN
      END SUBROUTINE Get_vparallel
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Floor_for_eflux ()
      IMPLICIT NONE
      INTEGER (iprec) :: ivalue_max, i, j
      REAL (rprec)    :: eflux_max
      DO j = 1, jsize
         eflux_max = eflux (isize, j, ie_el)
         ivalue_max = isize
         DO i = isize-1, imin_j(j), -1
            IF (eflux(i,j,ie_el) > eflux(i+1,j,ie_el)) THEN
               eflux_max  = eflux(i,j,ie_el)
               ivalue_max = i
            END IF
         END DO
         DO i = imin_j(j), ivalue_max - 1
            eflux(i,j,ie_el) = MAX (half*eflux_max, eflux(i,j,ie_el))
         END DO
      END DO
      !ion precipitation
      DO j = 1, jsize
         eflux_max = eflux (isize, j, ie_hd)
         ivalue_max = isize
         DO i = isize-1, imin_j(j), -1
            IF (eflux(i,j,ie_hd) > eflux(i+1,j,ie_hd)) THEN
               eflux_max  = eflux(i,j,ie_hd)
               ivalue_max = i
            END IF
         END DO
         DO i = imin_j(j), ivalue_max - 1
            eflux(i,j,ie_hd) = MAX (half*eflux_max, eflux(i,j,ie_hd))
         END DO
      END DO
      RETURN
      END SUBROUTINE Floor_for_eflux
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
    SUBROUTINE Move_plasma ( dt )
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: dt
!_____________________________________________________________________________
!
!  Time step subroutine to do simple euler time step                    
!                                                                       
!  Last update:
!   8-29-86                                                 
!   1-29-96 frt added boundary arrays and calls to bndy     
!   3-19-97 rws ibtime and nbf added as calling parameters  
!   10-02-98 sts fudge is sized as kcdim for electrons on grid
!   may 99 sts removed hardy coeffs--they are in module
!_____________________________________________________________________________
!
!
  
  IF (L_move_plasma_grid) THEN
    IF (i_advect == 1) THEN
       CALL Move_plasma_grid  (dt, 1_iprec, isize, j1, j2, 1_iprec)
       CALL Move_plasma_grid  (dt, 1_iprec, isize, j1, j2, 2_iprec)
    ELSE IF (i_advect == 2) THEN
!      CALL Move_plasma_grid (dt, 1, isize, j1, j2, 1)
       STOP 'This option is no longer available, aborting RCM'
    ELSE IF (i_advect == 3) THEN
        IF (L_doKaiClaw) THEN
          CALL Move_plasma_grid_KAIJU (dt)
        ELSE
          CALL Move_plasma_grid_new (dt)
        ENDIF
    ELSE
       STOP 'ILLEGAL I_ADVECT IN MOVING PLASMA'
    END IF
  END IF
!
    RETURN
    END SUBROUTINE Move_plasma
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Move_plasma_grid (dt, i_start, i_stop, j_start, j_stop,ie_ask)
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: dt
    INTEGER (iprec), INTENT (IN) :: ie_ask, i_start, i_stop, j_start, j_stop
!_____________________________________________________________________________
!   Subroutine to advance eta distribution for a time step
!   by a lifetime-based algorithm (raw doc dated 5/12/87)
!                                                                       
!   Last update: 05-11-88
!                01-29-96 ain ,min_j and calls to bndy added - frt
!   rws          06-05-97 etamov changed to reflect new use of
!                        eeta array in rcm697 version
!
!   CALLED FROM: TSTEP1
!_____________________________________________________________________________
!
!
    REAL (rprec) :: eeta2 (isize,jsize), veff  (isize,jsize), &
                    dvefdi(isize,jsize), dvefdj(isize,jsize), &
                    didt, djdt, biold, bjold, rate, mass_factor, a1, a2, fi, fj
    INTEGER (iprec) :: i, j, kc, ie, i_1, i_2, j_1, j_2
    LOGICAL :: pt_1_1, pt_1_2, pt_2_1, pt_2_2
!
real :: v_1_1, v_1_2, v_2_1, v_2_2
!
    DO kc = 1, kcsize
!
       IF (alamc(kc) < zero) THEN
          ie = 1
       ELSE
          ie = 2 ! but must change if O+ is added
       END IF
       IF (ie /= ie_ask) CYCLE
       mass_factor = SQRT (xmass(1) / xmass(ie))
       veff = v + vcorot - vpar  + alamc(kc)*vm
!
       dvefdi = Deriv_i (veff, imin_j)
       dvefdj = Deriv_j (veff, imin_j, j1, j2, 1.0E+25_rprec)
       WHERE (ABS(dvefdj) > 1.0E+24)
          dvefdj = 0.0
          dvefdi = 0.0
       END WHERE
!
!
       eeta2 (:,:) = eeta (:,:,kc)
       DO j = j1, j2
!      DO i = imin_j(j), i2
       DO i = i_start, i_stop
          IF (i <= imin_j(j)) CYCLE
          didt   =   dvefdj (i,j) / fac (i,j)
          djdt   = - dvefdi (i,j) / fac (i,j)
          biold  = REAL(i,rprec) - didt * dt 
          bjold  = Bjmod (REAL(j,rprec) - djdt * dt, jwrap, jsize )
          rate   = Ratefn (fudgec(kc), alamc(kc), sini (i,j), bir (i,j), &
                           vm (i,j), mass_factor)
!          IF (biold > Bndy(bndloc,bjold)) THEN 
          IF (biold > bndloc(j)) THEN 
!            Particle came from within the modeling region, find ETA_old by interp:
             i_1 = INT (biold)
             i_2 = i_1 + 1
             j_1 = INT (bjold)
             j_2 = j_1 + 1
!
!   (i_1,j_1) x--------------x (i_1,j_2)
!             |              |
!             |              |
!             |              |
!             |              |
!   (i_2,j_1) x--------------x (i_2,j_2)
!
             pt_1_1 = i_1 >= imin_j(j_1)
             pt_1_2 = i_1 >= imin_j(j_2)
             pt_2_1 = i_2 >= imin_j(j_1)
             pt_2_2 = i_2 >= imin_j(j_2)
!
             v_1_1 = eeta2 (i_1,j_1)
             v_1_2 = eeta2 (i_1,j_2)
             v_2_1 = eeta2 (i_2,j_1)
             v_2_2 = eeta2 (i_2,j_2)
!
             IF ((.NOT.pt_2_1) .OR. (.NOT. pt_2_2)) THEN
                 STOP 'ONE OF I_2 POINTS OUT OF BNDY'
             ELSE IF ((.NOT.pt_1_1) .AND. (.NOT.pt_1_2)) THEN
                 STOP 'BOTH I_1 POINTS OUT OF BNDY'
             ELSE IF (pt_1_1) THEN ! get 1,2 by interp.
                 
             ELSE ! get 1,1 by interp.
                 
             END IF
!            IF (i_1 < imin_j(j_1) .OR. i_2 < imin_j(j_1) .OR. &
!                i_1 < imin_j(j_2) .OR. i_2 < imin_j(j_2)) THEN
!              STOP 'PNT IS OUTSIDE MODELING REGION'
!            END IF
             fi = REAL(i_1)
             fj = REAL(j_1)
             a1 = (1.0-(biold-fi))*eeta2(i_1,j_1) + (biold-fi)*eeta2(i_2,j_1)
             a2 = (1.0-(biold-fi))*eeta2(i_1,j_2) + (biold-fi)*eeta2(i_2,j_2)
             eeta (i,j,kc) = ((1.0 - (bjold-fj)) * a1 + (bjold-fj) * a2)*EXP(-rate*dt)
          ELSE
!            Particle came from outside the modeling region, find ETA_old from b.c.:
             STOP 'PARTICLE IS FLOWING IN, NOT IMPLEMENTED'
          END IF
       END DO
       END DO
!
       CALL Circle (eeta(:,:,kc))
!
    END DO
!
    RETURN
!
    CONTAINS
!
    FUNCTION Ratefn (fudgx, alamx, sinix, birx, vmx, xmfact)
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: fudgx,alamx,sinix,birx,vmx,xmfact
    REAL (rprec)              :: Ratefn
!                                                                       
!   Function subprogram to compute precipitation rate
!   Last update:  04-04-88
!
    Ratefn = 0.0466_rprec*fudgx*SQRT(ABS(alamx))*(sinix/birx)*vmx**2
    Ratefn = xmfact * ratefn
    RETURN
    END FUNCTION Ratefn
    END SUBROUTINE Move_plasma_grid
!
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Circle_2d (r)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN OUT) :: r (:,:)
!
      INTEGER (iprec) :: jlast, i, j, imax, jmax
!
      imax = SIZE (r, DIM = 1)
      jmax = SIZE (r, DIM = 2)
      jlast = jmax - jwrap 
      DO i = 1, imax 
        DO  j = 1, jwrap - 1
          r (i, j) = r (i, jlast + j)
        END DO
        r (i, jmax) = r (i, jwrap) 
      END DO
      RETURN 
      END SUBROUTINE Circle_2d
!
!
!
!
!
      SUBROUTINE Circle_1d (r)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN OUT) :: r (:)
!
      INTEGER (iprec) :: jlast, j, jmax
!
      jmax = SIZE (r, DIM = 1)
      jlast = jmax - jwrap 
      DO  j = 1, jwrap - 1
        r (j) = r (jlast + j)
      END DO
      r (jmax) = r (jwrap)
      RETURN
      END SUBROUTINE Circle_1d
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    FUNCTION Eta_lambda_vgamma ( kbeg, kend, kcbeg, kcend, gamma)  RESULT (out)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: kbeg, kend, kcbeg, kcend
    REAL (rprec), INTENT (IN) :: gamma
    REAL (rprec) :: out (isize,jsize,iesize)
!______________________________________________________________________________
!
!   Subroutine computes quantity ETA*ABS(ALAM)*V**GAMMA at each grid point
!   for electrons and ions separately. KBEG, KEND, KCBEG, KCEND can be used
!   to restrict species to be included in the sum (use 1, ksize, 1, kcsize for
!   no restrictions). GAMMA is an input parameter:
!   ** if GAMMA = 0, then the computed sum is the adiabatic parameter PVGAMMA
!   ** if GAMMA = -5/3, compute energy density (or pressure without the 2/3 factor)
!   ** if GAMMA = -2/3, compute total energy of particles
!______________________________________________________________________________
!
    INTEGER (iprec) :: nbi, k, m, mbeg, mend, ipmax, ncount, i,j,n, ie,kc
    REAL (rprec) :: q, bimax, bicrss (100), charge
!
!
    out = zero
!
    DO ie = 1, iesize
!
        IF (ie == 1) THEN
           charge = - one
        ELSE
           charge = + one
        END IF
!
!       I. Compute sum for plasma on inner edges, electrons and ions separately :

!
!   II. Compute the sum for grid_based electrons or ions:
!
        DO j = 1, jsize
        DO i = 1, isize
!           IF (REAL(i,rprec) < Bndy(bndloc, REAL(j,rprec)) ) CYCLE
           IF (REAL(i,rprec) < bndloc(j) ) CYCLE
           DO kc = kcbeg, kcend
              q = alamc(kc) / charge
              IF (q > zero) THEN
                 out (i,j,ie) = out (i,j,ie) + ABS (alamc(kc) * eeta(i,j,kc))
              END IF
           END DO
        END DO
        END DO
!
    END DO
!
!
    DO ie = 1, iesize
       CALL Circle (out (:,:,ie))
    END DO
!
!
    DO ie = 1, iesize
       out (:,:,ie) = out(:,:,ie) * (vm (:,:)**((-three/two)*gamma))
    END DO
!
    RETURN
    END FUNCTION Eta_lambda_vgamma
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_qtcond ()
      IMPLICIT NONE
      INTEGER (iprec) :: n, i, j
      CHARACTER (LEN=80) :: form_string
      LOGICAL, SAVE :: called_already = .FALSE.
!
      IF (called_already) RETURN
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FORM = 'FORMATTED', &
            FILE = rcmdir//'rcmcond', ACTION = 'READ') 
!
        READ (LUN, '(I10.10)') n
        IF (n /= isize*jsize) STOP 'sizes do not match in qtcond'
        READ (LUN,'(A80)') form_string
        DO j = 1, jsize
        DO i = 1, isize
           READ (LUN,form_string) qtplam(i,j), qthall(i,j), qtped(i,j)
        END DO
        END DO
!
!
        READ (LUN, '(I10.10)') n
        IF (n /= jsize) STOP 'sizes do not match in qtcond'
        READ (LUN,'(A80)') form_string
        DO j = 1, jsize
           READ (LUN,form_string) ss(j)
        END DO
!
      CLOSE (LUN)
      called_already = .TRUE.
      RETURN
      END SUBROUTINE Read_qtcond
!
!
!
!
      SUBROUTINE Write_qtcond
      IMPLICIT NONE
      INTEGER (iprec) :: n, i, j
      CHARACTER (LEN=80) :: form_string
!
      OPEN (LUN, FILE = rcmdir//'rcmcond', FORM = 'FORMATTED', STATUS = 'REPLACE')
!
        form_string = '(3(TR2,ES23.15))'
        WRITE (LUN,'(I10.10)') SIZE(qtplam)
        WRITE (LUN,'(A80)') form_string
        DO j = 1, jsize
        DO i = 1, isize
           WRITE (LUN,form_string) qtplam(i,j), qthall(i,j), qtped(i,j)
        END DO
        END DO
!
        form_string = '(1(TR2,ES23.15))'
        WRITE (LUN, '(I10.10)') SIZE(ss)
        WRITE (LUN,'(A80)') form_string
        DO j = 1, jsize
           WRITE (LUN,form_string) ss(j)
        END DO
!
      CLOSE (LUN)
!
      RETURN
      END SUBROUTINE Write_qtcond
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_vdrop ()
      IMPLICIT NONE
!_____________________________________________________________________________
!
!     Subroutine to read cross polar cap potential drops
!     and place them in vinput array.  Times are stored in
!     ivtime array.  nvmax is actual number of potential drop
!     values.  These results are used in subroutine getv to
!     interpolate in time to get potential at any time of
!     interest.
!     rws  03-20-97
!_____________________________________________________________________________
!
      INTEGER (iprec) :: nv, nvmax
      LOGICAL         :: logical_flag
      LOGICAL, SAVE   :: called_already = .FALSE.
!
      IF (called_already) RETURN
!
      INQUIRE (FILE = rcmdir//'rcmpcp_inp', EXIST = logical_flag)
      IF (.NOT.logical_flag ) STOP 'READV: RCMPCP_INP not found'
      INQUIRE (UNIT = LUN, OPENED = logical_flag)
      IF (logical_flag) STOP 'READV: LUN is already open'
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FILE = rcmdir//'rcmpcp_inp')
      nvmax = 0
      DO
         READ (LUN,*, END = 19 )
         nvmax = nvmax + 1
      END DO
  19  CLOSE (UNIT = LUN)
!
      ALLOCATE (ivtime (nvmax), vinput (nvmax), vinput_phase(nvmax) )
!
      OPEN (UNIT = LUN, STATUS ='OLD', FILE = rcmdir//'rcmpcp_inp') 
      DO nv = 1, nvmax
         READ (LUN, *) ivtime (nv), vinput (nv), vinput_phase(nv)
      END DO
      CLOSE (UNIT = LUN)
      called_already = .TRUE.
!
      RETURN 
      END SUBROUTINE Read_vdrop
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    FUNCTION Get_imin_for_grid (j)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: j
    INTEGER (iprec) :: Get_imin_for_grid
!
    Get_imin_for_grid = MIN (imin_j(Bjmod(j,jwrap,jsize)), &
                             imin_j(Bjmod(j+1,jwrap,jsize)), &
                             imin_j(Bjmod(j-1,jwrap,jsize)) )
    RETURN
    END FUNCTION Get_imin_for_grid
!
!
!
!
     FUNCTION Lt_from_aloct (phi)
     IMPLICIT NONE
     REAL (rprec), INTENT (IN) :: phi
     REAL (rprec) :: Lt_from_aloct
!
!    Convert an RCM phi angle (aloct in ionosphere) to MLT
!    Output (result) is:   0.0  <=  result < 24.00
!
     IF (phi < zero .OR. phi > pi_two) THEN
        WRITE (*,*) 'IN LT_FROM_ALOCT, PHI IS OUT OF BOUNDS'
        STOP
     ELSE
        Lt_from_aloct = MODULO((phi-pi)*RTH, 24.0_rprec)
        Lt_from_aloct = MODULO(Lt_from_aloct, 24.0_rprec)
     END IF
     RETURN
     END FUNCTION Lt_from_aloct
!
!
!
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_grid ( )
      IMPLICIT NONE
      INTEGER (iprec) :: istat, isize_pr, jsize_pr, jwrap_pr
      CHARACTER (LEN=80) :: form_length
      OPEN (UNIT = LUN, STATUS = 'OLD', FORM = 'FORMATTED', &
            FILE = rcmdir//'rcmcrd11', IOSTAT = istat)
         IF (istat /= 0) STOP 'ERROR OPENING RCMCRD11'
         READ (LUN, '(A80)') form_length
         READ (UNIT = LUN, FMT = form_length) isize_pr, jsize_pr, jwrap_pr, dlam, dpsi, ri, re
      CLOSE (UNIT = LUN)
!
!
      IF (isize /= isize_pr .OR. jsize /= jsize_pr .OR. jwrap /= jwrap_pr) THEN
         WRITE (*,*) ' GRID SIZES IN rcmcrd11 DO NOT MATCH THOSE IN THE CODE'
         STOP
      END IF
!
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FORM = 'FORMATTED', &
            FILE = rcmdir//'rcmcrd21', IOSTAT = istat)
         IF (istat /= 0) STOP 'ERROR OPENING RCMCRD21'
         READ (LUN, '(A80)') form_length
         READ (UNIT=LUN, FMT = form_length) alpha
         READ (UNIT=LUN, FMT = form_length) beta
         READ (UNIT=LUN, FMT = form_length) colat
         READ (UNIT=LUN, FMT = form_length) aloct
         READ (UNIT=LUN, FMT = form_length) vcorot
         READ (UNIT=LUN, FMT = form_length) bir
         READ (UNIT=LUN, FMT = form_length) sini
      CLOSE (UNIT = LUN)
      RETURN
      END SUBROUTINE Read_grid
!
!
!
      SUBROUTINE Write_grid ( )
      IMPLICIT NONE
      INTEGER (iprec) :: istat, isize_pr, jsize_pr, jwrap_pr
      CHARACTER (LEN=80) :: form_string
      OPEN (UNIT = LUN, STATUS = 'REPLACE', FORM = 'FORMATTED', &
            FILE = rcmdir//'rcmcrd11', IOSTAT = istat)
         IF (istat /= 0) STOP 'ERROR OPENING RCMCRD11'
       form_string = '(3(TR2,I10),4(TR2,ES23.15))'
       WRITE (UNIT = LUN, FMT = '(A80)') form_string
       WRITE (UNIT = LUN, FMT=form_string) &
              isize, jsize, jwrap, dlam, dpsi, Re, Ri
      CLOSE (UNIT = LUN)
!
!
      OPEN (UNIT = LUN, STATUS = 'REPLACE', FORM = 'FORMATTED', &
            FILE = rcmdir//'rcmcrd21', IOSTAT = istat)
         IF (istat /= 0) STOP 'ERROR OPENING RCMCRD21'
       form_string = '(3(TR2,ES23.15))'
       WRITE (LUN, '(A80)') form_string
       WRITE (LUN, form_string) alpha
       WRITE (LUN, form_string) beta
       WRITE (LUN, form_string) colat
       WRITE (LUN, form_string) aloct
       WRITE (LUN, form_string) vcorot
       WRITE (LUN, form_string) bir
       WRITE (LUN, form_string) sini
      CLOSE (UNIT = LUN)
      RETURN
      END SUBROUTINE Write_grid
!
!
      SUBROUTINE Read_plasma_H5
        use ioh5
        use files
        implicit none
        logical :: doSP
        type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
        integer :: nvar


        doSP = .false.
        call ClearIO(IOVars) !Reset IO chain
        call AddInVar(IOVars,"alamc")
        call AddInVar(IOVars,"ikflavc")
        call AddInVar(IOVars,"fudgec")

        call ReadVars(IOVars,doSP,RCMGAMConfig)

        !Store data
        alamc(:)   = IOVars(1)%data
        ikflavc(:) = IOVars(2)%data
        fudgec(:)  = IOVars(3)%data
! reset to make sure species if ikflav ==1 alamc is set to negative, for electrons
        where(ikflavc==1)alamc = -abs(alamc)

      END SUBROUTINE Read_plasma_H5
!
!
!
      FUNCTION Dipole_Bfield (theta, phi, arg)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: theta, phi
      INTEGER (iprec), INTENT (IN) :: arg
      REAL (rprec) :: Dipole_Bfield
!
!       THETA and PHI are in radians, THETA measured from n.pole down,
!       PHI measured from noon to dusk to midnight etc.
!       Distances RMIN, XMIN, and YMIN are in units of RE
!       Since in the RCM, BESU is in [nT], and the factor of RE is ommited
!       from the formula, VM has units of (RE/nT)**(-2/3)
!
      REAL (rprec) :: rmin, xmin, ymin, bmin, vm
!
      rmin = one / SIN(theta)**2
      xmin = rmin * COS (phi)
      ymin = rmin * SIN (phi)
      bmin = besu * (one/rmin**3)
      vm   = (32.0_rprec/35.0_rprec* rmin**4 / Besu *  &
              SQRT(one-one/rmin)* &
              (one+half/rmin+three/eight/rmin**2+five/eight/two/rmin**3) &
              ) ** (-two/three)
!
      IF (arg == 1) THEN
         Dipole_Bfield = xmin
      ELSE IF (arg == 2) THEN
         Dipole_Bfield = ymin
      ELSE IF (arg == 3) THEN
         Dipole_Bfield = bmin
      ELSE IF (arg == 4) THEN
         Dipole_Bfield = vm
      ELSE
         STOP 'ILLEGAL ARGUMENT FOR DIPOLE_BFIELD'
      END IF
!
      RETURN
      END FUNCTION Dipole_Bfield
!
!
!
!
!
      FUNCTION Deriv_i (a, imin_j) RESULT (d_di)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: imin_j(:)
      REAL (rprec), INTENT (IN) :: a(:,:)
      REAL (rprec) :: d_di (SIZE(a,DIM=1), SIZE(a,DIM=2))
!
!_________________________________________________________________________
!     The idea is to use central differences for a second-order accuracy
!     of the first-order derivative, if we can. If one of the neighboring
!     points is outside the boundary, use forward or back difference
!     (first-order accuracy) if we can. If both points are outside, set
!     derivative to zero.
!_________________________________________________________________________
!
      INTEGER (iprec) :: i, j, j_size, i_size
!
      i_size = SIZE (a, DIM=1)
      j_size = SIZE (a, DIM=2)
!
      DO j = 1, j_size
!
         d_di (1:imin_j(j)-1,j) = 0.0_rprec
!
         i = imin_j(j)
         d_di (i,j) = -1.5_rprec*a(i,j) + 2.0_rprec*a(i+1,j) - 0.5_rprec*a(i+2,j)
!
         DO i = imin_j(j)+1, i_size - 1
            d_di (i,j) = 0.5_rprec*(a(i+1,j) - a(i-1,j))
         END DO
!
         i = i_size
         d_di (i,j) = +1.5_rprec*a(i,j) - 2.0_rprec*a(i-1,j) + 0.5_rprec*a(i-2,j)
!
      END DO
!
      CALL Circle (d_di)
!
      RETURN
      END FUNCTION Deriv_i
!
!
!
      FUNCTION Deriv_j (a, imin_j, j1, j2, error_value) RESULT (d_dj)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: a (:,:), error_value
      INTEGER (iprec), INTENT (IN) :: imin_j(:), j1, j2
      REAL (rprec) :: d_dj (SIZE(a,DIM=1), SIZE(a,DIM=2))
!___________________________________________________________________________
!
!     Take derivative of A with respect to J. The idea is to use
!     grid points that are inside the modeling region. Stanislav 10/27/2000.
!     3/6/2001: modified to be the same as Frank's version in rcm296.
!     Notice that pt_jpp is defined differently (>) than pt_jp(>=), and ditto
!     for pt_jmm and pt_jm. Took this from Frank's version.
!___________________________________________________________________________
!
      INTEGER (iprec) :: i_bnd, i, j, jm, jmm, jp, jpp, j_size
      LOGICAL :: pt_jp, pt_jpp, pt_jm, pt_jmm
!
!
      i_bnd = MAXVAL (imin_j)
      j_size = SIZE (imin_j)
!
      DO j = j1, j2
!
         jp = j + 1
         jm = j - 1
!
         d_dj (1:imin_j(j)-1,j) = 0.0_rprec
!
         DO i = imin_j(j), i_bnd-1
!
            pt_jp  = (i >= imin_j(jp))
            pt_jm  = (i >= imin_j(jm))
!
            IF (pt_jp .AND. pt_jm ) THEN
               d_dj (i,j) = half*(a (i,jp) - a (i,jm))
            ELSE IF (pt_jp ) THEN !j-1 is outside
               d_dj (i,j) = a(i,jp) - a(i,j)
            ELSE IF (pt_jm ) THEN  ! j+1 is outside
               d_dj (i,j) = (a (i,j) - a (i,jm))
            ELSE
               d_dj (i,j) = error_value
            END IF 
!
         END DO
!
         d_dj (i_bnd:,j) = 0.5_rprec*( a(i_bnd:,jp) - a(i_bnd:,jm))
      END DO
!
      CALL Circle (d_dj)
!
      RETURN
      END FUNCTION Deriv_j
!
!
!
      SUBROUTINE Rcm (itimei_in, itimef_in,&
                      idt_in, idt1_in, idt2_in, icontrol,stropt,nslcopt,iXML)
!---------------------------------------------
! notes
!     icontrol controls behaviour of the RCM
!       0:  initialize RCM size params and go back
!       (RCMINIT) 1:  initialize RCM grid, energy channels, quit:
!       2:  read in inputs, quit:
!       3: Set initial conditions on plasma (edges and grid-based)
!       4:  run RCM from itimei to itimef with time step idt, quit:
!       ICONWRITERESTART = 31337: write a restart record to RCM
!       ICONWRITEOUTPUT = ICONWRITERESTART + 1: write an HDF5 output
!       ICONRESTART = ICONWRITERESTART - 1:  Read HDF5 restart


      USE xml_input
      USE rcm_timing_module

      IMPLICIT NONE
!
      type(XML_Input_T), intent(in), optional :: iXML
      INTEGER (iprec), INTENT (IN) :: itimei_in, itimef_in, &
                                      idt_in, idt1_in,& 
                                      idt2_in, icontrol
      character(len=*), intent(in), optional :: stropt
      integer(iprec)  , intent(in), optional :: nslcopt
      CHARACTER(LEN=8) :: real_date
      CHARACTER (LEN=8) :: time_char
      CHARACTER(LEN=10) ::real_time
      CHARACTER(LEN=80) :: ST='', PS='', HD='', string_null=''
      LOGICAL :: FD, logical_flag
!                                                                       
!
      INTEGER (iprec), SAVE :: itimei, itimef, idt, idt1, idt2 
      INTEGER (iprec), SAVE :: itout1, itout2,  itcln, idebug, i_time, &
                         k, kc, n, i_avg
      REAL (rprec) :: dt


      CALL SYSTEM_CLOCK (timer_start(1), count_rate)


      IF (IsCoupledExternally) then
       itimef = itimef_in
       itimei = itimei_in
       idt    = idt_in
       idt1   = idt1_in
       idt2   = idt2_in
      END IF

      IF (icontrol == ICONWRITERESTART) then  ! write a restart record to RCM
         call WriteRCMH5(stropt,nslcopt,isRestart=.true.)
         return
      ENDIF
      IF (icontrol == ICONWRITEOUTPUT) then  ! write an HDF5 output
        call WriteRCMH5(stropt,nslcopt,isRestart=.false.)
        return
      ENDIF
      IF (icontrol == ICONRESTART) then
        !Read HDF5 restart
        call ReadRCMRestart(stropt,nslcopt)
        return
      ENDIF

   IF (icontrol == 0) then  ! initialize RCM size params and go back:


      !  Grid limit parameters
      i1   = 2          ! will reset anyway
      i2   = isize - 1
      j1   = jwrap
      j2   = jsize - 1
      iint = 1 
      jint = 1      

      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN

   END IF

   IF (icontrol == 1) then   ! initialize RCM grid, energy channels, quit:
      CALL Read_grid ()
      if (isGAMRCM) then
        CALL Read_plasma_H5()
      else
       STOP ' Wrong read plasma'
      endif
      
      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN
   END IF



   IF (icontrol == 2) then ! read in inputs, quit:
      !K: Splitting based on whether running coupled to Gamera
      if (isGAMRCM) then
        if(present(iXML)) then
          call RCM_Params_XML(iXML)
        else
          call RCM_Params_XML()
        endif
        CALL Read_dktime_H5(L_dktime)
      else

        OPEN (UNIT = LUN, FILE = rcmdir//'rcm.params', STATUS = 'OLD', &
                 ACTION = 'READ', FORM = 'FORMATTED')

  !
         READ (LUN, '(a80)') label%char! 5.  text label
         READ (LUN,*) idebug   ! 6.  0 <=> do disk printout
         READ (LUN,*) imin  !  5.  i-value of poleward bndy
         READ (LUN,*) ipot  !  6.  which potential solver to use
         READ (LUN,*) iwind !  9.  0 is no neutral winds
         READ (LUN,*) ibnd_type  ! 14.  type of bndy (1-eq.p, 2-iono)
         READ (LUN,*) ipcp_type  ! 14.  type of bndy (1-eq.p, 2-iono)
         READ (LUN,*) nsmthi! 15.  How much to smooth cond in I
         READ (LUN,*) nsmthj! 16.  How much to smooth cond in J
         READ (LUN,*) icond ! 17. 1 is active conductances, 2 is Hardy with kp  
         READ (LUN,*) ifloor! 18. if true, install a floor for EFLUX
         READ (LUN,*) icorrect! 19. if true, make lat. correction to EFLUX
  !
         READ (LUN,*) cmax    ! in rcm_mod_balgn
         READ (LUN,*) eeta_cutoff ! as a fraction
         READ (LUN,*) tol_gmres ! should be 1e-5
         READ (LUN,*) itype_bf  ! 1 is interpolate for HV, 2--MHD code, 3--receive through module
         READ (LUN,*) i_advect  ! 1-interpolate, 2rCLAWPACK/inter, 3-CLAWPACK
         READ (LUN,*) i_eta_bc  ! 1-time-dep. from file, 2-constant for run
         READ (LUN,*) kill_fudge ! .true. means no loss
            if (kill_fudge) then; fudgec = 0.0; end if
         READ (LUN,*) i_birk  ! birk calculation 1=default 3 = new
         READ (LUN,*) L_dktime
         READ (LUN,*) sunspot_number
         READ (LUN,*) L_move_plasma_grid


        ! now run parameters (bypassed via arguments if coupled):
        IF (.NOT.IsCoupledExternally) then
         READ (LUN,*) itimei   ! 1.  start time
         READ (LUN,*) itimef   ! 2.  end time
!         READ (LUN,*) irdr     ! 3.  record # to read in
!         READ (LUN,*) irdw     ! 4.  record # to write out
         READ (LUN,*) idt   !  1.  basic time step in program
         READ (LUN,*) idt1  !  2.  t-step for changing disk write rcds
         READ (LUN,*) idt2  !  3.  t-step for writing formatted output
        END IF

        CLOSE (UNIT = LUN)
          !  Read in other inputs, both constant and run-specific:
          
          IF (.NOT.IsCoupledExternally) CALL Read_vdrop  
          CALL Read_dktime (L_dktime)

      endif !isGAMRCM

      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN

   END IF


   IF (icontrol == 3) then   !-->  Set initial conditions on plasma (grid-based):
!
      ! Open file for formatted output and do initial print out :
      CALL Date_and_time (real_date, real_time)
      IF (itimei == 0) THEN
         ST = 'REPLACE'
         PS = 'APPEND'
         HD = 'BEGINNING NEW RUN'
      ELSE
         ST = 'OLD'
         PS = 'APPEND'
         HD = 'CONTINUE SAME RUN'
      END IF
      !K: Commenting out output
      !write(*,*) "L9604, rcm.printout", LUN_2, ST, PS

      OPEN  (LUN_2, FILE = rcmdir//'rcm.printout', STATUS = ST, POSITION = PS)
      OPEN  (LUN_3, FILE = rcmdir//'rcm.index',  STATUS = ST, POSITION = PS)
      CALL Initial_printout ()
      CLOSE (LUN_3)
      CLOSE (LUN_2)

      i1 = imin + 1


      IF (itimei == 0) THEN

         IF (.NOT.IsCoupledExternally) then
!            CALL Read_array (rcmdir//'rcmeeta_inp',   irdr, label, ARRAY_3D = eeta,ASCI=asci_flag)
         ELSE
            ! grid-based plasma should have been set up elsewhere and passed via module, do nothing here:
            !  CALL Read_array (rcmdir//'rcmeeta_inp',   irdr, label, ARRAY_3D = eeta,ASCI=asci_flag)
         END IF

      ELSE
        if (isGAMRCM) then
          imin_j = CEILING (bndloc)
          !If on first record, create fresh binary files
!          if (irdr == 1) CALL Disk_write_arrays ()
          if (itimei==0)then
              CALL Disk_write_arrays ()
              !call WriteRCMH5(stropt,nslcopt,isRestart=.false.)
          end if

        else  
!          CALL Read_array (rcmdir//'rcmbndloc', irdr, label, ARRAY_1D = bndloc)
          imin_j = CEILING (bndloc)
!          CALL Read_array (rcmdir//'rcmeeta',   irdr, label, ARRAY_3D = eeta)
!         IF (label%intg(6) /= itimei-idt )THEN
          IF (label%intg(6) /= itimei )THEN
            WRITE (*,*)' label%intg(6) =',label%intg(6),' itimei =',itimei
            !write(*,*) 'T in file /=  ITIMEI for EETA, RESTART IS CORRUPTED'
            STOP 'T in file /=  ITIMEI for EETA, RESTART IS CORRUPTED'
          END IF
        endif !isGAMRCM
      END IF

      ! IF hot restart, read V and check the time label:
      IF (.NOT.IsCoupledExternally) THEN
         IF (itimei /= 0) THEN
!         IF (irdr /= 1) THEN
            WRITE (*,'(A)', ADVANCE='NO') &
              'HOT restart, reading V from file to check time label...'
!            CALL READ_array (rcmdir//'rcmv', irdr, label, ARRAY_2D = v)
            IF (label%intg(6) /= itimei )THEN
              WRITE (*,*)' label%intg(6) =',label%intg(6),' itimei =',itimei
              STOP 'T in file /=  ITIMEI for V'
            ELSE
              WRITE (*,*) 'OK'
            END IF 
         END IF
      END IF


      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN

   END IF



   IF (icontrol == 4) then  ! run RCM from itimei to itimef with time step idt, quit:
      CALL SYSTEM_CLOCK (timer_start(2), count_rate)

      v_avg    = zero
      birk_avg = zero
      eeta_avg = zero
      i_avg    = 0
!
!                                                                       
!
!*******************  main time loop  *************************
!
      IF (idt1/idt*idt /= idt1) STOP 'RCM: idt1--idt'
      itout1 = itimei   ! next time to write to disk; set this to
!                         ITIMEI to write out initial configuration
      itout2 = itimei   ! next time to do formatted output
      itcln = itimei    ! next time to call ADD & ZAP 
!
      dt = REAL (idt)
!
      fac = 1.0E-3_rprec * bir * alpha * beta * dlam * dpsi * ri**2 * signbe
!
      DO i_time = itimei, itimef-idt, idt 

!
         CALL Comput (i_time, dt)
!
         v_avg    = v_avg    + v
         birk_avg = birk_avg + birk
         eeta_avg = eeta_avg + eeta
         i_avg    = i_avg    + 1
!
         IF (i_time == itout1) THEN
!
            birk_avg = birk_avg / REAL(i_avg)
            eeta_avg = eeta_avg / REAL(i_avg)
            v_avg    = v_avg    / REAL(i_avg)
            i_avg    = 0
!
            CALL Disk_write_arrays ()
            !call WriteRCMH5(stropt,nslcopt,isRestart=.false.)
!            call AddToList(i_time,rcm_timing)

            itout1 = MIN (itout1 + idt1, itimef-idt)
!
            ! this is a special case: if we wrote output and time is not
            ! last time (i.e., we are not exiting RCM), then reset average
            ! arrays. Otherwise, deal with them separately below:

            IF (i_time < itimef-idt) then
               birk_avg = zero
               eeta_avg = zero
               v_avg    = zero
               i_avg    = 0
            END IF
!
         END IF


         ! this will force RCM to stop with E-field and plasma
         ! in sync (at the same time). Also, since it is time
         ! to exit RCM, we average the arrays but not reset them to zero
         ! so that average arrays stay in memory:

         IF (i_time == itimef - idt) then
             IF (i_avg > 0) then
               birk_avg = birk_avg / REAL(i_avg)
               eeta_avg = eeta_avg / REAL(i_avg)
               v_avg    = v_avg    / REAL(i_avg)
               i_avg    = 0
             END IF
             CYCLE ! exit loop
         END IF
!
         CALL Move_plasma ( dt )
!
      END DO


      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      CALL SYSTEM_CLOCK (timer_stop(2), count_rate)      
      timer_values (2) = (timer_stop (2) - timer_start (2))/count_rate

      CALL Formatted_output ()

      RETURN

   END IF



   IF (icontrol == 5) then ! finalize RCM, quit:
      INQUIRE(FILE = rcmdir//'rcm.printout', EXIST=logical_flag)
      IF (.NOT.logical_flag) then
         ST='REPLACE'
      else
         ST='OLD'
      endif
      write(*,*) "L9786, rcm.printout",itimei
      OPEN  (LUN_2, FILE = rcmdir//'rcm.printout', STATUS=ST, POSITION = 'APPEND')
      WRITE (LUN_2,'(//A)') 'End RCM timing table'
      CLOSE (UNIT = LUN_2)
      CLOSE (UNIT = LUN_3)
      RETURN
   END IF

   
   WRITE (*,*) ' RCM was called with an invalid value of Icontrol, aborting ...'
   STOP


      CONTAINS
!
!
      SUBROUTINE Initial_printout ()
!
      WRITE (LUN_3,'(T2,A)',ADVANCE='NO') TRIM(HD)
      WRITE (LUN_3,'(A11,A4,A1,A2,A1,A2, A8,A2,A1,A2,A1,A2)') &
            '  TODAY IS ', real_date(1:4), '/', &
                           real_date(5:6), '/', &
                           real_date (7:8), &
            '  TIME: ', real_time(1:2), ':', &
                      real_time(3:4), ':', &
                      real_time(5:6)
      WRITE (LUN_3,902) 
      write(6,*)'here'
!
      WRITE (LUN_2,*) 'START OF RCM RUN:'
      WRITE (LUN_2,'(A,I6,A,I6)') 'WILL START AT ITIMEI=',itimei, &
                    '  AND STOP AT ITIMEF=',itimef
      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'time step =', idt 
      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'disk write time step=', idt1 
      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'printout time step=', idt2 
      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'imin =', imin 
      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'start at itimei=', itimei
      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'stop at itimef=', itimef 
!      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'read at itimei from REC ',irdr
!      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'start writing at REC ', irdw 
      WRITE (LUN_2,'(/T5,A)' ) 'SIZES PARAMETERS:'
         WRITE (LUN_2,'(T10,A,T20,I6)') 'isize=',isize
         WRITE (LUN_2,'(T10,A,T20,I6)') 'jsize=',jsize
         WRITE (LUN_2,'(T10,A,T20,I6)') 'ksize=',ksize
         WRITE (LUN_2,'(T10,A,T20,I6)') 'kcsize=',kcsize
         WRITE (LUN_2,'(T10,A,T20,I6)') 'iesize=',iesize
      WRITE (LUN_2,'(/T5,A)' ) 'GRID PARAMETERS:'
         WRITE (LUN_2,'(T10,A,T20,G9.2)') 'dlam=',dlam
         WRITE (LUN_2,'(T10,A,T20,G9.2)') 'dpsi=',dpsi
         WRITE (LUN_2,'(T10,A,T20,G9.2)') 're=',re
         WRITE (LUN_2,'(T10,A,T20,G9.2)') 're=',re
         WRITE (LUN_2,'(T10,A,T20,2G9.2)') 'xmass',xmass
      WRITE (LUN_2,'(/T5,A)' ) 'PLASMA EDGES PARAMETERS:'
      WRITE (LUN_2,'(/T5,A)' ) 'PLASMA GRID PARAMETERS:'
         DO kc = 1, kcsize
           WRITE (LUN_2,'(T10,A,I3,T20,A,G9.2,T45,A,ES9.2, T65,A,F5.2)') &
            'kc=', kc, 'alamc=',  alamc(kc), 'etac=', etac(kc), 'f=', fudgec(kc)
         END DO 
      WRITE (LUN_2,'(T2,A,T20,I2)')     'IPOT =',     ipot
      WRITE (LUN_2,'(T2,A,T20,I2)')     'ICOND = ',   icond
      WRITE (LUN_2,'(T2,A,T20,I2)')     'IBND = ',    ibnd_type
      WRITE (LUN_2,'(T2,A,T20,I2)')     'IPCP_TYPE=', ipcp_type


      WRITE (LUN_2,'(//A)') 'Begin RCM timing table'
      WRITE (LUN_2,'(T1,A,T16,A,T31,A,T46,A)')  'RCM itime [s]', '  record#', 'cpu_time [s]', 'sum_cpu_time [s]'

  902 FORMAT (T2,'TIME', T12,'ITIME' , T19,'REC#' ,&
              T26,'VDROP', T33,'FSTOFF',   &
              T46,'FMEB', T53, 'DST', T62,'FCLPS', T69,'VDROP_PHASE' )
      RETURN
      END SUBROUTINE Initial_printout
!
!
!
         SUBROUTINE Disk_write_arrays ()
!
!        Writing rcm arrays to files is done at each time step, but the
!        record number is changed only with the time step specified in
!        'rcm.params'. This ensures that if the model crashes, files
!        contain the most recent arrays.
!
!
!________We call OUTPUT subroutine with this flag. The policy is
!        that if we start at time=0, then we delete any old files
!        and start from scratch. Otherwise, continue to output to
!        existing files if they exist or create them if not.
!
         IF (i_time == 0) THEN
            FD = .TRUE.
         ELSE
            FD = .FALSE.
         END IF
!
         label%intg = 0
         label%real = zero
         label%char   = ''
!
!        UT TIME  = HH:MM:SS=ilabel(3):ilabel(4):ilabel(5)
!
         label%intg (2) = i_time ! UT in seconds
         label%intg (3) = (i_time) / 3600! hrs of UT
         label%intg (4) = MOD (label%intg (2), 3600) / 60 ! mints of UT
         label%intg (5) = MOD (label%intg (2), 60)  ! scs of UT time
         label%intg (6) = i_time       ! elapsed time in seconds
         label%intg (8) = isize
         label%intg (9) = jsize
         label%intg (10) = ksize
!        label%intg (12) used in OUTPUT and READ3D for kmax(=kdim)
!        label%intg (13) used in OUTPUT and READ3D for k-index
         label%intg (14) = - 1
!
!        label%real (1) = eb   !phoney loss
         label%real (2) = cmax
         label%real (12) = fmeb
         label%real (13) = fstoff
         label%real (14) = fdst
         label%real (15) = fclps
         label%real (16) = vdrop
!         label%real (17) = kp
!
         ST = 'OLD'
         PS = 'APPEND'
         OPEN  (LUN_3, FILE = rcmdir//'rcm.index',  STATUS = ST, POSITION = PS)
         WRITE (time_char,'(I2.2,A1,I2.2,A1,I2.2)') &
               label%intg(3), ':', label%intg(4), ':', label%intg(5)
         WRITE (*,'(T2,A21,I5.5,A10,TR4)') &
                'RCM:-->TIME_STEP, T=', i_time,'('//time_char//')'
!                                                                       
!        IF (i_time == itout1 .OR. i_time == itimef) THEN
            WRITE (lun_3,901) time_char, &
                     i_time, vdrop, fstoff, fmeb, fdst, fclps, vdrop_phase
            CLOSE (LUN_3)
  901       FORMAT (T2,A8, T12,I6,  T26,F5.1, &
                    T39,F5.2, T46,F5.2, T53,F7.1, T62,F5.1, T69, F6.2)
!        END IF
!
         IF (idebug == 0) THEN 

                 STOP 'idebug =0, should not be here'
!
!
!
         END IF
         RETURN
         END SUBROUTINE Disk_write_arrays
!
!
        !HDF5 Restart reader
        subroutine ReadRCMRestart(runid,nStp)
          use ioh5
          use files
          implicit none
          character(len=*), intent(in) :: runid
          integer(iprec), intent(in) :: nStp
          logical :: doSP !Do single precision
          character(len=strLen) :: H5File
          type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
          integer(iprec) :: nvar,nres

        !Prepare for reading
          doSP = .false. !Restarts are always double precision
          nres = nStp-1 !nStp holds number for *NEXT* restart output
          if (nres == -1) then
            !Use sym link
            H5File = trim(runid) // ".RCM.Res.XXXXX.h5"
          else
            !Use actual #
            write (H5File, '(A,A,I0.5,A)') trim(runid), ".RCM.Res.", nres, ".h5"
          endif
          write(*,*) 'Restarting RCM with file, ', trim(H5File)
          call ClearIO(IOVars) !Reset IO chain

        !List variables to read
          !Scalars (need to specify integers), order doesn't matter

          call AddInVar(IOVars,"itimei",vTypeO=IOINT )
          call AddInVar(IOVars,"isize" ,vTypeO=IOINT )
          call AddInVar(IOVars,"jsize" ,vTypeO=IOINT )
          call AddInVar(IOVars,"ksize" ,vTypeO=IOINT )
          call AddInVar(IOVars,"cmax"  ,vTypeO=IOREAL)
          call AddInVar(IOVars,"fmeb"  ,vTypeO=IOREAL)
          call AddInVar(IOVars,"fstoff",vTypeO=IOREAL)
          call AddInVar(IOVars,"fdst"  ,vTypeO=IOREAL)
          call AddInVar(IOVars,"fclps" ,vTypeO=IOREAL)
          call AddInVar(IOVars,"vdrop" ,vTypeO=IOREAL)
          call AddInVar(IOVars,"kp"    ,vTypeO=IOREAL)
          call AddInVar(IOVars,"i_avg" ,vTypeO=IOREAL)

          !Arrays
          call AddInVar(IOVars,"rcmxmin"  )
          call AddInVar(IOVars,"rcmymin"  )
          call AddInVar(IOVars,"rcmzmin"  )
          call AddInVar(IOVars,"rcmvm"    )
          call AddInVar(IOVars,"rcmbmin"  )
          call AddInVar(IOVars,"rcmbndloc")

          call AddInVar(IOVars,"rcmetac"   )
          call AddInVar(IOVars,"rcmeeta"   )
          call AddInVar(IOVars,"rcmeetaavg")

          call AddInVar(IOVars,"rcmpedlam" )
          call AddInVar(IOVars,"rcmpedpsi" )
          call AddInVar(IOVars,"rcmhall"   )
          call AddInVar(IOVars,"rcmeavg"   )
          call AddInVar(IOVars,"rcmeflux"  )
          call AddInVar(IOVars,"rcmbirk"   )
          call AddInVar(IOVars,"rcmbirkavg")

          call AddInVar(IOVars,"rcmetac")
          call AddInVar(IOVars,"rcmeeta")
          call AddInVar(IOVars,"rcmeetaavg")

          call AddInVar(IOVars,"rcmv")
          call AddInVar(IOVars,"rcmvavg")

          call AddInVar(IOVars,"alpha")
          call AddInVar(IOVars,"aloct")
          call AddInVar(IOVars,"colat")
          call AddInVar(IOVars,"beta")
          call AddInVar(IOVars,"bir")
          call AddInVar(IOVars,"sini")

          !Extra stuff (not in write arrays)
          call AddInVar(IOVars,"alamc")
          
        !Now do actual reading
          call ReadVars(IOVars,doSP,H5File)

        !Parse data and put it where it goes, need to do each variable
          !Scalars
          itimei = GetIOInt(IOVars,"itimei")
          
          cmax   = GetIOReal(IOVars,"cmax")
          fmeb   = GetIOReal(IOVars,"fmeb")
          fstoff = GetIOReal(IOVars,"fstoff")
          fdst   = GetIOReal(IOVars,"fdst")
          fclps  = GetIOReal(IOVars,"fclps")
          vdrop  = GetIOReal(IOVars,"vdrop")
          i_avg  = GetIOReal(IOVars,"i_avg")

          !Pull 2D arrays
          call IOArray2DFill(IOVars,"rcmxmin",xmin)
          call IOArray2DFill(IOVars,"rcmymin",ymin)
          call IOArray2DFill(IOVars,"rcmzmin",zmin)
          call IOArray2DFill(IOVars,"rcmbmin",bmin)
          
          call IOArray2DFill(IOVars,"rcmv",v)
          call IOArray2DFill(IOVars,"rcmvavg",v_avg)
          call IOArray2DFill(IOVars,"rcmvm",vm)

          call IOArray2DFill(IOVars,"rcmbirk",birk)
          call IOArray2DFill(IOVars,"rcmbirkavg",birk_avg)

          call IOArray2DFill(IOVars,"rcmhall",hall)
          call IOArray2DFill(IOVars,"rcmpedlam",pedlam)
          call IOArray2DFill(IOVars,"rcmpedpsi",pedpsi)

          call IOArray2DFill(IOVars,"alpha",alpha)
          call IOArray2DFill(IOVars,"aloct",aloct)
          call IOArray2DFill(IOVars,"colat",colat)
          call IOArray2DFill(IOVars,"beta",beta)
          call IOArray2DFill(IOVars,"bir",bir)
          call IOArray2DFill(IOVars,"sini",sini)


          !Pull 1D arrays
          call IOArray1DFill(IOVars,"rcmetac",etac)
          call IOArray1DFill(IOVars,"rcmbndloc",bndloc)
          call IOArray1DFill(IOVars,"alamc",alamc)

          !Pull 3D arrays
          call IOArray3DFill(IOVars,"rcmeavg",eavg)
          call IOArray3DFill(IOVars,"rcmeeta",eeta)
          call IOArray3DFill(IOVars,"rcmeetaavg",eeta_avg)
          call IOArray3DFill(IOVars,"rcmeflux",eflux)
          
        end subroutine ReadRCMRestart

        !HDF5 output routine
        !isRestart = Whether we're writing restart dump or regular output slice
        subroutine WriteRCMH5(runid,nStp,isRestart)
          use ioh5
          use files
          implicit none
          character(len=*), intent(in) :: runid
          integer(iprec), intent(in) :: nStp
          logical, intent(in) :: isRestart

          type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
          logical :: doSP !Do single precision output
          character(len=strLen) :: H5File,gStr,lnResF

        !Prepare for output
          !Reset IO chain
          call ClearIO(IOVars)
          !Distinguish output slices vs restarts
          if (isRestart) then
            doSP = .false. !Double precision restarts
            write (H5File, '(A,A,I0.5,A)') trim(runid), ".RCM.Res.", nStp, ".h5"
          else
            !Regular output
            doSP = .true.
            H5File = trim(runid) // ".rcm.h5"
            write (gStr, '(A,I0)') "Step#", nStp
          endif

        !Attributes
          call AddOutVar(IOVars,"time",1.0_rp*itimei)
          call AddOutVar(IOVars,"itimei",itimei)
          call AddOutVar(IOVars,"isize" ,isize )
          call AddOutVar(IOVars,"jsize" ,jsize )
          call AddOutVar(IOVars,"ksize" ,ksize )
          call AddOutVar(IOVars,"cmax"  ,cmax  )
          call AddOutVar(IOVars,"fmeb"  ,fmeb  )
          call AddOutVar(IOVars,"fstoff",fstoff)
          call AddOutVar(IOVars,"fdst"  ,fdst  )
          call AddOutVar(IOVars,"fclps" ,fclps )
          call AddOutVar(IOVars,"vdrop" ,vdrop )
          call AddOutVar(IOVars,"i_avg" ,i_avg )

        !Arrays
          call AddOutVar(IOVars,"rcmxmin",xmin)
          call AddOutVar(IOVars,"rcmymin",ymin)
          call AddOutVar(IOVars,"rcmzmin",zmin)
          call AddOutVar(IOVars,"rcmvm"  ,vm  )
          call AddOutVar(IOVars,"rcmbmin",bmin)
          call AddOutVar(IOVars,"rcmbndloc",bndloc)

          call AddOutVar(IOVars,"rcmetac"   ,etac)
          call AddOutVar(IOVars,"rcmeeta"   ,eeta)
          call AddOutVar(IOVars,"rcmeetaavg",eeta_avg)

          call AddOutVar(IOVars,"rcmpedlam" ,pedlam  )
          call AddOutVar(IOVars,"rcmpedpsi" ,pedpsi  )
          call AddOutVar(IOVars,"rcmhall"   ,hall    )
          call AddOutVar(IOVars,"rcmeavg"   ,eavg    )
          call AddOutVar(IOVars,"rcmeflux"  ,eflux   )
          call AddOutVar(IOVars,"rcmbirk"   ,birk    )
          call AddOutVar(IOVars,"rcmbirkavg",birk_avg)

          call AddOutVar(IOVars,"rcmv",v)
          call AddOutVar(IOVars,"rcmvavg",v_avg)

        !Extra stuff not in write_array
          call AddOutVar(IOVars,"alamc",alamc)
          call AddOutVar(IOVars,"aloct",aloct)
          call AddOutVar(IOVars,"colat",colat)
          call AddOutVar(IOVars,"alpha",alpha)
          call AddOutVar(IOVars,"beta" ,beta )
          call AddOutVar(IOVars,"bir"  ,bir  )
          call AddOutVar(IOVars,"sini" ,sini )
          
        !Done staging output, now let er rip
          if (isRestart) then
            call AddOutVar(IOVars,"nRes",nStp)
            call CheckAndKill(H5File) !Always overwrite restarts
            call WriteVars(IOVars,doSP,H5File)
            !Create link to latest restart
            write (lnResF, '(A,A,A,A)') trim(runid), ".RCM.Res.", "XXXXX", ".h5"
            call EXECUTE_COMMAND_LINE('ln -sf '//trim(H5File)//' '//trim(lnResF), wait=.false.)
          else
            call WriteVars(IOVars,doSP,H5File,gStr)
          endif
        end subroutine WriteRCMH5

        subroutine RCM_Params_XML(iXML)
          use xml_input
          use strings

          type(XML_Input_T), intent(in), optional :: iXML
          character(len=strLen) :: inpXML
          type(XML_Input_T) :: xmlInp

          if(present(iXML)) then
            call iXML%GetFileStr(inpXML)
          else
            !Find input deck filename
            call getIDeckStr(inpXML)
          endif
          
          !Create new XML reader w/ RCM as root
          xmlInp = New_XML_Input(trim(inpXML),'RCM',.true.)

          call xmlInp%Set_Val(label%char,"sim/runid","MHD code run")

          !Output
          call xmlInp%Set_Val(idebug,"output/idebug",1) ! 6.  0 <=> do disk printout

          !eflux
          call xmlInp%Set_Val(ifloor,"eflux/ifloor",.true.) ! 18. if true, install a floor for EFLUX
          call xmlInp%Set_Val(icorrect,"eflux/icorrect",.true.) ! 19. if true, make lat. correction to EFLUX

          !Grid
          call xmlInp%Set_Val(imin,"grid/imin",1)
          call xmlInp%Set_Val(ibnd_type,"grid/ibnd_type",4) ! 14.  type of bndy (1-eq.p, 2-iono)
          call xmlInp%Set_Val(ipcp_type,"grid/ipcp_type",13) ! 14.  type of bndy (1-eq.p, 2-iono)
          call xmlInp%Set_Val(nsmthi,"grid/nsmthi",0) ! 15.  How much to smooth cond in I
          call xmlInp%Set_Val(nsmthj,"grid/nsmthj",0) ! 16.  How much to smooth cond in J
          call xmlInp%Set_Val(L_move_plasma_grid,"grid/L_move_plasma_grid",.true.)

          !Catch-all params
          call xmlInp%Set_Val(ipot,"params/ipot",-1) !  6.  which potential solver to use
          call xmlInp%Set_Val(iwind,"params/iwind",0) !  9.  0 is no neutral winds
          call xmlInp%Set_Val(icond,"params/icond",3) ! 1 is active conductances, 2 is Hardy with kp, 3 is input

          call xmlInp%Set_Val(cmax,"params/cmax",3.0) ! in rcm_mod_balgn
          call xmlInp%Set_Val(eeta_cutoff,"params/eeta_cutoff",0.05) ! as a fraction

          !Charge exchange
          call xmlInp%Set_Val(kill_fudge,"chargex/kill_fudge",.false.) ! .true. means no loss
          if (kill_fudge) then
            fudgec = 0.0
          endif
          call xmlInp%Set_Val(L_dktime,"chargex/L_dktime",.true.)
          call xmlInp%Set_Val(sunspot_number,"chargex/sunspot_number",96.0)

          !Clawpack options
          call xmlInp%Set_Val(L_doKaiClaw,"clawpack/doKaiClaw",L_doKaiClaw)

          !Some values just setting
          tol_gmres = 1.0e-5
          itype_bf = 3 ! 1 is interpolate for HV, 2--MHD code, 3--receive through module
          i_advect = 3  ! 1-interpolate, 2rCLAWPACK/inter, 3-CLAWPACK
          i_eta_bc = 2! 1-time-dep. from file, 2-constant for run
          i_birk  = 1 ! birk calculation 1=default 3 = new
        end subroutine RCM_Params_XML

         SUBROUTINE Formatted_output ()
          !K: Suppressing output
          !write(*,*) "L10019, rcm.printout", LUN_2, ST, PS
         OPEN  (LUN_2, FILE = rcmdir//'rcm.printout', STATUS = 'OLD', POSITION = 'append')
         WRITE (LUN_2,'(T1,I10,T31,F10.2,T46,F10.2)')  &
        &  i_time, timer_values(2), timer_values(1) 
         close (lun_2)
         itout2 = itout2 + idt2
         RETURN
         END SUBROUTINE Formatted_output
!
!
      END SUBROUTINE Rcm
!
!
!
!
    SUBROUTINE Read_dktime (L_dktime)
    IMPLICIT NONE
    LOGICAL, INTENT (IN) :: L_dktime
!
!
    IF (L_dktime) THEN
       OPEN (LUN, FILE=rcmdir//'dktable', STATUS='OLD', ACTION='READ')
       READ (LUN,800) dktime
 800   FORMAT (8(E10.3))
       CLOSE (LUN)
    ELSE
        dktime = 0.0
    END IF
    RETURN
    END SUBROUTINE Read_dktime

    !K: HDF5 version of lifetime reader
    SUBROUTINE Read_dktime_H5(L_dktime)
      use ioh5
      use files
      IMPLICIT NONE
      LOGICAL, INTENT (IN) :: L_dktime
      logical :: doSP
      type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
      
      !real(rprec) :: dktime2(irdk,inrgdk,isodk, iondk)

      if (L_dktime) then
        !Read from HDF5
        doSP = .false.
        call ClearIO(IOVars) !Reset IO chain
        call AddInVar(IOVars,"dktable")
        call ReadVars(IOVars,doSP,RCMGAMConfig)

        dktime = reshape(IOVars(1)%data,[irdk,inrgdk,isodk, iondk])
        
        ! !Debugging
        ! dktime2 = reshape(IOVars(1)%data,[irdk,inrgdk,isodk, iondk])
        ! call Read_dktime(L_dktime)
        ! write(*,*) 'dktime1 = ', dktime
        ! write(*,*) 'dktime2 = ', dktime2
        ! write(*,*) 'Del = ', sum(abs(dktime-dktime2))
      endif
    END SUBROUTINE Read_dktime_H5

!
!
!
!
      FUNCTION Cexrat (isp,enrg,rloc,ssn,dktime,irdk,inrgdk,isoldk, &
                       iondk)
      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: isp, irdk, inrgdk, isoldk, iondk
      REAL(rprec), INTENT (IN) :: enrg, rloc, ssn, dktime (irdk,inrgdk,isoldk,iondk)
      REAL(rprec) :: Cexrat
!
!-------------------------------------------------------------------------
!  copyright rice university, 1993
!
!  version 1.00                                 05.09.90
!          2.00                                 02.04.90
!                                       msm delivery version
!          2.10                                 06.11.93
!               error output routed to unit 9
!
!  programmer: r. w. spiro
!
!  purpose:  function subprogram to return charge exchange loss rate
!          (sec**(-1)) for ions of species isp, energy enrg (ev) at
!          l=rloc (re) for sunspot number ssn.  this routine is based
!          on a table generated by james bishop of u. of michigan.
!
!  calling parameters
!        isp       species identifier
!                    isp=2 for h+ ions
!                    isp=3 for o+ ions
!        enrg      energy in ev
!        rloc      radial location (re)
!        ssn       sunspot number
!        dktime    table of ion decay times
!        irdk      radial dimension of dktime array
!        inrgdk    energy dimension of dktime array
!        isoldk    sunspot number dimension of dktime array
!        iondk     number of ion species in dktime array
!--------------------------------------------------------------------------------
!
      INTEGER(iprec), PARAMETER ::irsiz=18,inrgsz=13,isolsz=2,ionsiz=2
      REAL(rprec) ::  elgvec(inrgsz), rvec(irsiz),ssnvec(2), &
               enrglg, br, bnrg, ssnuse, bssn, decayt
      INTEGER(iprec) :: ispndx, ir, inrg
!
      DATA elgvec /2.50,2.75,3.00,3.25,3.50,3.75,4.00, &
                   4.25,4.50,4.75,5.00,5.25,5.50/
!
      DATA rvec /1.50,2.00,2.50,3.00, &
                 3.50,4.00,4.50,5.00, &
                 5.50,6.00,6.50,7.00, &
                 7.50,8.00,8.50,9.00, &
                 9.50,10.00/
!
      DATA ssnvec /0.0,100./
!
!
      IF (irsiz /= irdk .OR. inrgsz /= inrgdk .OR. &
          ionsiz /= iondk .OR. isolsz /= isoldk) THEN
         write(*,*) 'dimension error in function cexrat'
         write(*,*) 'irdk,inrgdk,iondk,isoldk',irdk,inrgdk,iondk,isoldk
         write(*,*) 'irsiz,inrgsz,ionsiz,isolsz',irsiz,inrgsz,ionsiz,isolsz
         write(*,*) 'stopping program in cexrat'
         STOP
      END IF
!
      enrglg = LOG10(enrg) !  work with log10 of particle energy
      ispndx=isp-1
!
      if_1: IF (rloc <= rvec(1)) THEN !  find br for interpolation
         br=1.0
      ELSE IF (rloc > rvec(irdk)) THEN
         br=irdk
      ELSE
         do_1: DO ir=1,irdk-1
            IF (rloc <= rvec(ir+1)) THEN
               br=ir+(rloc-rvec(ir))/(rvec(ir+1)-rvec(ir))
               EXIT do_1
            END IF
         END DO do_1
      END IF if_1
!
      if_2: IF (enrglg.le.elgvec(1)) THEN !  find bnrg for interpolation
         bnrg = 1.0
      ELSE IF (enrglg > elgvec(inrgdk)) THEN
         bnrg = inrgdk
      ELSE
         do_2: DO inrg=1,inrgdk-1
            IF (enrglg <= elgvec(inrg+1)) THEN
               bnrg=inrg+(enrglg-elgvec(inrg))/(elgvec(inrg+1)-elgvec(inrg))
               EXIT do_2
            END IF
         END DO do_2
      END IF if_2
!
!**********  change 9/30/91  *****************************************
!  if ssn.gt.ssnvec(2), then use ssnvec(2) for ssn
      ssnuse=ssn
      IF (ssnuse > ssnvec(2)) ssnuse=ssnvec(2)
!
!*********  end change  9/30/91  ************************************
!
!  find bssn for interpolation
      bssn=1.0+(ssnuse-ssnvec(1))/(ssnvec(2)-ssnvec(1))
!
!  decayt is decay time in seconds
      decayt = G3ntrp (dktime(1_iprec,1_iprec,1_iprec,ispndx),irdk,inrgdk,isoldk,br,bnrg,bssn)
!
      IF (ABS(decayt) < 1.0E-20) THEN
         write(*,*) 'decayt is less than 1.e-20 sec in cexrat'
         write(*,*) 'decayt=',decayt,'  br=',br,'  bnrg=',bnrg,'bssn=',bssn
         write(*,*) 'isp=',isp,'  enrg=',enrg,' rloc=',rloc,' ssn=',ssn
         write(*,*) 'ssnuse=',ssnuse
      END IF
!
!  to get charge exchange rate (sec**9-1)) cexrat, invert decayt
!
      cexrat=1.0/decayt
      RETURN
      END FUNCTION Cexrat
!
!
      FUNCTION G3ntrp (a,imax,jmax,kmax,bi,bj,bk)
      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: imax, jmax, kmax
      REAL(rprec), INTENT (IN) :: a(imax,jmax,kmax), bi, bj, bk
      REAL(rprec) :: G3ntrp
!
!---------------------------------------------------------------------------
!  copyright Rice University, 1993
!
!  VERSION 1.00                            DATE: 01.11.88
!          1.01A                                 02.02.89
!          2.00  MSM DELIVERY VERSION            01.28.93
!
!  PURPOSE: FUNCTION SUBPROGRAM TO PERFORM A GENERAL 3-D LINEAR
!           INTERPOLATION OF ARRAY A(I,J,K) AT PT(BV(1),BV(2),BV(3))
!
!  INPUT:
!       A          3-D ARRAY TO BE INTERPOLATED
!       IMAX       I DIMENSION OF ARRAY A
!       JMAX       J DIMENSION OF ARRAY A
!       KMAX       K DIMENSION OF ARRAY A
!       BI         FLOATING POINT VALUE TO INTERPOLATE IN I DIMENSION
!       BJ         FLOATING POINT VALUE TO INTERPOLATE IN J DIMENSION
!       BK         FLOATING POINT VALUE TO INTERPOLATE IN K DIMENSION
!
!
!  OUTPUT:
!       G3NTRP     INTERPOLATED VALUES OF ARRAY A
!----------------------------------------------------------------------
!
!
      INTEGER(iprec) ::  ndx(3),ndim(3), kstop, jstop, L, i, j, k
      REAL(rprec) ::  BV(3),COEF(3,2), fndx
!
      NDIM(1)=IMAX
      NDIM(2)=JMAX
      NDIM(3)=KMAX
      BV(1)=BI
      BV(2)=BJ
      BV(3)=BK
      DO L=1,3
         NDX(L)=BV(L)
         IF(NDX(L).LT.1) NDX(L)=1
         IF(NDX(L).GT.NDIM(L)-1) NDX(L)=NDIM(L)-1
         IF(NDX(L).LE.0) NDX(L)=1
         FNDX=REAL(NDX(L))
         COEF(L,1)=1.-BV(L)+FNDX
         COEF(L,2)=BV(L)-FNDX
      END DO
!
      G3NTRP=0.
      kstop = MIN(KMAX,2)
      jstop = MIN(JMAX,2)
      DO I=1,2
      DO J=1,jstop
      DO K=1,kstop
         G3ntrp=G3ntrp+ &
          coef(1,i)*coef(2,j)*coef(3,k)*a(ndx(1)+i-1,ndx(2)+j-1,ndx(3)+k-1)
      END DO
      END DO
      END DO
!
      RETURN
      END FUNCTION G3ntrp
!
!
!
!
!
!
!
!
!=========================================================================
!=========================================================================
!Attempt by K: to incorporate newer clawpack 04/20
SUBROUTINE Move_plasma_grid_KAIJU (dt)
  USE rcmclaw, only : claw2ez95
  
  IMPLICIT NONE

  real(rprec), intent(in) :: dt

  real(rprec), dimension(isize,jsize) :: eeta2,veff,dvefdi,dvefdj
  real(rprec), dimension(-1:isize+2,-1:jsize-1) :: loc_didt,loc_djdt,loc_Eta, loc_rate
  real(rprec) :: didt,djdt,mass_factor,r_dist
  integer(iprec) :: i,j,kc,ie,joff,icut,clawiter
  REAL (rprec) :: max_eeta, eps = 0.0 !sbao 07/2019

  joff=jwrap-1
  fac = 1.0E-3*signbe*bir*alpha*beta*dlam*dpsi*ri**2

  
  !$OMP PARALLEL DO default(NONE) &
  !$OMP schedule(dynamic) &
  !$OMP private (i,j,kc,ie,icut,clawiter) &
  !$OMP private (eeta2,veff,dvefdi,dvefdj,didt,djdt) &
  !$OMP private (mass_factor,loc_didt,loc_djdt) &
  !$OMP private (loc_Eta,loc_rate,r_dist,max_eeta) &
  !$OMP shared (alamc,eeta,v,vcorot,vpar,vm,imin_j,j1,j2,joff) &
  !$OMP shared (xmin,ymin,fac,fudgec,bir,sini,L_dktime,dktime,sunspot_number) &
  !$OMP shared (dt,eps)
  DO kc = 1, kcsize
    !If oxygen is to be added, must change this!
    IF (alamc(kc) <= 0.0) THEN
      ie = RCMELECTRON
    ELSE
      ie = RCMPROTON
    END IF

    IF (maxval(eeta(:,:,kc)) == 0.0) then
      cycle
    ENDIF
    mass_factor = SQRT (xmass(1)/xmass(ie))

    !K: Here we're adding corotation to total effective potential
    veff = v + vcorot - vpar + vm*alamc(kc)

    dvefdi = Deriv_i (veff, imin_j)
    dvefdj = Deriv_j (veff, imin_j, j1, j2, 1.0E+26_rprec)
    !K: Why only dvefdj and not dvefdi?
    WHERE (dvefdj > 1.0E+20)
      dvefdj = 0.0
    END WHERE

    loc_Eta  = zero
    loc_didt = zero
    loc_djdt = zero
    loc_rate = zero

    icut=0
    do j=j1,j2
      icut=max(icut,imin_j(j))
      do i=imin_j(j),isize-1
        if (eeta(i,j,kc) > 1.) then
          icut=max(icut,i)
        endif
      end do
    end do !j loop
    icut=icut+5

    DO j = j1, j2
      DO i = 2, isize-1
        loc_didt (i,j-joff) = + dvefdj (i-1,j) / fac(i-1,j)
        loc_djdt (i,j-joff) = - dvefdi (i,j-1) / fac(i-1,j)
        IF (i > icut) THEN
          loc_didt(i,j-joff) = 0.0
          loc_djdt(i,j-joff) = 0.0
        END IF
!
        IF (ie == RCMELECTRON) THEN
          loc_rate(i,j-joff) = Ratefn (fudgec(kc), alamc(kc), sini(i,j),&
                                       bir (i,j), vm(i,j), mass_factor)
        ELSE IF (ie == RCMPROTON) THEN
          IF (L_dktime .AND. i >= imin_j(j)) THEN
            r_dist = SQRT(xmin(i,j)**2+ymin(i,j)**2)
            loc_rate(i,j-joff) = Cexrat (ie, ABS(alamc(kc))*vm(i,j)    , &
                                          R_dist,sunspot_number, dktime, &
                                          irdk,inrgdk,isodk,iondk)
                                          
          ELSE
            loc_rate(i,j-joff) = 0.0
          END IF
        ELSE
          STOP 'UNKNOWN IE IN COMPUTING LOSS'
        END IF !ie = X

      END DO ! i loop

      loc_didt(isize,j-joff) = loc_didt(isize-1,j-joff)
      loc_djdt(isize,j-joff) = loc_djdt(isize-1,j-joff)
      loc_rate(isize,j-joff) = loc_rate(isize-1,j-joff)
    END DO ! j loop

    !Copy to local variables
    loc_Eta (1:isize, 1:jsize-jwrap) = eeta (1:isize, jwrap:jsize-1, kc)     

    !Call clawpack
    call claw2ez95(dt,loc_Eta(1:isize,1:jsize-jwrap),loc_didt(1:isize,1:jsize-jwrap),loc_djdt(1:isize,1:jsize-jwrap),loc_rate(1:isize,1:jsize-jwrap),clawiter)
    
    !Copy out
    DO j = j1, j2
      DO i = imin_j(j)+1, isize-1
        eeta (i, j, kc) = loc_Eta (i, j-joff)
      END DO
    END DO !j loop

    DO j = j1, j2
      IF (veff(imin_j(j+1),j+1)-veff(imin_j(j-1),j-1) < 0.0) THEN
        eeta (imin_j(j),j,kc) = loc_eta (imin_j(j),j-joff)
      END IF
    END DO !j loop

    !floor eeta 12/06 frt
    max_eeta = maxval(eeta(:,:,kc))
    eeta(:,:,kc) = MAX(eps*max_eeta,eeta(:,:,kc))
    CALL Circle (eeta(:,:,kc))    
  ENDDO !kc loop
  
END SUBROUTINE Move_plasma_grid_KAIJU

!=========================================================================
!
SUBROUTINE Move_plasma_grid_NEW (dt)
! USE Rcm_mod_subs
  IMPLICIT NONE
  REAL (rprec), INTENT (IN) :: dt
!_____________________________________________________________________________
!   Subroutine to advance eta distribution for a time step
!   by using new CLAWPACK advection routines
!                                                                       
!   Created:     12-05-00
!_____________________________________________________________________________
!
!
  REAL (rprec) :: eeta2 (isize,jsize), veff  (isize,jsize),  &
                  dvefdi(isize,jsize), dvefdj(isize,jsize), &
                  didt, djdt, mass_factor
  REAL (rprec) :: max_eeta, eps = 0.0 !sbao 07/2019
  INTEGER (iprec) :: i, j, kc, ie
!\\\
  integer(iprec) :: CLAWiter, joff, icut
!  real(rprec), dimension(-1:isize+2,-1:jsize-1):: loc_didt, loc_djdt
!  real(rprec), dimension(-1:isize+2,-1:jsize-1):: loc_Eta, loc_rate
  double precision, dimension(-1:isize+2,-1:jsize-1):: loc_didt, loc_djdt
  double precision, dimension(-1:isize+2,-1:jsize-1):: loc_Eta, loc_rate
  double precision, save :: xlower,xupper,ylower,yupper, T1,T2
  !real(rprec), save :: xlower,xupper,ylower,yupper, T1,T2
  logical, save :: FirstTime=.true.
  REAL(rprec) :: r_dist
  integer(iprec):: ii,istop
  
  joff=jwrap-1
  
  if (FirstTime) then
     T1=0.
  else
     T1=T2
  end if
  T2=T1+dt

  xlower = 1
  xupper = isize
  ylower = zero
  yupper = jsize-3
!///


 fac = 1.0E-3*signbe*bir*alpha*beta*dlam*dpsi*ri**2

!K: Experiment OMP binding on RCM calculation
!!$OMP PARALLEL DO schedule(dynamic) &
!!$PRIVATE (eeta2, veff, dvefdi, dvefdj, didt, djdt, &
!!$OMP                  & mass_factor, loc_didt, loc_djdt,loc_Eta,loc_rate, &
!!$OMP                  & ie, icut, j, i, r_dist, FirstTime, max_eeta) &
!!$OMP        & SHARED (alamc, eeta, v, vcorot, vpar, vm, imin_j, j1, j2, joff, &
!!$OMP                  xmin, ymin, fac, fudgec, bir, sini, L_dktime, dktime, sunspot_number, &                
!!$OMP                  T1, T2, xlower, ylower, xupper, yupper, CLAWiter, eps) &
!!$OMP        & DEFAULT (NONE)

  DO kc = 1, kcsize
!
!    If oxygen is to be added, must change this!
!
     IF (alamc(kc) <= 0.0) THEN
        ie = 1  ! electrons
     ELSE
        ie = 2  ! protons
     END IF
!
!
     IF (MAXVAL(eeta(:,:,kc)) == 0.0) CYCLE
!
     mass_factor = SQRT (xmass(1)/xmass(ie))
!
!    1. Compute the effective potential for the kc energy channel:
! 
    !K: Here we're adding corotation to total effective potential
     veff = v +vcorot - vpar + vm*alamc(kc)
!!!  CALL V_eff_polar_cap (veff)
!
!    2. Differentiate Veff with respect to I and J:
!
!!!  CALL Deriv_i_new (veff, isize, jsize, j1, j2, imin_J, dvefdi)
!!!  CALL Deriv_j_new (veff, isize, jsize, j1, j2, imin_J, dvefdj)
     
     dvefdi = Deriv_i (veff, imin_j)
     dvefdj = Deriv_j (veff, imin_j, j1, j2, 1.0E+26_rprec)
     WHERE (dvefdj > 1.0E+20)
       dvefdj = 0.0
     END WHERE
!
!
     loc_Eta  = zero
     loc_didt = zero
     loc_djdt = zero
     loc_rate = zero
!
!
      icut=0
      do j=j1,j2
         icut=max(icut,imin_j(j))
         do i=imin_j(j),isize-1
            if (eeta(i,j,kc) > 1.) icut=max(icut,i)
         end do
      end do
      icut=icut+5

     DO j = j1, j2
        DO i = 2, isize-1
           loc_didt (i,j-joff) = + dvefdj (i-1,j) / fac(i-1,j)
           loc_djdt (i,j-joff) = - dvefdi (i,j-1) / fac(i-1,j)
           IF (i > icut) THEN
              loc_didt(i,j-joff) = 0.0
              loc_djdt(i,j-joff) = 0.0
           END IF
!
            IF (ie == 1) THEN
               loc_rate(i,j-joff) = Ratefn (fudgec(kc), alamc(kc), sini(i,j),&
                                            bir (i,j), vm(i,j), mass_factor)
            ELSE IF (ie == 2) THEN
               IF (L_dktime .AND. i >= imin_j(j)) THEN
               r_dist = SQRT(xmin(i,j)**2+ymin(i,j)**2)
               loc_rate(i,j-joff) = Cexrat (ie, ABS(alamc(kc))*vm(i,j), &
                                            R_dist, &
                                            sunspot_number, dktime, &
                                            irdk,inrgdk,isodk,iondk)
               ELSE
               loc_rate(i,j-joff) = 0.0
               END IF
            ELSE
               STOP 'UNKNOWN IE IN COMPUTING LOSS'
            END IF
!
        END DO
!       eeta (1:imin_j(j)-1,j,kc) = etac (kc)
        loc_didt(isize,j-joff) = loc_didt(isize-1,j-joff)
        loc_djdt(isize,j-joff) = loc_djdt(isize-1,j-joff)
        loc_rate(isize,j-joff) = loc_rate(isize-1,j-joff)
     END DO
!
! boundary condition correction:
!    DO j = j1, j2
!       IF (loc_didt(imin_j(j),j-joff) < 0.0) THEN
!          eeta(imin_j(j),j,:) = eeta(imin_j(j)+1,j,:)
!       END IF
!    END DO
!
!
!Copy to local variables
     loc_Eta (1:isize, 1:jsize-jwrap) = eeta (1:isize, jwrap:jsize-1, kc)
!    !Set ghost cell values for clawpack solver
!    !  Pole
!    do i=1-2, 1-1
!       loc_Eta (i,j1-joff:j2-joff) = loc_Eta (1,j1-joff:j2-joff)
!       loc_didt(i,j1-joff:j2-joff) = loc_didt(1,j1-joff:j2-joff)
!       loc_djdt(i,j1-joff:j2-joff) = loc_djdt(1,j1-joff:j2-joff)
!    end do
!    !  Equator
!    do i=isize+1,isize+2
!       loc_Eta (i,j1-joff:j2-joff) = loc_Eta (isize,j1-joff:j2-joff)
!       loc_didt(i,j1-joff:j2-joff) = loc_didt(isize,j1-joff:j2-joff)
!       loc_djdt(i,j1-joff:j2-joff) = loc_djdt(isize,j1-joff:j2-joff)
!    end do
!    !  Periodic
!    loc_Eta (-1:isize+1,-1:0) = loc_Eta (-1:isize+1,jsize-4:jsize-3)
!    loc_didt(-1:isize+1,-1:0) = loc_didt(-1:isize+1,jsize-4:jsize-3)
!    loc_djdt(-1:isize+1,-1:0) = loc_djdt(-1:isize+1,jsize-4:jsize-3)
!    loc_Eta (-1:isize+1,jsize-joff:jsize-joff+1) = loc_Eta (-1:isize+1,1:2)
!    loc_didt(-1:isize+1,jsize-joff:jsize-joff+1) = loc_didt(-1:isize+1,1:2)
!    loc_djdt(-1:isize+1,jsize-joff:jsize-joff+1) = loc_djdt(-1:isize+1,1:2)
!    
     !Call clawpack
     FirstTime=.true.
!    print*,'calling clawpack with k=',kc
!    write(6,*)FirstTime,T1,T2,xlower,xupper,ylower,yupper,CLAWiter
!    write(6,*)isize,jsize
!    istop=-1
!    if(istop.eq.-1)stop
     CALL Claw2ez (FirstTime, T1,T2, xlower,xupper, ylower,yupper, &
                   CLAWiter, 2,isize-1+1,jsize-3, &
                   loc_Eta, loc_didt, loc_djdt, loc_rate)
     FirstTime=.false.
 
     !Copy out
     DO j = j1, j2
     DO i = imin_j(j)+1, isize-1
        eeta (i, j, kc) = loc_Eta (i, j-joff)
     END DO
     END DO
     DO j = j1, j2
        IF (veff(imin_j(j+1),j+1)-veff(imin_j(j-1),j-1) < 0.0) THEN
           eeta (imin_j(j),j,kc) = loc_eta (imin_j(j),j-joff)
        END IF
     END DO
!
! floor eeta 12/06 frt
     max_eeta = maxval(eeta(:,:,kc))
     eeta(:,:,kc) = MAX(eps*max_eeta,eeta(:,:,kc))
     CALL Circle (eeta(:,:,kc))
!
  END DO


  RETURN
!
! CONTAINS
! !
!   FUNCTION Ratefn (fudgx, alamx, sinix, birx, vmx, xmfact)
!     IMPLICIT NONE
!     REAL (rprec), INTENT (IN) :: fudgx,alamx,sinix,birx,vmx,xmfact
!     REAL (rprec)              :: Ratefn
! !                                                                       
! !   Function subprogram to compute precipitation rate
! !   Last update:  04-04-88
! !
!     Ratefn = 0.0466_rprec*fudgx*SQRT(ABS(alamx))*(sinix/birx)*vmx**2
!     Ratefn = xmfact * ratefn
!     RETURN
!   END FUNCTION Ratefn
END SUBROUTINE Move_plasma_grid_NEW

FUNCTION Ratefn (fudgx, alamx, sinix, birx, vmx, xmfact)
  IMPLICIT NONE
  REAL (rprec), INTENT (IN) :: fudgx,alamx,sinix,birx,vmx,xmfact
  REAL (rprec)              :: Ratefn
  !                                                                       
  !   Function subprogram to compute precipitation rate
  !   Last update:  04-04-88
  !
  Ratefn = 0.0466_rprec*fudgx*SQRT(ABS(alamx))*(sinix/birx)*vmx**2
  Ratefn = xmfact * ratefn
  RETURN
END FUNCTION Ratefn


  SUBROUTINE Deriv_i_NEW (array, isize, jsize, j1, j2, imin_j, derivi)
!   USE Rcm_mod_subs, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: isize, jsize, j1, j2, imin_j(jsize)
    REAL (rprec), INTENT (IN) :: array (isize,jsize)
    REAL (rprec), INTENT (OUT) :: Derivi (isize,jsize)
!
    INTEGER (iprec) :: i, j
!
    DO j = 1, jsize
    DO i = 1, isize
       IF (i == 1) THEN
          Derivi(i,j) = -1.5*array(i,j) + 2.0*array(i+1,j) - 0.5*array(i+2,j)
       ELSE IF (i == isize) THEN
          Derivi (i,j) =  +1.5*array(i,j) - 2.0*array(i-1,j) + 0.5*array(i-2,j)
       ELSE
          Derivi(i,j) = 0.5*(array(i+1,j)-array(i-1,j))
       END IF
    END DO
    END DO
    RETURN
  END SUBROUTINE Deriv_i_NEW
  
  SUBROUTINE Deriv_j_NEW (array, isize, jsize, j1, j2, imin_j, derivJ)
!   USE rcm_mod_subs, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: isize, jsize, j1, j2, imin_j(jsize)
    REAL (rprec), INTENT (IN) :: array (isize,jsize)
    REAL (rprec), INTENT (OUT) :: Derivj (isize,jsize)
!
    INTEGER (iprec) :: i, j
!
    DO j = j1, j2
    DO i = 1, isize
       Derivj (i,j) = (array(i,j+1)-array(i,j-1))*0.5
    END DO
    END DO
    CALL Circle (Derivj)
    RETURN
  END SUBROUTINE Deriv_j_NEW
!-------------------------------------
    FUNCTION Gntrp_2d_ang (array, bi, bj, ikind)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: ikind
    REAL (rprec), INTENT (IN)    :: array (:,:), bi, bj
    REAL (rprec)                 :: Gntrp_2d_ang
!
    INTEGER (iprec) :: ii, jj, ni, nj
    REAL (rprec)    :: fi,fj,a1,a2,v_1, v_2
!
!
!   Prepare indices for interpolation:
!
    ni = SIZE (array, 1)
    nj = SIZE (array, 2)
    IF (bj < 1.0 .OR. bj > nj) write(*,*) 'BJ out of range in gntrp_2d_ang'
    IF (bi < 1.0 .OR. bi > ni) write(*,*) 'BI out of range in gntrp_2d_ang'
!
    ii = INT (bi)
    fi = REAL (ii,rprec)
    jj = INT (bj)
    fj = REAL (jj,rprec)
    IF (ii == ni) ii = ii - 1
    IF (jj == nj) jj = jwrap   ! periodicity in J
!                                                                       
    v_1 = array (ii,jj)
    v_2 = array (ii+1,jj)
!   IF (ii < imin_j(jj))   v_1 = array(imin_j(jj),jj)
!   IF (ii+1 < imin_j(jj)) v_2 = array(imin_j(jj),jj)
    a1 = (1.0-(bi-fi))*v_1 + (bi-fi)*v_2
!
    v_1 = array (ii  , jj+1)
    v_2 = array (ii+1, jj+1)
!   IF (ii   < imin_j(jj+1)) v_1 = array (imin_j(jj+1),jj+1)
!   IF (ii+1 < imin_j(jj+1)) v_2 = array (imin_j(jj+1),jj+1)
    a2 = (1.0-(bi-fi))*v_1 + (bi-fi)*v_2
!                                                                       
    IF (ikind == 1) THEN
       IF (jj+1 == nj) a2 = a2 + pi_two
       IF (jj  == jwrap)      a1 = zero
       IF (jj+1 == jwrap)      a2 = a2 + pi_two      ! sts, feb 22
    END IF
!
    Gntrp_2d_ang = (1.0 - (bj-fj)) * a1 + (bj-fj) * a2
!
    RETURN
    END FUNCTION Gntrp_2d_ang
!
!
!
!
!
    FUNCTION Intrp_2d_grid (array, bi, bbj, jwrap, imin_j)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: jwrap, imin_j(:)
    REAL (rprec), INTENT (IN)    :: array (:,:), bi, bbj
    REAL (rprec)                 :: Intrp_2d_grid
!
    INTEGER (iprec) :: ii, jj, imax_array, jmax_array
    REAL (rprec)    :: bj, ca, cb, cc, cd, di, dj
!
!
    imax_array = SIZE (array, 1)
    jmax_array = SIZE (array, 2)
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    bj  = Bjmod ( bbj,jwrap,jmax_array)
    jj  = INT (bj)
!
!   (i,j)---------------------(i,j+1)
!   | A          |               B  |
!   |            |                  |
!   |            |di                |
!   |   dj       |       1-dj       |
!   |---------(bi,bj)---------------|
!   |            |                  |
!   |            |                  |
!   |            |1-di              |
!   |            |                  |
!   | C          |              D   |
!   (i+1,j)-----------------(i+1,j+1)
!
    di = bi - ii
    dj = bj - jj
    IF (ii < imin_j(jj)) THEN
       ca = 0.0
    ELSE
       ca = (one-di)*(one-dj) 
    END IF
    IF (ii+1 < imin_j(jj)) THEN
       cc = 0.0
    ELSE
       cc = di * (one-dj)
    END IF
    IF (ii < imin_j(jj+1)) THEN
       cb = 0.0
    ELSE
       cb = (one-di)*dj
    END IF
    IF (ii+1 < imin_j(jj+1)) THEN
       cd = 0.0
    ELSE
       cd = di*dj
    END IF
!                                                                       
    Intrp_2d_grid = ca*array(ii,jj)+cb*array(ii,jj+1)+ &
                      cc*array(ii+1,jj)+cd*array(ii+1,jj+1)
    Intrp_2d_grid = Intrp_2d_grid / (ca+cb+cc+cd)
!
    RETURN
    END FUNCTION Intrp_2d_grid
    FUNCTION Interp_1d (array, bi )
    IMPLICIT NONE
    REAL (rprec), INTENT (IN)    :: array (:), bi
    REAL (rprec)                 :: Interp_1d
!                                                                       
!   This function subprogram interpolates 1-dim. ARRAY to return
!   the value of ARRAY at the non-integer point BI.
!
!   Stanislav: if Bi < 1, then the array is extrapolated
!              linearly based on the values A(1,:) and A(2,:).
!              If Bi > imax, then array is linearly
!              extrapolated on the values A(imax-1,:) and A(imax,:).
!
    INTEGER (iprec) :: ii, imax_array
    REAL (rprec)    :: fi
!
    imax_array = SIZE (array)
!
    IF (bi < one .OR. bi > imax_array) STOP 'OUT OF BOUNDS IN INTERP_1D'
!
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
    Interp_1d = (one - (bi-fi) ) * array (ii) + (bi-fi) * array (ii+1)
    RETURN
    END FUNCTION Interp_1d



    FUNCTION Interp_2d (array, bi, bj)
    IMPLICIT NONE
    REAL (rprec), INTENT (IN)    :: array (:,:), bi, bj
    REAL (rprec)                 :: Interp_2d
!
    INTEGER (iprec) :: ii, jn, jj, jp1, imax_array, jmax_array
    REAL (rprec)    :: fi,fj,a1,a2
!
!
!   Prepare indices for interpolation:
!
    imax_array = SIZE (array, 1)
    jmax_array = SIZE (array, 2)
!
!
    IF (bi < one .OR. bi > imax_array .OR. bj < 1 .OR. bj > jmax_array) &
        STOP 'OUT OF BOUNDS IN INTERP_2D'
!
!
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
!   Decide which interpolation to perform and proceed:
!
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp. of 2-d array
!
      Interp_2d = (one-(bi-fi)) * array(ii,jn) + (bi-fi)*array(ii+1,jn)
!
    ELSE    !        2-d interpolation of 2-d array:
!
!         If jwrap <= bj < jmax, then jwrap-1 <= INT(bj) <= jmax-1
!         and jwrap <= INT(bj)+1 <= jmax
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (one-(bi-fi))*array(ii,jj)  + (bi-fi)*array(ii+1,jj)
       a2 = (one-(bi-fi))*array(ii,jp1) + (bi-fi)*array(ii+1,jp1)
!                                                                       
       Interp_2d = (one - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Interp_2d




    FUNCTION Interp_2d_of3d (array, bi, bj, index_3)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: index_3
    REAL (rprec), INTENT (IN)    :: array (:,:,:), bi, bj
    REAL (rprec)                 :: Interp_2d_of3d
!                                                                       
!   This is the same as Gntrp_2d but for a 3-dim array, see comments for
!   Gntrp_2d. A separate function is needed since if Gntrp_2d were used,
!   then we would need to pass array sections (the other option is Fortran
!   77 style of passing an offset array, but that should be avoided for
!   compiler checking and parallelization reasons).
!
    INTEGER (iprec) :: ii, jn, jj, jp1, imax_array, jmax_array, kmax_array
    REAL (rprec)    :: fi,fj,a1,a2
!
!
!   Prepare indices for interpolation:
!
    imax_array = SIZE (array, 1)
    jmax_array = SIZE (array, 2)
    kmax_array = SIZE (array, DIM = 3)
!
!
    IF (bi < one .OR. bi > imax_array .OR. bj < one .OR. bj > jmax_array) THEN
       WRITE (*,*) 'OUT OF BOUNDS IN INTERP_2D_OF_3D'
       STOP
    END IF
!
!
    IF (index_3 > kmax_array .OR. index_3 < 1) STOP 'INTRP_2D_OF3D: index_3 OUT OF RANGE'
    ii = MAX (1, MIN (INT (bi), imax_array-1))
    fi = REAL (ii,rprec)
!                                                                       
!                                                                       
!   Decide which interpolation to perform and proceed:
!
    jn    = NINT (bj)
    IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp. of 2-d array
!
      Interp_2d_of3d = (one-(bi-fi)) * array(ii,jn,index_3) + &
                  (bi-fi)*array(ii+1,jn,index_3)
!
    ELSE    !        2-d interpolation of 2-d array:
!
!         If jwrap <= bj < jmax, then jwrap-1 <= INT(bj) <= jmax-1
!         and jwrap <= INT(bj)+1 <= jmax
!
       jj  = INT (bj)
       fj  = REAL (jj,rprec)
       jp1 = jj + 1
!
       a1 = (one-(bi-fi))*array(ii,jj,index_3)  + (bi-fi)*array(ii+1,jj,index_3)
       a2 = (one-(bi-fi))*array(ii,jp1,index_3) + (bi-fi)*array(ii+1,jp1,index_3)
!                                                                       
       Interp_2d_of3d = (one - (bj-fj)) * a1 + (bj-fj) * a2
!
    END IF
    RETURN
    END FUNCTION Interp_2d_of3d
!
!   7. A few routines needed for circulariation and normalization in J:
!
    FUNCTION Bjmod_real (bj, jwrap, jsize)
    IMPLICIT NONE
    REAL (rprec),    INTENT (IN) :: bj
    INTEGER(iprec),  INTENT (IN) :: jwrap, jsize
    REAL (rprec)                 :: Bjmod_real
!_____________________________________________________________________________
!   last update: 11-28-84               by:rws
!                                                                       
!   this function subporgram returns bjmod with a value
!   between jwrap and jmax-1. In RCM, arrays in j (local time angle)
!   are dimensioned from 1 to jsize, but the grid wraps around and
!   overlaps such that array (jwrap) = array (jsize)
!   array (jwrap-1) = array (jsize-1), etc. In other words, only
!   elements from j=jwrap to j=jsize-1 are unique. This function takes
!   a non-integer j index, BJ, and makes sure that it is larger or
!   equal to jwrap but smaller than jsize-1. Then when array is interpolated
!   on two elements, j is never larger than jsize or smaller than jwrap-1.
!   For the case of jwrap = 3 and a 1-dim array, this looks like:
!
!   j-value:   1  2  3  4  5              jsize-2   jsize-1   jsize
!              x  x  x  x  x .................x        x     x
!              |  |  |                        |        |     |
!              |  |   --------->---------->---|--------|-----
!              |  --------------->------------|-->-----
!               ---------->----------->-------
!
!   Dependency:  none
!
    Bjmod_real = bj
!                                                                       
!   do_1: DO
!      IF (Bjmod_real < REAL (jsize - 1,rprec)) EXIT
!      Bjmod_real = Bjmod_real - REAL (jsize - jwrap,rprec)
!   END DO do_1
!
!   do_2: DO
!      IF (Bjmod_real >= REAL (jwrap,rprec)) EXIT
!      Bjmod_real = Bjmod_real + REAL(jsize - jwrap,rprec)
!   END DO do_2
!
bjmod_real = MODULO(bj-REAL(jwrap),REAL(jsize-jwrap-1)) + REAL(jwrap) 
    RETURN
    END FUNCTION Bjmod_real
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    FUNCTION Bjmod_int (bj, jwrap, jsize)
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: bj
    INTEGER (iprec), INTENT (IN) :: jwrap, jsize
    INTEGER (iprec)              :: Bjmod_int
!_____________________________________________________________________________
!   last update: 11-28-84               by:rws
!                                                                       
!   this function subporgram returns bjmod with a value
!   between jwrap and jmax-1. In RCM, arrays in j (local time angle)
!   are dimensioned from 1 to jsize, but the grid wraps around and
!   overlaps such that array (jwrap) = array (jsize)
!   array (jwrap-1) = array (jsize-1), etc. In other words, only
!   elements from j=jwrap to j=jsize-1 are unique. This function takes
!   a non-integer j index, BJ, and makes sure that it is larger or
!   equal to jwrap but smaller than jsize-1. Then when array is interpolated
!   on two elements, j is never larger than jsize or smaller than jwrap-1.
!   For the case of jwrap = 3 and a 1-dim array, this looks like:
!
!   j-value:   1  2  3  4  5              jsize-2   jsize-1   jsize
!              x  x  x  x  x .................x        x     x
!              |  |  |                        |        |     |
!              |  |   --------->---------->---|--------|-----
!              |  --------------->------------|-->-----
!               ---------->----------->-------
!
!   Dependency:  none
!
    Bjmod_int = bj
!                                                                       
    do_1: DO
       IF (Bjmod_int > jsize - 1) THEN
          Bjmod_int = Bjmod_int - (jsize - jwrap)
       ELSE
          EXIT do_1
       END IF
    END DO do_1
!
    do_2: DO
       IF (Bjmod_int < jwrap) THEN
           Bjmod_int = Bjmod_int + (jsize - jwrap)
       ELSE
           EXIT do_2
       END IF
    END DO do_2
!
    RETURN
    END FUNCTION Bjmod_int
!
!
    END MODULE Rcm_mod_subs
