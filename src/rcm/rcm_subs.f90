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
                    densrcm(isize,jsize)
    INTEGER (iprec) :: ipcp_type, ipot
!
!
!   Input PCP drop and its current value:
    INTEGER (iprec), ALLOCATABLE :: ivtime (:)
    REAL    (rprec), ALLOCATABLE :: vinput (:), vinput_phase(:)
    REAL    (rprec)              :: vdrop,      vdrop_phase
!
!
!   Input Kp values and its current value:
    INTEGER (iprec), ALLOCATABLE :: ikptime (:)
    REAL    (rprec), ALLOCATABLE :: Kpinput (:)
    REAL    (rprec)              :: Kp
!
!
!   Input ETAC values:
    INTEGER (iprec), ALLOCATABLE :: itime_etac (:)
    REAL (rprec),    ALLOCATABLE :: etac_inp(:,:)
!
    INCLUDE 'rcmdir.h'


    Logical :: IsCoupledExternally = .false.  ! flag to determine if RCM is standalone or not



    ! Variables for internal RCM timing:
    INTEGER(iprec) :: timer_start(10) = 0, timer_stop(10) = 0, count_rate
    REAL (rprec) :: timer_values (10)=0.0_rprec


    INTERFACE Read_array
       MODULE PROCEDURE Read_real_1d_array, Read_real_2d_array, Read_real_3d_array,&
                        Read_intg_1d_array, Read_intg_2d_array, Read_intg_3d_array
    END INTERFACE
!
    INTERFACE Write_array
       MODULE PROCEDURE Write_real_1d_array, Write_real_2d_array, Write_real_3d_array,&
                        Write_intg_1d_array, Write_intg_2d_array, Write_intg_3d_array
    END INTERFACE
!
    INTERFACE Gntrp
       MODULE PROCEDURE Gntrp_2d_ang
    END INTERFACE
!
    INTERFACE Interp
      MODULE PROCEDURE Interp_1d, Interp_2d, Interp_2d_of3d
    END INTERFACE
!
    INTERFACE Bjmod
       MODULE PROCEDURE Bjmod_int, Bjmod_real
    END INTERFACE
!
    INTERFACE Outp
       MODULE PROCEDURE Outp_real, Outp_integer, Outp_logical
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
      REAL (rprec) :: c (ncoeff, isize, jsize), c5w (isize, jsize)
!
!
      IF (IsCoupledExternally) then
         vdrop = (MAXVAL(vbnd) - MINVAL(vbnd))/1.0E+3
         vdrop_phase = 0.0
      ELSE
         vdrop = Get_vdrop    (ivtime,  vinput,  jtime)
         vdrop_phase = Get_vdrop_phase    (ivtime,  vinput_phase,  jtime)
      END IF
      Kp    = Get_kp       (ikptime, kpinput, jtime)
      CALL    Get_bfield   (ibtime,           jtime, itype_bf)
      CALL    Get_boundary (boundary, bndloc)
      IF (i_eta_bc == 1) THEN
         CALL    Get_eta_on_bndy (jtime)
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
      IF (icond == 1) THEN      ! active conductances:
!
         IF (ifloor)   CALL Floor_for_eflux ()
         IF (icorrect) CALL Correct_eflux   ()
         CALL Get_active_cond ( )
!
      ELSE IF (icond == 2) THEN ! use Hardy statistical model:
!
         CALL Get_hardy_cond ()
!
      ELSE IF (icond == 3) THEN
         pedlam = qtped
         pedpsi = qtplam
         hall   = qthall
      ELSE
            STOP 'COMPUT: icond not defined'
      END IF
!                                                                       
!
      IF (ibnd_type == 2) THEN
!
        IF (ipcp_type > 0 .AND. ipcp_type < 8) STOP 'IPCP_TYPE IS ONLY FOR HMR'
         deqdt = zero
         a (2) = boundary(2)%aa
         b (2) = boundary(2)%bb
         dx(2) = boundary(2)%xx
         dy(2) = boundary(2)%yy
         CALL Efield (ipcp_type, vdrop*1000, deqdt, a, b, dx, dy, colat,&
                     aloct, 0_iprec, -2_iprec, 1_iprec, v, vbnd)

      ELSE IF (ibnd_type == 1 .OR. ibnd_type == 3) THEN
!
         IF (ipcp_type < 11) STOP 'IPCP NOT RIGHT'
         vbnd = Get_v_on_boundary (ipcp_type)
         DO j = 1, jsize
            v (1:imin_j(j)-1, j) = vbnd (j)
         END DO
!
      ELSE IF (ibnd_type == 4) THEN
!        DO NOTHING, VBND IS ALREADY SET
         DO j = 1, jsize
            v (1:imin_j(j)-1,j) = vbnd (j)
         END DO
      ELSE
         STOP 'COMPUT: ibnd_type not implemented'
      END IF
!
!
      IF (ipot == -1)THEN ! use the mhd potential, ie do nothing
! the potential should have already have been loaded in in the rcm grid in torcm

      ELSE IF (ipot == 3) THEN        ! Old calling sequence for C:
!
         CALL Comput_coeff            (c)
         CALL Comput_c5_wind          (iwind, c5w )
         CALL Comput_c5_total         (c, c5w )
         CALL Comput_lowlat_boundary  (c)
         CALL Comput_highlat_boundary (c, vbnd)
         CALL Comput_v_Potnt3         (c, v)
!
      ELSE           !  Now new calling sequence :
! 
!
!        Note: new_cfive thinks c5w is the old c5w without d
!        denominator. Before activating winds, check this.  5/29/99
!
         CALL Comput_c5_wind          (iwind, c5w )
         CALL New_cfive               (c, c5w, birk)
         CALL New_coeff               (c, vbnd )
         CALL Comput_lowlat_boundary  (c )
         CALL Gmresm (bndloc, c, v)
!
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
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      FUNCTION Get_kp (ikptime, kpinput, jtime)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime, ikptime (:)
      REAL    (rprec), INTENT (IN) :: kpinput (:)
      REAL    (rprec)              :: Get_kp
!_____________________________________________________________________________
!     Subroutine to specify Kp index value at time jtime.
!     This is accomplished by interpolating kpinput in time.
!     rws     3/20/97; stanislav 5/28/99
!     If jtime <= ikptime(1) then kp = kpinput(1)
!     If jtime >  ikptime(nkpmax) then kp = kpinput(nkpmax)
!     all other cases--interpolated.
!_____________________________________________________________________________
!
      INTEGER (iprec) :: nkp, nkpmax
      REAL    (rprec) :: f
!
      nkpmax = SIZE (kpinput)
      DO nkp = 1, nkpmax
         IF (jtime <= ikptime (nkp) ) THEN
            IF (nkp == 1) THEN
               Get_kp = kpinput (1)
               RETURN 
            ELSE 
               f = REAL(jtime-ikptime(nkp-1), rprec)/ &
                   REAL(ikptime(nkp)-ikptime(nkp-1), rprec)
               Get_kp = (one - f) * kpinput (nkp - 1) + f * kpinput (nkp)
               RETURN 
            END IF 
         END IF 
      END DO 
      Get_kp = kpinput (nkpmax)
!
      RETURN 
      END FUNCTION Get_kp
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_bfield (ibtime, jtime, itype_bf )
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime, ibtime (:), itype_bf
!
!_____________________________________________________________________________
!
!     Subroutine to find magnetic field arrays r,p,be and vm 
!     at time itime by interpolating in time between 
!     precomputed magnetic field
!     models.  nbf is the number of precomputed bfield models
!     and ibtime is a vector giving the event times associated
!     with each of the models.
!
!     rws   3/19/97; stanislav 5/01/98
!    
!     Stanislav: nold_bf is initialized to 0; this is actually
!     enough to "save" its value, but just in case add SAVE 
!     attribute. nold_bf is 
!     incremented by 1 only when the appropriate sets of 
!     B-models are read in from the files for interpolation;
!     "appropriate" here is:
!      | first B-model if itime <= ibtime(1)
!      | last  B-model if itime >  ibtime(nbf)
!      | two B-models for ibtime(n-1) and ibtime(n) where n
!      |  is such that ibtime(n-1) < itime <= ibtime(n)
!_____________________________________________________________________________
!
!
      INTEGER (iprec), SAVE :: nold_bf = 0
      INTEGER (iprec) :: n, nn, Lrec, nbf, i, j
      REAL    (rprec) :: f, fstoff1, fstoff2, fdst1, fdst2, fmeb1, fmeb2, &
                         fclps1, fclps2
!
!
 IF (itype_bf == 2) THEN ! use friction code results
    Lrec = 1
    CALL Read_array (rcmdir//'rcmxmin_inp', LREC , label, ARRAY_2D=xmin, ASCI=asci_flag)
    CALL Read_array (rcmdir//'rcmymin_inp', LREC , label, ARRAY_2D=ymin, ASCI=asci_flag)
    CALL Read_array (rcmdir//'rcmvm_inp',   LREC , label, ARRAY_2D=vm  , ASCI=asci_flag)
    CALL Read_array (rcmdir//'rcmbmin_inp', LREC , label, ARRAY_2D=bmin, ASCI=asci_flag)
    RETURN
ELSE IF (itype_bf == 1) THEN ! interpolate HV
!
      nbf = SIZE (ibtime(:))
!
      IF (jtime <= ibtime(1)) THEN
!
         IF (nold_bf < 1) THEN
            Lrec = 1
            CALL Read_array (rcmdir//'rcmxmin_inp', LREC , label, ARRAY_2D = xmin, ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmymin_inp', LREC , label, ARRAY_2D = ymin, ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmvm_inp',   LREC , label, ARRAY_2D = vm  , ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmbmin_inp', LREC , label, ARRAY_2D = bmin, ASCI=asci_flag)
!
            fmeb   = label%real (12)
            fstoff = label%real (13)
            fdst   = label%real (14)
            fclps  = label%real (15)
            nold_bf   = 1
         END IF
!
      ELSE IF (jtime > ibtime(nbf)) THEN
!
         IF (nold_bf < nbf+1) THEN
!
            Lrec = nbf
            CALL Read_array (rcmdir//'rcmxmin_inp', LREC , label, ARRAY_2D = xmin, ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmymin_inp', LREC , label, ARRAY_2D = ymin, ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmvm_inp',   LREC , label, ARRAY_2D = vm  , ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmbmin_inp', LREC , label, ARRAY_2D = bmin, ASCI=asci_flag)
!
            fmeb   = label%real (12)
            fstoff = label%real (13)
            fdst   = label%real (14)
            fclps  = label%real (15)
            nold_bf   = nbf + 1
         END IF
!
      ELSE 
         nn = -999
         find_loop: DO n = 2, nbf
            IF (jtime <= ibtime(n)) THEN
               nn = n
               EXIT find_loop
            END IF
         END DO find_loop
         IF (nn == -999) STOP 'ibtime screwed up, stop in bfield'
!
         IF (nn /= nold_bf) THEN
!
            LREC = nn - 1
            CALL Read_array (rcmdir//'rcmxmin_inp', LREC , label, ARRAY_2D = x1, ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmymin_inp', LREC , label, ARRAY_2D = y1, ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmvm_inp',   LREC , label, ARRAY_2D = vm1, ASCI=asci_flag )
            CALL Read_array (rcmdir//'rcmbmin_inp', LREC , label, ARRAY_2D = b1, ASCI=asci_flag)
!
            fmeb1   = label%real (12)
            fstoff1 = label%real (13)
            fdst1   = label%real (14)
            fclps1  = label%real (15)
!
            LREC = nn
            CALL Read_array (rcmdir//'rcmxmin_inp', LREC , label, ARRAY_2D = x2, ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmymin_inp', LREC , label, ARRAY_2D = y2, ASCI=asci_flag)
            CALL Read_array (rcmdir//'rcmvm_inp',   LREC , label, ARRAY_2D = vm2, ASCI=asci_flag )
            CALL Read_array (rcmdir//'rcmbmin_inp', LREC , label, ARRAY_2D = b2, ASCI=asci_flag)
!
            fmeb2   = label%real (12)
            fstoff2 = label%real (13)
            fdst2   = label%real (14)
            fclps2  = label%real (15)
!
            nold_bf = nn
         END IF
!
         f = REAL(jtime-ibtime(nn-1), rprec) / &
             REAL(ibtime(nn)-ibtime(nn-1), rprec)
         xmin   = (one-f)*x1 + f*x2
         ymin   = (one-f)*y1 + f*y2
         bmin   = (one-f)*b1 + f*b2
         vm     = (one-f)*vm1 + f*vm2
         fstoff = (one-f)*fstoff1+f*fstoff2
         fmeb   = (one-f)*fmeb1+f*fmeb2
         fdst   = (one-f)*fdst1+f*fdst2
         fclps  = (one-f)*fclps1+f*fclps2
!
      END IF
!
      rmin = SQRT (xmin**2+ymin**2)
      pmin = ATAN2 (ymin, xmin)
!
      RETURN
ELSE IF (itype_bf == 3) THEN
!  DO NOTHING, IT IS PASSED THROUGH MODULE VARIABLES FROM SOMEWHERE ELSE
!  including rmin and pmin
      RETURN
ELSE
   STOP 'ILLEGAL BFIELD TYPE IN GET_BFIELD'
END IF
      END SUBROUTINE Get_bfield
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_eta_on_bndy (jtime )
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: jtime
!_____________________________________________________________________________
!
!_____________________________________________________________________________
!
!
      INTEGER (iprec), SAVE :: nold_t = 0
      INTEGER (iprec) :: n, nn, Lrec, n_t, kc, j, i, i_least
      REAL    (rprec) :: f
      REAL    (rprec), SAVE :: etac_1 (kcsize), etac_2 (kcsize)
!
!
      n_t = SIZE (itime_etac(:))
!
      IF (jtime <= itime_etac(1)) THEN
!
         etac = etac_inp (:,1)
!
      ELSE IF (jtime > itime_etac(n_t)) THEN
!
         etac = etac_inp (:,n_t)
!
      ELSE 
         nn = -999
         find_loop: DO n = 2, n_t
            IF (jtime <= itime_etac(n)) THEN
               nn = n
               EXIT find_loop
            END IF
         END DO find_loop
         IF (nn == -999) STOP 'ibtime screwed up, stop in get_eta_on_bndy'
!
         f = REAL(jtime-itime_etac(nn-1), rprec) / &
             REAL(itime_etac(nn)-itime_etac(nn-1), rprec)
         etac   = (one-f)*etac_inp(:,nn-1) + f*etac_inp(:,nn)
!
      END IF
!
!
      i_least = MIN(MINVAL(imin_j)-2,1)
!
      DO kc = 1, kcsize
         DO j = 1, jsize
            DO i = i_least, imin_j(j)
               eeta(i,j,kc) = etac(kc)
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE Get_eta_on_bndy
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Get_boundary (boundary, bndloc)
    IMPLICIT NONE
    TYPE (ellipse_def), DIMENSION (2), INTENT (OUT) :: boundary
    REAL (rprec), INTENT (IN OUT) :: bndloc (:)
!_____________________________________________________________________________
!   This function must be called AFTER calling BFIELD.
!
!   IBND_TYPE is an RCM control parameter, read from input file. It specifies
!             how to set the RCM high-latitude boundary.
!   IBND_TYPE = 1 is set boundary to ellipse in equatorial plane
!                following RAW document "Setup of First Maynard
!                Run" dated 3/5/97.
!            = 2 is set to ellipse in ionosphere for use with
!                Heppner-Maynard scaled model (EFIELD)
!            = 3 boundary at L=6.6 in the magn. equat. plane (10/24/200)
!            = 4 is to read directly from file (for friction code,
!                02/23/2001). No BOUNDARY is used.
!
!   Note that subroutine returns the parameters of ellipse,
!   but for different IMODE they have different meaning, namely:
!   IBND = 1:  AA, BB, XC, YC are in Re units, with XC, YC
!              being GSM coordinates, so XC>0 on dayside,
!              YC positive on duskside
!   IBND = 2:  AA, BB, XC, YC are in degrees measured from
!              north pole, displacements being >0 toward noon,
!              and dusk (same convention as in EFIELD).
!_____________________________________________________________________________
!
    INTEGER (iprec) :: j, idim, jdim, jind (3), n, ierr
    REAL (rprec) :: a, b, c, f_c, colat_bnd, theta (3), bndloc_tmp
    REAL (rprec), PARAMETER :: phi (3) = (/ 0.0_rprec, pi_by_two, pi /), &
!                                                  |     |        |
!                                                noon, dusk,    midnight
!
                       fdist (3) = (/0.95_rprec, 1.45_rprec, 1.9_rprec/)
!
    IF (ibnd_type /= 1 .AND. ibnd_type /= 2 .AND. &
        ibnd_type /= 3 .AND. ibnd_type /=4) STOP 'ILLEGAL VALUE FOR IBND_TYPE'
    idim = SIZE (xmin, DIM =1 )
    jdim = SIZE (xmin, DIM =2 )
!
!
!   I. Setting boundary in equatorial plane (also needed for IMODE=2).
!
!   1. Set ellipse in equatorial plane:
!
    IF (ibnd_type == 1) THEN
!
!      This (below) gives ellipse extending to 0.95*fstoff at noon,
!      2*fstoff at midnight, and 1.5*fstoff at dawn and dusk:
!
       boundary(1)%aa =  1.475_rprec * fstoff
       boundary(1)%bb =  1.500_rprec * fstoff
       boundary(1)%xx = -0.525_rprec * fstoff
       boundary(1)%yy =  zero
!
    ELSE IF (ibnd_type == 3) THEN
!
       boundary(1)%aa =  6.6_rprec
       boundary(1)%bb =  6.6_rprec
       boundary(1)%xx =  zero
       boundary(1)%yy =  zero
!
    ELSE IF (ibnd_type == 4) THEN
       ! receive it elsewhere via module association
       imin_j = CEILING (bndloc,rprec) 
       RETURN  
    ELSE
       STOP 'ILLEGAL VALUE OF IBND_TYPE' 
    END IF
!
!
    DO j = 1, jdim
!
!   2. Check if boundary is within the grid limits:
!
       IF (  (xmin(1,j)-boundary(1)%xx)**2 / boundary(1)%aa**2 + &
             (ymin(1,j)-boundary(1)%yy)**2 / boundary(1)%bb**2 &
             <= one .OR. &
             (xmin(idim,j)-boundary(1)%xx)**2 / boundary(1)%aa**2 + &
             (ymin(idim,j)-boundary(1)%yy)**2 / boundary(1)%bb**2 &
              >= one) THEN
          WRITE (*,'(A30,I3)') 'COULD NOT PLACE BOUNDARY AT J=',j
          STOP 'ABNORMAL STOP IN GETBND'
       END IF
!
!    3. Pinpoint the exact location of bndy by bisection:
!
       a = one
       b = REAL (idim, rprec)
       DO
          c   = half * (a+b)
          IF ( ABS ( Fequat_of_x (c, REAL(j,rprec)) ) < 100*machine_eps1) EXIT
          IF (       Fequat_of_x (c, REAL(j,rprec)) < zero) THEN
              b = c
          ELSE
              a = c
          END IF
       END DO
       bndloc_tmp = c
!
!
!     4. Adjust boundary if it is too close to a grid line:
!
!      IF (ABS(FLOOR(bndloc_tmp) - CEILING(bndloc_tmp)) < machine_eps2) THEN
!         bndloc(j) = bndloc_tmp
!         WRITE (*,*) 'BNDY ON INTEGER GRID LINE'
!      ELSE IF ( bndloc_tmp - FLOOR (bndloc_tmp) < machine_eps2 ) THEN
!         bndloc(j) = FLOOR (bndloc_tmp)
!         WRITE (*,*) 'adjusting bndy in GETBND, FLOOR'
!      ELSE IF ( CEILING(bndloc_tmp) - bndloc_tmp < machine_eps2) THEN
!         bndloc(j) = CEILING (bndloc_tmp)
!         WRITE (*,*) 'adjusting bndy in GETBND, CEILING'
!      ELSE
!         bndloc(j) = bndloc_tmp
!      END IF
    END DO
    imin_j = CEILING (bndloc) ! first grid point inside modeling region.
    IF (ibnd_type == 1 .OR. ibnd_type == 3) RETURN
!
!
!
!  II. Set up boundary in the ionosphere to an ellipse such that
!      it maps out to the equatorial plane to given locations at
!      noon, midnight and dusk.
!
!   1. Specify the locations in the equatorial plane in units of R_stoff:
!
!
!   2. Since the ellipse will be in the ionosphere, map the three points:
!
    DO n = 1, 3
       jind (n) = -1
       DO j = j1, j2
          IF ( ABS(aloct(1,j)-phi(n)) < machine_eps1) THEN
              jind (n) = j
          END IF
       END DO
       IF (jind(n) == -1) STOP 'UNABLE TO LOCATE ONE OF PNTS IN GETBND'
!
       a = REAL (idim,rprec)
       b = REAL (one, rprec)
       DO
         c = half * (a+b)
         f_c = Gntrp_2d_ang (rmin, c, REAL(jind(n),rprec),0_iprec)-fstoff*fdist(n)
         IF (ABS (f_c) < machine_eps2) EXIT
         IF (f_c < zero) THEN
            a = c
         ELSE
            b = c
         END IF
       END DO
       theta (n) = Gntrp_2d_ang (colat, c, REAL(j,rprec), 0_iprec)
    END DO
!
!   3. Compute parameters of ellipse:
!
    boundary(2)%bb = theta(2)
    boundary(2)%xx = half * (theta(1) - theta(3))
    boundary(2)%aa = theta(1) - boundary(2)%xx
    boundary(2)%yy = zero
!
!
!   4. From colatitudes of boundary points, estimate their I-values:
!
    DO j = 1, jdim
       colat_bnd = Thet (boundary(2)%aa, boundary(2)%bb, &
                         boundary(2)%xx, boundary(2)%yy, aloct(1,j) )
       a = one
       b = REAL (idim,rprec)
       DO
         c = half * (a+b)
         f_c = Gntrp_2d_ang (colat, c, REAL(j,rprec), 0_iprec) - colat_bnd
         IF (ABS (f_c) < machine_eps2) EXIT
         IF (f_c < zero) THEN
           a = c
         ELSE
           b = c
         END IF
       END DO
       bndloc_tmp = c
!
!     5. Adjust boundary if it is too close to a grid line:
!
!      IF (ABS(FLOOR(bndloc_tmp) - CEILING(bndloc_tmp)) < machine_eps2) THEN
!         bndloc(j) = bndloc_tmp
!         WRITE (*,*) 'BNDY ON INTEGER GRID LINE'
!      ELSE IF ( bndloc_tmp - FLOOR (bndloc_tmp) < machine_eps2 ) THEN
!         bndloc(j) = FLOOR (bndloc_tmp)
!         WRITE (*,*) 'adjusting bndy in GETBND, FLOOR'
!      ELSE IF ( CEILING(bndloc_tmp) - bndloc_tmp < machine_eps2) THEN
!         bndloc(j) = CEILING (bndloc_tmp)
!         WRITE (*,*) 'adjusting bndy in GETBND, CEILING'
!      ELSE
!         bndloc(j) = bndloc_tmp
!      END IF
    END DO
    imin_j = CEILING (bndloc) ! first grid point inside modeling region.
!
!
    boundary(2)%aa = boundary(2)%aa * RTD
    boundary(2)%bb = boundary(2)%bb * RTD
    boundary(2)%xx = boundary(2)%xx * RTD
    boundary(2)%yy = boundary(2)%yy * RTD
!
    RETURN
!
    CONTAINS
!
      FUNCTION Fequat_of_x (bi, bj)
      IMPLICIT NONE
      REAL (KIND=rprec), INTENT (IN) :: bi, bj
      REAL (KIND=rprec) :: Fequat_of_x
      REAL (KIND=rprec) :: xx, yy
      xx  = Gntrp_2d_ang (xmin, bi, bj, 0_iprec)
      yy  = Gntrp_2d_ang (ymin, bi, bj, 0_iprec)
      Fequat_of_x = (xx-boundary(1)%xx)**2 / boundary(1)%aa**2 + &
                    (yy-boundary(1)%yy)**2 / boundary(1)%bb**2 - one
      RETURN
      END FUNCTION Fequat_of_x
!
    END SUBROUTINE Get_boundary
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
      iedim_local = 1
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
               IF (sum1 (ie) > zero) THEN
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
               ELSE 
!                                                                       
!                 Case fudge=0: we want eflux=0 and eavg=0 for no precipitation.
!
                  eflux (i, j, ie) = zero
                  eavg  (i, j, ie) = zero
!
               END IF 
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
      RETURN
      END SUBROUTINE Floor_for_eflux
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Correct_eflux ()
      IMPLICIT NONE
      INTEGER (iprec) :: i, j, ikp_low, ikp_high
      REAL    (rprec) :: hardy_eflux_int, eflux_int, value_l, value_h, value,&
                         factor
      INTEGER (iprec), PARAMETER :: ie_ele = 1
!
!     IF (kp > 5.99) STOP ' kp too large in correct_eflux'
      ikp_low = INT (kp)
      ikp_high = ikp_low + 1
      IF (ikp_high > 6) THEN
         ikp_high = 6
         ikp_low  = 5
      END IF
      factor   = (kp - ikp_low)/(ikp_high-ikp_low)
      DO j = 1, jsize
!
!           1. Compute latitudinal integral of hardy's EFLUX:
!
            hardy_eflux_int = zero
            DO i = 2, isize-1
               CALL Elemod (ICASE = 1_iprec, IKP = ikp_low, &
                            GLAT = 90.0-colat(i,j)*RTD, &
                            AMLT = MODULO (12.0+aloct(i,j)*RTH,24.0), &
                            VALUE = value_l)
               CALL Elemod (ICASE = 1_iprec, IKP = ikp_high, &
                            GLAT = 90.0-colat(i,j)*RTD, &
                            AMLT = MODULO (12.0+aloct(i,j)*RTH,24.0), &
                            VALUE = value_h)
               value = value_l*(one-factor)+value_h*factor
               value = (10.0**value) *1.6E-09 * pi
               hardy_eflux_int = hardy_eflux_int + value*(colat(i,j)-colat(i-1,j))
            END DO
!
!
!           2. Compute latitudinal integral of uncorrected RCM's EFLUX:
!
            eflux_int  = zero
            DO i = imin_j(j), isize - 1
               eflux_int = eflux_int + eflux(i,j,1)*(colat(i,j)-colat(i-1,j))
            END DO
!
!
!           3. Make correction:
!
            IF (eflux_int > 0.0) THEN
               DO i = imin_j(j), isize
                  eflux (i,j,1) = eflux(i,j,1)*(0.5-0.3*SIN(aloct(1,j)))*&
                                  hardy_eflux_int / eflux_int
               END DO
            ELSE IF (eflux_int == 0.0) THEN
               IF (ANY(eflux(imin_j(j):isize,j,1) /=0.0)) THEN
                  STOP 'EFLUX_INT = ZERO BUT EFLUX IS NOT'
               ELSE
                  eflux (imin_j(j):isize,j,1) = eflux (imin_j(j):isize,j,1)
               END IF
            ELSE
               STOP 'EFLUX_INT IS NEGATIVE'
            END IF
!
      END DO
      DO i = 1, isize
      DO j = 1, jsize
         IF (eavg(i,j,ie_ele) < 1.0E-30) eflux(i,j,ie_ele) = 0.0
      END DO
      END DO

      RETURN
      END SUBROUTINE Correct_eflux
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_hardy_cond ()
      IMPLICIT NONE
      INTEGER (iprec) :: i,j, kp_low, kp_high
      REAL    (rprec) :: value, value_low, value_high, factor 
!
      kp_low = INT (kp)
      kp_high = kp_low + 1
      factor = (kp-kp_low)/(kp_high-kp_low)
      DO i = 1, isize
      DO j = 1, jsize
         CALL Elemod (ICASE = 3_iprec, IKP = kp_low, GLAT = 90.0-colat(i,j)*RTD, &
                      AMLT = MODULO(12.+aloct(i,j)*RTH,24.0), VALUE = value_low)
         CALL Elemod (ICASE = 3_iprec, IKP = kp_high, GLAT = 90.0-colat(i,j)*RTD, &
                      AMLT = MODULO(12.+aloct(i,j)*RTH,24.0), VALUE = value_high)
         value = value_low*(one-factor)+value_high*factor
         hall (i,j) = qthall(i,j) + two*value / sini(i,j)
         CALL Elemod (ICASE = 4_iprec, IKP = kp_low, GLAT = 90.0-colat(i,j)*RTD, &
                      AMLT = MODULO(12.+aloct(i,j)*RTH,24.0), VALUE = value_low)
         CALL Elemod (ICASE = 4_iprec, IKP = kp_high, GLAT = 90.0-colat(i,j)*RTD, &
                      AMLT = MODULO(12.+aloct(i,j)*RTH,24.0), VALUE = value_high)
         value = value_low*(one-factor)+value_high*factor
         pedpsi(i,j) = qtped(i,j) + two*value
         pedlam(i,j) = qtplam(i,j) + two*value/sini(i,j)**2
      END DO
      END DO
      RETURN
      END SUBROUTINE Get_hardy_cond
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Get_active_cond ( )
      IMPLICIT NONE
!
!______________________________________________________________________________
!
!  This subroutine calculates conductance enhancement
!  due to auroral electron precipitation and adds this 
!  to the quiet time conductances read in subroutine 
!  qtcond. This subroutine contains changes made for 
!  tjfr run and corrected formulas for conductances 
!  as put forth by robinson         
!  last update: 11-06-86                                                
!               02-07-96 frt - min_j array added                        
!                                                                       
!  Stanislav, april 14 1999: added arrays for precipitation
!             conductances so that they are smoothed and 
!             quiet-time conductances are not modified.
!             This was also accompanied by a change in BALGN4
!             that now runs from i = 1 not min_j.
!
!______________________________________________________________________________
!
      INTEGER (iprec) :: i, j
      REAL    (iprec) :: ezero, sigp, sigh
!
!
      pedpsi = zero
      pedlam = zero
      hall   = zero
!
!
      DO j = jwrap, jsize-1
!     DO i = Get_imin_for_grid(j)-1, isize
      DO i = imin_j(j), isize
        IF (eflux(i,j,ie_el) > 1.0E-6 .AND. eavg(i,j,ie_el) < 1.E-5) THEN
           WRITE (*,*) 'stopping in cond, see the code'
           WRITE (*,*) i,j,eflux(i,j,1),eavg(i,j,1), imin_j(j)
           STOP
        END IF 
        ezero = eavg (i, j, ie_el) / 1.0E3
        sigp  = SQRT(eflux(i,j,ie_el)) * 40.0 * ezero / (16.0 + ezero**2)
        sigh  = 0.45 * sigp * ezero**(0.85)
        pedpsi (i, j) = two * sigp
        pedlam (i, j) = two * sigp / (sini(i,j)**2)
        hall   (i, j) = two * sigh / sini (i, j)
      END DO
      END DO
!
      CALL Circle (pedpsi)
      CALL Circle (pedlam)
      CALL Circle (hall)
!
      CALL Smooth_j (pedpsi)
      CALL Smooth_j (pedlam)
      CALL Smooth_j (hall)
!                                                                       
      CALL Smooth_i (pedpsi)
      CALL Smooth_i (pedlam)
      CALL Smooth_i (hall)
!
      pedpsi (:,:) = pedpsi (:,:) + qtped (:,:)
      pedlam (:,:) = pedlam (:,:)+ qtplam (:,:)
      hall   (:,:) = hall   (:,:)+ qthall (:,:)
!
      RETURN
      CONTAINS
!
              SUBROUTINE Smooth_i (array)
              IMPLICIT NONE
              REAL(rprec), INTENT (IN OUT) :: array (:,:)
        !
              INTEGER (iprec) :: i,j,n, idim, jdim
              REAL (rprec), DIMENSION (SIZE(array,1),SIZE(array,2)) :: work
              idim = SIZE (array, DIM = 1)
              jdim = SIZE (array, DIM = 2)
        !
              DO n = 1, nsmthi
        !
                DO j = 1, jdim
!               DO  i = Get_imin_for_grid(j)+2, idim - 1
                DO  i = imin_j(j)+1, idim - 1
                   work(i, j) = (array(i-1,j) + four * array(i,j)+array(i+1,j))/six
                END DO
                work (imin_j(j), j) = array(imin_j(j), j)
                work (idim, j) = array (idim, j)
                END DO
        !
                DO j = 1, jdim
!               DO i = Get_imin_for_grid(j), idim
                DO i = imin_j(j), idim
                  array (i, j) = work (i, j)
                END DO
                END DO

        !
              END DO
        !
              CALL Circle (array)
        !
              RETURN
              END SUBROUTINE Smooth_i
        !
        !
        !
              SUBROUTINE Smooth_j (array)
              IMPLICIT NONE
              REAL(rprec), INTENT (IN OUT) :: array (:,:)
        !
              INTEGER (iprec) :: i,j,n, idim, jdim
              REAL (rprec) :: work (SIZE(array,1), SIZE (array,2))
        !
              idim = SIZE (array, DIM = 1)
              jdim = SIZE (array, DIM = 2)
        !
              DO n = 1, nsmthj
        !
                DO i = 1, idim
                DO j = j1, j2
                  work(i,j)=(array(i,j-1)+four*array(i,j)+array(i,j+1))/six
                END DO
                END DO
        !
                CALL Circle (work)
        !
                DO i = 1, idim
                DO j = 1, jdim
                  array (i, j) = work (i, j)
                END DO
                END DO
        !
              END DO
        !
              RETURN
              END SUBROUTINE Smooth_j
!
!
       
      END SUBROUTINE Get_active_cond
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      FUNCTION Get_V_on_boundary (ipcp_type)
      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: ipcp_type
      REAL (rprec) :: Get_v_on_boundary (jsize)
!
!
      INTEGER (iprec) :: j
      REAL    (rprec) :: r_eq, p_eq
!   
      IF (ipcp_type == 11) THEN
          DO j = 1, SIZE (aloct, DIM = 2)
             Get_v_on_boundary (j) = -vdrop * SIN(aloct(1,j)-vdrop_phase*HTR ) / two
          END DO
          Get_v_on_boundary = Get_v_on_boundary * 1.0E+3_rprec
!
      ELSE IF (ipcp_type == 13) THEN
!
!       Maynard and Chen [JGR, 1975]:
!
        DO j = 1, SIZE(aloct,DIM=2)
           r_eq = Gntrp_2d_ang (rmin, bndloc(j), REAL(j,rprec), 0_iprec)
           p_eq = Gntrp_2d_ang (pmin, bndloc(j), REAL(j,rprec), 1_iprec)
           Get_v_on_boundary (j) = (92.4_rprec / r_eq - &
              A_coeff_MC(Kp)*r_eq**2*SIN(p_eq)) * 1000.0_rprec
        END DO
!
      ELSE
         STOP 'VBOUND: IPCP_TYPE NOT IMPLEMENTED'
      END IF
      RETURN
!
      CONTAINS
      FUNCTION A_coeff_MC (Kp)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: Kp
      REAL (rprec) :: A_coeff_MC
      A_coeff_MC = 0.045_rprec / (one-0.159_rprec*Kp+0.0093_rprec*Kp**2)**3
      RETURN
      END FUNCTION A_coeff_MC
!
      END FUNCTION Get_v_on_boundary
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_Coeff (c)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN OUT) :: c (:,:,:)
!______________________________________________________________________________
!
! code based on subroutine coeff in spiro.agu83.fort         
! last update 06-25-85 by rws                                
!              02-07-96 frt                                 
!                - min_j replaced imin and min_j+1 replaces i1
!
! This subroutine computes the coefficients of the 
! discretized PDE of MI coupling, except for the  
! inhomogenious term. Formulas are from Jaggi and Wolf 1973
! JGR paper. The discretized PDE for the electrostatic 
! potential looks like:
!
! V(i,j) = C1(i,j)*V(i+1,j) + C2(i,j)*V(i-1,j)
!        + C3(i,j)*V(i,j+1) + C4(i,j)*V(i,j-1) + C5(i,j)
!         
! This subroutine computes the coefficients for grid
! points that are inside the high-latitude boundary
! as defined by MIN_J array; that is, for all points (i,j)
! such that i <= min_j(j). If the boundary is non-integer,
! then some coefficients will have to be modified. c1-c5
! will be recomputed in subroutine CASE3; however, that 
! subroutine requires that arrays A, B, and D be computed
! here and be available in CASE3 for all I values.
!
! STANISLAV: if the boundary coincides with an integer
!            I-value, then this subroutine computes 
!            coefficients for i = imin+1 to i=imax.
!            It does not compute the coefficients at
!            i = imin because we don't use the difference
!            equations there, but rather the Dirichlet
!            boundary condtion (value of V). 
!            On the row just inside of the boundary 
!            (i = imin + 1), we approximate first deri-
!            vatives in I by a 2-point forward difference
!            rather than 3-pt. This leads to a O(dlam)
!            accuracy approximation (compared to 
!            (O(dlam**2)) for other points, but this is
!            due to the conductivities changing sharply
!            (edge of auroral zone!), so the 2-pt diff
!            may simply be not representative of the 
!            derivative in question.
!
!-----------------------------------------------------------
!                                                             
!                                                            
!   this subroutine computes the coefficients c1,c2,c3 & c4.
!   these are coefficients of the elliptic magnetosphere-
!   ionosphere coupling  equation that is solved in potent.
!   computed values  of the coeffecients are stored in array c.  
!                                                             
!   this subroutine called from subroutine comput           
!
!______________________________________________________________________________
!
      INTEGER (iprec) :: i,j,k
      REAL (rprec) :: aa, bb, cc, dd, ee, ff, bmin, bc, hmin, hc, &
                      a (isize,jsize), b (isize,jsize), d (isize,jsize)
!                                                         
      c (:,:,:) = zero
!                                                          
      DO j = 1, jsize, jint
      DO i = 1, isize, iint
         a (i, j) = alpha (i, j) * pedpsi (i, j) / beta  (i, j) 
         b (i, j) = beta  (i, j) * pedlam (i, j) / alpha (i, j) 
         d (i, j) = two * ( b(i, j) / dlam**2 + a(i, j) / dpsi**2 )
      END DO
      END DO
      open(unit=2,file=rcmdir//'cond.dat',status='unknown')
       write(2,*)pedpsi,pedlam
      close(2)

!                                                        
      loop_30: DO  j = j1, j2, jint
         Loop_20: DO  i = imin_j(j), isize, iint
!
            IF (i < imin_j(j)) THEN
!
!              Definitely outside modeling region, skip point:
!
               CYCLE Loop_20
!
            ELSE

               IF (i < isize .AND. i > imin_j(j) + 1) THEN
!
!                 Strictly inside the modeling region,
!                 Use central differences for I-derivatives:
!
                  bb = b (i + iint, j) - b (i - iint, j)
                  ee = hall(i+iint,j)-hall(i-iint,j)
                  ee = signbe * ee
!
               ELSE IF (i == isize) THEN
!
!                 On the equatorial boundary,
!                 Use backward 3-pt difference for I-derivatives:
!
                  bb = three * b (i,j) - four * b (i-1,j) + b (i - 2, j)
                  ee = three*hall (i,j) - four * hall (i-1,j) + hall (i-2, j)
                  ee = signbe * ee
!
               ELSE
!
!                 On the second row of modeling region,
!                 Use forward 2-pt differences for I-derivatives:
!
                  bmin = two * b (i, j) - b (i + 1, j)
                  bc = half * b (i, j)
                  IF (bmin < bc) bmin = bc
                  bb = b (i + 1, j) - bmin
                  hmin = two * hall (i, j) - hall (i + 1, j)
                  hc = half * hall (i, j)
                  IF (ABS (hmin)  < ABS (hc) ) hmin = hc
                  ee = (hall (i + 1, j) - hmin)*signbe

               END IF
!
            END IF 
!                                                         
            cc = hall (i, j + jint) - hall (i, j - jint)
            cc = cc * signbe 
            dd = (bb - cc * dlam / dpsi) * qtr
            aa = a (i, j + jint) - a (i, j - jint) 
            ff = (aa + ee * dpsi / dlam) * qtr
            c (1, i, j) = (b (i, j) + dd) / (d (i, j) * dlam**2) 
            c (2, i, j) = (b (i, j) - dd) / (d (i, j) * dlam**2) 
            c (3, i, j) = (a (i, j) + ff) / (d (i, j) * dpsi**2) 
            c (4, i, j) = (a (i, j) - ff) / (d (i, j) * dpsi**2) 
!
         END DO loop_20
      END DO loop_30

      DO i = 1, isize
         DO k = 1, ncoeff
            c (k, i, j1 - 2) = c (k, i, j2 - 1) 
            c (k, i, j1 - 1) = c (k, i, j2) 
            c (k, i, j2 + 1) = c (k, i, j1) 
         END DO 
      END DO
!                                                        
!                                                           
      RETURN 
      END SUBROUTINE Comput_Coeff
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_C5_wind (iwind, c5w)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: iwind
      REAL (rprec), INTENT (IN OUT) :: c5w (:,:)
!                                                                       
      INTEGER (iprec) :: i, j
      REAL (rprec) :: th, th0, th1, th2, dth, dph, denom, dr0, dr1, dr2, &
                      djre, dp1, dp2, denom1, denom2, djph, dw (isize,jsize)
      REAL (rprec), PARAMETER :: bnorm = 1.0E-6_rprec

!
!
!     1. If iwind=0, c5w(i,j) is set equal to zero and                  
!        subroutine returns to calling program.       
!                                                                       
!
      IF (iwind == 0) THEN
         c5w (:,:) = zero
         RETURN
      ELSE
         STOP 'wind is not implemented yet'
      END IF
!
!                                                                       
!     2. If iwind.ne.0, calculation of c5w(i,j) in volts.               
!        Discretization of the rhs of div*(s*e)=birk-div*jw 
!        yields a term dw(i,j)*v(i,j). Then 
!        c5w(i,j)=-div*jw/dw(i,j).    
!        jw=part of ionospheric current due to action of
!        thermospheric winds.       
!                                                                       
!     2.1 calculation of dw(i,j)                                        
!                                                                       
      DO j = j1, j2, jint
      DO i = imin_j(j)+1, isize
         th = colat (i, j)
         IF (i == isize) THEN
            th0 = pi_by_two - colat (i, j)
            th1 = pi_by_two - colat (i - 1, j)
            th2 = pi_by_two - colat (i - 2, j)
            dth = - half*(one/COS(th2)**2-four/COS(th1)**2+three/COS(th0)**2)
         ELSE
            th1 = colat (i - 1, j)
            th2 = colat (i + 1, j)
            dth = half * (one/SIN(th1)**2 - one /SIN(th2)**2)
         END IF
         dph = half*(aloct(i,j+1) - aloct(i,j-1))
         IF (dph < zero) dph = dph + pi
!                                                                       
         dw (i, j) = two /ri*( &
            two*COS(th)*pedlam(i,j) / SIN(th)**2 / dth**2 + &
            SIN(th)**2*pedpsi(i,j) / (two*COS(th) ) / dph**2)
      END DO
      END DO
!                                                                       
!     2.2  calculation of -div*jw. meridional component.                
!        div*jw is multiplied by 1.e-6 to express bir in teslas.        
!
      DO j = j1, j2, jint
      DO i = imin_j(j)+1, isize
         IF (i == isize) THEN
!                                                                       
!        2.2.1 meridional part at i=imax. derivative is approximated
!              by a 3-point forward difference formula.
!                                                                       
            th0 = pi_by_two - colat (i, j)
            th1 = pi_by_two - colat (i - 1, j)
            th2 = pi_by_two - colat (i - 2, j)
            dth = - half*(one / COS(th2)**2 - four/COS(th1)**2 + &
                         three/COS(th0)**2)
            denom = two * dth
!
!           comment by yong on 7/26/90. introduce "signbe" for the
!           following lines of rcm where there is a hall:
            dr0 = cos (th0) * bnorm * bir (i, j) * ( &
                  pedlam (i,j) * pwe (i,j) + &
                  signbe * hall (i, j) * hwn (i, j) )
            dr1 = cos (th1) * bnorm * bir (i - 1, j) * ( &
                        pedlam (i-1,j) * pwe (i-1,j) + &
                        signbe * hall (i-1,j)*hwn(i-1,j) )
            dr2 = cos (th2) * bnorm * bir (i - 2, j) * ( &
                       pedlam(i-2,j)*pwe(i-2,j) + &
                       signbe * hall (i-2,j) * hwn (i-2,j))
!                                                                       
            djre = - (dr2 - four * dr1 + 3. * dr0) / denom
!                                                                       
         ELSE
!                                                                       
!           2.2.2 meridional part at i.lt.imax. derivative is
!                 approximated by central differences.
!                                                                       
            th1    = colat (i - 1, j)
            th2    = colat (i + 1, j)
            dth    = half*(one / sin (th1)**2 - one / sin (th2)**2)
            denom1 = two * dth
!                                                                       
            dr2    = bnorm * bir(i-1,j) * sin(th1) * &
                     ( pedlam (i-1,j) * pwe(i-1,j) + &
                       hall(i-1,j) * hwn(i-1,j) * signbe )
!                                                                       
            dr1 = bnorm * bir (i + 1, j) * sin (th2) * &
                  ( pedlam (i + 1, j) * pwe (i + 1, j) + &
                    signbe * hall (i + 1, j) * hwn (i + 1, j) )
!                                                                       
            djre = (dr2 - dr1) / denom1
         END IF
!
!        2.2.3 zonal part.derivative is approximated by
!              central differences.
!                                                                       
         th1 = colat (i, j - 1)
         th2 = colat (i, j + 1)
         dph = half * (aloct (i, j + 1) - aloct (i, j - 1) )
         IF (dph < zero) dph = dph + pi
         denom2 = two * dph
         dp2 = SIN(th2)**3 / (two * cos (th2) ) * &
               bnorm * bir (i, j + 1) * signbe * &
               (hall(i,j+1) * hwe(i,j+1) - pedpsi(i,j+1) * pwn(i,j+1))
         dp1 = (SIN(th1))**3 / (two * COS(th1) ) * &
               bnorm * bir (i,j-1) * signbe * &
               (hall(i,j-1) * hwe(i,j-1) - pedpsi(i,j-1) * pwn(i,j-1))
!
!        end of change for inserting "signbe" on 7/26/90
!                                                                       
         djph = (dp2 - dp1) / denom2
         c5w (i, j) = - (djre+djph) / dw (i, j)
      END DO
      END DO
!
      CALL Circle (c5w)
!
      RETURN 
      END SUBROUTINE Comput_C5_wind
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_C5_total (c, c5w)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: c5w (:,:)
      REAL (rprec), INTENT (IN OUT) :: c (:,:,:)
!
!
      INTEGER (iprec) :: i,j
      REAL (rprec)    :: d
!                                                                       
      DO j = 1, jsize
      DO i = imin_j(j), isize
         d = two * ( beta(i,j)  * pedlam(i,j) / (alpha(i,j) * dlam**2)  &
                   + alpha(i,j) * pedpsi(i,j) / (beta(i,j) *  dpsi**2) )
         IF (d <= 1.0e-30_rprec) THEN
            c (5, i, j) = zero
         ELSE 
            c (5, i, j) = alpha(i,j) * beta(i,j) * (ri**2) * birk(i,j) / d + &
                          c5w(i,j)
         END IF
      END DO
      END DO
!
      RETURN 
      END SUBROUTINE Comput_c5_total
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_lowlat_boundary (c)
      IMPLICIT NONE
      REAL (rprec), INTENT(IN OUT) :: c (:,:,:)
!                                                                       
!______________________________________________________________________________
!                                                                       
!  last update| 08-27-86        written by| g.a.mantjoukis     
!                                                                       
!  subroutine to compute equatorial bndy condition.based on             
!   mantjou.eqbndy.text                                                 
!                                                                       
!  current conservation requires that at i=imax,                        
!  v(imax,j)=c(1,imax,j)*v(imax+1,j)
!           +c(2,imax,j)*v(imax-1,j)            
!           +c(3,imax,j)*v(imax,j+1)
!           +c(4,imax,j)*v(imax,j-1)
!           +c(5,imax,j)
!  where v(imax+1,j) is outside the modeling region.                    
!                                                                       
!  the equatorial bndy condition gives an expression for
!  v(imax+1,j) in terms of quantities inside the modeling
!  region and is of the form    
!  v(imax+1,j)=ceq1*v(imax,j-1)+ceq2*v(imax,j)
!             +ceq3*v(imax,j+1)+ceq4*v(imax-1,j)+ceq5    
!  where ceq1 through ceq5 are calculated below.                        
!                                                                       
!  ss(j) is a cowling-type conductance (see mantjou.eqbndy.text)        
!       integrated over the cross-section of the equatorial band,at     
!       any given local time.                                           
!                                                                       
!  sw(j) is a wind-dependent quantity that contributes to ceq5          
!                                                                       
!                                                                       
!  to set bnd cond to no current across imax 
!  (ie., no eq electrojet) explicityly zero ss(j) for all j  
!                                                                       
!______________________________________________________________________________
!                                                                       
      INTEGER (iprec) :: i, j, n
      REAL (rprec) :: cf, ceq1, ceq2, ceq3, ceq4, ceq5, den
      REAL (rprec), PARAMETER :: bnorm = 1.0E-6_rprec
!
      i = isize
      DO j = j1, j2 
         cf = alpha (i,j) * dlam / beta (i,j) / dpsi / pedlam (i,j)
         ceq1 = cf * (signbe*hall (i, j) -  &
                half * (ss (j + 1) - four * ss (j) - ss (j - 1) ) / dpsi)
         ceq2 = - four * cf * ss (j) / dpsi
         ceq3 = cf * ( - signbe * hall (i, j) +  &
                half * (ss (j + 1) + four * ss (j) - ss (j - 1) ) / dpsi)
         ceq4 = one
         ceq5 = - two * ri * alpha (i, j) * dlam * bnorm * bir (i, j) &
                * (pwe(i,j) + signbe * hall(i,j) / pedlam (i,j) * hwn (i,j))
         ceq5 = ceq5 - cf * (sw (j + 1) - sw (j - 1) ) 
         den  = one - ceq2 * c (1, i, j)
!                                                                       
         c (5, i, j) = (c (5, i, j) + ceq5 * c (1, i, j) ) / den
         c (4, i, j) = (c (4, i, j) + ceq1 * c (1, i, j) ) / den
         c (3, i, j) = (c (3, i, j) + ceq3 * c (1, i, j) ) / den
         c (2, i, j) = (c (2, i, j) + ceq4 * c (1, i, j) ) / den
         c (1, i, j) = zero
      END DO 
!                                                                       
      DO n = 1, ncoeff 
         DO j = 1, jwrap - 1 
            c (n, isize, j) = c (n, isize, jsize - jwrap + j)
         END DO
         c (n, isize, jsize) = c (n, isize, jwrap)
      END DO 
!                                                                       
      RETURN 
      END SUBROUTINE Comput_lowlat_boundary
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!                                                                       
      SUBROUTINE Comput_highlat_boundary (c, vbnd)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: vbnd (:)
      REAL (rprec), INTENT (IN OUT) :: c (:,:,:)
!______________________________________________________________________________
!
! subroutine to determine coefficients to implement polar
!  boundary cond
! acm --- june 30,1995                                                  
!
! Stanislav: Jan 1999, corrected circularization of C array.
!
!---------------------------------------------------------
! please note that a,b,c,d,e now refer to the rcm code variables        
! alpha, beta are from rcm, my alpha -> alp1, beta -> bet1              
! will keep f as c5 (source)                                            
! will keep gamma                                                       
!                                                                       
!-----------------------------------------------------------
!
!   For each point inside the modeling region, there is a 
!   difference equation of the form
!
!   v(i,j) = c(1,i,j)*v(i+1,j) + c(2,i,j)*v(i-1,j) +
!          + c(3,i,j)*v(i,j+1) + c(4,i,j)*v(i,j-1) + c(5,i,j)
!
!   If a point (i,j) has all four neighbor points I+1,I-1,
!   J+1,J-1 inside the modeling region, then subroutine COEFF
!   (should have been called before) has already computed
!   coefficients c1-c5. Otherwise, this routine, APHI, 
!   will look for points that have at least one of the four
!   neighbor points outside the boundary, and recompute 
!   the coefficients for those points. Then, since values
!   of the potential on those points are not variables but
!   given boundary conditions, terms c(k,i,j)*v(i+/-1,j+/-1)
!   will be added to the constant term c(5,i,j), and c(k,i,j)
!   will be reset to zero. Subroutine CASE3 called from
!   this routine will compute c1-c4.
!
! 
      INTEGER (iprec) :: i, j, ip, l, im, ibig, k, min_j (jsize)
      REAL (rprec) :: alp1 (jsize), bet1 (jsize), gamma (jsize), vb (jsize), &
                      vc(jsize), f (isize, jsize), tol, rmin, rlf, rlb
      LOGICAL :: Lsflag (jsize)
!                                                                       
!
!
      tol = 1e-06 
!                                                                       
! 
!    This subroutine will assume (from now on) that on
!    input, min_j(j) > ain(j) (strictly). 
!
     f (:, :) = zero
     min_j (:) = imin_j(:)
!
! alp1 = separation between physical boundary and 
!        nearest grid line (above)
!        alp1 = min_j(j)-ain(j)   
! ain = `grid' distance, not an integer                                 
!                                                                       
!    bet1: if for given point (i,j), the boundary is "above"
!          (that is, ain(j) < i), and for point (i,j+1)
!          boundary is "below" (ain(j+1)>i+1), then line
!          segment connecting (ain(j),j) and (ain(j+1),j+1)
!          intersects grid line (i,j)-(i,j+1). Bet1 is the
!          distance from (i,j) to that point of intersection,
!          so that bet1 <=1. If bet1=1, then (i,j+1) point
!          is actually on the boundary (so it is inside). 
!          In cases when ain(j+1)< i+1, we set bet1=1.
!    gamma: same as bet1, but for j-1
!                                                                       
!
!    In this subroutine, we look at points that are close to
!    but still inside the physical boundary of the RCM 
!    modeling region. In the RCM notation, for each J-value,
!    boundary is at ain(j) non-integer I-value, and min_j(j)
!    holds the I-value of the first point inside the boundary
!    (so min_j(j) > ain(j) ).
!
!    For each line of constant J, we definitely need to 
!    modify the coefficients of the difference equation
!    at point (min_j(j),j)--first point inside the
!    modeling region. If the boundary is "steep" enough
!    (min_j(j+1)-min_j(j)>2 ? and same for j-1), then
!    possibly we need to modify coefficients for other
!    points along the line of J, until we reach I-value
!    such that for point (i,j), both neighboring points
!    (j+1,j-1) are entirely inside the modeling region.
!                              
!
!    I. Treat the first grid point inside the boundary,
!    (min_j(j),j).
!
!    I.1 Compute the distances:
!
      DO j = 2, jsize - 1
!
!      At given J-value, the boundary crosses J-line at 
!      ain(j). We assume that min_j(j)=INT(ain(j))+1, so
!      that     ain(j) < min_j(j) <= ain(j)+1.
!      ALP1 is the difference  min_j(j)-ain(j), therefore,
!      0 < alp1 <= 1. Also it is clear that
!      min_j(j)-1 < ain(j) < min_j(j). 
!
!
         rmin = REAL (min_j (j),rprec )
!
!
!       For given J, consider the first grid point inside the
!       modeling region, and the corresponding difference
!       equation on this point. 3 different and mutually exclusive 
!       cases are possible: 
!
         IF ( (bndloc(j) - INT (bndloc(j)) ) < tol) THEN
!
!          Case 1: effectively integer boundary:
!          ain(j) is very close to but larger than min_j(j)-1,
!          so we reset min_j(j) to min_j(j)-1. Now we treat point
!          (min_j(j),j) as on the boundary (ain(j),j), and
!          take VIN(J) as the value of V(min_j(j),j), so the
!          difference equation is v(min_j(j),j) = 
!          c(5,min_j(j),j)=vin(j) and we don't need
!          coefficients c1-c4.
!
            Lsflag (j) = .TRUE. 
            min_j (j) = INT(bndloc(j))
            alp1 (j)  = zero
            bet1 (j)  = one  ! this will not be needed
            gamma (j) = one  ! this will not be needed
!
         ELSE IF ( min_j(j) - bndloc(j) < tol) THEN
!
!          Case 2: effectively integer boundary:
!          ain(j) is barely less than min_j(j), so
!          we don't modify anything, but will treat
!          ain(j) as integer boundary at min_j(j). Take
!          value of V at ain(j) as value of V at 
!          (min_j(j),j), so the difference equation is
!          v(min_j(j),j) = c(5,min_j(j),j)=vin(j), and 
!          we don't need coefficients c1-c4.
!
            Lsflag (j) = .TRUE. 
            alp1 (j)   = zero
            bet1 (j)   = one  ! this will not be needed
            gamma (j)  = one  ! this will not be needed
!
         ELSE 
!                                                                       
!         Case 3: non-integer boundary.
!         Need to compute the distances from this point
!         (min_j(j),j) to the boundary and get coefficients
!         c1-c4.
!
!         Obviously, if we look from (min_j(j),j) to 
!         (min_j(j)-1,j), that point is outside modeling
!         region. Therefore, will use (ain(j),j) as 
!         "i-1,j" neighbor. Get the distance to that point:
!
          alp1 (j) = REAL(min_j(j),rprec) - bndloc (j)
!
!
!         Look from (min_j(j),j) point to (min_j(j),j+1)
!         point and see where it is with respect to the
!         boundary: 
!
          IF (min_j(j+1) <= min_j(j)) THEN   
!
!           Grid point (min_j(j),j+1) is inside the modeling region,
!           so we will use it as a neighbor in the difference
!           equation. Set the distance from (min_j(j),j) to 
!           (min_j(j),j+1) to one (grid unit):
! 
            bet1 (j) = one
            vb (j)   = one
!
          ELSE 
!
!           Point (min_j(j),j+1) lies outside the boundary,
!           and cannot be used in the difference equation.
!           We will replace it with the point where line
!           connecting (min_j(j),j) and (min_j(j),j+1) 
!           intersects the segment of boundary connecting
!           (ain(j),j) to (ain(j+1),j+1). Get distance from
!           (min_j(j),j) to that point (bet1) and interpolate V
!           at that point from values of V at (ain(j),j) and 
!           (ain(j+1),j+1) (vb):
!
             bet1 (j) = alp1(j) / (bndloc (j + 1) - bndloc (j) )
             rlf      = SQRT ( (one - bet1(j))**2 +     &
                               (bndloc(j+1)-min_j(j))**2)
             rlb      = SQRT (bet1(j)**2 + alp1(j)**2) 
             vb (j)   = (vbnd (j + 1) * rlb + vbnd (j) * rlf) / (rlf + rlb)
!
          END IF 
!
!
!         Look from (min_j(j),j) point to (min_j(j),j-1)
!         point and see where it is with respect to the
!         boundary: 
!
          IF (min_j(j-1) <= min_j(j)) THEN  
!
!           Grid point (min_j(j),j-1) is inside modeling region,
!           so we will use it in the difference equation.
!           equation by central differences. Set the distance
!           to that point to one grid unit:
! 
            gamma (j) = one
            vc (j)    = one
!
          ELSE 
!
!           Grid point (min_j(j),j-1) lies outside the boundary,
!           and cannot be used in the difference equation.
!           We will replace it with the point where line
!           connecting (min_j(j),j) and (min_j(j),j-1) 
!           intersects the segment of boundary connecting
!           (ain(j),j) to (ain(j-1),j-1). Get distance from
!           (min_j(j),j) to that point (gamma) and interpolate V
!           at that point from known values of V at (ain(j),j) and 
!           (ain(j-1),j-1) (vc):
!
            gamma (j) = alp1(j) / (bndloc (j - 1) - bndloc (j) )
            rlb       = SQRT ( (one - gamma (j) ) **2 + &
                               (bndloc (j - 1) - min_j(j)) **2)
            rlf       = SQRT (gamma(j)**2 + alp1(j)**2) 
            vc (j)    = (vbnd (j) * rlb + vbnd (j - 1) * rlf) / (rlf + rlb)
          END IF 
!
         END IF 
!
      END DO 
!                                                                       
!
!     Circularize arrays that have been modified:
!
      Lsflag (1)    = Lsflag (jsize - 2)
      Lsflag (jsize) = Lsflag (3)
!                                                                       
      min_j (1)     = min_j (jsize - 2)
      min_j (jsize)  = min_j (3)
!                                                                       
      bet1 (1)      = bet1 (jsize - 2)
      bet1 (jsize)   = bet1 (3)
!                                                                       
      vb (1)        = vb (jsize - 2)
      vb (jsize)     = vb (3)
!                                                                       
      gamma (1)     = gamma (jsize - 2)
      gamma (jsize)  = gamma (3)
!                                                                       
      vc (1)        = vc (jsize - 2)
      vc (jsize)     = vc (3)
!                                                                       
!
!    I.2 Get the coefficients for the first point in the
!    modeling region.
!    4 cases to consider: which reduces to 2            
!    (1) 1 leg crossing |                           
!    (2)                | all three should be covered by case3
!    (3)                |                          
!    (4) physical boundary pt close to grid point      
!                                                                       
!    In the cases where gamma or beta is less than one,
!    values of V that are multiplied by the appropriate
!    coefficients are taken from boundary conditions (so
!    they are known) and therefore a term like this is a 
!    constant one. We put the contribution of these points
!    into the source term and then set the coefficient to 0.
!                                                                       
      DO j = 2, jsize - 1
!
         i = min_j (j) 
!
         IF (Lsflag (j) ) THEN 
!
!           Point (min_j(j),j) is on the high-latitude
!           boundary. Then, neglect all coefficients 
!           simply use the boundary condition at this point:
!
            c (1:4, i, j) = zero
            f (i, j)      = f (i, j) + vbnd (j)
!
         ELSE 
!
!          Using computed distances to the neighboring 
!          points, get coefficients for the difference
!          equation:
!
            CALL Case3 (c, alp1(j), bet1(j), gamma(j),i,j)
!
!          And rearrange the terms:
!
            IF (gamma (j) < one) THEN
               f (i, j)    = f (i, j) + c (4, i, j) * vc (j) 
               c (4, i, j) = zero
            END IF 
!
            IF (bet1 (j) < one) THEN
               f (i, j)    = f (i, j) + c (3, i, j) * vb (j) 
               c (3, i, j) = zero
            END IF 
!
            f (i, j)    = f (i, j) + c (2, i, j) * vbnd (j)
            c (2, i, j) = zero
!
         END IF 
!                                                                       
      END DO 
!
!
!    Circularize arrays:
!                                                                       
      DO k = 1, 5 
         c (k, min_j (1), 1)       = c (k, min_j (jsize - 2), jsize - 2)
         f (   min_j (1), 1)       = f (   min_j (jsize - 2), jsize - 2)
         c (k, min_j (jsize), jsize) = c (k, min_j (3), 3)
         f (   min_j (jsize), jsize) = f (   min_j (3), 3)
      END DO 
!                                                                       
!
!    II. Second loop in J is to treat points with i>min_j(j)
!    that might need modification of their coefficients too.
!                                                                       
      DO j = 2, jsize - 1
!
         i = min_j (j) 
         ip = min_j (j + 1) 
         im = min_j (j - 1) 
         ibig = MAX (ip,im)        ! get the higher one    
         alp1 (j) = one
!
!       We have considered point (min_j(j),j) in the 
!       first J-loop above. Now consider points with
!       i>min_j(j). For point (min_j(j)+1,j), the 
!       point on the "right side" (j+1) with minimum
!       I-value still inside the boundary is (min_j(j+1),j+1).
!       For points on the "left side" (j-1), it is
!       (min_j(j-1),j-1). We only need to worry if at least
!       one of the 2 j-neighbors of (min_j(j)+1,j) is outside
!       the boundary, i.e., either min_j(j+1) > min_j(j) + 1 or
!       min_j(j-1) > min_j(j) + 1. This is equivalent to
!       MAX (min_j(j+1),min_j(j-1)) > min_j(j) + 1.
!       The loop in L below will only execute under
!       these conditions. It will run from min_j(j)+1 to
!       the last point (max I-value) that still does not
!       have two j-neighbors inside the boundary.
!
         DO L = i + 1, ibig - 1 
!
            IF (ip > L) THEN 
!
!             J+1 grid neighbor is outside the boundary.
               rmin     = REAL (L,rprec)
               bet1 (j) = (rmin - bndloc (j) ) / &
                           (bndloc (j + 1) - bndloc (j) )
               rlf      = SQRT ( (one - bet1 (j) ) **2 + &
                                 (bndloc (j + 1) - rmin) **2)
               rlb      = SQRT (bet1 (j) **2 + &
                                (rmin - bndloc (j) ) **2)
               vb (j)   = (vbnd (j + 1) * rlb + vbnd (j) * rlf)&
                          / (rlf + rlb) 
!
            ELSE 
!
!             J+1 grid neighbor is within the boundary, don't
!             need to modify the coefficient.
!
               bet1 (j) = one
               vb (j)   = one
!
            END IF
!
!
            IF (im > L) THEN 
!
!             The grid neighbor at J-1 is outside the
!             boundary, will modify the coefficient, so
!             need the distances:
!
               rmin      = REAL (L,rprec)
               gamma (j) = (rmin - bndloc (j) ) / &
                           (bndloc (j - 1) - bndloc (j) )
               rlb       = SQRT ( (one - gamma (j) ) **2 +  &
                                  (bndloc (j - 1) - rmin) **2)
               rlf       = SQRT (gamma (j) **2 +  &
                                 (rmin - bndloc (j) ) **2)
               vc (j)    = (vbnd (j) * rlb + vbnd (j - 1) * rlf)&
                           / (rlf + rlb) 
!
            ELSE 
!
!             Grid neighbor at J-1 is inside the
!             boundary, don't need to modify coefficient.
!
               gamma (j) = one
               vc (j)    = one
!
            END IF 
!
!          Compute modified coefficients:
!
            CALL Case3 (c, alp1(j), bet1(j), gamma(j), L, j)
!
!          Rearrange coefficients:
!
            IF (gamma (j) < one) THEN
               f (L, j)    = f (L, j) + c (4, L, j) * vc (j) 
               c (4, L, j) = zero
            END IF
!
            IF (bet1 (j) < one) THEN
               f (L, j)    = f (L, j) + c (3, L, j) * vb (j) 
               c (3, L, j) = zero
            END IF
!
            c (1, L, 1) = c (1, L, jsize - 2)
            c (2, L, 1) = c (2, L, jsize - 2)
            c (3, L, 1) = c (3, L, jsize - 2)
            c (4, L, 1) = c (4, L, jsize - 2)
            f (   L, 1) = f (   L, jsize - 2)
!                                                                       
            c (1, L, jsize) = c (1, L, 3)
            c (2, L, jsize) = c (2, L, 3)
            c (3, L, jsize) = c (3, L, 3)
            c (4, L, jsize) = c (4, L, 3)
            f (   L, jsize) = f (   L, 3)
!
         END DO 
!
      END DO  
!                                                                       
!
!     Add source term, f, to c5:       
!
      c (5,:,:) = f (:,:) + c (5,:,:)
!
      RETURN 
      END SUBROUTINE Comput_highlat_boundary
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!                                                                       
      SUBROUTINE Case3 (c, alp1, bet1, gamma, i, j)
      IMPLICIT NONE
      INTEGER (iprec) ,INTENT (IN) :: i, j
      REAL (rprec), INTENT (IN OUT) :: c (:,:,:), alp1, bet1, gamma
!                                                                       
!----------------------------------------------------------
! code based on subroutine coeff in rcm.f, copied from on
! 23-may-95     
! acm --- june 30,1995                                                  
!----------------------------------------------------------
!                                                                       
! This subroutine computes the coefficients c1,c2,c3 & c4.
!   these are   
! coefficients of the elliptic magnetosphere-ionosphere coupling        
! equation that is solved in potent.  computed values                   
! of the coeffecients are stored in array c.                            
!                                                                       
! This subroutine called from subroutine aphi                          
!                                                                       
! This computes the coefficients for the inner boundary
! with possible 3 crossing 
!
! Stanislav: this is exact translation of the F77 case3 
!            routine, no changed made.
!-----------------------------------------------------------
!                                                                       
      REAL  (rprec) ::  h, k, dpf, dpb, arcm, brcm, t1, aij, bij, cij, &
                        dij, eij, aa, bb, cc, ee, hmin,  &
                        hc, denom1, denom2, &
                        a (isize,jsize), b (isize,jsize), d (isize,jsize)
!                                                                       
      IF (alp1 == zero) RETURN
!                                                                       
! Modify the differentials for the cases where the crossings
! are not at lines:
!
      k   = dlam 
      h   = dlam * alp1 
      dpf = dpsi * bet1 
      dpb = dpsi * gamma 
!                                                                       
!
!    Since c1-c4 for (i,j) will be defined anew, reset them:
!
      c (1, i, j) = zero
      c (2, i, j) = zero
      c (3, i, j) = zero
      c (4, i, j) = zero
!                                                                       
      a (i-1:i+1, j-1:j+1) = alpha (i-1:i+1, j-1:j+1) * &
                             pedpsi (i-1:i+1, j-1:j+1) / beta (i-1:i+1, j-1:j+1)
      b (i-1:i+1, j-1:j+1) = beta (i-1:i+1, j-1:j+1) * &
                             pedlam (i-1:i+1, j-1:j+1) / alpha (i-1:i+1, j-1:j+1)
!                                                                       
!    d(i,j) depends on the differences which will not be the
!    same as for interior points ---- use smp generated eij  
!                                                                       
      arcm = a (i, j) 
      brcm = b (i, j) 
!                                                                       
!    Calculate the differences needed in the c's : 
!    aa,bb,cc,dd,ee,ff. They are used in approximating
!    coefficients of the MI-coupling PDE and arise from       
!    approximating derivatives of "known" functions.
!    These need to be modified on the boundaries due 
!    to points being no longer in the physical region   
!                                                                       
!    First, look at lambda derivatives:
!
      IF ( ABS(alp1-one) < machine_eps1 ) THEN
!
!       It is an interior point (distance from (i,j) to
!       (ain(j),j) is at least one grid cell) and we can
!       use a central difference formula:
!
         bb = b (i + iint, j) - b (i - iint, j) 
         ee = (hall(i + iint, j) - hall(i - iint, j))*signbe
!
      ELSE 
!
!       for alp1<1 need forward difference
!
         bb = b(i + 1, j) - MAX (two*b(i,j)-b(i+1,j),half*b(i,j))
         hmin = 2. * signbe * hall (i, j) - signbe * hall (i + 1, j) 
         hc = half * signbe * hall (i, j)
         IF (ABS (hmin) .lt.ABS (hc) ) hmin = hc 
         ee = signbe*hall(i+1,j) - hmin 
!
      END IF
!                                                                       
!
!    Now need to look at psi derivatives.    
!    If bet1<1 then need backward difference -- must x by 2 
!    because cent. doesn't have the half. 
!    If gamma<1 then need forward difference     
!    If both are < 1 then take derivative equal to 0  
!
      IF (bet1 < one .AND. gamma < one) THEN
            cc = zero
            aa = zero
      ELSE IF (bet1 < one .AND. ABS(gamma-one) < machine_eps1 ) THEN
            cc = (hall(i,j)-hall(i,j-jint))*two*signbe
            aa = (a (i, j) - a (i, j - jint) ) * two
      ELSE IF (ABS(bet1-one) < machine_eps1 .AND. gamma < one) THEN
         cc = signbe * (hall (i, j + jint) - hall (i, j) ) * two
         aa = (a (i, j + jint) - a (i, j) ) * two
      ELSE IF (ABS(bet1-one) < machine_eps1 .AND. &
               ABS(gamma-one) < machine_eps1)                  THEN
         cc = signbe * (hall (i, j + jint) - hall (i, j - jint)) 
         aa = a (i, j + jint) - a (i, j - jint) 
      ELSE
         STOP 'nonsense in CASE3'
      END IF
!                                                                       
!
!    Put in dx parts, and factor of 2:
!
      aa = aa / dpsi * half
      bb = bb / dlam * half
      cc = cc / dpsi * half
      ee = ee / dlam * half
!                                                                       
!
!    now calculate coeffnts c1-c4 using smp generated formulas:
!                                                                       
      aij = (two * brcm + bb * h + ( - 1) * cc * h) / (k * (h + k) )
      bij = (-one)* (((-two) * brcm + bb * k + (-1) * cc * k) / (h * (h + k) ) )
      cij = (two * arcm + aa * dpb + dpb * ee) / (dpf * (dpb + dpf) )
      dij = (-one)* (((-two) * arcm + aa * dpf + dpf * ee) / (dpb * (dpb + dpf)))
!                                                                       
!    d(i,j) depends on the differences which will not be
!    the same as for interior points -- use smp generated eij:
!                                                                       
      arcm = a (i, j) 
      brcm = b (i, j) 
!
      denom1 = (h + k) * h * k 
      denom2 = (dpb + dpf) * dpf * dpb 
!                                                                       
      t1 = ( - two) * ( (brcm * h) / denom1) &
         + ( - two) * ( (brcm * k) / denom1) &
         + (cc * h**2) / denom1  &
         + ( - one) * ( (cc * k**2) / denom1) &
         + ( - one) * ( (dpb**2 * ee) / denom2) &
         + (dpf**2 * ee) / denom2                
!                                                                       
      eij = ( - one) * ( (aa * dpb**2) / denom2) &
          + (aa * dpf**2) / denom2  &
          + ( - two) * ( (arcm * dpb) / denom2) &
          + ( - two) * ( (arcm * dpf) / denom2) &
          + ( - one) * ( (bb * h**2) / denom1)  &
          + (bb * k**2) / denom1 &
          + t1                                                     
!                                                                       
      d (i, j) = - eij 
!                                                                       
!
!    Have to divide by d (multiplier of v(i,j) [e] in n.r. SOR)
!
      c (1, i, j) = aij / d (i, j) 
      c (2, i, j) = bij / d (i, j) 
      c (3, i, j) = cij / d (i, j) 
      c (4, i, j) = dij / d (i, j) 
!                                                                       
      RETURN 
      END SUBROUTINE Case3                          
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE New_coeff (c, vbnd)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN)     :: vbnd (:)
      REAL (rprec), INTENT (IN OUT) :: c    (:,:,:)
!
!     ..Local variables:
!
      LOGICAL :: pt_loc (isize,jsize)
      LOGICAL :: a2_loc, a3_loc, a4_loc
!
      REAL(rprec) :: alp1, bet1, gamma, lam_f, lam_b, psi_f, psi_b
!
      REAL(rprec) :: a (isize,jsize), b (isize,jsize)
      REAL(rprec) :: bb, ee, aa, cc, denom
      REAL(rprec) :: g_lam, g_psi, f_lam, f_psi
      REAL(rprec) :: cl1_a1, cl1_a2, cl1_p
      REAL(rprec) :: cl2_a1, cl2_a2, cl2_p
      REAL(rprec) :: cp1_a3, cp1_a4, cp1_p
      REAL(rprec) :: cp2_a3, cp2_a4, cp2_p
      REAL(rprec) :: v_right, v_left, dis_r_1, dis_r_2
      REAL(rprec) :: dis_L_1, dis_L_2 , bmin, bc, hmin, hc
!
      INTEGER (iprec), SAVE :: N_calls = 0
      INTEGER (iprec):: i, j, kindex
!
!         j-1    j     j+1        
!
!   i-1..  x     A2     x
!                |
!                |              
!   i ...  A4----P----A3
!                |
!                |
!   i+1..  x     A1     x
!
!     In this subroutine we presume that high-lat. bndy is
!     adequately specified before calling it, that is, 
!     ain(j)+1 > min_j(j) >= ain(j) and |min(j)-ain(j)| > 1E-6
!     and min_j(j) >=2 (check all this outside?)
!
!     .. Executable statements:
!
      N_calls = N_calls + 1
!
      IF (isize < 3) THEN
         STOP 'idim must be > 2 in NEW_COEFF'
      END IF
!
!
!     1. Run over ALL grid points and flag them as being 
!     inside the modeling region (includes integer boundary
!     crossings) with pt_loc = 1 or outside the modeling
!     region with pt_loc = 0:
!
      DO j = 1, jsize
         DO i = 1, isize
            IF (i >= bndloc(j)) THEN
               pt_loc(i,j) = .TRUE.
            ELSE 
               pt_loc(i,j) = .FALSE.
            END IF
         END DO
      END DO
! 
!
!-------------------------------------------------------
!     Preliminary calculations:
!
      DO i = 1, isize
         DO j = 1, jsize
!
            a (i,j) = alpha(i,j)*pedpsi(i,j)/beta(i,j)
            b (i,j) = beta(i,j)*pedlam(i,j)/alpha(i,j)
!
         END DO
      END DO
!------------------------------------------------------- 
!
!     2. Run over grid points inside boundaries. For each
!     point in the modeling region, generate coefficients
!     of the difference equation. Our PDE has the form:
!
!     -b*D2V/Dlambda^2 -a*D2V/Dpsi^2 + 
!     + (-Db/Dlambda + Dhall/Dpsi)*DV/Dlambda +
!     + (-Da/Dpsi - Dhall/Dlambda)*DV/Dpsi = RHS
!     or
!     g_lambda*D2V/Dlambda^2 + g_psi*D2V/Dpsi^2 +
!     f_lambda*DV/Dlambda    + f_psi*DV/Dpsi     = RHS
!
!     We need: to approximate derivatives of V by finite
!     differences, and approximate "coefficients" f_lambda,
!     f_psi, g_lambda, g_psi by finite differences. For 
!     approximating coefficients, we can only use neighbors
!     that are integer grid points since conductances are 
!     evaluated on the grid. For approximating derivatives of
!     V, we will use values of V on non-integer grid neighbors
!     computed from the boundary condition on V.
!
      Loop_j: DO j = jwrap, jsize-1
         Loop_i:  DO i = imin_j(j), isize
!
            IF (i < 2) THEN
               STOP 'I < 2 IN NEW_COEFF'
            END IF
!           IF (i < imin_j(j)) CYCLE
!
!            IF (.NOT.pt_loc(i,j)) CYCLE
!
!           For each grid point (i,j)=P in the modeling
!           region, we will need four neighbors.
!           Determine how many of the grid neighbors are
!           inside the modeling region, and flag them. 
!
            a2_loc = pt_loc(i - 1,j)
            a3_loc = pt_loc(i,j + 1)
            a4_loc = pt_loc(i,j - 1)
!
!
!           Determine distances from (i,j) to its 3 neighbors
!           (A_1 is always one grid unit away), in grid units:
!
            IF (a2_loc) THEN      ! A_2 is also grid point
               alp1 = one
            ELSE                  ! A_2 is a boundary point
               alp1 = REAL(i,rprec) - bndloc(j)
            END IF
!
            IF (alp1 < 1E-4) THEN 
!
!               This is a special case when the point (i,j) is
!               on the boundary (integer boundary crossing).
!               Handle as Dirichlet boundary condition:
!
                c(1:4,i,j) = zero
                c(5,i,j)   = vbnd(j)
                CYCLE Loop_i
            END IF
!
            IF (a3_loc) THEN      ! A_3 is a grid point
               bet1 = one
            ELSE                  ! A_3 is on the boundary
               bet1    = (i-bndloc(j)) / (bndloc(j+1)-bndloc(j))
               dis_r_1 = SQRT( (i-bndloc(j))**2 + bet1**2)
               dis_r_2 = SQRT( (bndloc(j+1)-i)**2 + (one-bet1)**2)
               v_right = (dis_r_1*vbnd(j+1)+dis_r_2*vbnd(j)) &
                         / (dis_r_1+dis_r_2)
            END IF
!
            IF (a4_loc) THEN      ! A_4 is a grid point
               gamma = one
            ELSE                  ! A_4 is on the boundary
               gamma = (i-bndloc(j))/(bndloc(j-1)-bndloc(j))
               dis_L_1 = SQRT( (i-bndloc(j))**2 + gamma**2)
               dis_L_2 = SQRT( (bndloc(j-1)-i)**2 + (one-gamma)**2)
               v_left  = (dis_L_1*vbnd(j-1)+dis_L_2*vbnd(j)) &
                         / (dis_L_1+dis_L_2)
            END IF
!
!
!           Approximate coefficients with lambda-derivatives:
!
            IF (ABS(alp1-one) < machine_eps1 ) THEN

               IF (i < isize) THEN
!           
!                 (i,j) is an interior grid point, and we
!                 can use central differences formula for
!                 lambda-derivatives:
! 
                  bb = b(i+1,j) - b(i-1,j)
                  ee = (hall(i+1,j) - hall(i-1,j))*signbe

               ELSE
!
!                 (i,j) in on low-latitude boundary, need
!                 to use backward differences for deriv:
!
                  bb = three * b(i,j) - four * b(i-1,j)+b(i-2,j)
                  ee = three * hall(i,j)-four * hall(i-1,j)+hall(i-2,j)
                  ee = ee * signbe

               END IF
!
            ELSE
!         
!              alp1 < 1, so "i-1,j" grid point is outside
!              the boundary, and we need forward difference:
!
!!!!           bb = b(i+1,j) - &
!!!!                MAX (0.5*b(i,j),2.*b(i,j)-b(i+1,j))
!!!!           hmin = 2.*hall(i,j) - hall(i+1,j)
!!!!           hc   = 0.5*hall(i,j)
!!!!           IF (ABS(hmin) < ABS(hc)) hmin = hc
!!!!!!         ee   = signbe*(hall(i+1,j)-hmin)
!st               bb = -three*b(i,j)+ four*b(i+1,j)-b(i+2,j)
!st               ee = -three*hall(i,j)+four*hall(i+1,j)-hall(i+2,j)
!st               ee = ee * signbe
              bmin = 2. * b (i, j) - b (i + 1, j)
              bc = 0.5 * b (i, j)
              IF (bmin < bc) bmin = bc
              bb = b (i + 1, j) - bmin
              hmin = 2. * hall (i, j) - hall (i + 1, j)
              hc = 0.5 * hall (i, j)
              IF (ABS (hmin)  < ABS (hc) ) hmin = hc
              ee = (hall (i + 1, j) - hmin)*signbe

!
            END IF
!
!
!           Approximate coefficients with psi-derivatives:
!
            IF (ABS(bet1-one) < machine_eps1 .AND. &
                ABS(gamma-one) < machine_eps1)            THEN
!
!               (i,j) is an inner point, can use central
!               differences:
!
                cc = signbe*(hall(i,j+1)-hall(i,j-1))
                aa = a(i,j+1) - a(i,j-1)
!
            ELSE IF (bet1 < one .AND. ABS(gamma-one) < machine_eps1) THEN
!
!               use backward difference, mult. by 2:
!
                cc = two*signbe*(hall(i,j)-hall(i,j-1))
                aa = two*(a(i,j)-a(i,j-1))
!
            ELSE IF (ABS(bet1-one) < machine_eps1 .AND. gamma < one) THEN
!
!               use forward difference, mult. by 2:
!
                cc = two * signbe*(hall(i,j+1)-hall(i,j))
                aa = two * (a(i,j+1)-a(i,j))
!
            ELSE
!
!               gamma and bet1 are < 1, set derivs to zero:
!
                aa = zero
                cc = zero
!
            END IF
!               
!
            g_lam       = - b(i,j)
            g_psi       = - a(i,j)
            f_lam       = - bb/dlam/2. + cc/dpsi/2. 
            f_psi       = - ee/dlam/2. - aa/dpsi/2. 
!
!
!           Approximate partial derivatives of V in the PDE
!           DV/Dlambda, DV/Dpsi, D2V/Dlambda^2, D2V/Dpsi^2. Lambda
!           derivatives will be linear combinations of V(P), V(A_1),
!           and V(A_2); psi derivatives will be linear combinations
!           of V(P), V(A_3), and V(A_4); here P=(i,j). We use
!           notation
!           DV/Dlambda =    cl1_a1*V(A_1)+cl1_a2*V(A_2)+cl1_p*V(P)
!           D2V/Dlambda^2 = cl2_a1*V(A_1)+cl2_a2*V(A_2)+cl2_p*V(P)
!           DV/Dpsi =       cp1_a3*V(A_3)+cp1_a4*V(A_4)+cp1_p*V(P)
!           D2V/Dpsi^2 =    cp2_a3*V(A_3)+cp2_a4*V(A_4)+cp2_p*V(P)
!
!
!           Compute the distances to the 4 neighbors:
!
            lam_f = dlam
            lam_b = dlam * alp1
            psi_f = dpsi * bet1
            psi_b = dpsi * gamma
!
            cl1_a1 = + lam_b / lam_f / (lam_f+lam_b)
            cl1_a2 = - lam_f / lam_b / (lam_f+lam_b)
            cl1_p  = + (lam_f-lam_b) / (lam_f*lam_b)
!
            cl2_a1 = + 2. / lam_f / (lam_f+lam_b)
            cl2_a2 = + 2. / lam_b / (lam_f+lam_b)
            cl2_p  = - 2. / (lam_f*lam_b)
!
            cp1_a3 = + psi_b / psi_f / (psi_f+psi_b)
            cp1_a4 = - psi_f / psi_b / (psi_f+psi_b)
            cp1_p  = + (psi_f-psi_b) / (psi_f*psi_b)
!
            cp2_a3 = + 2. / psi_f / (psi_f+psi_b)
            cp2_a4 = + 2. / psi_b / (psi_f+psi_b)
            cp2_p  = - 2. / (psi_f*psi_b)
!
            denom     = g_lam*cl2_p  + g_psi*cp2_p &
                      + f_lam*cl1_p  + f_psi*cp1_p
!
            c (1,i,j) = g_lam*cl2_a1 + f_lam*cl1_a1
            c (2,i,j) = g_lam*cl2_a2 + f_lam*cl1_a2
            c (3,i,j) = g_psi*cp2_a3 + f_psi*cp1_a3
            c (4,i,j) = g_psi*cp2_a4 + f_psi*cp1_a4
!
            c (1,i,j) = - c(1,i,j)/denom
            c (2,i,j) = - c(2,i,j)/denom
            c (3,i,j) = - c(3,i,j)/denom
            c (4,i,j) = - c(4,i,j)/denom
            c (5,i,j) = + c(5,i,j)/denom
!
            IF (.NOT.a2_loc) THEN
               c (5,i,j) = c(5,i,j) + c(2,i,j)*vbnd(j)
               c (2,i,j) = 0.
            END IF
!
            IF (.NOT.a3_loc) THEN
               c (5,i,j) = c(5,i,j) + c(3,i,j)*v_right
               c (3,i,j) = 0.
            END IF
!
            IF (.NOT.a4_loc) THEN
               c (5,i,j) = c(5,i,j) + c(4,i,j)*v_left
               c (4,i,j) = 0.
            END IF
!
         END DO Loop_i
!        don't put anything between these two ENDDOs
      END DO Loop_j
!
!
      DO i = 1, isize
      DO kindex = 1, ncoeff
         c(kindex, i, j1 - 2) = c (kindex, i, j2 - 1)
         c(kindex, i, j1 - 1) = c (kindex, i, j2    )
         c(kindex, i, j2 + 1) = c (kindex, i, j1    )
      END DO
      END DO
!
      RETURN
      END SUBROUTINE New_coeff
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE New_cfive ( c, c5w, birk)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: birk (:,:), c5w (:,:)
      REAL (rprec), INTENT (IN OUT) :: c (:,:,:)
!
      INTEGER (iprec) :: i,j
!
      DO j = 1, jsize
      DO i = imin_j(j), isize
         c (5,i,j) = alpha(i,j)*beta(i,j)*(Ri**2)*birk(i,j) + c5w(i,j)
      END DO
      END DO
!
      RETURN 
      END SUBROUTINE New_cfive                          
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Comput_V_Potnt3 (c, u)
      IMPLICIT NONE
      REAL (rprec) , INTENT (IN)   :: c (:,:,:)
      REAL (rprec), INTENT(IN OUT) :: u (:,:)
!
      INTEGER (iprec), PARAMETER :: maxits = 50000
      REAL (rprec) :: cut,rjac, anorm, resid, omega
      INTEGER  (iprec) :: n,i,j, imax, jmax, min_j(jsize)
!
!
      min_j = imin_j(:)
      imax = SIZE (v, DIM = 1)
      jmax = SIZE (v, DIM = 2)
      cut = 10.0_rprec        ! cut for sum of residuals in volts
      rjac = one - pi**2 /two / (jmax**2) ! what is it for a non-square grid?
!                                                                       
      omega = one 
!                                                                       
!
!
      iterate_loop: DO n = 1, maxits 
!
         anorm = zero
!                                                                       
!        Inner boundary using coefficients   
!
         DO j = 2, jmax - 1 
            i = min_j (j) 
            u (i, j) = c (5, i, j) + c (1, i, j) * u (i + 1, j) &
                                   + c (2, i, j) * u (i - 1, j) &
                                   + c (3, i, j) * u (i, j + 1) &
                                   + c (4, i, j) * u (i, j - 1)
         END DO 
!                                                                       
!        Outer boundary using coefficients             
!        took out c1 because imax+1 is out of bounds  
!
         DO j = 2, jmax - 1 
            i = imax 
            u (i, j) = c (5, i, j) + c (2, i, j) * u (i - 1, j) &
                                   + c (3, i, j) * u (i, j + 1) &
                                   + c (4, i, j) * u (i, j - 1) 
         END DO 
!
!        Use periodicity to get other points:            
!                                                                       
         u (min_j (1), 1) = u (min_j (jmax - 2), jmax - 2) 
         u (min_j (jmax), jmax) = u (min_j (3), 3) 
!                                                                       
         u (imax, 1) = u (imax, jmax - 2) 
         u (imax, jmax) = u (imax, 3) 
!                                                                       
         DO j = 2, jmax-2
         DO i = min_j(j)+1, imax - 1 
!
           u (i, 1)        = u (i,jmax-2)
           u (i, jmax - 1) = u (i, 2)
           u (i,jmax)      = u (i,3)
!
           IF (MOD (i + j, 2) == MOD (n, 2) ) THEN 
              resid = - c (1, i, j) * u (i + 1, j) &
                      - c (2, i, j) * u (i - 1, j) &
                      - c (3, i, j) * u (i, j + 1) &
                      - c (4, i, j) * u (i, j - 1) &
                      + u (i, j) - c (5, i, j) 
!
              anorm = anorm + ABS (resid) 
              u (i, j) = u (i, j) - omega * resid 
!
           END IF 
!
         END DO 
         END DO 
!
         IF (n == 1) THEN 
            omega = one / (one-half * rjac**2) 
         ELSE 
            omega = one / (one-qtr * rjac**2 * omega) 
         END IF 
!                                                                       
         IF (n >  1 .AND. anorm < cut) THEN
             WRITE (*,*) 'in potnt3',n,anorm
             EXIT iterate_loop
         END IF
!                                                                       
      END DO iterate_loop
!
      IF (n >= maxits) STOP 'maxits exceeded'
!
      RETURN
      END SUBROUTINE Comput_V_Potnt3
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!   Potential solver package:
!
!
    SUBROUTINE Gmresm (bndloc, c, v)
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: c (:,:,:), bndloc (:)
    REAL (rprec), INTENT (IN OUT) :: v (:,:)
!
!
!     In this subroutine we will solve difference equations on a part of
!     the RCM grid. The part of the grid where V is found is defined as
!     j1 <= j <= j2,  MINVAL(min_j) <= i <= idim, where j1 and j2 are
!     RCM variables from common block NDICES. On this part of the grid,
!     we will renumber the grid points and re-write the RCM difference
!     equations as A*X=B, where X is the vector of unknowns, A is a matrix
!     and B is a vector both formed from C1-C5 coefficients.
!
!     Since A is sparse, it is stored in the CRS format (compressed-row).
!     Linear system A*X=B is solved iteratively by GMRES(m) algorithm
!     (Generalized Minimized Residuals with restarts).
!
!
!     Local variables:
!       nmax  : max # of grid points to treat
!       M     : max # of Krylov vectors (when exceeded, goes into restart)
!       itermx: max # of restarts allowed
!       ilow  : min value of min_j over all J
!       ni    : I-size (# of pts) of the square part of grid treated here.
!       nj    : J-size of the square part of grid treated
!       nij   : total number of grid points over which method loops.
!       nnz   : the total number of nonzero elements of matrix A
!       nzmax : max # of non-zero coefficients of A.
!       amatrx: non-zero elements of A ordered into a linear sequence
!       colind: vector holding column numbers for each element of AMATRX
!       ROWPTR: vector with locations of 1st element of A in each row in AMATRX
!       DIAPTR: vector with locations of diagonal elements of A in AMATRX
!       PIVOTS: vector with inverses of diagonal of the preconditioner matrix
!       TOL   : relative error. Empirically, 1e-5 is equiv. to
!               sum of resids < 1 Volt
!       H     : Hessenberg matrix holding dot products
!       B     : holds right-hand side vector of linear system to be solved
!       CS, SN: vectors holding Givens rotations parameters
!       X0    : vector holding initial approx. on entry and solution on exit
!       RESID : vector of residuals (b-A*x)
!       AX, W : work vectors
!       Y     : vector holding coefficients of the solution expanded in the
!               Krylov vectors
!       S     : initially unit vector, then rotated by Givens rotations
!
!       This algorithm uses preconditioning based on incomplete LU
!       factorization. Namely, it uses D-ILU type.
!

    INTEGER (iprec), PARAMETER :: M = 300, itermx = 300
    INTEGER (iprec)            :: nzmax, nij, nnz, imin_here
    INTEGER (iprec), ALLOCATABLE   :: row_ptr (:), i_column (:), diag_ptr(:)
    REAL    (rprec), ALLOCATABLE   :: a_mtrx(:), b_mtrx(:), pivots(:), &
                                      x0(:), x(:,:), resid(:), w (:), ax (:)
    REAL (rprec), POINTER :: window (:)
    REAL (rprec), TARGET  :: h (m+1,m), s (m+1)
    REAL (rprec)          :: sn (m), cs (m), y (m)
    INTEGER (iprec) :: nj, ni, krow, ii, jj, i, jkryl, irstrt
!
      REAL (rprec) :: bnorm, rnorm, relerr
!
!     REAL :: error_s, error_1
!
!  1. Arrange difference equations into a matrix form A*x=b. This call returns
!     A as (AMATRX, COLIND, ROWPTR), and also B (RHS) and X0:
!
      CALL Gmresm_define_matrix ( )
!
!
! 2.  Compute preconditioner:
!
      CALL Gmresm_compute_DILU ()
!
!
!
!  3. Now begin GMRES(m) algorithm to solve A*x=b with X0 as init. approx.
!
      bnorm = SQRT ( DOT_PRODUCT (b_mtrx, b_mtrx))
      IF (ABS(bnorm) < machine_tiny) bnorm = one
!
!
!     Begin GMRES(m) iterations (restarts):
!
      Restart_loop: DO irstrt = 1, itermx
!
!
!        ... Compute the norm of initial residuals:
!
         ax     = b_mtrx - Gmresm_Mtrx_times_vect (X0)
         resid  = Gmresm_Msolve (ax )
         rnorm  = SQRT (DOT_PRODUCT (resid, resid ))
         relerr = rnorm / bnorm
         IF (relerr < TOL_gmres ) RETURN     ! V already holds solution
!
!        .. Set 1st Krylov vector to R/||R||:
!
         x (:, 1) = resid / rnorm
!
!        .. Set up unit vector E1 of length RNORM:
!
         s (1) = rnorm
         s (2:m+1) = zero
!
!        .. Loop to generate orthonormal vectors in Krylov subspace:
!
         iterate_loop: DO jkryl = 1, M
!
!            ... Compute A*X(Jkryl) and solve M*w=A*X(kryl) for w:
!
             ax = Gmresm_Mtrx_times_vect (x (:,jkryl) )
             w  = Gmresm_Msolve ( ax )
!
!            ... Form J-th column of H-matrix and X (Jkryl+1)
!                (modified Gramm-Schmidt process):
!
            DO i = 1, jkryl
               H (i,jkryl)   = DOT_PRODUCT ( w , x (:,i) )
               w             = w  - h (i,jkryl) * x (:,i)
            END DO
            h (jkryl+1,jkryl)  = SQRT (DOT_PRODUCT (w, w))
            x (:, jkryl+1)     = w / h(jkryl+1,jkryl)
!
!
!           .. Update QR-factorization of H. For that, 
!           .... first, apply 1, ..., (Jkryl-1)th rotations
!                to the new (Jkryl-th) column of H:
!
            DO i = 1, Jkryl-1
               window => h (i:i+1, jkryl)
               window = Gmresm_Rotate_vector ( window, cs (i), sn (i) )
            END DO
!
!           .... second, compute the Jkryl-th rotation that
!                will zero H (jkryl+1,jkryl):
!
            window => h (jkryl:jkryl+1, jkryl)
            CALL Gmresm_Get_rotation ( window, cs (jkryl), sn (jkryl) )
!
!           .... third, apply Jkryl-th rotation to Jkryl-th column of H
!                and to S (rhs):
!
            window => h (jkryl:jkryl+1, jkryl)
            window = Gmresm_Rotate_vector ( window, cs (jkryl), sn (jkryl) )
            h (jkryl+1,jkryl) = zero
!
            window => s (jkryl : jkryl+1)
            window = Gmresm_Rotate_vector ( window, cs (jkryl), sn (jkryl) )
!
!
!           .. Approximate the norm of current residual:
!
            relerr = ABS (s (jkryl+1)) / bnorm
!
            IF (relerr < TOL_gmres) THEN
!
!              .. Residual is small, compute solution, exit:
!
               y(1:Jkryl) = Gmresm_Solve_upper_triang (A = h, B_RHS = s, N = Jkryl)
               DO i = 1, Jkryl
                  x0  = x0 + y(i)* X(:,i)
               END DO
               EXIT restart_loop
!
            END IF
!
         END DO iterate_loop
!
!
!        We got here because after a maximum number of Krylov vectors
!        was reached, approximated norm of residual was not small enough.
!        However, need to compute approx solution and check the actual norm
!        of residual (because the approx. norm may not be accurate due to
!        round offs):
!
         y(1:m) = Gmresm_Solve_upper_triang (A = h, B_RHS = s, N = m)
         DO i = 1, m
            x0 = x0 + y(i)* X(:,i)
         END DO
! 
         resid = b_mtrx - Gmresm_Mtrx_times_vect ( x0 )
         rnorm   = SQRT (DOT_PRODUCT (resid, resid))
         relerr = rnorm / bnorm
!
!        .. If the actual norm of residual is indeed small, exit:
! 
         IF (relerr < TOL_gmres) EXIT restart_loop
!
!        .. If not, continue by restarting...
!
      END DO restart_loop
!
!
!     Finished GMRES(m) loop. We get here either because
!     the solution was found, or because maximum number of
!     iterations was exceeded. Check for this:
!     
      IF (relerr >= TOL_gmres) THEN
         STOP 'convergence in GMRES(m) not achieved, stopping'
      END IF
!
!
!  4. Solution was found. The final step is to decode solution and put it
!     back into V array.
!
      DO jj = j1, j1+nj-1
         DO ii = imin_here, isize
            krow = ni*(jj-j1) + (ii-imin_here+1)
            IF (ii >= CEILING(bndloc(jj))) v (ii,jj) = X0 (krow)
         END DO
      END DO
!
      CALL Circle ( v )
!
!
! ************* Residual check ********************
!  error_s = 0.
!  DO jj = j1, j2
!     DO ii = min_j(jj), idim -1
!        error_1 = v (ii,jj) - c(1,ii,jj)*v(ii+1,jj) &
!                  -c(2,ii,jj)*v(ii-1,jj)&
!                  -c(3,ii,jj)*v(ii,jj+1)&
!                  -c(4,ii,jj)*v(ii,jj-1) - c(5,ii,jj)
!        error_s = error_s + ABS (error_1)
!     END DO
!     ii = idim
!     error_1 = v(ii,jj) - c(2,ii,jj)*v(ii-1,jj)&
!               -c(3,ii,jj)*v(ii,jj+1)-c(4,ii,jj)*v(ii,jj-1)&
!               -c(5,ii,jj)
!     error_s = error_s + ABS (error_1)
!  END DO
!  WRITE (*, &
!  &'(A11,ES9.2, 2X, A7,ES9.2, 2X, A14,ES9.2, 2X, A5,I3, 2X,&
!  & A2,I3)')  &
!        'SUM(resid)=', error_s, &
!        'RNORM2=', rnorm, &
!        'RNORM2/BNORM2=',rnorm/bnorm, &
!        'ITER=', irstrt,  &
!        'J=', jkryl
!

!
      RETURN
!
!
    CONTAINS  !----------------------------------------------------
!
!
!
      SUBROUTINE Gmresm_Define_matrix ( )
      IMPLICIT NONE
!_____________________________________________________________________________
!     Subroutine returns matrix A stored in 3 vectors AMATRX, COLIND, ROWPTR
!     also NNZ is the number of non-zero elements of A, and DIAPTR vector
!     holds locations of the diagonal elements of A in AMATRX.
!
!     This subroutine will compute:
!     -- nij, size of smallest rect. grid area enclosing the modeling region,
!     -- nzmax, upper limit on the number of non-zero coeffs of matrix A,
!     -- nnz, actual number of non-zero elements of A,
!     -- b_mtrx(1:nij), right-hand side of system of linear equations,
!     -- a_mtrx(1:nnz)--matrix A of linear system encoded in CRS,
!     -- row_ptr(1:nij+1), i_column(1:nnz), diag_ptr (1:nij), encoding of A in CRS,
!
!     -- also allocate pivots (1:nzmax) for the pre-conditioner.
!_____________________________________________________________________________
!
!     Local variables:
!
      INTEGER (iprec):: i, j, L, krow
!
!
!     In this subroutine we take the coefficients of the RCM difference
!     equations approximating the MI-coupling PDE, and reformulate these
!     equations as to cast them into a linear system A*X=B, where A is
!     an NxN square matrix, X is the unknown vector whose elements are
!     unknown values of the potential on grid points V(i,j), and B is the
!     right-hand-side vector.
!     Apparently, such reformulation requires: (1) to number all grid points
!     in the modeling region sequentially into a 1-dimensional sequence, and
!     then (2) to form A from c1-c4 and B from c5 RCM coefficients.
!     As A is going to be sparse, an additional task is to store (encode) A
!     in the Compressed-Row-Storage format for using in the potential solver.
!
!     Matrix A is stored in one-dim REAL array AMATRX and two INTEGER 1-dim
!     arrays COLIND and ROWPTR. We simply go along each row of A starting with
!     the 1st row, then 2nd, etc, and for each non-zero element a(p,q), we
!     write AMATRX(L)=a(p,q), COLIND(L)=q, and L-index numbers those non-zero
!     elements sequentially. ROWPTR(p) has the L-index of where p-th row
!     starts in AMATRX.
!
!  1. Numbering grid points into a 1-dim. sequence.
!     Imagine the RCM grid as extending vertically in I from I=1 (highest lat,
!     top) to I=IDIM (lowest lat., bottom) and horizontally in J from J=j1
!     (noon, left) to J=J2 (last point before noon, right). If MIN_J(j) gives
!     the first I-point inside the modeling region, then we will consider all
!     grid points (i,j) such that ilow <= i <= idim, j1 <= j <= j2, where
!     ilow is MINVAL(min_j(:)). This will result in inclusion of some points
!     that are outside the modeling region, but we will define the difference
!     equations for them such that they don't matter.
!
!     The rectangular region of the grid we treat has the size NI by NJ, with
!     total of NIJ points in it:
!
      imin_here = MINVAL (CEILING(bndloc))
      nj  = j2 - j1 + 1
      ni  = isize - imin_here + 1
      nij = ni*nj
      nzmax = nij * ncoeff
      CALL Gmresm_Make_storage ( )
!
!
!     Number all points and store only non-zero elements of A. Order grid
!     points along J-lines from ILOW to IDIM, occasionally including points
!     outside the modeling region. Each point (i,j)
!     has number KROW (so that coefficients of the difference equation on that
!     point are on the krow-th row of A).
!
      L = 0
      DO j = j1, j2
      DO i = imin_here, isize
!
         krow          = ni * (j - j1) + (i-imin_here + 1)
         b_mtrx (krow) = c (5,i,j)   ! this will be RHS vector
         X0 (krow)     = v (i,j) ! initial approximation taken from prev solution
         row_ptr(krow) = L + 1
!
         IF (i < CEILING(bndloc(j))) THEN
!
!           we are outside the main modeling region. In this case, the value
!           of the potential at this point is irrelevant, but we need to define
!           coefficients so as to proceed with this point as efficiently as possible.
!           Therefore, make V(i,j) a solution; that is, add a difference equation
!           V (i,j) = V (i,j)_on_input; the row of matrix A is then:
!
!           . . . . . .   1 . . . . . . . . . .
!           and RHS is just V(i,j)
!
            L               = L + 1
            a_mtrx (L)      = one
            i_column (L)    = krow
            b_mtrx (krow)   = v (i,j)
            diag_ptr (krow) = L
!
         ELSE IF (j == j1) THEN   ! first NI rows of A:
!
!           we are on the left (J=J1) side of the grid, and j-1 neighbors
!           wrap around from the other side (j=j2)--periodic boundary
!           condition in J. So for these first NI rows of the matrix A,
!           c4 is the last coefficient. The first NI matrix rows look like:
!
!            1 -c1   0 .............-c3........ . . .      -c4...........0
!          -c2  1  -c1 ............. 0  -c3.... . . .       0  -c4.......0
!            0  -c2  1 -c1.................-c3. . . .       ......-c4....0
!           . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!                   -c2 1 ....................... -c3  . .  ..........  -c4
!
            IF (i /= imin_here) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (2,i,j)
               i_column (L) =  krow - 1
            END IF
!
            L               = L + 1
            a_mtrx (L)      = one
            i_column (L)    = krow
            diag_ptr (krow) = L
!
            IF (i /= isize) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (1,i,j)
               i_column (L) =  krow + 1
            END IF
!
            L            =  L + 1
            a_mtrx (L)   = -c (3,i,j)
            i_column (L) =  krow + NI
!
            L            =  L + 1
            a_mtrx (L)   = -c (4,i,j)
            i_column (L) =  nij - NI + krow
!
         ELSE IF (j == j2) THEN  ! last NI rows of A:
!
!           we are on the right (J=J2) side of the grid, periodic boundary
!           condition means c3-neighbors will come from J=J1, so that c3
!           coefficients will be first, c4-last on each row of matrix A. The
!           last NI rows of matrix A look like:
!
!          -c3  0 ..... . . .     -c4............  1  -c1 ................
!           0  -c3..... . . .      0 -c4......... -c2  1  -c1 ............
!           ......-c3.. . . .      .....-c4......  0  -c2  1  -c1 ........
!           . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
!           ...............-c3. . . ............-c4 .................-c2  1
!
            L            =  L + 1
            a_mtrx (L)   = -c (3,i,j)
            i_column (L) =  i - imin_here + 1
!
            L            =  L + 1
            a_mtrx (L)   = -c (4,i,j)
            I_column (L) =  krow - NI
!
            IF (i /= imin_here) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (2,i,j)
               i_column (L) =  krow - 1
            END IF
!
            L               = L + 1
            a_mtrx (L)      = one
            i_column (L)    = krow
            diag_ptr (krow) = L
!
            IF (i /= isize) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (1,i,j)
               i_column (L) =  krow + 1
            END IF
!
         ELSE          ! the rest of A:
!
!           we are neither on the left (j=j1) nor on the right (j=j2) side
!           of the grid, so don't need to worry about periodic boundary
!           conditions in J. RCM difference equation looks like:
!           V(i,j)=c(1,i,j)*V(i+1,j)+c(2,i,j)*V(i-1,j)+c(3,i,j)*V(i,j+1)+
!                  c(4,i,j)*V(i,j-1)+c(5,i,j);
!           except that if i=idim, then there is no C(1) term, and
!           if i=min_j(j), there is no c(2) term.
!           For such a grid point, the row of A looks like:
!
!           . . . . . -c4 . . . . . . -c2  1  -c1 . . . . . . -c3 . . . .
!
!           except that there is no c2 if i=ilow and no c1 if i=idim.
!
            L            =  L + 1
            a_mtrx (L)   = -c (4,i,j)
            i_column (L) =  krow - ni
!
            IF (i /= imin_here) THEN
               L           =  L + 1
               a_mtrx(L)   = -c (2,i,j)
               i_column(L) =  krow - 1
            END IF
!
            L               = L + 1
            a_mtrx (L)      = one
            i_column (L)    = krow
            diag_ptr (krow) = L
!
            IF (i /= isize) THEN
               L            =  L + 1
               a_mtrx (L)   = -c (1,i,j)
               i_column (L) =  krow + 1
            END IF
!
            L            =  L + 1
            a_mtrx (L)   = -c (3,i,j)
            i_column (L) =  krow + ni
!
         END IF
      END DO
      END DO
!
      nnz = L ! number of non-zero elements of A
      row_ptr (nij+1) = nnz + 1    ! by definition of CRS format.
!
      RETURN
      END SUBROUTINE Gmresm_Define_matrix
!
!
!*****************************************************************************
!
!
      SUBROUTINE Gmresm_Make_storage ( )
      IMPLICIT NONE
      IF (ALLOCATED (a_mtrx)) DEALLOCATE (a_mtrx)
      IF (ALLOCATED (b_mtrx)) DEALLOCATE (b_mtrx)
      IF (ALLOCATED (i_column)) DEALLOCATE (i_column)
      IF (ALLOCATED (diag_ptr)) DEALLOCATE (diag_ptr)
      IF (ALLOCATED (row_ptr)) DEALLOCATE (row_ptr)
      IF (ALLOCATED (pivots)) DEALLOCATE (pivots)
      IF (ALLOCATED (ax)) DEALLOCATE (ax)
      IF (ALLOCATED (w)) DEALLOCATE (w)
      IF (ALLOCATED (resid)) DEALLOCATE (resid)
      IF (ALLOCATED (x0)) DEALLOCATE (x0)
      IF (ALLOCATED (x)) DEALLOCATE (x)
      ALLOCATE ( a_mtrx (nzmax), b_mtrx (nij), pivots (nij), &
                 i_column (nzmax), row_ptr (nij+1), diag_ptr (nij), &
                 ax (nij), w (nij), resid (nij), x0(nij), x (nij,m+1) )
      RETURN
      END SUBROUTINE Gmresm_Make_storage
!
!
!*****************************************************************************
!
!
      SUBROUTINE Gmresm_Compute_DILU ()
      IMPLICIT NONE
!_____________________________________________________________________________
!     Compute the preconditioner M. Matrix A is split as A = L_a + D_a + U_a
!     (strictly-lower triangular, diagonal and strictly-upper triangular).
!     Then M = L * U = (D + L_a) * D^(-1) * (D + U_a), so only need to find
!     and store D (one diagonal). In fact, PIVOTS holds inverses of D since
!     will divide by them later. D-ILU preconditioner M is kept in PIVOTS.
!     All structures are accessed from the host subroutine.
!
!     This subroutine only modifies (computes) PIVOTS (1:nij).
!_____________________________________________________________________________
!
      INTEGER (iprec) :: irow, jcol, krow
      REAL    (rprec) :: element
      LOGICAL         :: found
!
      pivots  = a_mtrx (diag_ptr )
      DO irow = 1, nij
         pivots (irow) = one / pivots (irow)
         DO jcol = diag_ptr(irow)+1, row_ptr(irow+1)-1
            found = .FALSE.
            DO krow = row_ptr (i_column (jcol)), diag_ptr (i_column (jcol)) - 1
               IF (i_column (krow) == irow) THEN
                  found = .TRUE.
                  element = a_mtrx (krow)
               END IF
            END DO
            IF (found) THEN
               pivots (i_column (jcol)) = &
                  pivots (i_column (jcol)) - element*pivots(irow)*a_mtrx(jcol)
            END IF
         END DO
      END DO
      RETURN
      END SUBROUTINE Gmresm_Compute_DILU
!
!
!*****************************************************************************
!
!
      FUNCTION Gmresm_Mtrx_times_vect (x)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN)         :: x (:)
      REAL (rprec), DIMENSION (SIZE(x)) :: Gmresm_Mtrx_times_vect
!____________________________________________________________________________
!     subroutine to form matrix-vector product. Matrix A of size NxN
!     is assumed to be sparse with NNZ non-zero elements and is stored
!     in the compressed-row (CRS) format.
!     We compute y = A*x, where y and x are both vectors of length NNZ.
!     On entry, X holds x, and on exit Y is the result.
!____________________________________________________________________________
!
      INTEGER (iprec) :: i, j
!
      DO i = 1, SIZE (x)
         Gmresm_Mtrx_times_vect (i) = zero
         DO j = row_ptr (i), row_ptr (i+1)-1
            Gmresm_Mtrx_times_vect (i) = &
               Gmresm_Mtrx_times_vect (i) + a_mtrx(j) * x (i_column(j) )
         END DO
      END DO
      RETURN
      END FUNCTION Gmresm_Mtrx_times_vect
!
!
!*****************************************************************************
!
!
      SUBROUTINE Gmresm_Get_rotation ( vector_in, cos_theta, sin_theta )
      IMPLICIT NONE
      REAL (rprec), INTENT (IN)  :: vector_in (2)
      REAL (rprec), INTENT (OUT) :: cos_theta, sin_theta
!_____________________________________________________________________________
!     Compute a Givens (plane) rotation that will act on 2 elements of a
!     vector, A and B, and will zero B. Returns cosine and sine of THETA, the
!     angle of rotation. The transformation is
!
!        A_prime = A * COS(theta) - B * SIN(theta)
!        B_prime = A * CIN(theta) - B * COS(theta)
!
!     In matrix-vector terms,
!        X = (...... A ..... B .....)^T,
!        T =
!        X_prime = T * X = (..... A_prime ...... 0 .....)^T,
!        only 2 elements of X are changed by the rotation.
!_____________________________________________________________________________
!
      REAL (rprec) :: temp
!
      IF ( ABS(vector_in (2)) < machine_tiny ) THEN
         cos_theta = one
         sin_theta = zero
      ELSE IF ( ABS ( vector_in(2) ) > ABS ( vector_in (1) ) ) THEN
         temp = -vector_in (1) / vector_in (2)
         sin_theta = one / SQRT( one + temp**2 )
         cos_theta = temp * sin_theta
      ELSE
         temp = -vector_in (2) / vector_in (1)
         cos_theta = one / SQRT( one + temp**2 )
         sin_theta = temp * cos_theta
      END IF
!
      RETURN
      END SUBROUTINE Gmresm_Get_rotation
!
!
!*****************************************************************************
!
!
      FUNCTION Gmresm_Rotate_vector ( vec_in, cos_theta, sin_theta )
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: vec_in (2), cos_theta, sin_theta
      REAL (rprec)              :: Gmresm_Rotate_vector (2)
!___________________________________________________________________________
!     Apply a plane (Givens) rotation with cos_theta, sin_theta) to a vector.
!     Rotation acts on only two elements of the vector, X and Y.
!___________________________________________________________________________
!
      REAL (rprec):: temp
!
      temp                = cos_theta * vec_in (1) - sin_theta * vec_in (2)
      Gmresm_Rotate_vector(2)    = sin_theta * vec_in (1) + cos_theta * vec_in (2)
      Gmresm_Rotate_vector(1)    = temp
!
      RETURN
      END FUNCTION Gmresm_Rotate_vector
!
!
!*****************************************************************************
!
!
      FUNCTION Gmresm_Solve_upper_triang (a, b_rhs, n)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: n
      REAL (rprec), INTENT (IN) :: a (:,:), b_rhs (:)
      REAL (rprec) :: Gmresm_Solve_upper_triang (n)
!
!     Given an upper triangular matrix A and right-hand side vector B(n),
!     solves linear system A*x=b. A is a Hessenberg matrix (nmax+1 by nmax),
!     but only n by n section is used.
!
      INTEGER (iprec) :: j
!
      IF (UBOUND(a,DIM=1) /= UBOUND (a,DIM=2) + 1) STOP 'PROBLEM 1 IN SOLVETR'
      IF (n > UBOUND (a,DIM = 2) .OR. n < 1) STOP 'PROBLEM 2 IN SOLVETR'
!
      Gmresm_Solve_upper_triang = b_rhs(1:n)
      DO j = N, 1, -1
         Gmresm_Solve_upper_triang (j)     = Gmresm_Solve_upper_triang (j) / a (j,j)
         Gmresm_Solve_upper_triang (1:j-1) = Gmresm_Solve_upper_triang (1:j-1) - &
                                      Gmresm_Solve_upper_triang (j) * a(1:j-1,j)
      END DO
      RETURN
      END FUNCTION Gmresm_Solve_upper_triang
!
!
!*****************************************************************************
!
!
      FUNCTION Gmresm_Msolve (x)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: x (:)
      REAL (rprec), DIMENSION (SIZE(x)) :: Gmresm_Msolve
!_____________________________________________________________________________
!     This subroutine solves the system L *  U * y = x, where
!     M = L * U = (D + L_a) * (I + D^(-1) * U_a) is the D-ILU preconditioner.
!     Matrices L_a and U_a are strictly lower and strictly upper triangular,
!     so that A = D_a + L_a + U_a, and D comes from incomplete LU factorization
!     when computing preconditioner. Solution proceeds in the regular way by
!     forward- and then back-substition (solving L*z=x, then U*y=z).
!     A is in the compressed-row-storage format (A_MTRX, I_COLUMN, ROW_PTR).
!     Diagonal matrix D (in fact, its inverse) is stored in the PIVOTS:
!     PIVOTS(i)=1/D(i,i), and DIAPTR vector holds locations of d_i_i in amatrx.
!     Since A_MTRX and PIVOTS do not change in the potential solver once is
!     has been called, we access them by host association from GMRESM. Only
!     vector X changes from invocation to invocaton of MSOLVE, and we pass it
!     as an argument.
!**** NOTE: book by Barrett et al ("templates ...") has an error in the back-
!           substitution algorithm (p.73). Here I do it correctly.
!_____________________________________________________________________________
!
      INTEGER (iprec) :: i, j
      REAL (rprec)    :: tmp
!
      DO i = 1, SIZE (x)
         tmp = zero
         DO j = row_ptr (i), diag_ptr (i) - 1
            tmp = tmp + a_mtrx (j) * Gmresm_msolve (i_column (j))
         END DO
         Gmresm_Msolve (i) = pivots (i) * (x(i) - tmp)
      END DO
      DO i = SIZE(x), 1, -1
         tmp = zero
         DO j = diag_ptr (i) + 1, row_ptr (i+1)-1
            tmp = tmp + a_mtrx(j) * Gmresm_Msolve (i_column(j))
         END DO
         Gmresm_Msolve (i) = Gmresm_Msolve (i) - pivots (i) * tmp
      END DO
      RETURN
      END FUNCTION Gmresm_Msolve
!
!
    END SUBROUTINE Gmresm
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
       CALL Move_plasma_grid_new (dt)
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
          IF (biold > Bndy(bndloc,bjold)) THEN 
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
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
      SUBROUTINE Read_bfield ()
      IMPLICIT NONE
!
!_____________________________________________________________________________
!
!     Subroutine to read all of the labels of the input
!     bfield models to determine the event times for each
!     of the records.  This information is placed in the
!     1-d array ibtime. These mark times are set in the
!     program creating the B-field arrays.
!     rws    3/19/97
!_____________________________________________________________________________
!
      INTEGER (iprec) :: n, nbf
      LOGICAL :: error_flag
      LOGICAL, SAVE :: called_already = .FALSE.
!
      IF (called_already) RETURN
!
IF (itype_bf == 2 .OR. itype_bf == 3) THEN
  nbf = 1
  ALLOCATE (ibtime (nbf))
  ibtime = -999
ELSE IF(itype_bf == 1) THEN 
      n = 1
      DO
         CALL Read_array ('rcmxmin_inp', n, label, ARRAY_2D = xmin, &
                           ERROR_FLAG = error_flag, ASCI = asci_flag)
         IF (error_flag) EXIT
         n = n + 1
      END DO
      nbf = n - 1

      ALLOCATE (ibtime (nbf))

      DO n = 1, nbf
        CALL Read_array ('rcmxmin_inp', n, label, ARRAY_2D = xmin, ASCI=asci_flag)
        ibtime (n) = label%intg (6)
      END DO 
ELSE
   STOP 'ILLEGAL ITYPE_BF IN READ_BFIELD'
END IF
called_already = .TRUE.
      RETURN
      END SUBROUTINE Read_bfield
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
      SUBROUTINE Read_eta_on_bndy ()
      IMPLICIT NONE
!_____________________________________________________________________________
!
!_____________________________________________________________________________
!
      INTEGER (iprec) :: n, n_t
      LOGICAL :: error_flag
      LOGICAL, SAVE :: called_already = .FALSE.
!
      IF (called_already) RETURN
IF (i_eta_bc == 1) THEN
      n = 1
      DO
         CALL Read_array ('rcmetac_inp', n, label, ARRAY_1D = etac, &
                           ERROR_FLAG = error_flag, ASCI = asci_flag)
         IF (error_flag) EXIT
         n = n + 1
      END DO
      n_t = n - 1

      ALLOCATE (itime_etac (n_t), etac_inp(kcsize,n_t))

      DO n = 1, n_t
        CALL Read_array ('rcmetac_inp', n, label, ARRAY_1D = etac, ASCI=asci_flag)
        itime_etac (n) = label%intg (6)
        etac_inp (:,n) = etac
      END DO 
ELSE IF (i_eta_bc == 2) THEN
      n_t = 1
      ALLOCATE (itime_etac (n_t), etac_inp(kcsize,n_t))
      itime_etac = 0
      etac_inp = 0.0
ELSE 
      STOP 'ILLEGAL VALUE OF I_eta_bd ON INPUT'
END IF
      called_already = .TRUE.
      RETURN
      END SUBROUTINE Read_eta_on_bndy
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
           IF (REAL(i,rprec) < Bndy(bndloc, REAL(j,rprec)) ) CYCLE
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
      SUBROUTINE Read_kp ()
      IMPLICIT NONE
!_____________________________________________________________________________
!
!     Subroutine to read Kp index values
!     and place them in kpinput array.  Times are stored in
!     ikptime array.  nkpmax is actual number of Kp
!     values.  These results are used in subroutine getkp to
!     interpolate in time to get Kp at any time of interest.
!     rws  03-20-97 stanislav 05-28-99
!_____________________________________________________________________________
!
      INTEGER (iprec) :: nkp, nkpmax
      LOGICAL :: logical_flag
      LOGICAL, SAVE :: called_already = .FALSE.
!                                                                       
      IF (called_already) RETURN
!
      IF (IsCoupledExternally) then
         nkpmax = 1
         ALLOCATE (ikptime(nkpmax), kpinput(nkpmax))
         ikptime(1) = 9999999
         kpinput(1) = 3.0
         called_already = .TRUE.
         RETURN
      END IF
         
      INQUIRE (FILE = rcmdir//'rcmkp_inp', EXIST = logical_flag)
      IF (.NOT.logical_flag ) STOP 'READKP: RCMKP_INP not found'
      INQUIRE (UNIT = LUN, OPENED = logical_flag)
      IF (logical_flag) STOP 'READKP: LUN is already open'
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FILE = rcmdir//'rcmkp_inp')
      nkpmax = 0
      DO
         READ (LUN, *, END = 19)
         nkpmax = nkpmax + 1
      END DO
   19 CLOSE (LUN)
!
      ALLOCATE (ikptime (nkpmax), kpinput (nkpmax) )
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FILE = rcmdir//'rcmkp_inp')
      DO nkp = 1, nkpmax
         READ (LUN, *) ikptime (nkp), kpinput (nkp)
      END DO
      CLOSE (UNIT = LUN)
      called_already = .TRUE.
!
      RETURN 
      END SUBROUTINE Read_kp
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Read_winds (iwind)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: iwind
!
!_____________________________________________________________________________
!                                                                       
!     last update 08-27-86        written by g.a.mantjoukis
!                                                                       
!     subroutine to compute pedersen and hall winds on
!     the rcm-specified grid
!                                                                       
!     iwind=0    no wind                                                
!     iwind=1    tarpley-type wind (calculatoin)                        
!     iwind=2    roble sq winds                                         
!                                                                       
!-------------------------------------------------------------
!                                                                       
      INTEGER (iprec) :: i, j, is, n
      CHARACTER (LEN=80) :: form_string
      REAL    (rprec) :: s, sm, vnorm, sv, su, phi, v, u, ath, c8, c6, c4, c2, c1,&
                         v8, v6, v4, v2, v0, u6, u4, u2, u0
      LOGICAL, SAVE :: called_already = .FALSE.
!
      IF (called_already) RETURN
!
      IF (iwind == 0) THEN
!        sw   = zero
!        pwe  = zero
!        hwn  = zero
!        hwe  = zero
!        pwn  = zero
         OPEN (LUN, FILE = rcmdir//'rcmwind', STATUS = 'OLD', FORM = 'FORMATTED')
!
            READ (LUN,'(I10.10)') n
            IF (n /= isize*jsize) STOP 'size mismatch in read_winds'
            READ (LUN,'(A80)') form_string
            DO j = 1, jsize
            DO i = 1, isize
               READ (LUN, form_string) pwe(i,j), hwn(i,j), hwe(i,j), pwn(i,j)
            END DO
            END DO
!
            READ (LUN,'(I10.10)') n
            IF (n /= jsize) STOP 'size mismatch in read_winds'
            READ (LUN,'(A80)') form_string
            DO j = 1, jsize
               READ (LUN, form_string) sw(j)
            END DO
!
         CLOSE (LUN)
!
IF (ANY(pwn /= 0.0) .or. ANY(pwe /= 0.0) .or. &
    ANY(hwn /= 0.0) .or. ANY(hwe /= 0.0)) STOP 'WINDS NON-ZERO'
!
      ELSE IF (iwind == 1) THEN
!
          pwn = zero
          pwe = zero
!
          u0 =   0.4225615_rprec
          u2 = - 1.3858640_rprec
          u4 = - 1.1390120_rprec
          u6 = - 0.4121196_rprec
!                                                                       
          v0 = - 0.4549289_rprec
          v2 = + 3.0388490_rprec
          v4 = - 3.6561400_rprec
          v6 = - 4.5478570_rprec
          v8 = - 1.9232250_rprec
    !
          DO j = 1, jsize
          DO i = 1, isize
             IF (colat (i, j) > 1.029744 .AND. &
                 colat (i, j) < 1.064651         ) CYCLE
             c1 = COS (colat (i, j) )
             c2 = c1**2
             c4 = c1**4
             c6 = c1**6
             c8 = c1**8
             ath = 130.0 / (qtr - c2)
    !
             u = ath * three * c1 * (u0 + u2 * c2 + u4 * c4 + u6 * c6)
             v = ath * (v0 + v2 * c2 + v4 * c4 + v6 * c6 + v8 * c8)
    !
             phi = aloct (i, j) + pi
             su = sin (phi + 250.0 * pi / 180.0)
             sv = sin (phi + 340.0 * pi / 180.0)
    !
    !        there is a minus sign in unorm since sin(phi)=-sin(phi+pi)
             vnorm = - 13.97155_rprec
             pwn (i, j) = - u * su / vnorm
             pwe (i, j) =   v * sv / vnorm
          END DO
          END DO
    !
    !     fixing singularity at colat=60.if for any j,more than
    !     one i's have 59<colat(i,j)<61 the fixing is no good.
    !     in latter case the program will stop.
    !
          DO_60: DO j = 1, jsize
             is = 0
             DO_50: DO i = 1, isize
             IF (pwn (i, j) /= zero) CYCLE DO_50
             is = is + 1
             IF (is /= 1) STOP 'singularity in winds needs fixing'
             sm = colat (i + 1, j) - colat (i - 1, j)
             s = colat (i, j) - (colat (i + 1, j) + colat (i - 1, j) )*half
             pwn (i, j) = (half + s / sm) * pwn (i + 1, j) &
                        + (half - s / sm) * pwn (i - 1, j)
             pwe (i, j) = (half + s / sm) * pwe (i + 1, j) &
                        + (half - s / sm) * pwe (i - 1, j)
             END DO DO_50
          END DO DO_60
!
          hwn = pwn
          hwe = pwe
!
      END IF
      called_already = .TRUE.
      RETURN
      END SUBROUTINE Read_winds
!
!
!
      SUBROUTINE Write_winds ()
      IMPLICIT NONE
      INTEGER (iprec) :: i,j
      CHARACTER (LEN=80) :: form_string
!
      OPEN (LUN, FILE = rcmdir//'rcmwind', STATUS = 'REPLACE', FORM = 'FORMATTED')
!
      form_string = '(4(TR2,ES23.15))'
      WRITE (LUN, '(I10.10)') SIZE(pwe)
      WRITE (LUN,'(A80)') form_string
      DO j = 1, jsize
      DO i = 1, isize
         WRITE (LUN,form_string) pwe(i,j), hwn(i,j), hwe(i,j), pwn(i,j)
      END DO
      END DO
!
      form_string = '(1(TR2,ES23.15))'
      WRITE (LUN, '(I10.10)') SIZE(sw)
      WRITE (LUN,'(A80)') form_string
      DO j = 1, jsize
         WRITE (LUN,form_string) sw(j)
      END DO
!
      CLOSE (LUN)
!
       RETURN
       END SUBROUTINE Write_winds
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_real_1d_array (filename, rec_num, label, &
                                   array_1d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    REAL (rprec),      INTENT (OUT) :: array_1d (:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: length, istat
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_1d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_1d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_1d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_REAL_1D_ARRAY'
       WRITE (*,'(T2,A)') 'FILE NAME: '//filename
       STOP
    END IF
!
    OPEN (UNIT = LUN, FILE = filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, &
             ERR = 1, FMT = form_string)        label, array_1d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1                   ) label, array_1d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename,' rec=',rec_num
       STOP
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_real_1d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_real_2d_array (filename, rec_num, label, &
                                   array_2d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    REAL (rprec),      INTENT (OUT) :: array_2d (:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    INTEGER (iprec)    :: length, istat
    CHARACTER (LEN=80) :: form_string
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_2d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_2d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_2d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_REAL_2D_ARRAY'
       STOP
    END IF
!
    OPEN (UNIT = LUN, FILE = filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1, FMT = form_string) label, array_2d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1                   ) label, array_2d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename,istat,rec_num
       STOP
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_real_2d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_real_3d_array (filename, rec_num, label, array_3d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    REAL (rprec),      INTENT (OUT) :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: k, rec_num_mod
    LOGICAL            :: flag, asci_format
!
!
    DO k = 1, SIZE (array_3d, DIM = 3)
       rec_num_mod = (rec_num - 1) * SIZE(array_3d, DIM = 3) + k
       CALL Read_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                        error_flag, asci)
    END DO
!
    RETURN
    END SUBROUTINE Read_real_3d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_intg_1d_array (filename, rec_num, label, &
                                   array_1d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    INTEGER (iprec),   INTENT (OUT) :: array_1d (:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: length, istat
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_1d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_1d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,I10))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_1d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_INTG_1D_ARRAY'
       STOP
    END IF
!
    OPEN (UNIT = LUN, FILE = filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1, FMT = form_string) label, array_1d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1                   ) label, array_1d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename
       STOP
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_intg_1d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_intg_2d_array (filename, rec_num, label, &
                                   array_2d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    INTEGER (iprec),   INTENT (OUT) :: array_2d (:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: length, istat
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_2d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_2d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,I10))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_2d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_INTG_2D_ARRAY'
       STOP
    END IF
!
    OPEN (UNIT = LUN, FILE = filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1, FMT = form_string) label, array_2d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1                   ) label, array_2d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename
       STOP
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_intg_2d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Read_intg_3d_array (filename, rec_num, label, &
                                   array_3d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    INTEGER (iprec),   INTENT (OUT) :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    INTEGER (iprec)    :: k, rec_num_mod
!
!
    DO k = 1, SIZE (array_3d, DIM = 3)
       rec_num_mod = (rec_num - 1) * SIZE(array_3d, DIM = 3) + k
       CALL Read_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                        error_flag, asci)
    END DO
!
    RETURN
    END SUBROUTINE Read_intg_3d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_real_1d_array (filename, rec_num, label, array_1d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    REAL (rprec),      INTENT (IN)      :: array_1d (:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    CHARACTER (LEN=11):: form_type_char
    CHARACTER (LEN=80) :: form_string
    LOGICAL           :: flag, asci_format
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       STOP 'UNIT ALREADY OPEN'
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_1d
!
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_1d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxxxx(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_1d)
       form_type_char = 'FORMATTED  '
    END IF

    OPEN (UNIT=LUN, FILE = filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = form_type_char,&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF

    IF (asci_format) THEN
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat, FMT = form_string) label, array_1d
    ELSE
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat                   ) label, array_1d
    END IF
    IF (istat /= 0) STOP 'ERROR READING FILE'

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) STOP 'ERROR CLOSING FILE'

    RETURN
    END SUBROUTINE Write_real_1d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_real_2d_array (filename, rec_num, label, array_2d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    REAL (rprec),      INTENT (IN)      :: array_2d (:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    CHARACTER (LEN=11):: form_type_char
    CHARACTER (LEN=80) :: form_string
    LOGICAL           :: flag, asci_format
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       STOP 'UNIT ALREADY OPEN'
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_2d
!
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_2d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxxxx(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_2d)
       form_type_char = 'FORMATTED  '
    END IF

    OPEN (UNIT=LUN, FILE = filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = form_type_char,&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF

    IF (asci_format) THEN
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat, FMT = form_string) label, array_2d
    ELSE
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat                   ) label, array_2d
    END IF
    IF (istat /= 0) STOP 'ERROR READING FILE'

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) STOP 'ERROR CLOSING FILE'

    RETURN
    END SUBROUTINE Write_real_2d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_real_3d_array (filename, rec_num, label, array_3d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    REAL (rprec),      INTENT (IN)      :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: k, rec_num_mod
    CHARACTER (LEN=7) :: status_char
!
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
!
    DO k = 1, SIZE (array_3d, DIM = 3)
       rec_num_mod = (rec_num - 1) * SIZE (array_3d, DIM = 3) + k
       IF (k == 1) THEN
          CALL Write_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                            SETUP = setup, ASCI = asci)
       ELSE
          CALL Write_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                            ASCI = asci)
       END IF
    END DO
!
    RETURN
    END SUBROUTINE Write_real_3d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_intg_1d_array (filename, rec_num, label, array_1d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    INTEGER (iprec),   INTENT (IN)      :: array_1d (:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    CHARACTER (LEN=11):: form_type_char
    CHARACTER (LEN=80) :: form_string
    LOGICAL           :: flag, asci_format
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       STOP 'UNIT ALREADY OPEN'
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_1d
!
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_1d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxxxx(TR2,I10))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_1d)
       form_type_char = 'FORMATTED  '
    END IF

    OPEN (UNIT=LUN, FILE = filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = form_type_char,&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF

    IF (asci_format) THEN
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat, FMT = form_string) label, array_1d
    ELSE
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat                   ) label, array_1d
    END IF
    IF (istat /= 0) STOP 'ERROR READING FILE'

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) STOP 'ERROR CLOSING FILE'

    RETURN
    END SUBROUTINE Write_intg_1d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_intg_2d_array (filename, rec_num, label, array_2d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    INTEGER (iprec),   INTENT (IN)      :: array_2d (:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    CHARACTER (LEN=11):: form_type_char
    CHARACTER (LEN=80) :: form_string
    LOGICAL           :: flag, asci_format
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       STOP 'UNIT ALREADY OPEN'
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
    INQUIRE (IOLENGTH = length ) label%intg, label%real, label%char, array_2d
!
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_2d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxxxx(TR2,I10))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_2d)
       form_type_char = 'FORMATTED  '
    END IF

    OPEN (UNIT=LUN, FILE = filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = form_type_char,&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF

    IF (asci_format) THEN
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat, FMT = form_string) label, array_2d
    ELSE
       WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat                   ) label, array_2d
    END IF
    IF (istat /= 0) STOP 'ERROR READING FILE'

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) STOP 'ERROR CLOSING FILE'

    RETURN
    END SUBROUTINE Write_intg_2d_array
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Write_intg_3d_array (filename, rec_num, label, array_3d, setup, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    INTEGER (iprec),   INTENT (IN)      :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup, asci
!
    INTEGER (iprec)   :: k, rec_num_mod
    CHARACTER (LEN=7) :: status_char
!
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
!
    DO k = 1, SIZE (array_3d, DIM = 3)
       rec_num_mod = (rec_num - 1) * SIZE (array_3d, DIM = 3) + k
       IF (k == 1) THEN
          CALL Write_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                            SETUP = setup, ASCI = asci)
       ELSE
          CALL Write_array (filename, rec_num_mod, label, array_3d (:,:,k), &
                            ASCI = asci)
       END IF
    END DO
!
    RETURN
    END SUBROUTINE Write_intg_3d_array
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Efield (ipatt, pcp, deqdt, a, b, dx, dy, colat, aloct, &
                         mode, icntrl, itop, v, vbndry)
!
!     this subroutine returns values of the potential v at all grid
!     points, specified by colat(i,j) and aloct(i,j).
!
!
!    general call parameters for subroutine efield
!
!       ipatt = heppner-maynard pattern no.  must be specified as an
!               integer between 1 and 7.
!               ipatt=1 means pattern a, for bz south, by=0.
!               ipatt=2 means pattern bc, for bz south, bynorth>0
!               ipatt=3 means pattern de, for bz south,bynorth<0
!               ipatt=4 means pattern bcp,twisted bc, for bz>0
!               ipatt=5 means pattern bcpp,twisted bc,bz strong>0
!               ipatt=6 means patt. dep, twisted de,bz weak >0
!               ipatt=7 means patt. depp, twisted de,bz strong>0
!
!       pcp = polar-cap potential drop in volts.  must be specified
!             unless using unscaled heppner-maynard potential(icntrl=-1)
!
!       deqdt = estimated d/dt of latitude of equatorward edge of
!               auroral zone at local midnight, in degrees/hour.  this
!               is a parameter in low-latitude e-field model.  not
!               needed if icntrl.le.0.
!
!      vectors a(3), b(3), dx(3), dy(3) describe the ellipses that form
!      the boundaries between regions 1, 2, and 3.  boundary 1 is the
!      equatorward edge of the polar cap.  boundary 2 is the equatorward
!      edge of region 1 or the equatorward edge of the field-reversal
!      region.  boundary 3 is the shielding layer.
!       a(l) = radius of ellipse measured in x(sunward) direction.
!       b(l) = radius of ellipse measured in y(duskward) direction.
!       dx(l) = sunward displacement of coord.system center from pole.
!       dy(l) = duskward displ. of coord. system center from pole.
!      to use efield with rcm, choose a(2) =.5*(colat(imin,noon)+
!        colat(imin,midnt)), b(2)=.5*(colat(imin,dusk)+colat(imin,dawn))
!      dx(2) = -offset*180./pi, dy(2) = 0.  ellipse parameters for
!      boundaries 1 and 3 need not be specified. All a's, b's, dx's and
!      dy's are in degrees.
!
!      colat and aloct are the usual magnetic colatitude and magnetic
!      local time angles in radians.
!
!      mode is a dummy parameter, not currently used.
!
!      imax, jmax, jwrap have their usual meaning, same as in rcm.
!
!      itmdim = number of time labels in the v-matrix.  for use with
!               rcm, it should be set equal to 1.
!
!      itmcur = time-label number of v-matrix currently being computed.
!               for use with rcm, it should be set equal to 1.
!
!      icntrl = 1 means different formulas are used in regions 1,2,3.
!             = 0 means heppner-maynard formula is used for all
!                 latitudes, but scaled to externally specified pcp and
!                 ellipse parameters.
!             = -1 means heppner-maynard formula is used for all
!                  latitudes, unscaled.
!             = -2 is the mode used for use with rcm.  heppner-maynard
!                  is used poleward of ellipse 2, and on ellipse 2,
!                  scaled to fit specified potential drop and dimensions
!                  of ellipse 2.
!
!      itop = 1 means that a(1), b(1), dx(1), dy(1) are computed
!               internally in efield.  otherwise they must be passed to
!               efield from outside.  itop is normally = 1 for use in
!
!               msm.  itop should be set to 1 for use with rcm.
!      v = potential matrix computed in efield.  it is dimensioned at
!          latdim x ltdim x itmdim, where latdim and ltdim are specified
!          in parameter statement in efield and itmdim is passed.
!
!
!
!
!     subroutine efield calls subroutines epot(heppner-maynard),
!     low(low-latitudes), aurora(region 2), reg1(region 1).
!
!
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: icntrl, itop, ipatt, mode
      REAL (rprec), INTENT (IN) :: pcp, deqdt, colat (:,:), aloct (:,:)
      REAL (rprec), INTENT (IN OUT) :: a(3), b(3), dx(3), dy(3), v (:,:)
      REAL (rprec), INTENT (OUT) :: vbndry (:)
!
!
      REAL (rprec):: thetaa (SIZE(vbndry)),thetab (SIZE(vbndry)), &
!             thetac (SIZE(vbndry)),&
              coef (18,18), xco (18,18), ahm (2,7), &
              bhm (2,7), dxhm(2,7), dyhm(2,7), vmin(7), vmax(7),&
              fg, al, gmlat, glong, vtemp
      INTEGER (iprec) :: ieff, l, j, i, nnmax, idim, jdim
!
!     parameter values
!
!       ntape = file number from which maynard-rich coefficients are
!               read.  dataset name is 'rawolf.efcoef.data'.
!
!
!     Stanislav: since now the grid is such that aloct(i,j) is
!     the same for all i and given J, ieff parameter is not
!     relevant. But watch out if this changes later. Also, the
!     boundary is assumed to be an ellipse (since we are
!     calling this subroutine!), but it does not have to 
!     coincide with the integer grid line. Location of bndy
!     is given fully by the parameters of the ellipse.
!
      ieff   = 1
      idim   = SIZE (v, DIM = 1)
      jdim   = SIZE (v, DIM = 2)
!
      IF (icntrl /= -2.OR.itop /= 1) STOP 'EFIELD IS FOR RCM ONLY'
!
      IF (ABS(deqdt) > zero .OR. mode /= 0) STOP 'EFIELD ERROR'
      nnmax = 16
!
!
      CALL Input (coef, ipatt, xco, ahm, bhm, dxhm, dyhm, vmin, vmax)
!
!
!   fg = scaling factor that converts hm potential to
!         present situation
!
      IF (icntrl /= -1) THEN
         fg = pcp/(vmax(ipatt)-vmin(ipatt))
      ELSE
         fg = one
      END IF
      IF (ABS(fg - one) > machine_eps1) STOP 'ERROR IN EFIELD'
!
!
!     compute a(1),b(1),dx(1),dy(1) if itop = 1.
!
      IF (itop == 1) THEN
          a(1) = a(2) * ahm(1,ipatt) / ahm(2,ipatt)
          b(1) = b(2) * bhm(1,ipatt) / bhm(2,ipatt)
         dx(1) = dx(2) + a(2) * (dxhm(1,ipatt)-dxhm(2,ipatt)) /&
                                ahm(2,ipatt)
         dy(1) = dy(2) + b(2) * (dyhm(1,ipatt)-dyhm(2,ipatt)) /&
                                bhm(2,ipatt)
      END IF
!
!!    vmaxx = vmax(ipatt)*fg
!!    vminn = vmin(ipatt)*fg
!!    sa=sin((a(3)-dx(3))*pi/180.)
!!    vpenet = -13112.*g*deqdt*sa
!!    vbar = 0.6*vmaxx + 0.4*vminn -0.22*vpenet
!
!
!    computation of boundary locations as fcns of j.
!
!   guard against irrelevant a's or b's being zero or negative.
!
      DO  l=1,3
         if(a(l).lt.0.001) a(l)=0.001
         if(b(l).lt.0.001) b(l) = 0.001
      END DO
!
!
      DO j = 1, jdim
         al        = aloct (ieff,j)
         thetaa(j) = thet(a(1),b(1),dx(1),dy(1),al)*DTR
         thetab(j) = thet(a(2),b(2),dx(2),dy(2),al)*DTR
!        thetac(j) = thet(a(3),b(3),dx(3),dy(3),al)*DTR
      END DO
!
!     main do loop over grid
!
!
!
      DO j = 1, jdim
      DO i = 1, idim
!
         IF (colat(i,j) <= thetaa(j)) THEN
!
!          Region zero: polar cap:
!
            glong = MOD (aloct(i,j)+pi, pi_two)
            gmlat = 90.0 - colat(i,j)*RTD
            CALL Epot (coef, gmlat, glong, vtemp, &
                       ipatt, nnmax, a, b, dx, dy, xco, &
                       ahm, bhm, dxhm, dyhm, vmin, vmax, pcp, icntrl)
            v (i,j) = vtemp
!
         ELSE IF (colat(i,j) < thetab(j)+5.0E-4 .AND. &
                  colat(i,j) > thetaa(j) .AND. &
                  icntrl == -2) THEN
!
!          use heppner-maynard directly in region 1 for 
!          use with rcm.
!
            glong = MOD (aloct(i,j) + pi, pi_two)
            gmlat = 90.0 - colat(i,j)*RTD
            CALL Epot (coef, gmlat, glong, vtemp, &
                       ipatt, nnmax, a, b, dx, dy, xco, &
                       ahm, bhm, dxhm, dyhm, vmin, vmax, pcp, icntrl)
            v ( i,j ) = vtemp
         END IF
!
!
      END DO
      END DO
!
!
!   segment for defining boundary potentials and e-fields
!
!     potential at b
!
      DO j = 1, jdim
         glong = MOD (aloct(ieff,j)+pi, pi_two)
         gmlat = 90.0 - thetab(j)*RTD
         CALL Epot (coef, gmlat, glong, vbndry(j), ipatt,&
                    nnmax, a, b, dx, dy, xco, ahm, bhm, &
                    dxhm, dyhm, vmin, vmax, pcp, icntrl)
      END DO
!
!
      RETURN
      END SUBROUTINE Efield
!
!
!
      SUBROUTINE Input (coef, ipatt, xco, ahm, bhm, dxhm, dyhm, vmin, vmax)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: ipatt
      REAL (rprec), INTENT (OUT) :: coef (18,18), dxhm(2,7), &
                            dyhm(2,7), xco(18,18), &
                            ahm(2,7), bhm(2,7), &
                            vmin(7), vmax(7)
!
!
      INTEGER (iprec) :: ip, n, i, j, ia, mmm, nnn
      REAL (rprec):: ccf
!     CHARACTER (LEN=50) :: mlbl
!
!  Subroutine that provides coefficients for Heppner-Maynard
!  polar cap convection model for use with RCM. yco arrays
!  holds Legendre coefficients (same for all HM models),
!  these come from the HM model.
! 
!  AAHM, BBHM, DDXHM, DDYHM arrays hold parameters of two
!  ellipses for each HM model. First ellipse (with index 1) 
!  is the poleward boundary of field-reversal region, and
!  second ellipse is the equatorward boundary of that region.
!  Second ellipse corresponds to the main RCM bounday.
!  VVMAX and VVMIN hold maximun (dawn) and minimun (dusk) 
!  potentials for each HM model.  All these arrays (as I
!  understand) were obtained by Bob Spiro.
!
!
      REAL (rprec), PARAMETER :: yco (18,18) = RESHAPE ( (/ &
 .282095e+00, .488603e+00, .109255e+01, .228523e+01, .468333e+01, .951188e+01, .192265e+02, .387523e+02, .779645e+02, &
 .156658e+03, .314501e+03, .630964e+03, .126523e+04, .253611e+04, .508196e+04, .101809e+05, .203918e+05, .408366e+05,&
!
 .488603e+00, .488603e+00, .546274e+00, .144531e+01, .331161e+01, .719031e+01, .151999e+02, .316411e+02, .652298e+02,&
 .133599e+03, .272366e+03, .553392e+03, .112151e+04, .226837e+04, .458082e+04, .923904e+04, .186151e+05, .374743e+05,&
!
 .946175e+00, .109255e+01, .546274e+00, .590044e+00, .177013e+01, .440314e+01, .101333e+02, .223736e+02, .481754e+02,&
 .102038e+03, .213661e+03, .443701e+03, .915709e+03, .188083e+04, .384866e+04, .785168e+04, .159791e+05, .324537e+05,&
!
 .186588e+01, .228523e+01, .144531e+01, .590044e+00, .625836e+00, .207566e+01, .555021e+01, .134918e+02, .310971e+02,&
 .693209e+02, .151081e+03, .324033e+03, .686782e+03, .144253e+04, .300864e+04, .623988e+04, .128827e+05, .264983e+05,&
!
 .370249e+01, .468333e+01, .331161e+01, .177013e+01, .625836e+00, .656382e+00, .236662e+01, .674590e+01, .172496e+02,&
 .414272e+02, .955522e+02, .214328e+03, .471128e+03, .102002e+04, .218269e+04, .462762e+04, .973844e+04, .203694e+05, &
!
 .736787e+01, .951188e+01, .719031e+01, .440314e+01, .207566e+01, .656382e+00, .683184e+00, .264596e+01, .798499e+01,&
 .213929e+02, .534153e+02, .127330e+03, .293800e+03, .661878e+03, .146420e+04, .319336e+04, .688612e+04, .147131e+05,&
!
 .146845e+02, .192265e+02, .151999e+02, .101333e+02, .555021e+01, .236662e+01, .683184e+00, .707163e+00, .291571e+01,&
 .926339e+01, .259102e+02, .671087e+02, .165101e+03, .391572e+03, .903721e+03, .204248e+04, .454057e+04, .996084e+04,&
!
 .292940e+02, .387523e+02, .316411e+02, .223736e+02, .134918e+02, .674590e+01, .264596e+01, .707163e+00, .728927e+00,&
 .317732e+01, .105778e+02, .307916e+02, .825507e+02, .209304e+03, .509767e+03, .120459e+04, .278052e+04, .629979e+04,&
!
 .584734e+02, .779645e+02, .652298e+02, .481754e+02, .310971e+02, .172496e+02, .798499e+01, .291571e+01, .728927e+00,&
!  in line above, I changed the first constant from .30... to 
!    .31... since Heppner-Maynard official release has it
!    this way, in constract to Bob Spiro's verion. Stanislav
!    March 8 1999.
!
 .748901e+00, .343190e+01, .119255e+02, .360281e+02, .997819e+02, .260366e+03, .650553e+03, .157290e+04, .370647e+04,  &
!
 .116766e+03, .156658e+03, .133599e+03, .102038e+03, .693209e+02, .414272e+02, .213929e+02, .926339e+01, .317732e+01,&
 .748901e+00, .767395e+00, .368030e+01, .133043e+02, .416119e+02, .118840e+03, .318704e+03, .816138e+03,  .201755e+04, &
!
 .233240e+03, .314501e+03, .272366e+03, .213661e+03, .151081e+03, .955522e+02, .534153e+02, .259102e+02, .105778e+02, &
.343190e+01, .767395e+00, .784642e+00, .392321e+01, .147120e+02, .475361e+02, .139761e+03, .384731e+03,   .100877e+04, &
!
 .465998e+03, .630964e+03, .553392e+03, .443701e+03, .324033e+03, .214328e+03, .127330e+03, .671087e+02, .307916e+02, &
 .119255e+02, .368030e+01, .784642e+00, .800822e+00, .416119e+01, .161472e+02, .537941e+02, .162579e+03, .458849e+03, &
!
 .931187e+03, .126523e+04, .112151e+04, .915709e+03, .686782e+03, .471128e+03, .293800e+03, .165101e+03, .825507e+02, &
 .360281e+02, .133043e+02, .392321e+01, .800822e+00, .816077e+00, .439471e+01, .176082e+02, .603802e+02, .187325e+03, &
!
 .186100e+04, .253611e+04, .226837e+04, .188083e+04, .144253e+04, .102002e+04, .661878e+03, .391572e+03, .209304e+03,&
 .997819e+02, .416119e+02, .147120e+02, .416119e+01, .816077e+00, .830522e+00, .462415e+01, .190939e+02, .672889e+02, &
!
 .371962e+04, .508196e+04, .458082e+04, .384866e+04, .300864e+04, .218269e+04, .146420e+04, .903721e+03, .509767e+03,&
 .260366e+03, .118840e+03, .475361e+02, .161472e+02, .439471e+01, .830522e+00, .844251e+00, .484985e+01, .206029e+02,&
!
 .743510e+04, .101809e+05, .923904e+04, .785168e+04, .623988e+04, .462762e+04, .319336e+04, .204248e+04, .120459e+04, &
 .650553e+03, .318704e+03, .139761e+03, .537941e+02, .176082e+02, .462415e+01, .844251e+00, .857341e+00, .507210e+01, &
!
 .148629e+05, .203918e+05, .186151e+05, .159791e+05, .128827e+05, .973844e+04, .688612e+04, .454057e+04, .278052e+04,&
 .157290e+04, .816138e+03, .384731e+03, .162579e+03, .603802e+02, .190939e+02, .484985e+01, .857341e+00, .869857e+00,&
!
 .297130e+05, .408366e+05, .374743e+05, .324537e+05, .264983e+05, .203694e+05, .147131e+05, .996084e+04, .629979e+04,&
 .370647e+04, .201755e+04, .100877e+04, .458849e+03, .187325e+03, .672889e+02, .206029e+02, .507210e+01,   .869857e+00&
  /), (/18,18/), ORDER = (/2,1/)  )
!
      REAL, PARAMETER :: aahm(2,7) = RESHAPE ( (/ &
            17.45,12.13,16.06,13.72,14.79,12.82,15., &
            20.00,16.17,18.88,19.31,18.3,16.6,17.93 /), &
            (/2,7/), ORDER = (/2,1/) )
      REAL, PARAMETER :: bbhm(2,7) = RESHAPE ( (/ &
            14.26,13.78,14.31,13.78,14.73,11.70,12.07, &
            16.97,16.7,17.07,17.23,17.93,17.5,16.6 /), &
            (/2,7/), ORDER = (/2,1/) )
      REAL, PARAMETER :: ddxhm(2,7) = RESHAPE ( (/ &
            -2.66,-5.53,-3.19,-.11,-2.45,-1.65,-2.23, &
            -3.30,-7.45,-4.1,-1.76,-2.77,-2.98,-3.14/), &
            (/2,7/), ORDER = (/2,1/) )
      REAL, PARAMETER :: ddyhm(2,7) = RESHAPE ( (/ &
            1.60,.8,1.12,.05,.48,3.4,1.22, &
            1.33,.32,1.54,.85,.69,1.86,.21 /), &
            (/2,7/), ORDER = (/2,1/) )
      REAL, PARAMETER :: vvmax(7)= (/ &
            34007.,55354.,14390.,11287.,9329.,13249., 13221. /)
      REAL, PARAMETER :: vvmin(7)= (/ &
            -42280.,-16003.,-60935.,-16250.,-12947., -14428.,-13460./)
!
!
      DO ip = 1,7
         DO  n = 1, 2
            ahm(n,ip) = aahm(n,ip)    
            bhm(n,ip) = bbhm(n,ip)   
            dxhm(n,ip) = ddxhm(n,ip)
            dyhm(n,ip) = ddyhm(n,ip)
         END DO
         vmax(ip) = vvmax(ip)
         vmin(ip) = vvmin(ip)
       END DO
!
!
!   Initialize coef from unit ntape:
!
      DO i = 1, 18
         DO j = 1, 18
            coef (i,j) = 0.
            xco  (i,j) = yco (i,j)
         END DO
      END DO
!
      OPEN (UNIT = LUN, FILE = rcmdir//'HMRCOEF', STATUS = 'OLD')
      DO ia = 1, ipatt
         READ (LUN, FMT = '(A)')
         DO
            READ (LUN,*) nnn,mmm,i,j,ccf
            IF (nnn == -1) EXIT
            IF (i > 18 .OR. j > 18) STOP 'out of bounds in INPUT'
            IF (mmm < 0 ) STOP 'MMM IN EPOT'
            coef (i,j) = ccf
         END DO
      END DO
      CLOSE (UNIT = LUN)
!
      RETURN
      END SUBROUTINE Input
!
!
!
      SUBROUTINE Epot (coef, tlat, tlon, value,  ipatt, nnmax, a, b, &
                       dx, dy, xco, ahm, bhm, dxhm, dyhm, vmin, vmax, &
                       pcp, icntrl)
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: ahm(2,7), bhm(2,7), dxhm(2,7), dyhm(2,7),&
                                   vmin(7), vmax(7), a(3), b(3), dx(3), dy(3),&
                                   coef(18,18), xco (18,18)
      REAL (rprec), INTENT (IN) :: tlat, tlon, pcp
      REAL (rprec), INTENT (OUT) :: value
      INTEGER (iprec), INTENT (IN) :: icntrl, ipatt, nnmax
!
      REAL (rprec) :: dp(18,18), sp(18), cp(18),  &
!                     fn(18), fm(18),  &
                      p (18,18) = zero, const(18,18)= zero
!
!
      INTEGER (iprec) :: nmax, n, m
      REAL (rprec)   :: cph, sph, st, ct, beta, alpha, xlon, xl, xlat, &
                        xcol, yy, xx, tlong, tcol, pol
      REAL (rprec), PARAMETER :: cmin=50.0_rprec, cmax = 90.0_rprec
!
!
      nmax = nnmax
!
      value = -1.0E-9_rprec
      IF (tlat <= cmin) THEN
         WRITE (*,*) 'latitude ',tlat,' out of model range'
         RETURN
      END IF
!
!
!   segment for scaling grid-pt. location to fit h-m ellipse.
!
      tcol = 90.0_rprec - tlat
      tlong = tlon - pi
      IF (icntrl /= -1) THEN
         xx = dxhm(1,ipatt)+ahm(1,ipatt)*(tcol*COS(tlong)-dx(1))/a(1)
         yy = dyhm(1,ipatt)+bhm(1,ipatt)*(tcol*SIN(tlong)-dy(1))/b(1)
      ELSE
         xx = tcol*COS(tlong)
         yy = tcol*SIN(tlong)
      END IF
      xcol  = SQRT(xx**2 + yy**2)
      xlat  = 90.0_rprec - xcol
      xl    = ATAN2(yy,xx)
      xlon  = xl+pi
      alpha = two / (cmax-cmin)
      beta  = one - alpha*cmax
      ct    = xlat*alpha+beta
      st    = SQRT(one - ct*ct)
      sph   = SIN(xlon)
      cph   = COS(xlon)
!
      IF (ABS(p(1,1)) < machine_tiny) THEN
         p(1,1)  = one
         dp(1,1) = zero
         sp(1)   = zero
         cp(1)   = one
         DO n = 2, 18
!            fn(n) =n
            DO m = 1, n
!               fm(m) = m-1
               const(n,m) = REAL((n-2)**2-(m-1)**2, rprec)/&
                            REAL((2*n-3)*(2*n-5), rprec)
            END DO
         END DO
      END IF
!
      sp(2)  = sph
      p(1,1) = one
      cp(2)  = cph
      do m = 3, nmax
         sp(m)=sp(2)*cp(m-1)+cp(2)*sp(m-1)
         cp(m)=cp(2)*cp(m-1)-sp(2)*sp(m-1)
      END DO
!
      value=coef(1,1)
!
      do n = 2, nmax
      do m = 1, n
         IF (n == m) THEN
            p(n,n)=st*p(n-1,n-1)
           dp(n,n)=st*dp(n-1,n-1)+ct*p(n-1,n-1)
         ELSE
            IF (n /= 2) p(n,m)=ct*p(n-1,m)-const(n,m)*p(n-2,m)
            IF (n == 2) p(n,m)=ct*p(n-1,m)
            dp(n,m)=ct*dp(n-1,m)-st*p(n-1,m)-const(n,m)*dp(n-2,m)
         END IF
      END DO
      END DO
!
!
!
      p(1,1) = p(1,1)*xco(1,1)
!
      value = coef(1,1)*p(1,1)
      DO n = 2, nmax
      DO m = 1, n
         IF (m /= 1) THEN
            pol      = p(n,m)*xco(n,m)
            p(m-1,n) = cp(m)*pol
            p(n,m)   = sp(m)*pol
            value    = value+p(m-1,n)*coef(m-1,n)+p(n,m)*coef(n,m)
         ELSE
            p(n,m) = p(n,m)*xco(n,m)
            value  = value+p(n,m)*coef(n,m)
         END IF
      END DO
      END DO
      IF (icntrl /= -1) THEN
         value = (pcp/(vmax(ipatt)-vmin(ipatt)))*value*1000.
      ELSE
         value = 1000.0_rprec*value
      END IF
!
!
      RETURN
      END SUBROUTINE Epot
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      SUBROUTINE Elemod (icase, ikp, glat, amlt, value)
      IMPLICIT NONE
      INTEGER (iprec), INTENT (IN) :: icase, ikp
      REAL (rprec),    INTENT (IN) :: glat, amlt
      REAL (rprec),    INTENT (OUT):: value
!
!______________________________________________________________________________
!
!   Stanislav Tue May 25 20:40:19 MDT 1999: adopt the hardy
!             code to use with RCM. This subroutine communicates
!             with RCM through the arguments, but also needs
!             module rcm_mod_hardy to store coeffs.
!
!   SUBROUTINE TO EVALUATE THE HARDY AVERAGE AURORAL MODEL
!     DESCRIBED IN:  HARDY, ET AL., "J. GEOPHYS. RES.",
!     VOL. 92, PAGE 12,275-12,294 (1987).
!
!
!   INPUTS:
!
!     ICASE=1   ENERGY FLUX
!           2   NUMBER FLUX
!           3   HALL CONDUCTIVITY
!           4   PEDERSON CONDUCTIVITY
!
!     IKP       KP DIVISION 0 TO 6
!
!     GLAT      GEOMAGNETIC LATITUDE
!
!     AMLT      MAGNETIC LOCAL TIME
!
!   OUTPUTS:
!
!     VALUE     LOG10 ENERGY FLUX IN KEV/CM**2-S-SR (ICASE=1)
!               LOG10 NUMBER FLUX IN ELECTRONS/CM**2-S-SR (ICASE=2)
!               CONDUCTIVITY IN MHOS (ICASE=3,4)
!  
!   INTERNAL VARIABLES
!
!     CRD       COEFFICIENTS FOR MAXIMUM FUNCTION VALUE
!     CHAT      COEFFICIENTS FOR LATITUDE OF MAXIMUM VALUE
!     CS1       COEFFICIENTS FOR UP-SLOPE
!     CS2       COEFFICIENTS FOR DOWN-SLOPE
!     CUTL      LOW LATITUDE CUT-OFF VALUE
!     CUTH      HIGH LATITUDE CUT-OFF VALUE
!
!   FILES:
!
!     THE FILE ELECOEF.DAT MUST BE PRESENT IN THE DEFAULT
!     DIRECTORY.
!
!   NOTES:
!
!     THIS VERSION OPERATES ON VAX/VMS OR IBM-PC MACHINES.
!
!   21 JUNE 1993 -- WJM
!______________________________________________________________________________
!
      REAL (rprec), DIMENSION(4), PARAMETER :: &
                    cutl=(/6.,6.,0.,0./), cuth=(/7.,7.,.55,.55/)
!      CHARACTER(LEN=80) :: aline
      INTEGER (iprec) :: j, jcase, jco, jkp, kp, ipc, ips
      REAL (rprec) :: xarg, rd, hat, s1, s2, xa, c, s
      INTEGER (iprec), SAVE :: iread_hardy_first = 0
      REAL (rprec), SAVE :: crd(13,7,4), chat(13,7,4), cs1(13,7,4), cs2(13,7,4)
!
!
      IF (ikp > 6) STOP 'ELEMOD: KP > 6 !'  ! stanislav
!
!
      IF (iread_hardy_first == 0) THEN
        iread_hardy_first = 1
!
        OPEN (UNIT = LUN, FILE = rcmdir//'elecoef.dat', STATUS = 'OLD')
        READ (LUN,'(A80)')  ! aline
        DO jcase = 1, 4
           DO jkp = 1, 7
              DO jco = 1, 13
                 READ (LUN, '(26X,4F12.7)') crd (jco,jkp,jcase), &
                                            chat (jco,jkp,jcase),&
                                            cs1 (jco,jkp,jcase), &
                                            cs2 (jco,jkp,jcase)
              END DO
           END DO
        END DO
        CLOSE (LUN)
      END IF
!
!
      IF (glat < 50.0) RETURN
!
      kp   = ikp + 1
      xarg = amlt*3.14159265/12.
      rd   = crd(1,kp,icase)
      hat  = chat(1,kp,icase)
      s1   = cs1(1,kp,icase)
      s2   = cs2(1,kp,icase)
!
      DO j = 1, 6
         xa  = j*xarg
         c   = COS(xa)
         s   = SIN(xa)
         ipc = j+1
         ips = j+7
!
         rd = rd+c*crd(ipc,kp,icase)+s*crd(ips,kp,icase)
         hat= hat+c*chat(ipc,kp,icase)+s*chat(ips,kp,icase)
         s1 = s1+c*cs1(ipc,kp,icase)+s*cs1(ips,kp,icase)
         s2 = s2+c*cs2(ipc,kp,icase)+s*cs2(ips,kp,icase)
!
      END DO
!
      value = Epst (glat,rd,hat,s1,s2,cutl(icase),cuth(icase))
!
      RETURN
!
      CONTAINS
!
          FUNCTION Epst (clat, rd, hat, s1, s2, xmin, xmax)
          IMPLICIT NONE
          REAL (rprec), INTENT (IN) :: clat, rd, hat, s1, s2, xmin, xmax
          REAL (rprec):: Epst
    !
    !
          REAL (rprec) :: d, ex, xl, ep
    !
          d  = clat-hat
          ex = EXP (d)
          xl = (one - s1 / s2 * ex) / (one - s1 / s2)
          xl = LOG (xl)
          ep = rd+s1*d+(s2-s1)*xl
    !
          IF (clat < hat .AND. ep < xmin) ep=xmin
          IF (clat > hat .AND. ep < xmax) ep=xmax
    !
          epst = ep
          RETURN
          END FUNCTION epst
!
      END SUBROUTINE Elemod
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      FUNCTION Thet (aa, bb, xc, yc, phi)
!
!     this function gives solution to equation for ellipse in
!     in flat polar coordinates.  it specifies colatitude as a
!     function of local time.
!
      IMPLICIT NONE
      REAL (rprec), INTENT (IN) :: aa, bb, xc, yc, phi
      REAL (rprec) :: thet
!
      REAL (rprec) :: cphi, sphi, ca, cb, cc
!
!      
         cphi = COS(phi)
         sphi = SIN(phi)
         ca = (cphi/aa)**2 + (sphi/bb)**2
         cb = -xc*cphi/aa**2 - yc*sphi/bb**2
         cc = (xc/aa)**2 + (yc/bb)**2 - 1.
         thet = (-cb+SQRT(cb**2-ca*cc))/ca
      RETURN
      END FUNCTION Thet
!
!
!
!   6. A set of facilities to interpolate 1-D, 2-D arrays and 2-D sections
!      of 3-D arrays, real numbers:


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
    IF (bj < 1.0 .OR. bj > nj) pause 'BJ out of range in gntrp_2d_ang'
    IF (bi < 1.0 .OR. bi > ni) pause 'BI out of range in gntrp_2d_ang'
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
!
!
!
!
!   FUNCTION Gntrp_2d (array, bi, bbj)
!   IMPLICIT NONE
!   REAL (rprec), INTENT (IN)    :: array (:,:), bi, bbj
!   REAL (rprec)                 :: Gntrp_2d
!
!**********************************************************
! jwrap is always assumed to be the same!!!!!!!!!!!!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  Author: r.w. spiro
!
!  Last update: 5-08-88 by rws         
!              10-05-98 by stanislav: 
!               corrected a minor bug by moving
!               line "fi=int(ii)" after correction for ii.
!              02-22-99 stanislav:
!               added a line for case of Ikind = 1 (see code)
!
!  If A is a 2-d array, then one  dimensional interpolation
!  is used if bj is within .0001 units of an integer j value.
!  Otherwise, a full 2-d interpolation is performed.
!  If a is an angular variable (eg., p or aloct) set ikind=1;
!  otherwise, set ikind=0
!  Description of input parameters:                        
!     a    = array to be interpolated, either 1-d or 2-d  
!     imax = leading dimension of a                      
!     jmax = second dinension of a.  jmax=1 if a is 1-d. 
!     jwrap= wrap around size for 2-d array.  jwrap=0 for 1-d.          
!     bi   = i location of interpolated point           
!     bbj  = j location of interpolated point. bj=1.0 for 1-d.
!     ikind=0 for all ano-angular arrays                                
!          =1 for angular arrays (eg., p and aloct)                     
!
!  Stanislav: if Bi < 1, then the array is extrapolated 
!             linearly based on the values A(1,:) and A(2,:).
!             If Bi > imax, then array is linearly 
!             extrapolated on the values A(imax-1,:) and 
!             A(imax,:).
!  Stanislav: if this function is to be used with IKIND=1
!             (for aloct), then the array must be 
!             circularized according to the RCM
!             rules (j=3,50, and then j+51 = j_3, 
!             j_1 = j_49, j_2 = j_50).
!                                                                       
!  Dependency:  BJMOD
!
!   INTEGER (iprec) :: ii, jn, jj, jp1, imax_array, jmax_array
!   REAL (rprec)    :: fi,fj,a1,a2,bj
!
!
!   Prepare indices for interpolation:
!
!   imax_array = SIZE (array, 1)
!   jmax_array = SIZE (array, 2)
!   ii = MAX (1, MIN (INT (bi), imax_array-1))
!   fi = REAL (ii,rprec)
!                                                                       
!                                                                       
!   Decide which interpolation to perform and proceed:
!
!   bj    = Bjmod ( bbj,jwrap,jmax_array)
!   jn    = NINT (bj)
!   IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp. of 2-d array
!
!     IF ( ii < imin_j(jn) .OR. ii+1 < imin_j(jn)) STOP 'BNDY PROBLEM IN GNTRP_N, 1'
!     Gntrp_2d = (one-(bi-fi)) * array(ii,jn) + (bi-fi)*array(ii+1,jn)
!
!   ELSE    !        2-d interpolation of 2-d array:
!
!         If jwrap <= bj < jmax, then jwrap-1 <= INT(bj) <= jmax-1
!         and jwrap <= INT(bj)+1 <= jmax
!
!      jj  = INT (bj)
!      fj  = REAL (jj,rprec)
!      jp1 = jj + 1
!
!     IF ( ii < imin_j(jj) .OR. ii+1 < imin_j(jj)) STOP 'BNDY PROBLEM IN GNTRP_N, 2'
!     IF ( ii < imin_j(jp1) .OR. ii+1 < imin_j(jp1)) STOP 'BNDY PROBLEM IN GNTRP_N, 2'
!      a1 = (one-(bi-fi))*array(ii,jj)  + (bi-fi)*array(ii+1,jj)
!      a2 = (one-(bi-fi))*array(ii,jp1) + (bi-fi)*array(ii+1,jp1)
!                                                                       
!      Gntrp_2d = (one - (bj-fj)) * a1 + (bj-fj) * a2
!
!   END IF
!   RETURN
!   END FUNCTION Gntrp_2d
!
!
!
!   FUNCTION Gntrp_2d_of3d (array, bi, bbj, ikind, index_3)
!   IMPLICIT NONE
!   INTEGER (iprec), INTENT (IN) :: ikind, index_3
!   REAL (rprec), INTENT (IN)    :: array (:,:,:), bi, bbj
!   REAL (rprec)                 :: Gntrp_2d_of3d
!                                                                       
!   This is the same as Gntrp_2d but for a 3-dim array, see comments for
!   Gntrp_2d. A separate function is needed since if Gntrp_2d were used,
!   then we would need to pass array sections (the other option is Fortran
!   77 style of passing an offset array, but that should be avoided for
!   compiler checking and parallelization reasons).
!
!   Dependency:  BJMOD
!
!   INTEGER (iprec) :: ii, jn, jj, jp1, imax_array, jmax_array, kmax_array
!   REAL (rprec)    :: fi,fj,a1,a2,bj
!
!
!   Prepare indices for interpolation:
!
!   imax_array = SIZE (array, 1)
!   jmax_array = SIZE (array, 2)
!   kmax_array = SIZE (array, DIM = 3)
!   IF (index_3 > kmax_array) STOP 'GNTRP_2D_OF3D: index_3 OUT OF RANGE'
!   ii = MAX (1, MIN (INT (bi), imax_array-1))
!   fi = REAL (ii,rprec)
!                                                                       
!                                                                       
!   Decide which interpolation to perform and proceed:
!
!   bj    = Bjmod ( bbj, jwrap, jmax_array)
!   jn    = NINT (bj)
!   IF (ABS (bj-REAL(jn,rprec)) < 1.0E-4_rprec) THEN  ! 1-d interp. of 2-d array
!
!     Gntrp_2d_of3d = (one-(bi-fi)) * array(ii,jn,index_3) + &
!                 (bi-fi)*array(ii+1,jn,index_3)
!
!   ELSE    !        2-d interpolation of 2-d array:
!
!         If jwrap <= bj < jmax, then jwrap-1 <= INT(bj) <= jmax-1
!         and jwrap <= INT(bj)+1 <= jmax
!
!      jj  = INT (bj)
!      fj  = REAL (jj,rprec)
!      jp1 = jj + 1
!
!      a1 = (one-(bi-fi))*array(ii,jj,index_3)  + (bi-fi)*array(ii+1,jj,index_3)
!      a2 = (one-(bi-fi))*array(ii,jp1,index_3) + (bi-fi)*array(ii+1,jp1,index_3)
!                                                                       
!      IF (ikind == 1) THEN
!         IF (jp1 == jmax_array) a2 = a2 + pi_two
!         IF (jj  == jwrap)      a1 = zero
!         IF (jp1 == jwrap)      a2 = a2 + pi_two      ! sts, feb 22
!      END IF
!
!      Gntrp_2d_of3d = (one - (bj-fj)) * a1 + (bj-fj) * a2
!
!   END IF
!   RETURN
!   END FUNCTION Gntrp_2d_of3d
!
!
!
!   6.2 A set of facilities to interpolate 1-D, 2-D arrays and 2-D sections
!      of 3-D arrays, real numbers, but these arrays are not circularized in J:


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
!
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
!
!
!   Routines to handle the high-lat. boundary (specify its location,
!   determine location by interpolation, etc.):
!
    FUNCTION Bndy (bndloc, bj)
    IMPLICIT NONE
    REAL (rprec), INTENT (IN) :: bndloc(:), bj
    REAL (rprec)              :: Bndy
!
!   Function to determine whether you are outside RCM
!   boundary as defined by array AIN; returns boundary 
!   location (non-integer I-value) for given bj
!   Written 1/25/96 frt                                   
!                                                                       
!   Dependency: BJMOD
!
    REAL (rprec)    :: bip, bim, bj_point, bjp, bjm
    INTEGER (iprec) :: jm, jp
!                                                                       
!   Call to BJMOD returns bj_point in [jwrap,jsize), then jm is in
!   [jwrap,jsize) and jp1 is in [jwrap+1,jsize], so all three indices are
!   within the required range.
!
    bj_point = Bjmod (bj, jwrap, jsize)
    jm       = INT (bj_point)
    jp       = jm + 1
    bjp = REAL (jp,rprec)
    bjm = REAL (jm,rprec)
!                                                                       
    bim = bndloc (jm)
    bip = bndloc (jp)
!                                                                       
    bndy = bim * (bjp - bj_point) + bip * (bj_point - bjm)
!                                                                       
    RETURN
    END FUNCTION Bndy
!
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
    SUBROUTINE Outp_real (r, ibeg, iend, iinc, jbeg, jend, jinc, &
                          xscale, ilabel, title, ntp, ncol)
    IMPLICIT NONE
    INTEGER (iprec), INTENT(IN) :: ibeg, iend, iinc, jbeg, jend, jinc, &
                                   ntp, ncol, ilabel(20)
    CHARACTER (LEN = * ), INTENT (IN) :: title
    REAL (rprec), INTENT (IN) :: r (:,:), xscale
!                                                                       
!------------------------------------------------------------
!
! author: r.w. spiro        last modified: 06-25-85                     
!                                                                       
! DESCRIPTION OF INPUT PARAMETERS:                          
!      r(isize,jsize)=array to be output                                
!      this subroutine can output selected elements of array     
!      ibeg= initial i value to be output                               
!      iend= final i value to be output                                 
!      iinc= i value increment for output                               
!      (note: (iend-ibeg)/iinc should be an integer)                    
!      jbeg,jend, and jinc are defined similarly                        
!      scale=scale factor                                               
!         if(scale.eq.0.) scale is calculated to give best 
!            display      
!         all elements or array are divided by scale before
!            being output
!      ilabel= vector of length 20 that gives the label                 
!      title=character string that identifies the array 
!             being output    
!      ntp= output unit number                                          
!      ncol= number of columns for output device (80 or 132)            
!                                                                       
!------------------------------------------------------------
!                                                                       
    INTEGER (iprec), PARAMETER :: jcol = 16
    REAL (rprec) ::  y (jcol), scale_tmp, sum0, sum1, sum2, test, ave, sd
    INTEGER (iprec) :: isize, jsize, mxline, jjcol, isclmx, i, j, itest, &
                       ipower, jstart, istart, ifinal, iutt, lcount, &
                       jfinal, jj, ii, jcount
!
!   initialization
!
    isize = SIZE (r, DIM = 1)
    jsize = SIZE (r, DIM = 2)
    scale_tmp = xscale
    isclmx = - 12
    sum0 = zero
    sum1 = zero
    sum2 = zero
    mxline = 57
    jjcol = (ncol - 4) / 8
    IF (ibeg > isize .OR. iend > isize .OR. &
        jbeg > jsize .OR. jend > jsize .OR. &
        ibeg < 1 .OR. iend < 1 .OR. jbeg < 1 .OR. jend < 1) THEN
        STOP 'INDICES WRONG IN OUTP'
    END IF
!                                                                       
!   if scale=0. then compute auto scale factor
!                                                                       
    IF (ABS(scale_tmp) < machine_tiny) THEN
!
!      Determine maximum # of places to left of decimal pt
       DO i = ibeg, iend, iinc
       DO j = jbeg, jend, jinc
          IF (ABS(r(i,j)) > machine_tiny) THEN
             test = LOG10 (ABS (r (i, j) ) )
             IF (test >= zero) then
                test = test + .000001
             ELSE
                test = test - .999999
             END IF
             itest = test
          ELSE
             itest = 0
          END IF
          IF (r (i, j) < zero) itest = itest + 1
          isclmx = MAX (isclmx, itest)
       END DO
       END DO
!
!      determine scale factor such that max # of places
!      to left of decimal point is 3:
!
       ipower = isclmx - 2
       scale_tmp = 10.0_rprec**ipower
    ELSE
       ipower = NINT (LOG10 (scale_tmp) )
    END IF
!                                                                       
!   Compute mean and std dev. of outputted array elements:
!
    DO i = 1, isize
    DO j = 1, jsize
       sum0 = sum0 + one
       sum1 = sum1 + r (i,j)
!       sum2 = sum2 + r (i,j) **2
    END DO
    END DO
    ave = sum1 / sum0
!    sd  = SQRT ((sum2-sum0*ave**2)/(sum0-one))
    sd = sum2
    sd = zero
!                                                                       
!   Output data:
!                                                                       
    jstart = jbeg
    istart = ibeg
    ifinal = iend
!
!   Below ilabel(6) is ITIME variable, ilabel(1) is IRDW  variable,
!                              
    WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6), ilabel(3),ilabel(4), &
                     ilabel(5)
    iutt = ilabel (2)
!
    WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
    lcount = 0
!                                                                       
!   Begin output loop:
!
    main_loop: DO
         jfinal = (jjcol - 1) * jinc + jstart 
         IF (jfinal > jend) jfinal = jend 
         IF (jstart > jend) EXIT main_loop
         IF (mxline-lcount <= 5) then 
!                                 start a new page                      
!
            WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                            ,ilabel(3),ilabel(4),ilabel(5)
            iutt = ilabel (2) 
            WRITE (ntp, 810) TRIM (title), ipower, iutt, ave, sd 
            lcount = 0 
         END IF
!
         WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc) 
         lcount = lcount + 2
         DO ii = istart, ifinal, iinc 
            lcount = lcount + 1 
            IF (lcount > mxline) then 
!                                    start a new page                   
!
               WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                               ,ilabel(3),ilabel(4),ilabel(5)
               iutt = ilabel (2) 
               WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd 
               WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc) 
               lcount = 2 
            END IF 
            jcount = 0 
            DO jj = jstart, jfinal, jinc 
               jcount = jcount + 1 
               y (jcount) = r (ii, jj) / scale_tmp
            END DO 
            WRITE (ntp, 840) ii, (y (jj), jj = 1, jcount) 
         END DO
         jstart = jfinal + jinc 
         istart = ibeg 
    END DO main_loop
!
    RETURN
800 FORMAT (//,TR3,'# ', I6.6, TR2, 'ut=', I6.6, TR2,     &
            'ITIME=',I5.5,'(', I2.2,':',I2.2,':',  I2.2,')' )
810 FORMAT (/,T4, A, '/(1.E', I2, ')', TR4, 'ut=', I6.6, TR2, &
            'ave=', ES11.4, TR2, 'sd=', ES11.4)
830 FORMAT ( / ,T3, 16(TR5,I3))
840 FORMAT (T2, I3, 16(F8.3))
    END SUBROUTINE Outp_real
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Outp_integer (r, ibeg, iend, iinc, jbeg, jend, jinc, &
                             xscale, ilabel, title, ntp, ncol)
    IMPLICIT NONE
    INTEGER (iprec), INTENT(IN) :: ibeg, iend, iinc, jbeg, jend, jinc, &
                                   ntp, ncol, ilabel(20)
    CHARACTER (LEN = * ), INTENT (IN) :: title
    REAL (rprec),INTENT (IN) :: xscale
    INTEGER (iprec), INTENT (IN) :: r (:,:)
!                                                                       
!------------------------------------------------------------
!
! author: r.w. spiro        last modified: 06-25-85                     
!                                                                       
! DESCRIPTION OF INPUT PARAMETERS:                          
!      r(isize,jsize)=array to be output                                
!      this subroutine can output selected elements of array     
!      ibeg= initial i value to be output                               
!      iend= final i value to be output                                 
!      iinc= i value increment for output                               
!      (note: (iend-ibeg)/iinc should be an integer)                    
!      jbeg,jend, and jinc are defined similarly                        
!      scale=scale factor                                               
!         if(scale.eq.0.) scale is calculated to give best 
!            display      
!         all elements or array are divided by scale before
!            being output
!      ilabel= vector of length 20 that gives the label                 
!      title=character string that identifies the array 
!             being output    
!      ntp= output unit number                                          
!      ncol= number of columns for output device (80 or 132)            
!                                                                       
!------------------------------------------------------------
!                                                                       
    INTEGER (iprec), PARAMETER :: jcol = 16
    REAL (rprec)::  y (jcol), scale_tmp, sum0, sum1, sum2, test, ave, sd
    INTEGER (iprec) :: mxline, jjcol, isclmx, i, j, itest, ipower, &
                       jstart, istart, ifinal, iutt, lcount, jfinal, &
                       jj, ii, jcount, isize, jsize
!
!   initialization
!
    isize = SIZE (r, DIM = 1)
    jsize = SIZE (r, DIM = 2)
    scale_tmp = xscale
    isclmx = - 12
    sum0 = zero
    sum1 = zero
    sum2 = zero
    mxline = 57
    jjcol = (ncol - 4) / 8
    IF (ibeg > isize .OR. iend > isize .OR. &
        jbeg > jsize .OR. jend > jsize .OR. &
        ibeg < 1 .OR. iend < 1 .OR. jbeg < 1 .OR. jend < 1) THEN
        STOP 'INDICES WRONG IN OUTP'
    END IF
!
!   if scale=0. then compute auto scale factor
!                                                                       
    IF (ABS(scale_tmp) < machine_tiny) then
!
!      Determine maximum # of places to left of decimal pt
       DO i = ibeg, iend, iinc
       DO j = jbeg, jend, jinc
          IF (ABS(r(i,j)) > machine_tiny) THEN
             test = LOG10 (ABS (REAL(r (i, j)) ) )
             IF (test >= zero) THEN
                test = test + .000001
             ELSE
                test = test - .999999
             END IF
             itest = test
          ELSE
             itest = 0
          END IF
          IF (r (i, j) < zero) itest = itest + 1
          isclmx = MAX (isclmx, itest)
       END DO
       END DO
!
!      determine scale factor such that max # of places
!      to left of decimal point is 3:
!
       ipower = isclmx - 2
       scale_tmp = 10.0_rprec**ipower
    ELSE
       ipower = NINT (LOG10 (scale_tmp) )
    END IF
!                                                                       
!   Compute mean and std dev. of outputted array elements:
!
    DO i = 1, isize
    DO j = 1, jsize
       sum0 = sum0 + one
       sum1 = sum1 + REAL (r (i,j))
!       sum2 = sum2 + REAL (r (i,j))**2
    END DO
    END DO
    ave = SUM1 / sum0
!    sd  = SQRT ((sum2-sum0*ave**2)/(sum0-one))
    sd = sum2
    sd = zero
!                                                                       
!   Output data:
!                                                                       
    jstart = jbeg
    istart = ibeg
    ifinal = iend
!
!                   Below ilabel(6) is ITIME variable,
!                         ilabel(1) is IRDW  variable,
!                              
    WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6), &
                     ilabel(3),ilabel(4),ilabel(5)
    iutt = ilabel (2)
!
    WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
    lcount = 0
!                                                                       
!   Begin output loop:
!
    main_loop: DO
       jfinal = (jjcol - 1) * jinc + jstart
       IF (jfinal > jend) jfinal = jend
       IF (jstart > jend) EXIT main_loop
       IF (mxline-lcount <= 5) then
!                               start a new page
!
          WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                          ,ilabel(3),ilabel(4),ilabel(5)
          iutt = ilabel (2)
          WRITE (ntp, 810) TRIM (title), ipower, iutt, ave, sd
          lcount = 0
       END IF
!
       WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc)
       lcount = lcount + 2
       DO ii = istart, ifinal, iinc
          lcount = lcount + 1
          IF (lcount > mxline) then
!                                  start a new page
!
             WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                             ,ilabel(3),ilabel(4),ilabel(5)
             iutt = ilabel (2)
             WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
             WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc)
             lcount = 2
          END IF
          jcount = 0
          DO jj = jstart, jfinal, jinc
             jcount = jcount + 1
             y (jcount) = r (ii, jj) / scale_tmp
          END DO
          WRITE (ntp, 840) ii, (y (jj), jj = 1, jcount)
       END DO
       jstart = jfinal + jinc
       istart = ibeg
    END DO main_loop
!
800 FORMAT (//, T3, '# ', I6.6, TR2, 'ut=', I6.6, TR2,     &
            'ITIME=', I5.5, '(', I2.2, ':', I2.2, ':', I2.2, ')' )
810 FORMAT (/, T4, A, '/(1.E', I2, ')', TR4, 'ut=', I6.6, TR2, &
            'ave=', ES11.4, TR2, 'sd=', ES11.4)
830 FORMAT ( / , T3, 16(TR5,i3))
840 FORMAT (T2, I3, 16(TR2,I6))
    RETURN
    END SUBROUTINE Outp_integer
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
    SUBROUTINE Outp_logical (r, ibeg, iend, iinc, jbeg, jend, jinc, &
                             xscale, ilabel, title, ntp, ncol)
    IMPLICIT NONE
    INTEGER (iprec), INTENT(IN) :: ibeg, iend, iinc, jbeg, jend, jinc, &
                                   ntp, ncol, ilabel(20)
    CHARACTER (LEN = * ), INTENT (IN) :: title
    REAL (rprec),INTENT(IN) :: xscale
    LOGICAL, INTENT(IN) :: r (:,:)
!                                                                       
!------------------------------------------------------------
!
! author: r.w. spiro        last modified: 06-25-85                     
!                                                                       
! DESCRIPTION OF INPUT PARAMETERS:                          
!      r(isize,jsize)=array to be output                                
!      this subroutine can output selected elements of array     
!      ibeg= initial i value to be output                               
!      iend= final i value to be output                                 
!      iinc= i value increment for output                               
!      (note: (iend-ibeg)/iinc should be an integer)                    
!      jbeg,jend, and jinc are defined similarly                        
!      scale=scale factor                                               
!         if(scale.eq.0.) scale is calculated to give best 
!            display      
!         all elements or array are divided by scale before
!            being output
!      ilabel= vector of length 20 that gives the label                 
!      title=character string that identifies the array 
!             being output    
!      ntp= output unit number                                          
!      ncol= number of columns for output device (80 or 132)            
!                                                                       
!------------------------------------------------------------
!                                                                       
    INTEGER (iprec), PARAMETER :: jcol = 16
    LOGICAL ::  y (jcol)
!
    REAL (rprec) :: scale_tmp, ave, sd
    INTEGER (iprec) :: mxline, jjcol, ipower, &
                       jstart, istart, ifinal, iutt, lcount, jfinal, &
                       jj, ii, jcount, isize, jsize
!
!
    isize = SIZE (r, DIM = 1)
    jsize = SIZE (r, DIM = 2)
    mxline = 57
    jjcol = (ncol - 4) / 8
    iutt   = 0
    ave    = 0
    sd     = 0
    scale_tmp = xscale
    ipower = 0
    IF (ABS(scale_tmp) > machine_tiny) THEN
       STOP 'IF CALL OUTP_LOGICAL, SET SCALE = 0'
    END IF
    IF (ibeg > isize .OR. iend > isize .OR. &
        jbeg > jsize .OR. jend > jsize .OR. &
        ibeg < 1 .OR. iend < 1 .OR. jbeg < 1 .OR. jend < 1) THEN
        STOP 'INDICES WRONG IN OUTP'
    END IF
!
!
!   Output data:
!                                                                       
    jstart = jbeg
    istart = ibeg
    ifinal = iend
!
!                 Below ilabel(6) is ITIME variable, ilabel(1) is IRDW
!                              
    WRITE (ntp,800) &
          ilabel(1),ilabel(2),ilabel(6), ilabel(3),ilabel(4),ilabel(5)
    iutt = ilabel (2)
!
    WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
    lcount = 0
!                                                                       
!   Begin output loop:
!
    main_loop: DO
       jfinal = (jjcol - 1) * jinc + jstart
       IF (jfinal > jend) jfinal = jend
       IF (jstart > jend) EXIT main_loop
       IF (mxline-lcount <= 5) then
!                               start a new page
!
          WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                          ,ilabel(3),ilabel(4),ilabel(5)
          iutt = ilabel (2)
          WRITE (ntp, 810) TRIM (title), ipower, iutt, ave, sd
          lcount = 0
       END IF
!
       WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc)
       lcount = lcount + 2
       DO ii = istart, ifinal, iinc
          lcount = lcount + 1
          IF (lcount > mxline) then
!                                  start a new page
!
             WRITE (ntp, 800) ilabel(1),ilabel(2),ilabel(6)&
                             ,ilabel(3),ilabel(4),ilabel(5)
             iutt = ilabel (2)
             WRITE (ntp, 810) TRIM(title), ipower, iutt, ave, sd
             WRITE (ntp, 830) (jj, jj = jstart, jfinal, jinc)
             lcount = 2
          END IF
          jcount = 0
          DO jj = jstart, jfinal, jinc
             jcount = jcount + 1
             y (jcount) = r (ii, jj)
          END DO
          WRITE (ntp, 840) ii, (y (jj), jj = 1, jcount)
       END DO
       jstart = jfinal + jinc
       istart = ibeg
    END DO main_loop
!
!
800 FORMAT (//,T3,'# ', I6.6, TR2, 'ut=', I6.6, TR2,     &
            'ITIME=', I5.5, '(', I2.2, ':', I2.2,':', I2.2, ')' )
810 FORMAT (/,T4, A, '/(1.E', I2, ')', TR4,'ut=', I6.6, TR2, &
            'ave=', ES11.4, TR2, 'sd=', ES11.4)
830 FORMAT ( /, T3, 16(TR5,I3))
840 FORMAT (T2, I3, 16L8)
!
    RETURN
    END SUBROUTINE Outp_logical
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
      FUNCTION Check_logical_units ()
      IMPLICIT NONE
      LOGICAL :: Check_logical_units, L1, L2
        INQUIRE (UNIT = LUN, EXIST = L1, OPENED = L2)
        IF (.NOT.L1) THEN
           Check_logical_units = .FALSE.
           RETURN
        ELSE IF (L2) THEN
           Check_logical_units = .FALSE.
           RETURN
        END IF
        INQUIRE (UNIT = LUN_2, EXIST = L1, OPENED = L2)
        IF (.NOT.L1) THEN
           Check_logical_units = .FALSE.
           RETURN
        ELSE IF (L2) THEN
           Check_logical_units = .FALSE.
           RETURN
        END IF
        INQUIRE (UNIT = LUN_3, EXIST = L1, OPENED = L2)      
        IF (.NOT.L1) THEN
           Check_logical_units = .FALSE.
           RETURN
        ELSE IF (L2) THEN
           Check_logical_units = .FALSE.
           RETURN
        END IF
        Check_logical_units = .TRUE.
      RETURN
      END FUNCTION Check_logical_units
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
!
      SUBROUTINE Read_plasma ()
      INTEGER (iprec) :: n,k
      CHARACTER (LEN=80) :: form_string
!
      OPEN (UNIT = LUN, STATUS = 'OLD', FORM = 'FORMATTED', FILE = rcmdir//'rcmlas1')
!
!
        form_string = '(2(TR2,ES23.15), (TR2,I10.10), (TR2,ES23.15))'
         READ (LUN,'(TR2,I10.10)') n
         IF (n /= kcsize) STOP 'problem with rcmlas1, grid-based'
         READ (LUN,'(A80)') form_string
         DO k = 1, n
            READ (LUN, form_string) alamc(k), etac(k), ikflavc(k), fudgec(k)
         END DO
!
      CLOSE (UNIT = LUN)
!
      RETURN
      END SUBROUTINE Read_plasma
!
!
!
      SUBROUTINE Write_plasma ()
      INTEGER (iprec) :: n,k
      CHARACTER (LEN=80) :: form_string
!
!
      OPEN (LUN, FILE = rcmdir//'rcmlas1', FORM = 'FORMATTED', STATUS='REPLACE')
!
        form_string = '(2(TR2,ES23.15), (TR2,I10.10), (TR2,ES23.15))'
        WRITE (LUN, '(I5.5)') SIZE (alamc)
        WRITE (LUN, '(A80)') form_string
        DO k = 1, SIZE (alamc)
           WRITE (LUN, form_string) alamc(k), etac(k), ikflavc(k), fudgec(k)
        END DO
!
      CLOSE (LUN)
!
      RETURN
      END SUBROUTINE Write_plasma
!
!
!
!
      FUNCTION Get_time_char_string (label) RESULT (time_string)
      IMPLICIT NONE
      TYPE (label_def), INTENT (IN) :: label
      CHARACTER (LEN=8) :: time_string
      WRITE (time_string,'(I2.2,A1,I2.2,A1,I2.2)') &
             label%intg(3),':',label%intg(4),':',label%intg(5)
      RETURN
      END FUNCTION Get_time_char_string
!
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
    SUBROUTINE Write_real_3d_array_old (filename, rec_num, label, array_3d, setup)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)      :: rec_num
    REAL (rprec),      INTENT (IN)      :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)      :: filename
    TYPE (label_def),  INTENT (IN)      :: label
    LOGICAL, OPTIONAL, INTENT (IN)      :: setup
!
    INTEGER (iprec)   :: length, istat
    CHARACTER (LEN=7) :: status_char
    LOGICAL           :: flag
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       STOP 'UNIT ALREADY OPEN'
    END IF
    status_char = 'OLD    '
    IF (PRESENT(setup) ) THEN
       IF (setup) status_char = 'REPLACE'
    END IF
    INQUIRE (IOLENGTH = length ) label, array_3d

    OPEN (UNIT=LUN, FILE = filename, STATUS = status_char, &
          ACCESS = 'DIRECT',  RECL = length, FORM = 'UNFORMATTED',&
          IOSTAT = istat, ACTION = 'WRITE')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF

    WRITE (UNIT = LUN, REC=rec_num, IOSTAT=istat) label, array_3d
    IF (istat /= 0) STOP 'ERROR READING FILE'

    CLOSE (UNIT=LUN, IOSTAT = istat)
    IF (istat /= 0) STOP 'ERROR CLOSING FILE'

    RETURN
    END SUBROUTINE Write_real_3d_array_old
!
    SUBROUTINE Read_real_3d_array_old (filename, rec_num, label, &
                                   array_3d, error_flag, asci)
    IMPLICIT NONE
    INTEGER (iprec),   INTENT (IN)  :: rec_num
    REAL (rprec),      INTENT (OUT) :: array_3d (:,:,:)
    CHARACTER (LEN=*), INTENT (IN)  :: filename
    TYPE (label_def),  INTENT (OUT) :: label
    LOGICAL, OPTIONAL, INTENT (OUT) :: error_flag
    LOGICAL, OPTIONAL, INTENT (IN)  ::  asci
!
    CHARACTER (LEN=11) :: form_type_char
    CHARACTER (LEN=80) :: form_string
    INTEGER (iprec)    :: length, istat
    LOGICAL            :: flag, asci_format
!
!
    asci_format = .FALSE.
    IF (PRESENT(asci)) THEN
       IF (asci) asci_format = .TRUE.
    END IF
!
    INQUIRE (IOLENGTH = length ) label, array_3d
    IF (PRESENT (error_flag)) error_flag = .FALSE.
    form_type_char = 'UNFORMATTED'
    IF (asci_format) THEN
       length = 20*(2+8)+20*(2+23)+(2+80) + SIZE(array_3d)*(2+23)
       form_string = '(20(TR2,I8),20(TR2,ES23.15),(TR2,A80),xxxxXXX(TR2,ES23.15))'
       WRITE (form_string(39:45),'(I7.7)') SIZE(array_3d)
       form_type_char = 'FORMATTED  '
    END IF
!
    INQUIRE (UNIT = LUN, OPENED = flag)
    IF (flag) THEN
       WRITE (*,'(T2,A)') 'UNIT ALREADY OPEN, IN READ_REAL_3D_ARRAY'
       STOP
    END IF
!
    OPEN (UNIT = LUN, FILE = filename, STATUS = 'OLD', ACCESS = 'DIRECT',&
          RECL = length, FORM = form_type_char, IOSTAT = istat, ACTION = 'READ')
    IF (istat /= 0) THEN
       WRITE (*,*) 'ERROR OPENING FILE: ',filename
       STOP
    END IF
!
    IF (asci_format) THEN
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1, FMT = form_string) label, array_3d
       IF (istat /= 0) GO TO 1
    ELSE
       READ (LUN, REC=rec_num, IOSTAT=istat, ERR = 1                   ) label, array_3d
       IF (istat /= 0) GO TO 1
    END IF
    CLOSE (UNIT = LUN)
    RETURN
 1  IF (PRESENT (error_flag)) THEN
       error_flag = .TRUE.
    ELSE
       WRITE (*,*) 'ERROR READING, file is: ', filename
       STOP
    END IF
!
    CLOSE (UNIT=LUN)
!
    RETURN
    END SUBROUTINE Read_real_3d_array_old
!
!
!
!
      SUBROUTINE Rcm (itimei_in, itimef_in, irdr_in, irdw_in, &
                      idt_in, idt1_in, idt2_in, icontrol,stropt,nslcopt)
      USE rcm_timing_module
      IMPLICIT NONE
!
      INTEGER (iprec), INTENT (IN) :: itimei_in, itimef_in, &
                                      irdr_in, irdw_in, idt_in, idt1_in,& 
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
      INTEGER (iprec), SAVE :: itimei, itimef, idt, idt1, idt2, irdw, irdr
      INTEGER (iprec), SAVE :: itout1, itout2,  itcln, idebug, i_time, &
                         k, kc, n, i_avg
      REAL (rprec) :: dt


      CALL SYSTEM_CLOCK (timer_start(1), count_rate)


      IF (IsCoupledExternally) then
       itimef = itimef_in
       itimei = itimei_in
       irdr   = irdr_in
       irdw   = irdw_in
       idt    = idt_in
       idt1   = idt1_in
       idt2   = idt2_in
      END IF

      IF (icontrol == 31337) then  ! write a restart record to RCM
         call WriteRCMH5(stropt,nslcopt,isRestart=.true.)
         return
      ENDIF
      IF (icontrol == 31338) then  ! write an HDF5 output
        call WriteRCMH5(stropt,nslcopt,isRestart=.false.)
        return
      ENDIF
      IF (icontrol == 31336) then
        !Read HDF5 restart
        call ReadRCMRestart(stropt,nslcopt)
        return
      ENDIF
   IF (icontrol == 0) then  ! initialize RCM size params and go back:

      IF (.NOT.Check_logical_units ( ) ) STOP 'LUNs NOT AVAILABLE'


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
      CALL Read_plasma ()

      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN
   END IF



   IF (icontrol == 2) then ! read in inputs, quit:
      !K: Replace RCM params file with XML data
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
       READ (LUN,*) irdr     ! 3.  record # to read in
       READ (LUN,*) irdw     ! 4.  record # to write out
       READ (LUN,*) idt   !  1.  basic time step in program
       READ (LUN,*) idt1  !  2.  t-step for changing disk write rcds
       READ (LUN,*) idt2  !  3.  t-step for writing formatted output
      END IF

      CLOSE (UNIT = LUN)
      

      !  Read in other inputs, both constant and run-specific:

      IF (.NOT.IsCoupledExternally) CALL Read_vdrop  
      CALL Read_kp   
      CALL Read_bfield
      IF (.NOT.IsCoupledExternally) CALL Read_eta_on_bndy
      IF (.NOT.IsCoupledExternally) CALL Read_qtcond
      IF (.NOT.IsCoupledExternally) CALL Read_winds (iwind)
      CALL Read_dktime (L_dktime)


      CALL SYSTEM_CLOCK (timer_stop(1), count_rate)      
      timer_values (1) = (timer_stop (1) - timer_start (1))/count_rate + timer_values(1)

      RETURN

   END IF


   IF (icontrol == 3) then   !-->  Set initial conditions on plasma (edges and grid-based):
!
      ! Open file for formatted output and do initial print out :
      CALL Date_and_time (real_date, real_time)
!     IF (itimei == 0) THEN
      IF (irdr == 1) THEN
         ST = 'REPLACE'
         PS = 'APPEND'
         HD = 'BEGINNING NEW RUN'
      ELSE
         ST = 'OLD'
         PS = 'APPEND'
         HD = 'CONTINUE SAME RUN'
      END IF
      write(*,*) "L9604, rcm.printout", LUN_2, ST, PS
      OPEN  (LUN_2, FILE = rcmdir//'rcm.printout', STATUS = ST, POSITION = PS)
      OPEN  (LUN_3, FILE = rcmdir//'rcm.index',  STATUS = ST, POSITION = PS)
      CALL Initial_printout ()
      CLOSE (LUN_3)
      CLOSE (LUN_2)

      i1 = imin + 1


      IF (itimei == 0) THEN

         IF (.NOT.IsCoupledExternally) then
            CALL Read_array (rcmdir//'rcmeeta_inp',   irdr, label, ARRAY_3D = eeta,ASCI=asci_flag)
         ELSE
            ! grid-based plasma should have been set up elsewhere and passed via module, do nothing here:
            !  CALL Read_array (rcmdir//'rcmeeta_inp',   irdr, label, ARRAY_3D = eeta,ASCI=asci_flag)
         END IF

      ELSE

         CALL Read_array (rcmdir//'rcmbndloc', irdr, label, ARRAY_1D = bndloc)
         imin_j = CEILING (bndloc)
         CALL Read_array (rcmdir//'rcmeeta',   irdr, label, ARRAY_3D = eeta)
!         IF (label%intg(6) /= itimei-idt )THEN
         IF (label%intg(6) /= itimei )THEN
!            WRITE (*,*)' label%intg(6) =',label%intg(6),' itimei-idt =',itimei-idt
            WRITE (*,*)' label%intg(6) =',label%intg(6),' itimei =',itimei
            STOP 'T in file /=  ITIMEI for EETA, RESTART IS CORRUPTED'
         END IF
      END IF





      ! IF hot restart, read V and check the time label:

      IF (.NOT.IsCoupledExternally) THEN
!        IF (itimei /= 0) THEN
         IF (irdr /= 1) THEN
            WRITE (*,'(A)', ADVANCE='NO') &
              'HOT restart, reading V from file to check time label...'
            CALL READ_array (rcmdir//'rcmv', irdr, label, ARRAY_2D = v)
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

            call AddToList(i_time,rcm_timing)

            itout1 = MIN (itout1 + idt1, itimef-idt)
            irdw   = irdw + 1
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
      write(*,*) "L9786, rcm.printout",irdr,itimei
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
      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'read at itimei from REC ',irdr
      WRITE (LUN_2,'(T5,A,T35,I5.5)') 'start writing at REC ', irdw 
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
      WRITE (LUN_2,'(/T5,A)' ) 'PCP INPUT PARAMETERS:'
         WRITE (LUN_2,'(T5,A,T20,I5)') 'nvmax =', SIZE(ivtime)
         DO n = 1, SIZE(ivtime)
            WRITE (LUN_2,'(T5,A,I5,T20,A,G9.2)') &
                  'ivtime=',ivtime(n), 'vinput=',vinput(n)
         END DO
      WRITE (LUN_2,'(/T5,A)' ) 'BFIELD MARKTIMES:'
         WRITE (LUN_2,'(T5,A,T20,I6)') 'nbf =', SIZE(ibtime)
         DO n = 1, SIZE(ibtime)
            WRITE (LUN_2,'(T5,A,I5)') 'ibtime=',ibtime(n)
         END DO
      WRITE (LUN_2,'(/T5,A)' ) 'KP MARKTIMES:'
         WRITE (LUN_2,'(T5,A,T20,I6)') 'nkpmax =', SIZE(ikptime)
         DO n = 1, SIZE(ikptime)
            WRITE (LUN_2,'(T5,A,I9)') 'ikptime=',ikptime(n)
         END DO
      WRITE (LUN_2,'(T2,A,T20,I2)')     'IPOT =',     ipot
      WRITE (LUN_2,'(T2,A,T20,I2)')     'ICOND = ',   icond
      WRITE (LUN_2,'(T2,A,T20,I2)')     'IBND = ',    ibnd_type
      WRITE (LUN_2,'(T2,A,T20,I2)')     'IPCP_TYPE=', ipcp_type


      WRITE (LUN_2,'(//A)') 'Begin RCM timing table'
      WRITE (LUN_2,'(T1,A,T16,A,T31,A,T46,A)')  'RCM itime [s]', '  record#', 'cpu_time [s]', 'sum_cpu_time [s]'

  902 FORMAT (T2,'TIME', T12,'ITIME' , T19,'REC#' ,&
              T26,'VDROP', T33,'KP', T39,'FSTOFF',   &
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
         label%intg (1) = irdw  ! record # for disk printout
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
         label%real (17) = kp
!
         ST = 'OLD'
         PS = 'APPEND'
         OPEN  (LUN_3, FILE = rcmdir//'rcm.index',  STATUS = ST, POSITION = PS)
         WRITE (time_char,'(I2.2,A1,I2.2,A1,I2.2)') &
               label%intg(3), ':', label%intg(4), ':', label%intg(5)
         WRITE (*,'(T2,A21,I5.5,A10,TR4)') &
                '-->TIME_STEP, T=', i_time,'('//time_char//')'
!                                                                       
!        IF (i_time == itout1 .OR. i_time == itimef) THEN
            WRITE (lun_3,901) time_char, &
                     i_time, irdw, vdrop, kp, fstoff, fmeb, fdst, fclps, vdrop_phase
            CLOSE (LUN_3)
  901       FORMAT (T2,A8, T12,I6, T19,I5, T26,F5.1, T33,F5.2,&
                    T39,F5.2, T46,F5.2, T53,F7.1, T62,F5.1, T69, F6.2)
!        END IF
!
         IF (idebug == 0) THEN 
!
!
!        1. Write magnetic field and bndy location:
!
            CALL Write_array (rcmdir//'rcmxmin', irdw, label, ARRAY_2D = xmin, SETUP = FD)
            CALL Write_array (rcmdir//'rcmymin', irdw, label, ARRAY_2D = ymin, SETUP = FD)
            CALL Write_array (rcmdir//'rcmzmin', irdw, label, ARRAY_2D = zmin, SETUP = FD)
            CALL Write_array (rcmdir//'rcmvm',   irdw, label, ARRAY_2D = vm,   SETUP = FD)
            CALL Write_array (rcmdir//'rcmbmin', irdw, label, ARRAY_2D = bmin, SETUP = FD)
            write(*,*)'writing rcmbndloc at rec=',irdw
            CALL Write_array (rcmdir//'rcmbndloc', irdw, label, ARRAY_1D = bndloc, SETUP = FD)
!
!
!         2. Write out plasma info:
!
            CALL Write_array (rcmdir//'rcmetac', irdw, label, ARRAY_1D = etac, SETUP=FD)
            CALL Write_array (rcmdir//'rcmeeta',   irdw, label, ARRAY_3D = eeta, SETUP=FD )
            CALL Write_array (rcmdir//'rcmeetaavg',irdw,label,ARRAY_3D=eeta_avg,SETUP=FD )
!
!
!         3. Write ionospheric quantities:
!
            CALL Write_array (rcmdir//'rcmpedlam', irdw, label, ARRAY_2D = pedlam, SETUP=FD )
            CALL Write_array (rcmdir//'rcmpedpsi', irdw, label, ARRAY_2D = pedpsi, SETUP=FD )
            CALL Write_array (rcmdir//'rcmhall',   irdw, label, ARRAY_2D = hall  , SETUP=FD )
            CALL Write_array (rcmdir//'rcmeavg',   irdw, label, ARRAY_3D = eavg  , SETUP=FD )
            CALL Write_array (rcmdir//'rcmeflux',  irdw, label, ARRAY_3D = eflux , SETUP=FD )
            CALL Write_array (rcmdir//'rcmbirk',  irdw, label, ARRAY_2D = birk , SETUP=FD )
            CALL Write_array (rcmdir//'rcmbirkavg',irdw,label,ARRAY_2D=birk_avg,SETUP=FD )
!
!
!           4. Write out magnetospheric quantities:
!
            CALL Write_array (rcmdir//'rcmv',    irdw, label, ARRAY_2D = v, SETUP=FD )
            CALL Write_array (rcmdir//'rcmvavg', irdw, label, ARRAY_2D = v_avg, SETUP=FD )
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
          type(IOVAR_T), dimension(50) :: IOVars !Lazy hard-coding max variables
          integer(iprec) :: nvar

        !Prepare for reading
          doSP = .false. !Restarts are always double precision
          write (H5File, '(A,A,I0.5,A)') trim(runid), ".RCM.Res.", nStp, ".h5"
          call ClearIO(IOVars) !Reset IO chain

        !List variables to read
          !Scalars
          !(Need to specify integers)
          call AddInVar(IOVars,"irdw",vTypeO=IOINT)
          call AddInVar(IOVars,"itimei",vTypeO=IOINT)
          call AddInVar(IOVars,"isize",vTypeO=IOINT)
          call AddInVar(IOVars,"jsize",vTypeO=IOINT)
          !(Add rest of scalars here)

          !Arrays
          call AddInVar(IOVars,"rcmetac")
          call AddInVar(IOVars,"rcmeeta")
          call AddInVar(IOVars,"rcmeetaavg")
          !(Add rest of arrays here)

        !Now do actual reading
          call ReadVars(IOVars,doSP,H5File)

        !Parse data and put it where it goes
          !Data is stored as 1D array, for multi-d arrays need to be reshaped

          !etac (1D)
          nvar = FindIO(IOVars,"rcmetac")
          etac = IOVars(nvar)%data

          !eeta (3D), need to reshape 1D data into 3D array
          nvar = FindIO(IOVars,"rcmeeta")
          eeta = reshape(IOVars(nvar)%data,[IOVars(nvar)%dims(1),IOVars(nvar)%dims(2),IOVars(nvar)%dims(3)])

          !rcmv (2D)
          nvar = FindIO(IOVars,"rcmv")
          v = reshape(IOVars(nvar)%data,[IOVars(nvar)%dims(1),IOVars(nvar)%dims(2)])

          !(Add rest of arrays here)

          write(*,*) 'RCM restart not finished ...'
        end subroutine ReadRCMRestart

        !HDF5 output routine 
        subroutine WriteRCMH5(runid,nStp,isRestart)
          use ioh5
          use files
          implicit none
          character(len=*), intent(in) :: runid
          integer(iprec), intent(in) :: nStp
          logical, intent(in) :: isRestart

          type(IOVAR_T), dimension(50) :: IOVars !Lazy hard-coding max variables
          logical :: doSP !Do single precision output
          character(len=strLen) :: H5File,gStr,lnResF

        !Prepare for output
          !Reset IO chain
          call ClearIO(IOVars)
          !Distinguish output slices vs restarts
          if (isRestart) then
            doSP = .false. !Double precision restarts
            write (H5File, '(A,A,I0.5,A)') trim(runid), ".RCM.Res.", nStp, ".h5"
            !write(*,*) 'RCM: Writing restart to ', trim(H5File)
          else
            !Regular output
            doSP = .true.
            H5File = trim(runid) // ".rcm.h5"
            write (gStr, '(A,I0)') "Step#", nStp
            !write(*,*) 'RCM: Writing output to ', trim(H5File), '/',trim(gStr)
          endif


        !Attributes
          call AddOutVar(IOVars,"time",1.0_rp*itimei)
          call AddOutVar(IOVars,"itimei",itimei)
          call AddOutVar(IOVars,"irdw"  ,irdw  )
          call AddOutVar(IOVars,"isize" ,isize )
          call AddOutVar(IOVars,"jsize" ,jsize )
          call AddOutVar(IOVars,"ksize" ,ksize )
          call AddOutVar(IOVars,"cmax"  ,cmax  )
          call AddOutVar(IOVars,"fmeb"  ,fmeb  )
          call AddOutVar(IOVars,"fstoff",fstoff)
          call AddOutVar(IOVars,"fdst"  ,fdst  )
          call AddOutVar(IOVars,"fclps" ,fclps )
          call AddOutVar(IOVars,"vdrop" ,vdrop )
          call AddOutVar(IOVars,"kp"    ,kp    )

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

        !Done staging output, now let er rip
          if (isRestart) then
            call CheckAndKill(H5File) !Always overwrite restarts
            call WriteVars(IOVars,doSP,H5File)
            !Create link to latest restart
            write (lnResF, '(A,A,A,A)') trim(runid), ".RCM.Res.", "XXXXX", ".h5"
            call EXECUTE_COMMAND_LINE('ln -sf '//trim(H5File)//' '//trim(lnResF), wait=.false.)
          else
            call WriteVars(IOVars,doSP,H5File,gStr)
          endif
        end subroutine WriteRCMH5

        subroutine RCM_Params_XML()
          use xml_input
          use strings

          character(len=strLen) :: inpXML
          type(XML_Input_T) :: xmlInp

          !Find input deck filename
          call getIDeckStr(inpXML)

          !Create XML reader
          xmlInp = New_XML_Input(trim(inpXML),'RCM',.true.)

        end subroutine RCM_Params_XML

         SUBROUTINE Formatted_output ()
      write(*,*) "L10019, rcm.printout", LUN_2, ST, PS
         OPEN  (LUN_2, FILE = rcmdir//'rcm.printout', STATUS = 'OLD', POSITION = 'append')
         WRITE (LUN_2,'(T1,I10,T16,I10,T31,F10.2,T46,F10.2)')  &
        &  i_time, irdw, timer_values(2), timer_values(1) 
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
!
      SUBROUTINE V_eff_polar_cap (veff)
      IMPLICIT NONE
      REAL(rprec), INTENT (IN OUT) :: veff (:,:)
!
!---------------------------------------------------------------------
!     Subroutine to define V_eff in the polar cap (beyound the RCM
!     modeling region) such that it matches V_eff and and its 
!     derivative on the boundary and becomes uniform dawn-dusk Efield
!     in the polar cap. Experssion from Dick Wolf, Aug. 2001.     
!---------------------------------------------------------------------
!
      INTEGER(iprec) :: i_b, j_18, j_06, i, j, j_24, jm
      REAL(rprec) :: delta_i = 3.0, veff_b, dveff_di_b, d_pcp, pcp
!
      j_18 = (jsize-jwrap)/4 +  jwrap
      j_06 = j_18 + (jsize-jwrap)/2
      j_24 = jwrap + (jsize-jwrap)/2
!print*,'veff:',j_18, j_06
!pause
      d_pcp = 0.0
      pcp   = 0.0
      DO j = jwrap+1, j_24-1
         jm = jsize - (j-jwrap)
         IF (colat(imin_j(j),j) + colat(imin_j(jm),jm) > d_pcp) &
             d_pcp = colat(imin_j(j),j) + colat(imin_j(jm),jm)
         IF (ABS(v(imin_j(j),j) - v(imin_j(jm),jm)) > pcp) &
             pcp = ABS(v(imin_j(j),j)-v(imin_j(jm),jm))
! print*,'finding pcp:', j, jm, d_pcp, pcp
      END DO
      d_pcp = 2.0*colat(imin_j(j_18),j_18)
!     pcp   = vdrop *1.0E+3
!print*,'veff: pcp=',pcp,d_pcp
!pause
!
      DO  j = 1, jsize
         i_b = imin_j (j)
         veff_b = veff (i_b,j)
         dveff_di_b = (veff(i_b+1,j)-veff(i_b,j))
         DO i = 1, i_b - 1
            veff(i,j) = (veff_b - dveff_di_b*(i_b-i)) / &
                         (1.0+(i_b-i)**2/delta_i**2) - &
                        ((i_b-i)/delta_i)**2 / (1.0+((i_b-i)/delta_i)**2) * &
                        pcp/d_pcp*sin(colat(i,j))*sin(aloct(i,j))
         END DO
      END DO
      RETURN
      END SUBROUTINE V_eff_polar_cap
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
  REAL (rprec) :: eps=1.0e-8, max_eeta
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
!$OMP PARALLEL PRIVATE (eeta2, veff, dvefdi, dvefdj, didt, djdt, &
!$OMP                  & mass_factor, loc_didt, loc_djdt,loc_Eta,loc_rate, &
!$OMP                  & ie, icut, j, i, r_dist, FirstTime, max_eeta) &
!$OMP        & SHARED (alamc, eeta, v, vcorot, vpar, vm, imin_j, j1, j2, joff, &
!$OMP                  xmin, ymin, fac, fudgec, bir, sini, L_dktime, dktime, sunspot_number, &                
!$OMP                  T1, T2, xlower, ylower, xupper, yupper, CLAWiter, eps) &
!$OMP        & DEFAULT (NONE)
!$OMP DO SCHEDULE (dynamic)

  DO kc = 1, kcsize
!
!    If oxygen is to be added, must change this!
!
     IF (alamc(kc) < 0.0) THEN
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
!$OMP END DO
!$OMP END PARALLEL


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
END SUBROUTINE Move_plasma_grid_NEW
 

  SUBROUTINE Deriv_i_NEW (array, isize, jsize, j1, j2, imin_j, derivi)
!   USE Rcm_mod_subs, ONLY : iprec, rprec
    IMPLICIT NONE
    INTEGER (iprec), INTENT (IN) :: isize, jsize, j1, j2, imin_j(jsize)
    REAL (rprec), INTENT (IN) :: array (isize,jsize)
    REAL (rprec), INTENT (OUT) :: Derivi (isize,jsize)
!
    INTEGER (iprec) :: i, j, idim, jdim
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
  
!
!
    END MODULE Rcm_mod_subs
