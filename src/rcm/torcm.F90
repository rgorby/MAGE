
MODULE torcm_mod
  USE rcm_precision
  USE constants, only: big_vm,tiote,nt,ev,boltz
  USE rice_housekeeping_module, ONLY: use_plasmasphere,LowLatMHD,L_write_vars_debug
  USE rcm_mhd_interfaces
  USE rcmdefs, ONLY : RCMTOPCLOSED,RCMTOPNULL,RCMTOPOPEN,DenPP0
  USE kdefs, ONLY : TINY,qp
  USE math, ONLY : RampDown
  USE etautils
  implicit none

  logical, parameter :: doSmoothEta = .true. !Whether to smooth eeta at boundary
  integer(iprec), private, parameter :: NumG = 4 !How many buffer cells to require

  contains
!==================================================================      
      SUBROUTINE Torcm (RM, itimei, ierr, icontrol) 

      USE Rcm_mod_subs, ONLY : isize,jsize, jwrap, kcsize, iesize, &
                               vm, bmin, xmin, ymin, pmin, rmin,v, & 
                               alamc, etac, ikflavc, fudgec, eeta, &
                               imin_j, bndloc, vbnd,               &
                               colat, aloct, bir, sini,            &
                               ibnd_type,rcmdir
      USE conversion_module
      USE earthhelper, ONLY : GallagherXY
      
! NOTE: This version fixes the rcm boundary condition at rec=1
!===================================================================
!
! purpose: to setup rcm for a run, converts mhd information to 
!          rcm information
!
! inputs:
!     rec_in   = record to be updated with new values
!               rec == 1: resets the rcm and initializes the files
!                rec > 1: only boundary conditions are updated
!    itimei = time to write to rcm records
!
! outputs:
!       ierr = error flag if there is a problem
!
!===================================================================
!
!  last update:  05.18.90
!                04.13.95
!                      96
!                20.01.00 frt mc version
!                    8.03 frt computes integral average pressure
!                06.01.04 frt cleaned up version
!                09.10.12 frt added code to handle tilt
!                19.05.20 frt removed use of records in rcm 
!                         frt remove use of idim,jdim,kdim
!
!===================================================================

      IMPLICIT NONE
      type(rcm_mhd_T),intent(inout) :: RM
      REAL(rprec), INTENT (IN) :: itimei
      INTEGER(iprec), INTENT (IN) :: icontrol
      INTEGER(iprec), INTENT (IN OUT) :: ierr

      INTEGER(iprec) :: kin,jmid,ibnd, inew, iold
      INTEGER(iprec) :: jm,jp,ip,itmax
      INTEGER(iprec) :: min0,i,j,k,n,ns,klow,n_smooth
      LOGICAL,PARAMETER :: use_ellipse = .true.
      LOGICAL, SAVE :: doReadALAM = .true.
      REAL(rprec) :: dpp,wMHD,wRCM

      !K: 8/20, rewritten to try to better incorporate immersed boundary BCs
      !iopen: -1 (RCMTOPCLOSED), CLOSED & inside RCM ellipse
      !       +1 (RCMTOPOPEN)  , OPEN & definitely outside RCM domain
      !        0 (RCMTOPNULL)  , CLOSED & outside RCM domain

      ierr = 0

      !Set density/pressure factors
      call SetFactors(RM%planet_radius)

      !Start by reading alam channels if they're not yet set
      !Rewriting this bit to not read_alam every call, K: 8/20
      if (doReadALAM) then
        CALL Read_alam (kcsize, alamc, ikflavc, fudgec, almdel, almmax, almmin, iesize, ierr)
        doReadALAM = .false.
        IF (ierr < 0) RETURN
        !Go ahead and do some other init stuff while you're here
        LowLatMHD = RM%llBC
      endif

      !Set lowest RC channel
      if (use_plasmasphere) then
        klow = 2
      else
        klow = 1
      endif

      !If T>0, save certain arrays before modifying anything
      IF (icontrol==RCMADVANCE.or.icontrol==RCMRESTART) then  
        bndloc_old = bndloc
        imin_j_old = imin_j
        if(minval(bndloc) <= 0)then
          write(6,*) ' TORCM: boundary problem'
          write(6,*)' bndloc:',bndloc
          ierr = -1
          return
          
        end if
      END IF  

      ! get B-field lines starting from RCM ionospheric grid
      ! points, compute flux-tube volume of field lines that are
      ! closed, and mark open ones with a mask array OPEN
      ! (1=open, -1=closed). For open field lines, FTV is set to big_vm:

      CALL Calc_ftv (RM,big_vm,ierr)
      IF (ierr < 0) THEN
        write(6,*) 'calc_ftv problem'
        RETURN
      ENDIF

    !-----
    !Domain calculation

      ! Now, set RCM high-latitude grid boundary. Initially,
      ! for each MLT, we find the grid point with highest
      ! latitude and still on closed field lines; the point
      ! next to it (down in latitude) is set to be boundary.
      ! Result of this subsection is to populate arrays BNDLOC
      ! and IMIN_J with values:
      do j=1,jsize
        bndloc(j) = 2  ! if everything else fails, can use this...
        do i=isize,2,-1
          if (iopen(i,j) >= 0) then !null or open
            bndloc(j) = i + 2 !Adding buffer cell/s here
            exit
          endif
        end do
      end do
      ! reset imin_j
      imin_j = ceiling(bndloc)
      IF (L_write_vars_debug) then
        write(6,*)' bndy ',bndloc(j),j,vm(imin_j(j),j)
      END IF

      ! fits an ellipse to the boundary
      !K: 8/20, changing ellipse so that it doesn't reset vm unless open
      if (use_ellipse) then
        CALL Set_ellipse(isize,jsize,rmin,pmin,vm,big_vm,bndloc,iopen)
      end if

      call reset_rcm_vm(isize,jsize,bndloc,big_vm,imin_j,vm,iopen)

      !Smooth boundary location
      !NOTE: This can only shrink RCM domain, not increase
      !K: 8/20, changing reset_rcm_vm to only reset VM on open fields
      !Calculate number of smoothing iterations based on longitudinal cell size
      !Approx. 1hr MLT
      !n_smooth = 5
      n_smooth = nint( 15.0/(360.0/jsize) )

      do ns=1,n_smooth
        call smooth_boundary_location(isize,jsize,jwrap,bndloc)
        call reset_rcm_vm(isize,jsize,bndloc,big_vm,imin_j,vm,iopen) ! adjust Imin_j
      enddo
      

    !-----
    !MHD thermodynamics
      IF (icontrol==RCMCOLDSTART .and. use_plasmasphere) THEN
        eeta(:,:,1) = 0.0 !Make sure no plasmasphere for first cold start calculation
      ENDIF

      !---->Set new EETA from MHD code pressure. 
      !     On open field lines, values of ETA will be zero:
      CALL Gettemp (ierr)
      IF (ierr < 0) RETURN
      CALL Press2eta()       ! this populates EETA_NEW array
      
      if ( (maxval(eeta_new) <=0) .or. any(isnan(eeta_new)) ) then
        write(6,*)' something is wrong in eeta_new'
        write(6,*) 'maxval = ', maxval(eeta_new)
        write(6,*) 'Num nans = ', count(isnan(eeta_new))
        ierr = -1
        return
      end if

      !Handle cold start, must happen after eeta_new is calculated (temp/press2eta)
      IF (icontrol==RCMCOLDSTART) THEN
        write(6,*)' TORCM: initializing the RCM arrays at t=',itimei
        bndloc_old = bndloc
        imin_j_old = imin_j
        !Now set RCM domain values to MHD state values
        do j=1,jsize
          do i=imin_j(j),isize
            !These will all be closed cells
            !Don't worry about plasmasphere, that channel will get reset anyways
            eeta(i,j,:) = eeta_new(i,j,:)
          enddo
        enddo

      ENDIF

      ! just in case:
      imin_j     = CEILING(bndloc)
      imin_j_old = CEILING(bndloc_old)

      ! initialize the dynamic plasmasphere, reset static part if necessary
      !sbao 03282020
      if (use_plasmasphere) then
        call set_plasmasphere(icontrol,isize,jsize,kcsize,xmin,ymin,rmin,vm,eeta,imin_j)
        
        !Adding some smoothing to the plasmapause
        call SmoothPPause(eeta(:,:,1),vm,iopen,imin_j)
      endif


      !Incorporate MHD-produced values into newly-acquired RCM cells
      DO j=1,jsize
        inew = imin_j(j)
        iold = imin_j_old(j)
        if (inew <= iold) then
          eeta(inew:iold,j,klow:) = eeta_new(inew:iold,j,klow:)
        endif
      ENDDO


    !-----
    !Fill in grid ghosts
      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j,dpp)
      do j=1,jsize
        do i=1,isize
          if (iopen(i,j) == RCMTOPOPEN) then
            !Zap everything here
            eeta(i,j,:) = 0.0
            vm  (i,j)   = big_vm
            !Reset mapping on open lines
            rmin(i,j) = 0.0
            pmin(i,j) = 0.0
            xmin(i,j) = 0.0
            ymin(i,j) = 0.0

          else if (iopen(i,j) == RCMTOPNULL) then
            !This is closed field region but outside RCM domain
            !Use MHD information for RC channels
            eeta(i,j,klow:) = eeta_new(i,j,klow:)
            
          endif !iopen

          if (use_plasmasphere .and. iopen(i,j) /= RCMTOPOPEN) then
          	!Check for below plasmasphere cutoff
          	dpp = GetDensityFactor()*1.0*eeta(i,j,1)*vm(i,j)**1.5
          	if (dpp < DenPP0*1.0e+6) eeta(i,j,1) = 0.0
          endif !plasmasphere

        enddo !i
      enddo !j

      !Do some reverse-blending near RCM outer boundary
      n = NumG
      do j=1,jsize
        do i=0,n
          ip = imin_j(j)+i !i cell
          if (ip<isize) then
            wRCM = 1.0/(1.0 + 5.0*beta_average(ip,j)/6.0)
            wMHD = (1-wRCM)/(2.0**i)
            wRCM = (1-wMHD)

            eeta(ip,j,klow:) = wRCM*eeta(ip,j,klow:) + wMHD*eeta_new(ip,j,klow:)
          endif
        enddo !i loop
      enddo !j loop

      if (doSmoothEta) then
        ! smooth eeta at the boundary
        CALL Smooth_eta_at_boundary(isize,jsize,kcsize,jwrap,eeta,iopen,imin_j)
      endif

    !-----
    !Finish up and get out of here

      ! import ionosphere
      call Ionosphere_toRCM(RM)

      !Do sanity check
      DO j = 1, jsize
        DO i = imin_j(j),isize
          IF (vm(i,j) <= 0.0) THEN
            write(6,*) 'vm problem in TORCM'
            ierr = -1
            return
          ENDIF
        END DO
      END DO

      !Return updates to topology, CLOSED=>NULL in buffer region to RM object
      RM%iopen   = iopen   (:,jwrap:jsize)

      if (any(isnan(eeta))) then
        write(6,*) 'Bad eeta at end of torcm!'
        ierr = -1
        return
      endif

      RETURN
      END SUBROUTINE Torcm

      !Smooth ragged edges at the plasmapause
!----------------------------------------------------------
      SUBROUTINE SmoothPPause(etapp,vm,iopen,imin_j)
        USE rcmdefs,ONLY : isize,jsize
        IMPLICIT NONE
        REAL(rprec)   , INTENT(INOUT) :: etapp(isize,jsize)
        REAL(rprec)   , INTENT(IN)    :: vm(isize,jsize)
        INTEGER(iprec), INTENT(IN)    :: iopen(isize,jsize)
        INTEGER(iprec), INTENT(IN)    :: imin_j(jsize)

        REAL(rprec) :: bndlocpp(jsize),bndlocpp_new(jsize)
        INTEGER(iprec) :: imin_jpp(jsize)
        INTEGER(iprec) :: i,j,n,n_smooth,di
        REAL(rprec) :: dpp,Ac(3)
        INTEGER(iprec), parameter :: NumI = NumG

        !Find plasmapause
        do j=1,jsize
          !For fixed j start from far and go near, find first value above DenPP0
          do i = 2,isize
            dpp = GetDensityFactor()*1.0*etapp(i,j)*vm(i,j)**1.5
            if (dpp >= DenPP0*1.0e+6) then
              bndlocpp(j) = i-1
              imin_jpp(j) = i-1
              exit
            endif !Above plasmapause cutoff
          enddo !i loop
        enddo !j loop

        !Now do smoothing on real-valued boundary
        n_smooth = nint( 5.0/(360.0/jsize) ) !Arbitrary 5deg window

        do n=1,n_smooth
          call SmoothJBnd(bndlocpp,bndlocpp_new)
          !Store back w/ trap for only earthward
          bndlocpp = max(bndlocpp,bndlocpp_new)
        enddo

        !Convert to indices
        imin_jpp = ceiling(bndlocpp)
        
        !Smooth around estimated plasmapause
        Ac = [0.25,0.5,1.0] !Smoothing coefficients
        do di = -NumI/2,+NumI/2
          call SmoothJEta(etapp,iopen,imin_jpp,di,Ac)
        enddo

        ! !Zero out below plasmapause
        ! do j=1,jsize
        !   etapp(1:imin_jpp(j),j) = 0.0 !Zero out below this
        ! enddo

      END SUBROUTINE SmoothPPause

      !Smooth real-valued bndloc, store in bndloc_new
      SUBROUTINE SmoothJBnd(bndloc,bndloc_new)
        USE rcmdefs,ONLY : jsize
        IMPLICIT NONE
        REAL(rprec), INTENT(IN ) :: bndloc(jsize)
        REAL(rprec), INTENT(OUT) :: bndloc_new(jsize)

        REAL(rprec), PARAMETER :: am = 1, a0 = 2, ap = 1
        INTEGER(iprec) :: j,jm,jp

        bndloc_new = 0.0

        do j=1,jsize
          ! 1 <=> jdim -jwrap +1
          ! jdim <=> jwrap
      
          jm  = j - 1
          if(jm < 1) jm = jsize - jwrap 
          jp  = j + 1
          if(jp > jsize) jp = jwrap + 1
          !Use 3-pt stencil
          bndloc_new(j) = (am*bndloc(jm) + a0*bndloc(j) + ap*bndloc(jp))/(am+a0+ap)

        enddo

      END SUBROUTINE SmoothJBnd

      !Smooth eta k-slice over ibnd_j + di using coefficients Ac (--,-,0)
      SUBROUTINE SmoothJEta(etak,iopen,ibnd_j,di,Ac)
        USE rcmdefs,ONLY : isize,jsize,jwrap
        IMPLICIT NONE

        REAL(rprec), INTENT(INOUT) :: etak(isize,jsize)
        INTEGER(iprec), INTENT(IN) :: iopen(isize,jsize)
        INTEGER(iprec), INTENT(IN) :: ibnd_j(jsize)
        INTEGER(iprec), INTENT(IN) :: di
        REAL(rprec), INTENT(IN)    :: Ac(3)

        logical :: isOpen(5)
        REAL(rprec) :: etaks(jsize)
        REAL(rprec) :: a1,a2,a3,a4,a5
        INTEGER(iprec) :: jmm,jm,j,jp,jpp
        a1 = Ac(1) ; a2 = Ac(2) ; a3 = Ac(3) ; a4 = a2 ; a5 = a1

        do j=1,jsize
        !Get indices
          ! 1 <=> jdim -jwrap +1
          ! jdim <=> jwrap
          jmm = j - 2
          if(jmm < 1)jmm = jsize - jwrap - 1
          jm  = j - 1
          if(jm < 1) jm = jsize - jwrap
          jpp = j + 2
          if(jpp > jsize)jpp = jwrap + 2
          jp  = j + 1
          if(jp > jsize) jp = jwrap + 1

          !Check topology
          isOpen(1) = (iopen(ibnd_j(jmm)+di,jmm) == RCMTOPOPEN)
          isOpen(2) = (iopen(ibnd_j(jm )+di,jm ) == RCMTOPOPEN)
          isOpen(3) = (iopen(ibnd_j(j  )+di,j  ) == RCMTOPOPEN)
          isOpen(4) = (iopen(ibnd_j(jp )+di,jp ) == RCMTOPOPEN)
          isOpen(5) = (iopen(ibnd_j(jpp)+di,jpp) == RCMTOPOPEN)

          if ( any(isOpen) ) then
            !Keep old values b/c too close to OCB
            etaks(j) = etak(ibnd_j(j)+di,j)
          else
            !Only smooth if all closed/null cells
            etaks(j) = ( a1*etak(ibnd_j(jmm)+di,jmm) + &
                         a2*etak(ibnd_j(jm )+di,jm ) + &
                         a3*etak(ibnd_j(j  )+di,j  ) + &
                         a4*etak(ibnd_j(jp )+di,jp ) + &
                         a5*etak(ibnd_j(jpp)+di,jpp) )/(a1+a2+a3+a4+a5)
          endif !isOpen
        enddo !j loop

        !Now go back and reset values
        do j=1,jsize
          etak(ibnd_j(j)+di,j) = etaks(j)
        enddo

      END SUBROUTINE SmoothJEta
!----------------------------------------------------------
      SUBROUTINE Calc_ftv (RM,big_vm,ierr) 
      USE conversion_module
      USE RCM_mod_subs,ONLY : isize,jsize,kcsize,bmin,vm,rmin,pmin,&
                              xmin,ymin,zmin,vbnd,jwrap,radcurv,losscone
      
      IMPLICIT NONE
      type(rcm_mhd_T), intent(in) :: RM
      REAL(rprec), INTENT (IN) :: big_vm
      INTEGER(iprec), INTENT (OUT) :: ierr
!
!===================================================================
!
! calc_ftv: gets flux tube volumes, mapping parameters (r,p)
!            be, sini, open/closed flag, and pressure on each 2D
!
! inputs: 
!      isize = dimension in latitude
!      jsize = dimension in local time
!     kcsize = channel dimensions (not used here)
!  x0_sm,y0_sm,z0_sm = location of grid points in ionosphere in sm coordinates
!   big_vm  = value to set open field lines
!
!  outputs:
!        be = magnitude of b field in the equatorial plane in nT
!        vm = flux tube volume**(-2/3)
!             note: vm for open fieldlines is flagged with a
!             large value (big_vm) (input)
!      sini = dip angle in the ionosphere from ground
!      rmin = equatorial r location of mapping in Re
!      pmin = equatorial lt location of mapping in radians
!    iopen  = 2d array flagging open/undefined (0/1) or closed (-1) pnts
!             on rcm grid
!    press  = pressure on field line, data from the mhd code in Pa
!    den    = density on field line, data from the mhd code in ple/cc 
!
!===================================================================
!
!  this routine maps either to the magnetopause        
!  or to the conjugate hemisphere and generates        
!  a potential distribution on a grid in the polar cap 
!  it does this by tracing field lines                 
!  using the trac tracer.                             
!  uses gehavo for the field line tracing              
!  created 03/90 last mod 08/95                       
!  01/19 - frt -this version gets all the values from the MHD code
!

      integer(iprec) :: i,j
      REAL(rprec) :: vm_rcmmhd(isize,jsize-jwrap+1) !Holder for vm calc

      !Pull RCM-MHD variables into RCM arrays and scale/wrap
      !RCM-MHD => RCM (sizes)
      call EmbiggenWrap(RM%Pave,press)
      call EmbiggenWrap(RM%Nave,den  )

      call EmbiggenWrap(RM%x_bmin(:,:,1)/RM%planet_radius,xmin)
      call EmbiggenWrap(RM%x_bmin(:,:,2)/RM%planet_radius,ymin)
      call EmbiggenWrap(RM%x_bmin(:,:,3)/RM%planet_radius,zmin)

      call EmbiggenWrap(RM%bmin/nt,bmin)
      call EmbiggenWrap(RM%beta_average,beta_average)

      call EmbiggenWrap(RM%radcurv,radcurv)
      call EmbiggenWrap(RM%losscone,losscone)

      call EmbiggenWrapI(RM%iopen,iopen)
      call EmbiggenWrap (RM%wImag,wImag)

    ! compute vm and find boundaries
      !Calculate vm on RCM/MHD-sized grid
      where (RM%vol > 0.0)
        vm_rcmmhd = 1.0/(RM%vol*nt)**(2.0/3.0) ! (nt/re)^0.667
      elsewhere
        vm_rcmmhd = big_vm
      end where

      !Convert to RCM-sized grid
      call EmbiggenWrap(vm_rcmmhd,vm)
      
      !Make sure topology is right
      where (vm<0)
        iopen = RCMTOPOPEN
      endwhere

      !Find boundary
      vbnd(:) = 2
      do j=jwrap,jsize
        do i=isize,2,-1
          if ( iopen(i,j) /= RCMTOPCLOSED ) then
            vbnd(j) = i
            exit
          endif
        enddo
      enddo

      ! wrap vbnd
      do j=1,jwrap-1
        vbnd(j)   = vbnd(jsize-jwrap+j)
      end do

      ! compute rmin,pmin
      rmin = sqrt(xmin**2 + ymin**2 + zmin**2)
      pmin = atan2(ymin,xmin)

      if (any(isnan(vm))) then
        write(6,*) 'RCM: NaN in Calc_FTV'
        write(6,*) 'Num NaNs = ', count(isnan(vm))
        ierr = -1
        return
      endif

      ierr = 0

      END SUBROUTINE Calc_ftv

!--------------------------------------------------
!
      SUBROUTINE Gettemp (ierr)
!
! routine to compute an estimate for temperature given
! the pressure assume that alam*vm = energy = kt
! 2/00 frt
! bug fix to fac 2/19 frt
! K: 8/20 - fix to temperature to remove plasmasphere contribution
! inputs:
! isize,jsize - rcm grid dimensions
!       r,p - equatorial mapping location of a grid point in Re and rad
!     press - press in Pa
!        vm - flux tube volume^(-2/3) in (nt/Re/^(-2/3)   
!      open - open/closed flag (-1 = closed)
!       den - ple density in ple/m^3
! outputs:
!      ti,te - temperaure if ions and electrons in kelvin
!      ierr  - error flag
!
!-------------------------------------------------
      USE Rcm_mod_subs, ONLY : isize,jsize,eeta,vm
      use conversion_module
      IMPLICIT NONE
!
      REAL(rprec), PARAMETER :: fac = tiote/(1.+tiote)
      REAL(rprec) :: dmhd,dpp
!
      INTEGER(iprec) :: i,j,ierr
      
!
!    set the temperature:

      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j,dmhd,dpp)
      DO j=1,jsize
        DO i=1,isize
          !Get corrected density from MHD
          call PartFluid(i,j,dmhd,dpp)
          IF ( (iopen(i,j) /= RCMTOPOPEN) .and. (dmhd>TINY) ) THEN
            ti(i,j) = fac*press(i,j)/dmhd/boltz
            te(i,j) = ti(i,j)/tiote
          ELSE
            ti(i,j) = 0.0
            te(i,j) = 0.0
          ENDIF

        ENDDO
      ENDDO

      ierr = 0

      if (maxval(press) <=0.) then
        write(6,*)' maxval pressure < 0 in gettemp'
        ierr = -1
      end if
!
      RETURN
      END SUBROUTINE Gettemp

!
!===================================================================
      !Attempt to separate hot/cold components of MHD fluid
      !TODO: Rewrite this messy code
      SUBROUTINE PartFluid(i,j,drc,dpp)
      USE conversion_module
      USE RCM_mod_subs, ONLY : ikflavc,vm,alamc,isize,jsize,kcsize,eeta
    
      IMPLICIT NONE
      integer(iprec), intent(in) :: i,j
      real(rprec), intent(out) :: drc,dpp

      real(rprec) :: dmhd,dpp_rcm,drc_rcm,dtot_rcm,prc_rcm
      real(rprec) :: drc1,drc2
      logical :: doM1

        drc = 0.0
        dpp = 0.0
      !Trap for boring cases
        if (iopen(i,j) == RCMTOPOPEN) then
          return
        endif

        !Get MHD bulk density
        dmhd = den(i,j)
        
        if (.not. use_plasmasphere) then
          dpp = 0.0
          drc = dmhd !All goes to RC fluid
          return
        endif

        if (dmhd <= TINY) then
          drc = 0.0
          dpp = dpp_rcm
          return
        endif
        
      !Get RCM info
        call eta2DP(eeta(i,j,:),vm(i,j),drc_rcm,dpp_rcm,prc_rcm)
        if (dpp_rcm <= TINY) then
          drc = dmhd
          dpp = 0.0
          return
        endif

      !If still here either null or closed and have some work to do
        !Try two ways:
        !1: drc = dmhd-dpp
        !2: Use RCM RC/PP ratio to split dmhd
        drc1 = 0.0
        drc2 = 0.0

        !Method 1
        if (dpp_rcm <= dmhd) then
          !Subtract plasmasphere from RCM from MHD density
          drc1 = dmhd-dpp_rcm
        else
          drc1 = min(dmhd,abs(dpp_rcm-dmhd))
        endif
        if (drc1 < TINY) drc1 = 0.0

        !Method 2
        if ( (drc_rcm > TINY) .and. (dpp_rcm > TINY) ) then
          dtot_rcm = drc_rcm + dpp_rcm
          drc2 = dmhd*(drc_rcm/dtot_rcm)
        else
          drc2 = 0.0
        endif

      !Pick which split to use
        if ( (drc1 <= TINY) .and. (drc2 <= TINY) ) then
          !Both methods failed, just return unsplit density
          drc = dmhd
          dpp = 0.0
          return
        endif

        !If still here at least one good method
        doM1 = .false.
        
        if ( (drc1 > TINY) .and. (drc2 > TINY) ) then
          !Both methods successful, choose min
          if (drc1 < drc2) then
            doM1 = .true.
          else
            !Second method
            doM1 = .false.
          endif
        else if (drc1 > TINY) then
          !Only method 1 worked
          doM1 = .true.
        else
          !Only method 2 worked
          doM1 = .false.
        endif

        if (doM1) then
          !Method 1
          drc = drc1
          dpp = dpp_rcm
        else
          !Second method
          drc = drc2
          dpp = dmhd*(dpp_rcm/dtot_rcm)
        endif

        !write(*,*) 'M: dmhd / drc / dpp = ', doM1,1.0e-6*dmhd,1.0e-6*drc,1.0e-6*dpp
      END SUBROUTINE PartFluid
!
!===================================================================      
      SUBROUTINE Press2eta() 
      USE conversion_module      
      USE RCM_mod_subs, ONLY : ikflavc,vm,alamc,isize,jsize,kcsize,eeta
      USE rcm_precision
      IMPLICIT NONE
      
      real(rprec) :: drc,dpp,prc
      integer(iprec) :: i,j


      !$OMP PARALLEL DO default(shared) &
      !$OMP schedule(dynamic) &
      !$OMP private(i,j,drc,dpp,prc)
      do j=1,jsize
        do i=1,isize
          !Get corrected density from MHD
          call PartFluid(i,j,drc,dpp)

          if ( (iopen(i,j) /= RCMTOPOPEN) .and. (drc>TINY) .and. (ti(i,j)>TINY) ) then
            !Good stuff, let's go
            prc = press(i,j)
            call DP2eta(drc,prc,vm(i,j),eeta_new(i,j,:))
        !Not good MHD
          else
            eeta_new(i,j,:) = 0.0
          endif
        enddo
      enddo !J loop

      END SUBROUTINE Press2eta
      

!===================================================================
    SUBROUTINE Set_ellipse(idim,jdim,rmin,pmin,vm,big_vm,bndloc,iopen)

! routine that fits an ellipse and resets the modeling
! boundary to be inside the ellipse
! 6/03 frt
! inputs:
!	idim,jdim - size of the 2d rcm arrays (lat, long)
!	xe,ye - equatorial mapping point of field line from rcm grid
!	vm - computed flux tube volume ^(-2/3) set to big_vm if open
!            or outside the ellipse boundary (output)
! 	bndloc/imin_j - boundary location (output)
!	a1 - dayside end of ellipse
!	a2 - nightside location of the ellipse
!	b - semi minor axis (y) of ellipse
      USE rice_housekeeping_module, ONLY : ellBdry
      
      implicit none
      integer(iprec) :: idim,jdim
      real(rprec) :: rmin(idim,jdim), pmin(idim,jdim)
      real(rprec) :: xe(idim,jdim), ye(idim,jdim)
      real(rprec) :: vm(idim,jdim)
      integer(iprec) :: iopen(idim,jdim)
      real(rprec) :: bndloc(jdim)
      real(rprec) :: big_vm,a1,a2,a,bP,bM,x0,ell
      real(rprec) :: xP,xM,yMaxP,yMaxM
      integer(iprec) :: i,j,jm,jp,newi,oldi
      logical :: isGood(idim,jdim)

!  x0 = (a1 + a2)/2
!   a = (a1 - a2)/2
!
!                   b
!           |       |
!           |       |
! x<a1------0-------x0-------------a2
!           |       |
!           |       |
!                   b
!

      xe = rmin * cos(pmin)
      ye = rmin * sin(pmin)

!K: Replacing these hard-coded values with ellipse type set by XML file
      a1 = ellBdry%xSun
      a2 = ellBdry%xTail
      !Separate positive/negative y-axis bounds
      bP = ellBdry%yDD
      bM = ellBdry%yDD

      if (ellBdry%isDynamic) then
        !Tune to current equatorial bounds
        xP   = maxval(xe     ,mask=iopen<0)
        xM   = minval(xe     ,mask=iopen<0)
        !Test dusk/dawn y extents
        yMaxP = maxval( ye,mask=iopen<0) !Dusk extent
        yMaxM = maxval(-ye,mask=iopen<0) !Dawn extent
        !Enforce max's from XML ellipse
        a1 = min(a1,xP)
        a2 = max(a2,xM)
        bP = min(bP,yMaxP)
        bM = min(bM,yMaxM)
      endif
      
      x0 = (a1 + a2)/2.
      a  = (a1 - a2)/2.

      !Fill isGood array
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,ell,jm,jp)
      do j=1,jdim
        do i=1,idim
          if (iopen(i,j) == RCMTOPOPEN) then
            isGood(i,j) = .false.
          else
          !This is closed field line region
            !Check ellipse, separate dawn/dusk extent
            if (ye(i,j)>=0) then
              ell = ((xe(i,j)-x0)/a)**2+(ye(i,j)/bP)**2
            else
              !Dawn
              ell = ((xe(i,j)-x0)/a)**2+(ye(i,j)/bM)**2
            endif !Dusk/dawn

            !Also check longitudinal neighbors
            jm = j-1
            if (jm < 1) jm = jdim - jwrap
            jp = j+1
            if(jp > jdim) jp = jwrap + 1

            isGood(i,j) = (ell<=1.0) .and. (iopen(i,jm) /= RCMTOPOPEN) .and. (iopen(i,jp) /= RCMTOPOPEN)
          endif
        enddo
      enddo !j loop

      !Now set boundary based on isGood and number of required buffer cells
      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j)
      do j=1,jdim
        do i=idim,1,-1
          !Loop from the top and find first bad cell
          if (.not. isGood(i,j)) exit
        enddo !i loop
        !Cell i,j is bad or got to i=1
        newi = min(i+1+NumG,idim)
        !newbndloc(j) = newi
        oldi = bndloc(j)
        !Throttle change in bndloc to NumG if possible
        if (newi > oldi+NumG) then
          !Boundary is moving inwards, limit inwards motion if the cells are good
          if ( all(isGood(oldi:newi,j)) ) then
            !All the cells are good, so limit retreat to NumG cells
            newi = min(oldi+NumG,idim)
          endif
        else if (newi < oldi-NumG) then
          !Trying to expand by too many cells
          newi = max(oldi-NumG,1)
        endif
        bndloc(j) = newi
      enddo !j loop

    end subroutine Set_ellipse


!------------------------------------

      SUBROUTINE Smooth_eta_at_boundary(idim,jdim,kdim,jwrap,eeta,iopen,imin_j)
! this routine attempts to smooth out high frequency noise at the boundary
! of the rcm 
! written 2/06 frt
      IMPLICIT NONE
      INTEGER(iprec) :: idim,jdim,kdim,jwrap
      INTEGER(iprec) :: imin_j(jdim)
      INTEGER(iprec) :: i,j,k,jm,jmm,jp,jpp
      REAL(rprec) :: eeta(idim,jdim,kdim)
      integer(iprec), intent(in) :: iopen(idim,jdim)
      REAL(rprec) :: eetas2d(jdim,kdim)
! these are the smoothing weights
      REAL(rprec) :: a1,a2,a3,a4,a5
      ! REAL(rprec), PARAMETER :: a1 = 1.0  
      ! REAL(rprec), PARAMETER :: a2 = 1.0  
      ! REAL(rprec), PARAMETER :: a3 = 2.0  
      ! REAL(rprec), PARAMETER :: a4 = 1.0  
      ! REAL(rprec), PARAMETER :: a5 = 1.0
      integer(iprec) :: klow,di
      integer(iprec), parameter :: NumI = NumG

      logical :: isOpen(5)
! now do the smoothing


      !Smooth RC channels          
      if (use_plasmasphere) then
        klow = 2
      else
        klow = 1
      endif

      !Loop over di levels at boundary and do j-smoothing
      do di=0,NumI
        !Use more centered weights as we move into the domain
        a3 = 1.00
        a2 = 0.50/(2.0**di)
        a1 = 0.25/(2.0**di)
        a4 = a2
        a5 = a1

        do j=1,jdim

          ! 1 <=> jdim -jwrap +1
          ! jdim <=> jwrap
          jmm = j - 2
          if(jmm < 1)jmm = jdim - jwrap - 1       
          jm  = j - 1
          if(jm < 1) jm = jdim - jwrap       
          jpp = j + 2
          if(jpp > jdim)jpp = jwrap + 2
          jp  = j + 1
          if(jp > jdim) jp = jwrap + 1

          isOpen(1) = (iopen(imin_j(jmm)+di,jmm) == RCMTOPOPEN)
          isOpen(2) = (iopen(imin_j(jm )+di,jm ) == RCMTOPOPEN)
          isOpen(3) = (iopen(imin_j(j  )+di,j  ) == RCMTOPOPEN)
          isOpen(4) = (iopen(imin_j(jp )+di,jp ) == RCMTOPOPEN)
          isOpen(5) = (iopen(imin_j(jpp)+di,jpp) == RCMTOPOPEN)

          if ( any(isOpen) ) then
            !Keep old values b/c too close to OCB
            eetas2d(j,klow:kdim) = eeta(imin_j(j  )+di,j  ,klow:kdim)
          else
            !Only smooth if all closed/null cells
            !Only smooth RC, plasmasphere would diffuse too much
            do k=klow,kdim
              eetas2d(j,k) = ( a1*eeta(imin_j(jmm)+di,jmm,k) + &
                               a2*eeta(imin_j(jm )+di,jm ,k) + &
                               a3*eeta(imin_j(j  )+di,j  ,k) + &
                               a4*eeta(imin_j(jp )+di,jp ,k) + &
                               a5*eeta(imin_j(jpp)+di,jpp,k) )/(a1+a2+a3+a4+a5)
            enddo !k loop
          endif !OCB

        enddo !j loop

        !Now go back and reset values
        do j=1,jdim
          eeta(imin_j(j)+di,j,klow:kdim) = eetas2d(j,klow:kdim)
        enddo

      enddo !di loop
      
      
      END SUBROUTINE Smooth_eta_at_boundary

!------------------------------------
      SUBROUTINE smooth_boundary_location(idim,jdim,jwrap,bndloc)
      USE rice_housekeeping_module, ONLY : L_write_vars_debug
      IMPLICIT NONE
      INTEGER(iprec), INTENT(IN) :: idim,jdim,jwrap
      REAL(rprec), INTENT(IN OUT) :: bndloc(jdim)
     
      REAL(rprec) :: bndloc_new(jdim)

      if (L_write_vars_debug) then
       write(6,*)' smooth_boundary_location, old boundary'
       write(6,*)bndloc
      end if

      call SmoothJBnd(bndloc,bndloc_new)

      !Store back w/ trap for only earthward
      bndloc = max(bndloc_new,bndloc)

      if (L_write_vars_debug) then
       write(6,*)' smooth_boundary_location, new boundary'
       write(6,*) bndloc
      end if

      return
    
      END SUBROUTINE smooth_boundary_location

!
      subroutine set_plasmasphere(icontrol,idim,jdim,kdim,xmin,ymin,rmin,vm,eeta,imin_j)
! subroutine set_plasmasphere(idim,jdim,kdim,rmin,pmin,vm,eeta,imin_j)
! crude routine to set a plasmasphere model in the rcm
! alam(1) should be set to a small value (0.01)
! 2/07 frt
! Use the gallagher model for initial condition, update plasmaspheric eeta in each RCM call
! alam(1) is set to be 0
! sbao 03/25

      USE earthhelper, ONLY : GallagherXY
      USE rice_housekeeping_module, ONLY: InitKp, staticR
      USE kdefs, ONLY : TINY
      
      IMPLICIT NONE

      integer(iprec) :: idim,jdim,kdim,icontrol
      real(rprec) :: dens_gal = 0.0
      integer(iprec) :: imin_j(jdim)
      real(rprec) :: vm(idim,jdim),xmin(idim,jdim),ymin(idim,jdim),rmin(idim,jdim)
      real(rprec) :: eeta(idim,jdim,kdim)

      integer(iprec) :: i,j,k

      if (icontrol == RCMCOLDSTART) then
        do j=1,jdim
          do i=imin_j(j),idim
            if(vm(i,j) > 0.0)then
              dens_gal = GallagherXY(xmin(i,j),ymin(i,j),InitKp)*1.0e6
              eeta(i,j,1) = dens_gal/(GetDensityFactor()*vm(i,j)**1.5)
            end if
          end do
        end do
      else
        ! reset the static part of the plasmasphere sbao 07292020
        !Tweak by K: 8/7/20
        if (staticR > TINY) then
          !$OMP PARALLEL DO default(shared) &
          !$OMP schedule(dynamic) &
          !$OMP private(i,j,dens_gal)
          do j=1,jdim
            do i=imin_j(j),idim
              if(rmin(i,j) <= staticR .and. vm(i,j) > 0.0)then
                !eeta (i,j,1) = eeta_pls0 (i,j)
                dens_gal = GallagherXY(xmin(i,j),ymin(i,j),InitKp)*1.0e6
                eeta(i,j,1) = dens_gal/(GetDensityFactor()*vm(i,j)**1.5)
              end if
            end do
          end do
        end if !staticR
      endif !RCMCOLDSTART

      return

      end subroutine set_plasmasphere

!------------------------------------------      
      subroutine reset_rcm_vm(idim,jdim,bndloc,big_vm,imin_j,vm,iopen)
! this routine resets imin_j, vm, and open based on a newly set bndloc      
      implicit none
      integer(iprec), intent(in) :: idim,jdim
      integer(iprec),intent(inout) :: imin_j(jdim),iopen(idim,jdim)
      real(rprec), intent(inout) :: bndloc(jdim)
      real(rprec), intent(in) :: big_vm
      real(rprec), intent(inout) :: vm(idim,jdim)

      integer(iprec) :: i,j,iC

      imin_j = CEILING(bndloc)

      !Loop through grid and reset closed cells below boundary to NULL and poison vm in open cells

      !$OMP PARALLEL DO default(shared) &
      !$OMP private(i,j,iC)
      do j=1,jdim
        iC = imin_j(j)
        do i=1,idim
          if ( (iopen(i,j) == RCMTOPCLOSED) .and. (i<iC) ) then
            !Closed field region outside RCM domain
            !Make this a buffer cell but keep physical vm
            iopen(i,j) = RCMTOPNULL
          endif

          if (iopen(i,j) == RCMTOPOPEN) then
            !Open cell, poison it
            vm(i,j) = big_vm
          endif

        enddo !i loop
      enddo !j loop

      !Finish up by resetting boundary
      do j=1,jsize
        bndloc(j) = 2  ! if everything else fails, can use this...
        do i=isize,2,-1
          if (iopen(i,j) >= 0) then !null or open
            bndloc(j) = i + 1
            exit
          endif
        end do
      end do
      ! reset imin_j
      imin_j = ceiling(bndloc)

      end subroutine reset_rcm_vm


!
!======================================
      subroutine allocate_conversion_arrays(isz,jsz,kcsz)
! used to allocate memory for the exchange arrays      
! 7/09 frt
      use conversion_module
      implicit none
      integer(iprec),intent(in) :: isz,jsz,kcsz
      integer(iprec) :: idim,jdim,kdim
! if the arrays are allocated, then return
      if(allocated(x0))return

      idim = isz
      jdim = jsz
      kdim = kcsz

      write(*,*)' Allocating conversion arrays'

      ! 1d arrays
      allocate(imin_j_old(jdim))
      allocate(bndloc_old(jdim))
      allocate(almmin(kdim))
      allocate(almmax(kdim))
      allocate(almdel(kdim))
      ! 2d arrays
      allocate(x0_sm(idim,jdim))
      allocate(y0_sm(idim,jdim))
      allocate(z0_sm(idim,jdim))
      allocate(x0(idim,jdim))
      allocate(y0(idim,jdim))
      allocate(z0(idim,jdim))

      allocate(den(idim,jdim))
      allocate(press(idim,jdim))
      allocate(deno(idim,jdim))
      allocate(presso(idim,jdim))
      allocate(te(idim,jdim))
      allocate(ti(idim,jdim))
      allocate(to(idim,jdim))
      allocate(beta_average(idim,jdim))
      allocate(iopen(idim,jdim))
      allocate(wImag(idim,jdim))

      ! 3d arrays
      allocate(eeta_new(idim,jdim,kdim))
    
     return

     end subroutine allocate_conversion_arrays 

END MODULE torcm_mod
