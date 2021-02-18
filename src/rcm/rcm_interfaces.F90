! File to hold RCM interfaces to other codes

   SUBROUTINE Grid_torcm (high_lat, low_lat, offseti, Re_external, Ri_external,doStretch)

      USE Rcm_mod_subs, imin_rcm=>imin
      USE conversion_module, ONLY: x0,y0,z0
      USE rice_housekeeping_module, ONLY: rcm_tilted
      USE constants, ONLY : RCMCorot
      IMPLICIT NONE

      REAL(rprec), INTENT (IN) :: high_lat, low_lat, offseti, Re_external, Ri_external
      LOGICAL, INTENT(IN) :: doStretch

      ! This routine will generate an ionospheric 2-D grid for the RCM
      !
      ! INPUTS: 
      !  high_lat    is starting latitude of grid (highest latitude, in degrees)
      !  low_lat     is ending latitude of grid   (lowest latitude,  in degrees)
      !  offseti     is offset location (NOTE: not used)
      !  Re_external is Earth radius in meters
      !  Ri_external is radius of Earth's ionospheric shell in meters
      !                  (hopefully Ri_external > Re_external)
      !    
      ! OUTPUTS/RESULTS:
      !  an ascii file with grid arrays and parameters
      ! 
      !  x0,y0,z0 = GSM locations of ionospheric grid points to be used for 
      !             tracing magnetic field lines later
      !    [x0,y0,z0 arrays are assigned values in the module where
      !     they are declared, not saved in a file]
      ! 5/20/20 - removed assumption of jwarp=3, frt


      INTEGER (iprec), PARAMETER :: imin = 1
      INTEGER(iprec) :: i, j
      REAL(rprec) :: start, end, offset
      REAL(rprec) :: glat (isize), teta (isize), phi (jsize)
      real(rprec) :: x,dir
      CHARACTER (LEN=15) :: grid_file='grid.dat'

      write(6,*)' setting up RCM grid'

      Re   = Re_external / 1.0E+3 ! in km
      Ri   = Ri_external / 1.0E+3 ! in km

      i1   = isize + 1
      i2   = isize - 1
      j1 = jwrap
      j2 = jsize - 1
      dlam = 1.0 / REAL(isize-1,rprec)
      dpsi = (2.0*pi) / REAL(jsize-j1,rprec)   

 
      start  = high_lat*DTR
      end    = low_lat*DTR
      offset = offseti*DTR

      if(jsize==0)then
        write(*,*)'grid_torcm:jsize =',jsize
        stop
      end if

      DO i = imin, isize
        if (doStretch) then
          glat(i) = start + (end-start) * Fun(REAL(i-imin)/REAL(isize-imin,rprec))
        else
          !Do uniform spacing
          glat(i) = start + (end-start) * REAL(i-imin)/REAL(isize-imin,rprec)
        endif
        teta(i) = pi/2 - glat(i)
      END DO

      DO j = j1,j2 
          phi(j) = (REAL(j,rprec)-REAL(j1,rprec))*pi /((REAL(jsize,rprec)-REAL(j1,rprec))/2.)
      END DO
 
      DO j=1,j1-1
       phi(j) = phi(jsize-j1+j)
      END DO
      phi (jsize)   = phi(j1)

      DO i = imin, isize
      DO j = j1, j2
         colat(i,j) = ACOS(COS(teta(i))*COS(offset) + &
                           SIN(teta(i))*SIN(offset)*COS(phi(j)))
         aloct(i,j) = ATAN2(SIN(teta(i))*SIN(phi(j)), &
                            SIN(teta(i))*COS(phi(j))*COS(offset)- &
                            COS(teta(i))*SIN(offset))
         IF (aloct(i,j).lt.0.0) aloct(i,j) = aloct(i,j)+2.0*pi
         IF (aloct(i,j).gt.2.0*pi) aloct(i,j)=aloct(i,j)-2.0*pi
      END DO
       DO j=1,j1-1
         colat(i,j)=colat(i,jsize-j1+j)
         aloct(i,j)=aloct(i,jsize-j1+j)
       END DO
         colat(i,jsize)=colat(i,j1)
         aloct(i,jsize)=aloct(i,j1)
      END DO

      DO j = 1, jsize
        alpha (1,j) = (teta(2)-teta(1))/dlam
        alpha (isize,j) = (teta(isize)-teta(isize-1))/dlam
        DO i = imin+1,isize-1
           alpha(i,j) = 0.5*(teta(i+1)-teta(i-1))/dlam
        END DO
        do i=imin,isize
           beta(i,j) = SIN(teta(i))
        END DO
      END DO
!
! only add corotation in the non-tilted world
! all the other quantities have to be calculated

!K: Notes for tilting
!vcorot = 0.0 (or RCMCorot=0)
!sini = two*COS(colat)/SQRT(one+three*COS(colat)**2)
!bir = two*(Re / Ri)**3*besu*COS(colat)

      if(rcm_tilted)then
       vcorot = 0.0
       sini   = 0.0
       bir    = 0.0
      else
       vcorot = -RCMCorot*(Re / Ri)*SIN(colat)**2
       sini   = two*COS(colat)/SQRT(one+three*COS(colat)**2)
       bir    = two*(Re / Ri)**3*besu*COS(colat)
      end if


      RETURN
      CONTAINS
!
         FUNCTION Fun (x)
         USE Rcm_mod_subs, ONLY : rprec
         IMPLICIT NONE
         REAL(rprec), INTENT (IN) :: x
         REAL(rprec) :: Fun
!
         REAL(rprec), PARAMETER :: a = 50.0, b=1.0, xm=0.001
!
         IF (x <= xm) THEN
           Fun = (b-a)/xm*x**2/2.0+a*x
         ELSE
           Fun = ((a-b)*x**2/2.0+(b-a*xm)*x+xm/2.0*(a-b))/(1.0-xm)
         END IF
         Fun = Fun/(a+b)*2.0
         RETURN
         END FUNCTION Fun

   END SUBROUTINE Grid_torcm

   SUBROUTINE Ionosphere_torcm(RM)

      ! This is where we transfer ionospheric 
      ! quantities onto the RCM ionospheric grid.

      USE Rcm_mod_subs
      USE Ionosphere_exchange
      USE rice_housekeeping_module
      use rcm_mhd_interfaces

      IMPLICIT NONE
      type(rcm_mhd_T),intent(inout) :: RM
     
      INTEGER (iprec) :: i,j

     ! Import ionosphere variables (Potential, Pedersen & Hall
     ! conductances) from MIX coupler/solver and store results in variables
     ! in ionosphere_intermediate module.
     ! also transfers average pressure, density, flux tube volume, bmin, and xmin(x,y,z).
     ! Note:  The following calls MUST occur before this:
     !         - intermediate_grid::setupIg()


!      CALL ImportIonosphere
      v (:,jwrap:jsize) = RM%pot (:,      :)
      do j=1,jwrap-1
       v (:,          j) = v   (:,jsize-jwrap+j)
      end do

      qtplam (:,jwrap:jsize) = RM%sigmap (:,      :)
      do j=1,jwrap-1
       qtplam (:,          j) = qtplam   (:,jsize-jwrap+j)
      end do

      qtped  (:,jwrap:jsize) = RM%sigmap (:,      :)
      do j=1,jwrap-1
       qtped (:,          j) = qtped   (:,jsize-jwrap+j)
      end do

      qthall (:,jwrap:jsize) = RM%sigmah (:,      :)
      do j=1,jwrap-1
       qthall (:,          j) = qthall   (:,jsize-jwrap+j)
      end do

      ! Double conductances to account for two hemispheres
      ! and correct for magnetic field inclination

      qtplam = qtplam * 2.0 * sinI 
      qtped  = qtped  * 2.0 / sinI
      qthall = qthall * 2.0 
      
 
! writeout max and minvalues of ionosphere values
      if(L_write_vars_debug)then 
         write(6,*)' Ionosphere torcm:'
         write(6,*)' Max ionospheric potential =',maxval(v)
         write(6,*)' Min ionospheric potential =',minval(v)
         write(6,*)'    Max ionospheric qtplam =',maxval(qtplam)
         write(6,*)'    Min ionospheric qtplam =',minval(qtplam)
         write(6,*)'     Max ionospheric qtped =',maxval(qtped)
         write(6,*)'     Min ionospheric qtped =',minval(qtped)
         write(6,*)'    Max ionospheric qthall =',maxval(qthall)
         write(6,*)'    Min ionospheric qthall =',minval(qthall)
      end if

      ! In case RCM will solve for its own electric field, we
      ! should set the boundary potential array. Here we assume
      ! that RCM's high-latitude boundary has been already set
      ! (call to torcm was made), and that the boundary is along
      ! RCM grid points. If RCM uses potential from MIX, this 
      ! will go unused.

      DO j = 1, jsize
         vbnd (j) = v (imin_j(j), j)
      END DO

      RETURN
   END SUBROUTINE Ionosphere_torcm
