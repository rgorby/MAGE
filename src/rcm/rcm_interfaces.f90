! File to hold RCM interfaces to other codes

   SUBROUTINE Grid_torcm (high_lat, low_lat, offseti, Re_external, Ri_external)

      USE Rcm_mod_subs, imin_rcm=>imin
      USE conversion_module, ONLY: idim,jdim,kdim,x0,y0,z0
      USE rice_housekeeping_module, ONLY: rcm_tilted
      USE constants, ONLY : RCMCorot
      IMPLICIT NONE

      REAL(rprec), INTENT (IN) :: high_lat, low_lat, offseti, Re_external, Ri_external


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


      INTEGER (iprec), PARAMETER :: imin = 1
      INTEGER(iprec) :: i, j, jmin
      REAL(rprec) :: start, end, offset
      REAL(rprec) :: glat (isize), teta (isize), phi (jsize)
      real(rprec) :: x,dir
      CHARACTER (LEN=15) :: grid_file='grid.dat'


      Re   = Re_external / 1.0E+3 ! in km
      Ri   = Ri_external / 1.0E+3 ! in km

      idim = isize
      jdim = jsize

      i1   = imin + 1
      i2   = idim - 1
      jmin = 1 
      j1   = jmin + 2 
      j2   = jdim - 1
      dlam = 1.0 / REAL(idim-1,rprec)
      dpsi = (2.0*pi) / REAL(jdim-jwrap,rprec)   

 
      start  = high_lat*DTR
      end    = low_lat*DTR
      offset = offseti*DTR

      if(jdim==0)then
        write(*,*)'jdim =',jdim
        stop
      end if

      DO i = imin, idim
          glat(i) = start + (end-start) * Fun(REAL(i-imin)/REAL(idim-imin,rprec))
          teta(i) = pi/2 - glat(i)
      END DO

      DO j = j1, j2
          phi(j) = (REAL(j,rprec)-REAL(j1,rprec))*pi /((REAL(jdim,rprec)-REAL(j1,rprec))/2.)
      END DO
      phi (jmin)   = phi(j2-1)
      phi (jmin+1) = phi(j2)
      phi (jdim)   = phi(j1)

      DO i = imin, idim
      DO j = j1, j2
         colat(i,j) = ACOS(COS(teta(i))*COS(offset) + &
                           SIN(teta(i))*SIN(offset)*COS(phi(j)))
         aloct(i,j) = ATAN2(SIN(teta(i))*SIN(phi(j)), &
                            SIN(teta(i))*COS(phi(j))*COS(offset)- &
                            COS(teta(i))*SIN(offset))
         IF (aloct(i,j).lt.0.0) aloct(i,j) = aloct(i,j)+2.0*pi
         IF (aloct(i,j).gt.2.0*pi) aloct(i,j)=aloct(i,j)-2.0*pi
      END DO
         colat(i,jmin)=colat(i,j2-1)
         colat(i,jmin+1)=colat(i,j2)
         colat(i,jdim)=colat(i,j1)
         aloct(i,jmin)=aloct(i,j2-1)
         aloct(i,jmin+1)=aloct(i,j2)
         aloct(i,jdim)=aloct(i,j1)
      END DO

      DO j = 1, jdim
        alpha (1,j) = (teta(2)-teta(1))/dlam
        alpha (idim,j) = (teta(idim)-teta(idim-1))/dlam
        DO i = imin+1,idim-1
           alpha(i,j) = 0.5*(teta(i+1)-teta(i-1))/dlam
        END DO
        do i=imin,idim
           beta(i,j) = SIN(teta(i))
        END DO
      END DO
!
! only add corotation in the non-tilted world
! all the other quantities have to be calculated
      if(rcm_tilted)then
       vcorot = 0.0
       sini   = 0.0
       bir    = 0.0
      else
       vcorot = -RCMCorot*(Re / Ri)*SIN(colat)**2
       sini   = two*COS(colat)/SQRT(one+three*COS(colat)**2)
       bir    = two*(Re / Ri)**3*besu*COS(colat)
      end if


      ! this writes a file for plotting (I think... /ss):

      OPEN (LUN, FILE = rcmdir//grid_file, status='UNKNOWN')
      DO j = 1, jdim
      do i = imin, idim
         x0(i,j) = (Ri/Re)*SIN(colat(i,j))*COS(aloct(i,j))
         y0(i,j) = (Ri/Re)*SIN(colat(i,j))*SIN(aloct(i,j))
         z0(i,j) = (Ri/Re)*COS(colat(i,j))
         WRITE (LUN,'(2(1x,i3),4(1x,f10.4))') &
                 i,j,x0(i,j),y0(i,j),z0(i,j)
200      FORMAT (TR2,a2,TR7,a5,TR10,a10,TR4,a5)
201      FORMAT (TR1,i3,TR5,f12.8,TR3,f12.8) 
      END DO
      END DO
      CLOSE (LUN)


      CALL Write_grid  ! this writes standard RCM grid files


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
! FIXME: assumes jwrap=3
      v (:,jwrap:jsize) = RM%pot (:,      :)
      v (:,          1) = v   (:,jsize-2)
      v (:,          2) = v   (:,jsize-1)

      qtplam (:,jwrap:jsize) = RM%sigmap (:,      :)
      qtplam (:,          1) = qtplam (:,jsize-2)
      qtplam (:,          2) = qtplam (:,jsize-1)

      qtped  (:,jwrap:jsize) = RM%sigmap (:,      :)
      qtped  (:,          1) = qtplam (:,jsize-2)
      qtped  (:,          2) = qtplam (:,jsize-1)

      qthall (:,jwrap:jsize) = RM%sigmah (:,      :)
      qthall (:,          1) = qthall (:,jsize-2)
      qthall (:,          2) = qthall (:,jsize-1)

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
