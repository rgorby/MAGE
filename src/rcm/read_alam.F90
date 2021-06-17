      SUBROUTINE Read_alam (kdim, alam, iflav,fudge, almdel, &
                            almmax, almmin, iedim, ierr)
      use rcmdefs
      USE rcm_precision
      use ioh5
      use files

      IMPLICIT NONE
      INTEGER(iprec), INTENT (IN) :: kdim, iedim
      INTEGER(iprec), INTENT (OUT) :: ierr
      REAL(rprec), INTENT (IN OUT) :: alam(kdim), almdel(kdim), &
                             almmax(kdim), almmin(kdim), fudge(kdim)
      INTEGER(iprec), INTENT (IN OUT), DIMENSION (kdim) :: iflav
!
!     This routine was re-designed so that it communicates with 
!     the outside world only through its list of arguments.
!     Thus, can use it to get information for grid-based energy
!     channels, or for edges. (Stanislav)
!
      INTEGER(iprec) :: NIn,k, iflavin, ie, lun = 1
      INTEGER(iprec) :: num_chan (iedim), k_start, k_stop
      REAL(rprec)  :: alamin, amin, amax,fudgein
      LOGICAL :: lflag_1, lflag_2
! 
      type(IOVAR_T), dimension(RCMIOVARS) :: IOVars !Lazy hard-coding max variables
      logical :: doSP

      INCLUDE 'rcmdir.h'
! 
      if (isGAMRCM) then

        !Use Gamera HDF5 stuff
        doSP = .false.
        call ClearIO(IOVars) !Reset IO chain
        call AddInVar(IOVars,"ikflavc")
        call AddInVar(IOVars,"alamc")
        call AddInVar(IOVars,"fudgec")
        call ReadVars(IOVars,doSP,RCMGAMConfig)

        !Test all read variables to have the right size
        do k=1,3
          NIn = IOVars(k)%N
          if (NIn /= kdim) then
            write(*,*) 'RCM configuration error, mismatched k sizes ...'
            stop
          endif
        enddo
        
        !Lazily replicating loop from below
        DO k = 1, kdim
          iflavin = IOVars(1)%data(k)
          alamin  = IOVars(2)%data(k)
          fudgein  = IOVars(3)%data(k)
          
          IF (iflavin == 1) THEN
             alam (k) = alamin
             IF (alam(k) > 0.0) alam(k) = - alam(k)
!             fudge (k) = 0.3333
             fudge (k) = fudgein
             iflav (k) = 1
          ELSE IF (iflavin == 2) THEN
             alam (k) = alamin
!             fudge (k) = 0.0000
             fudge (k) = fudgein
             iflav (k) = 2
          ELSE
             STOP 'ILLEGAL TYPE OF SPECIES'
          END IF

        ENDDO !k

      else
        write(*,*) "This ain't gonna work!"
        stop
  !       !Do old-style binary records
  !       INQUIRE (UNIT = lun, EXIST = lflag_1, OPENED = lflag_2)
  !       IF (.not.lflag_1 .OR. lflag_2) THEN
  !          STOP 'UNABLE TO USE LUN IN READ_ALAM'
  !       END IF
  !       OPEN (LUN, file = rcmdir//'enchan.dat', STATUS = 'old', IOSTAT = k)
  !       IF (k /= 0 ) THEN
  !          WRITE (*,'(T2,A)') 'ERROR OPENING FILE IN READ_ALAM'
  !          ierr = -1
  !          RETURN
  !       END IF
  ! !
  !       READ (LUN,*) k
  !       IF (k /= kdim)then
  !         write(*,*)' k = ',k,' kdim =',kdim
  !         STOP 'kcsizes inconsistent in READ_ALAM'
  !       END IF
  ! !
  !       DO k = 1, kdim
  !         READ (LUN,*) iflavin, alamin
  !         IF (iflavin == 1) THEN
  !            alam (k) = alamin
  !            IF (alam(k) > 0.0) alam(k) = - alam(k)
  !            fudge (k) = 0.3333
  !            iflav (k) = 1
  !         ELSE IF (iflavin == 2) THEN
  !            alam (k) = alamin
  !            fudge (k) = 0.0000
  !            iflav (k) = 2
  !         ELSE
  !            STOP 'ILLEGAL TYPE OF SPECIES'
  !         END IF
  !       END DO
  !       CLOSE (LUN)
      endif !isGAMRCM
!
!
!     check to see how many different types of speces there are
!
      num_chan(:) = 0
      DO k = 1, kdim
        ie = iflav (k)
        num_chan (ie) = num_chan (ie) + 1
        IF (num_chan(ie) > 1) THEN
           IF (ABS(alam(k)) < ABS(alam(k-1))) THEN
              STOP ' ERROR: enchan channels not in increasing order'
           END IF
           IF (ie > iedim) STOP 'ILLEGAL SPECIES IN READ_ALAM'
        END IF
      END DO
!
!
      DO k = 1, kdim
        ie = iflav (k)
        IF (num_chan(ie) == 1) THEN
           almmax (k) = 2.0*alam(k)
           almmin (k) = 0.0
        ELSE 
          IF (ie == 1) THEN
            k_start = 1
            k_stop  = num_chan (ie)
          ELSE IF (ie == 2) THEN
            k_start = num_chan(ie-1)+1
            k_stop  = num_chan(ie-1)+num_chan(ie)
          ELSE
            STOP 'WRONG IE'
          END IF

          IF (k == k_start) THEN
            almmax(k) = 0.5*(ABS(alam(k))+ABS(alam(k+1)))
            almmin(k) = 0.0
          ELSE IF (k == k_stop) THEN
            almmax(k) = 1.5*ABS(alam(k)) - 0.5*ABS(alam(k-1))
            almmin(k) = 0.5*(ABS(alam(k-1))+ABS(alam(k)))
          ELSE
            almmax(k) = 0.5*(ABS(alam(k))+ABS(alam(k+1)))
            almmin(k) = 0.5*(ABS(alam(k-1))+ABS(alam(k)))
          END IF
        END IF
        almdel(k)= ABS(almmax(k))-ABS(almmin(k))
      END DO 
!
      ierr = 0
      RETURN
 100  WRITE (*,'(T2,A)') 'error on opening enchan.dat'
      ierr = -1

      RETURN
      END SUBROUTINE Read_alam
