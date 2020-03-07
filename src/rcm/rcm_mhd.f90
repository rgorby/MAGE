subroutine rcm_mhd(mhdtime,mhdtimedt,RM,iflag,iXML)
! version to couple to gamera
! units are assumed to mks, except for distances which are in Re.
! iflag = 0 - setup arrays, read in parameters
! iflag = 1 - run rcm
! iflag = 2 - Resart RCM
! iflag = -1 - stop, write out timing information
! iflag = -2 - Write restart (icontrol = 31337)
! iflag = -3 - Write H5 output (icontrol = 31338)

! 2/19 frt
  use rcm_precision
  use Rcm_mod_subs
  use rcm_mhd_interfaces
  use ionosphere_exchange, only : setupIon, tearDownIon
  use constants, ONLY : radius_earth_m, radius_iono_m
  use rice_housekeeping_module
  use rcm_timing_module
  use files

  implicit none
  type(XML_Input_T), intent(in), optional :: iXML
  type(rcm_mhd_T),intent(inout) :: RM
  real(rprec), intent(in) :: mhdtime,mhdtimedt
  integer(iprec) :: dayOfYear 
!  integer(iprec) :: idim,jdim
!  real(rprec), intent(in) :: v(:,:),pave(:,:),xion(:,:),yion(:,:),zion(:,:)
!  real(rprec), intent(in) :: x_bmin(:,:), y_bmin(:,:),z_bmin(:,:) bmin(:,:),mask(:,:)
!  real(rprec), intent(in) :: pot(:,:),gasgamma
!  real(rprec), intent(out) :: Prcm(:,:), Nrcm(:,:)
  integer(iprec), intent(in) :: iflag
  CHARACTER(LEN=3) :: months(12)=(/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT', 'NOV','DEC'/)

  integer(iprec) :: ierr   !> Error code...
  integer(iprec) :: system !> This code makes a call to the
                           !! non-standard SYSTEM function.  Must
                           !! define as integer to avoid compile error.

!  double precision :: time  !> Modified Julian Date + decimal fraction of
!                            !! day of current sim time (according to LFM)
!  double precision :: time0 !> Modfieid Julian Date + decimal fraction of
!                            !! day when coupling starts
!  double precision :: dt    !> Interval between RCM exchanges in seconds

!  integer(iprec) :: timercm = 0 !> Time (in seconds) 

  integer(iprec) :: itimei !> RCM(...) param:  start time
  integer(iprec) :: itimef !> RCM(...) param:  end time
  real(rprec) :: time0 = 0 ! coupling start time
  integer(iprec) :: ircm_dt
  integer(iprec) :: itimef_old = -1
  integer(iprec) :: irdr = 1 !> RCM(...) param:  record # to read in
  integer(iprec) :: irdw = 1 !> RCM(...) param:  record # to write out
  integer(iprec) :: idt   !> RCM(...) param:  basic time step in program
  integer(iprec) :: idt1  !> RCM(...) param:  time step for
                          !! changing disk & write records
  integer(iprec) :: idt2  !> RCM(...) param:  time step for
                          !! writting formatted output

  integer(iprec) :: rec = 1 !> Record number that torcm/tomhd use to

  real(rprec) :: t1, t2  !> Used for performance timing

  !> Model coupling variables
  integer(iprec) :: exchangeNum = 0
  logical :: isFirstExchange
  logical :: isLastExchange

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! bypass for now
     time0 = 0 ! FIXME set for now

  IsCoupledExternally = .TRUE.  ! switch RCM to "coupled" mode before doing anything else

  if (doRCMVerbose) write (*,'(TR1,A,L7)') 'Welcome to the RCM, IsCoupledExternally=', IsCoupledExternally
! setup rcm,time in integer format
  itimei = nint(mhdtime-time0,iprec)
  itimef = nint(mhdtime + mhdtimedt-time0,iprec)
  ircm_dt = itimef - itimei
  
! finish up
  if(iflag==RCMWRITETIMING)then
    call write_rcm_timing(rcm_timing)
    return
  end if

! Write restart file
   if (iflag==RCMWRITERESTART) then
      CALL Rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2,icontrol=ICONWRITERESTART,stropt=RM%rcm_runid,nslcopt=RM%RCM_nRes)
      return
   endif
! Write output slice
   if (iflag==RCMWRITEOUTPUT) then
      CALL Rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2,icontrol=ICONWRITEOUTPUT,stropt=RM%rcm_runid,nslcopt=RM%RCM_nOut)
      return
   endif

! initialize
  if( (iflag == RCMINIT) .or. (iflag == RCMRESTART) ) then !Do this for initialization and restart?
    !Make sure RCM directory exists
    CALL CheckDirOrMake(Rcmdir)

    !Read RCM/MHD params from XML
    if(present(iXML)) then
      CALL RCM_MHD_Params_XML(iXML)
    else
      CALL RCM_MHD_Params_XML
    endif

    ! setup rcm
    CALL Rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2,icontrol=0_iprec)

    call allocate_conversion_arrays (isize,jsize,kcsize)

    ! Set up RCM ionospheric grid:
    !call Grid_torcm (75.0_rprec, 15.0_rprec, 0.0_rprec, radius_earth_m, radius_iono_m)  ! set up RCM ionospheric grid here
    call Grid_torcm (HighLatBD,LowLatBD, 0.0_rprec, radius_earth_m, radius_iono_m)  ! set up RCM ionospheric grid here

   ! Setup Ionosphere intermediate Grid by equating it to the RCM grid, without angular overlap:
    call setupIon(RM)
  
    CALL Rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2, icontrol=1_iprec)

    ! icontrol of 2 also needs the input xml file
    CALL Rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2, icontrol=2_iprec, iXML=iXML)
    

    ! restart
    if (iflag == RCMRESTART) then
      !HDF5 RESTART
      
      !Check if timing file exists
      if ( CheckFile(Rcmdir//"rcm_timing.dat") ) then
        !Reset timing for record business
        call read_rcm_timing(rcm_timing)
        call find_record(itimei,rcm_timing,rec)
      else
        !Create basic setup
        rec = 1
        call  AddToList(itimei,rcm_timing )
        !Write null timing data
        call write_rcm_timing(rcm_timing)
      endif
      
      irdr = rec
      irdw = rec

      !Read in HDF5 restart data
      CALL Rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2,icontrol=ICONRESTART,stropt=RM%rcm_runid,nslcopt=RM%RCM_nRes)

      return

    endif

    if(itimei>0)then
      !OLD RESTART CODE
      write(*,*) 'Ths code should not get called ...'
      stop

      write(*,*) 'I am here in itime>0!!!'
      write(*,*)' RCM RESTART at t=',itimei

      call read_rcm_timing(rcm_timing)
      call find_record(itimei,rcm_timing,rec)

      irdr = rec
      irdw = rec

      CALL Rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2, icontrol=3_iprec)

    end if

    return

  end if !RCMINIT or RCMRESTART



     if(iflag==RCMADVANCE)then ! run the rcm

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Determine exchange times...
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! set the record number; Note: nint seems to fix some rounding errors.
     ! 2 records per timestep 
      if(itimei==0)then ! first call to rcm
!       call  AddToList(itimei,rcm_timing )
       exchangeNum = 0

      else

      !K: uncommented to try to help with making restarts work
      !call find_record(itimei,rcm_timing,rec)
      !rec = nint(2.0*(mhdtime-time0)/mhdtimedt + 1,iprec)

     ! now move to the next record
       rec = rec + 1
      end if

    isFirstExchange = (exchangeNum==0)

    if (doRCMVerbose) then
      write(*,*)'-----------rcm_mhd: rec=',rec

      write(*,*) 'itimei = ', itimei
      write(*,*) 'exchangeNum = ', exchangeNum
      WRITE (*,'(//)')
      write (*,'(a,i4,a,g12.4,a,i4)') 'RCM: time=',itimei,'  time0=',time0, '  Delta_t[s]=',ircm_dt
      write (*,'(a,i4,a,i3,a,g12.4)') 'RCM: ____T_rcm[s] =', itimei,'  rec=',rec,'  T_MHD=',mhdtime
      WRITE (*,'(//)')

    endif
     
    idt = Idt_overwrite ! RCM internal time step in seconds
    ! Frequency (in seconds) to change disk & write records
    idt1 = itimef - itimei
    ! Frequency (in seconds) to write formatted output
    idt2 = idt1

    ! now round to to fit the correct number rcm timesteps
    itimef = itimei + idt *((itimef-itimei)/idt)
    itimef_old = itimef

     IF (isFirstExchange .AND. itimei == 0) then 
        !special case. Tell RCM to start from rec=1, but most initial
        ! conditions will actually be done in TORCM
        irdr = 1
     ELSE
        irdr = rec - 1
     END IF

     irdw = rec

     if (isFirstExchange) then ! Set RCM initial conditions on plasma:
        call rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2, icontrol=3_iprec)
     end if


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Import data from MIX
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call cpu_time(t1)
    if (doRCMVerbose) then
      write(6,'(2(a,i4))')'RCM: calling torcm with rec =',rec,' itimei=',itimei
      call print_date_time(6_iprec)
    endif

     call torcm(RM,rec,itimei,ierr)
     if (ierr > 0 )then
        stop 'RCM: error in torcm '
     endif
     exchangeNum = exchangeNum + 1
     call cpu_time(t2)

     if (doRCMVerbose) write(*,'(a,g14.4,a)')'RCM: torcm cpu time= ',t2-t1,' seconds'

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Advance RCM
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call cpu_time(t1)
    if (doRCMVerbose) then 
      write(6,'(a,i5,a,i5,a,i5,a)')'RCM: call rcm at itimei =',itimei,' to itimef =',itimef,' dt=',ircm_dt, ' sec'
      call print_date_time(6_iprec)
    endif

     ! now run the rcm
     call rcm (itimei, itimef, irdr, irdw, idt, idt1, idt2, icontrol=4_iprec)
     rec = rec + 1 ! update record after rcm has run

     call cpu_time(t2)
    if (doRCMVerbose) then
      write(*,'(a,g14.4,a)')'RCM_MHD:   rcm cpu time= ',t2-t1,' seconds'
      call print_date_time(6_iprec)
    endif

    ! Do not export data if this is both the first & last exchange.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Export data to MHD code
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (doRCMVerbose) write(6,'(a,i5.5)')'RCM_MHD: calling tomhd with rec =',rec
     call cpu_time(t1)
     call Tomhd (RM,rec, ierr)

     call cpu_time(t2)
    if (doRCMVerbose) then
      call print_date_time(6_iprec)
      write(*,*)'RCM: tomhd cpu time= ',t2-t1,' seconds'
   endif

     if (ierr > 0 )then
        stop 'RCM: error in tomhd '
     end if
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! end do couplingTimeLoop
  end if

  if(iflag==2)then ! stop
  call rcm (itimei,itimef,irdr,irdw,idt,idt1,idt2,icontrol=5_iprec)
!  call Finalize()    ! Matches Initialize() above
  call tearDownIon(RM) ! Matches setupIon() above
  end if

  return

end subroutine rcm_mhd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine min2hr(time,hours)
! returns a string of hr:min:sec from an input in minutes
! 2/04 frt
!      USE Rcm_mod_subs, ONLY : iprec,rprec
      USE rcm_precision
      implicit none
      real(rprec) :: time
      real(rprec) :: time_hours
      real(rprec) :: time_minutes
      real(rprec) :: time_seconds
      character (len=*) :: hours
      character (len=3) :: char_time_hours
      character (len=2) :: char_time_minutes,char_time_seconds
! time is in minutes

      time_hours = floor(time/60.)
      time_minutes = floor(time - 60*time_hours)
      time_seconds = 60*(time -floor(time))

      write(char_time_hours,'(i3.3)')int(time_hours)
      write(*,*)'RCM: char_time_hours =',char_time_hours
!     if(int(time_hours) < 10)then
!             char_time_hours = '00'//adjustr(char_time_hours)
!     endif
!     write(*,*)' char_time_hours =',char_time_hours
      write(char_time_minutes,'(i2.2)')int(time_minutes)
      write(char_time_seconds,'(i2.2)')int(time_seconds)
!     write(*,*)' char_time_seconds =',char_time_seconds
       
      hours = adjustr(char_time_hours) //':'//char_time_minutes//':' &
      //adjustl(char_time_seconds)

      return
      end subroutine min2hr

!-----------------------------------------------

      subroutine print_date_time(LUN)
! prints date and time, 6/08, frt (from rcm)      
!      USE Rcm_mod_subs, ONLY : iprec
      USE rcm_precision, ONLY : iprec
      IMPLICIT NONE
      CHARACTER(LEN=8) :: real_date
      CHARACTER(LEN=10) ::real_time
      INTEGER(iprec), INTENT(IN) :: LUN

      CALL Date_and_time (real_date, real_time)
      WRITE (LUN,'(A11,A4,A1,A2,A1,A2, A8,A2,A1,A2,A1,A2)') &
              '  Today is ', real_date(1:4), '/', &
                             real_date(5:6), '/', &
                             real_date(7:8), & 
                  ' Time: ', real_time(1:2), ':', &
                             real_time(3:4), ':', &
                             real_time(5:6)
      return
      end subroutine print_date_time



