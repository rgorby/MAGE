module rcm_mhd_mod

    use rcm_precision
    use Rcm_mod_subs
    use rcm_mhd_interfaces
    use torcm_mod
    use tomhd_mod
    use rcm_mhd_io
    use ionosphere_exchange, only : setupIon, tearDownIon
    use constants, ONLY : radius_earth_m, radius_iono_m
    use rice_housekeeping_module
    use rcm_timing_module
    use files
    use clocks

    implicit none

    contains

    subroutine rcm_mhd(mhdtime,mhdtimedt,RM,iflag,iXML)
    ! version to couple to gamera
    ! units are assumed to mks, except for distances which are in Re.
    ! iflag = 0 - setup arrays, read in parameters (RCMINIT)
    ! iflag = 1 - run rcm (RCMADVANCE)
    ! iflag = 2 - Restart RCM (RCMWRITERESTART)
    ! iflag = 10 - start from a cold start (RCMCOLDSTART) 
    ! iflag = -1 - stop, write out timing information (RCMWRITETIMING)
    ! iflag = -2 - Write restart (icontrol = 31337) (RCMWRITERESTART)
    ! iflag = -3 - Write H5 output (icontrol = 31338) (RCMWRITEOUTPUT)


    ! 2/19 frt

        implicit none
        type(XML_Input_T), intent(in), optional :: iXML
        type(rcm_mhd_T),intent(inout) :: RM
        real(rprec), intent(in) :: mhdtime,mhdtimedt
        integer(iprec), intent(in) :: iflag

        integer(iprec) :: ierr   !> Error code...
        integer(iprec) :: system !> This code makes a call to the
                               !! non-standard SYSTEM function.  Must
                               !! define as integer to avoid compile error.


        real(rprec) :: itimei !> RCM(...) param:  start time   sbaotime
        real(rprec) :: itimef !> RCM(...) param:  end time
        real(rprec) :: time0 = 0. ! coupling start time
        real(rprec) :: ircm_dt
        real(rprec) :: itimef_old = -1
        integer(iprec) :: nstep   !> RCM(...) param:  number of sub-time steps in program
        !integer(iprec) :: idt   !> RCM(...) param:  basic time step in program
        !real(rprec) :: idt1  !> RCM(...) param:  time step for
                              !! changing disk & write records
        !real(rprec) :: idt2  !> RCM(...) param:  time step for
                              !! writting formatted output

        real(rprec) :: t1, t2  !> Used for performance timing

        !> Model coupling variables
        integer(iprec) :: exchangeNum = 0
        logical :: isFirstExchange
        logical :: isLastExchange

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! bypass for now
        time0 = 0. ! FIXME set for now

!        IsCoupledExternally = .TRUE.  ! switch RCM to "coupled" mode before doing anything else

        !if (doRCMVerbose) write (*,'(TR1,A,L7)') 'Welcome to the RCM, IsCoupledExternally=', IsCoupledExternally
        
        ! setup rcm,time in integer format
        itimei = mhdtime   !floor(mhdtime-time0,iprec)
        itimef = mhdtime + mhdtimedt !floor(mhdtime + mhdtimedt-time0,iprec)
        ircm_dt = itimef - itimei
  
    ! finish up
        if(iflag==RCMWRITETIMING)then
          return ! do nothing
        end if

    ! Write restart file
        if (iflag==RCMWRITERESTART) then
            CALL Rcm (itimei, itimef, nstep, icontrol=ICONWRITERESTART,stropt=RM%rcm_runid,nslcopt=RM%RCM_nRes)
            return
        endif

    ! Write output slice
        if (iflag==RCMWRITEOUTPUT) then
            CALL Rcm (itimei, itimef, nstep, icontrol=ICONWRITEOUTPUT,stropt=RM%rcm_runid,nslcopt=RM%RCM_nOut)
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
            CALL Rcm (itimei, itimef, nstep, icontrol=0_iprec)

            call allocate_conversion_arrays (isize,jsize,kcsize)
            
            call Grid_torcm (HighLatBD,LowLatBD, 0.0_rprec, RM%planet_radius, RM%iono_radius,doLatStretch)  ! set up RCM ionospheric grid here
            ! Setup Ionosphere intermediate Grid by equating it to the RCM grid, without angular overlap:
            call setupIon(RM)

            CALL Rcm (itimei, itimef, nstep, icontrol=1_iprec)

            ! icontrol of 2 also needs the input xml file
            CALL Rcm (itimei, itimef, nstep, icontrol=2_iprec, iXML=iXML)

            if (iflag == RCMINIT) then
                exchangeNum = 0
            endif

        ! restart
            if (iflag == RCMRESTART) then

                !Read in HDF5 restart data
                CALL Rcm (itimei, itimef, nstep, icontrol=ICONRESTART,stropt=RM%rcm_runid,nslcopt=RM%RCM_nRes)
                exchangeNum = floor(itimef/(itimef-itimei)) ! Need to find another way of calculating exchangeNum

                return
            endif !restart

            return !We're done here

        end if !RCMINIT or RCMRESTART



        if(iflag==RCMADVANCE.or.iflag==RCMCOLDSTART) then ! run the rcm

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Determine exchange times...
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            isFirstExchange = (exchangeNum==0)

            if (doRCMVerbose) then
!                write(*,*)'-----------rcm_mhd: rec=',rec

                write(*,*) 'itimei = ', itimei
                write(*,*) 'exchangeNum = ', exchangeNum
                WRITE (*,'(//)')
                write (*,'(a,f12.4,a,f12.4,a,i4)') 'RCM: time=',itimei,'  time0=',time0, '  Delta_t[s]=',ircm_dt
                write (*,'(a,i6,a,f12.4)') 'RCM: _T_rcm[s] =', itimei, ' T_MHD=',mhdtime
                WRITE (*,'(//)')
            endif
         
            !idt = real(Idt_overwrite) ! RCM internal time step in seconds
            ! Frequency (in seconds) to change disk & write records
            ! idt1 = itimef - itimei

            !!!Ensure no problem w/ RCM's integer time
            !!idt must divide advance time
            !if ( (mod(idt1,idt)) /= 0) then
            !    write(*,*) 'RCM Integer Time Divisibility Error ...'
            !    stop
            !endif

            ! Frequency (in seconds) to write formatted output
            !idt2 = idt1

            ! now round to to fit the correct number rcm timesteps
            !itimef = itimei + idt *((itimef-itimei)/idt)
            nstep = nStep 
            itimef_old = itimef


            if (isFirstExchange) then ! Set RCM initial conditions on plasma:
                call rcm (itimei, itimef, nstep, icontrol=3_iprec)
            end if


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Import data from MIX
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call cpu_time(t1)
            if (doRCMVerbose) then
                write(6,'(a,f12.4,a,i6)')'RCM: calling torcm with  itimei=',itimei,' iflag=',iflag
                call print_date_time(6_iprec)
            endif

            call Tic("TORCM")
            call torcm(RM,itimei,ierr,iflag)
            call Toc("TORCM")

            if (ierr < 0 ) then
                write(*,*) 'RCM: error in torcm '
                call BlackBoxRCM(RM)
            endif
            exchangeNum = exchangeNum + 1
            call cpu_time(t2)

            if (doRCMVerbose) write(*,'(a,g14.4,a)')'RCM: torcm cpu time= ',t2-t1,' seconds'

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Advance RCM
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            call cpu_time(t1)
            if (doRCMVerbose) then 
                write(6,'(a,f12.4,a,f12.4,a,i5,a)')'RCM: call rcm at itimei =',itimei,' to itimef =',itimef,' dt=',ircm_dt, ' sec'
                call print_date_time(6_iprec)
            endif

            ! now run the rcm
            call Tic("xRCMx")
            call rcm (itimei, itimef, nstep, icontrol=4_iprec,stropt=RM%rcm_runid,nslcopt=RM%RCM_nOut)
            call Toc("xRCMx")
!            rec = rec + 1 ! update record after rcm has run

            call cpu_time(t2)
            if (doRCMVerbose) then
                write(*,'(a,f12.4,a)')'RCM_MHD:   rcm cpu time= ',t2-t1,' seconds'
                call print_date_time(6_iprec)
            endif

            ! Do not export data if this is both the first & last exchange.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Export data to MHD code
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (doRCMVerbose) write(6,'(a)')'RCM_MHD: calling tomhd '
            
            call cpu_time(t1)
            call Tic("TOMHD")
            call Tomhd (RM, ierr)
            call Toc("TOMHD")

            call cpu_time(t2)
            if (doRCMVerbose) then
                call print_date_time(6_iprec)
                write(*,*)'RCM: tomhd cpu time= ',t2-t1,' seconds'
            endif

            if (ierr < 0 ) then
                write(*,*) 'RCM: error in tomhd '
                call BlackBoxRCM(RM)
            endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        end if

        if (iflag==RCMRESTART)then ! stop
            call rcm (itimei,itimef,nstep,icontrol=5_iprec)
            !  call Finalize()    ! Matches Initialize() above
            call tearDownIon(RM) ! Matches setupIon() above
        end if

        return

    end subroutine rcm_mhd

    !Write black box and die
    subroutine BlackBoxRCM(RCMApp)
        type(rcm_mhd_t), intent(inout) :: RCMApp

    !TODO: Can add some console output here

    !Output last words
        RCMApp%rcm_runid = "CRASH" // trim(RCMApp%rcm_runid)
        call initRCMIO(RCMApp,isResO=.false.)
        call WriteRCM(RCMApp,0,0.0_rp,0.0_rp)
        RCMApp%rcm_nOut = 0
        call rcm_mhd(0.0_rp,TINY,RCMApp,RCMWRITEOUTPUT)

    !Die with dignity
        write(*,*) 'RCM Commiting suicide in 300s ...'
        call sleep(300) !Sleep before blowing up
        write(*,*) 'Goodbye cruel world'
        stop !Self destruct
    end subroutine BlackBoxRCM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine print_date_time(LUN)
        ! prints date and time, 6/08, frt (from rcm)
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

end module rcm_mhd_mod

