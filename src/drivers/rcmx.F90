! dummy mhd test program to drive the rcm_mhd code
! assumes dipole inputs
! 1/19 frt
program rcmx
    use kdefs
    use files
    use rcm_mhd_interfaces
    use rcm_mhd_mod, ONLY: rcm_mhd
    USE rcm_precision
    use rcm_mhd_io
    use strings
    use xml_input
    use rcm_mod_subs, ONLY: vm,bmin,xmin,ymin,zmin,v

    implicit none

    character(len=strLen) :: RunID
    character(len=strLen) :: XMLStr
    type(XML_Input_T) :: inpXML
    logical :: doRestart
    integer(iprec) :: i,j,n
    real(rprec) :: Lvalue
    real(rprec),parameter ::mdipole = 3.0e-5 ! dipole moment in T
    real(rprec) :: mhdtime
    integer(iprec) :: rcmbndy 
    real(rprec), parameter :: Lmax = 4.0 ! location of Pressure max
    real(rprec), parameter :: pmax = 5.0e-8 ! pressure max in Pa
    real(rprec), parameter :: pmin = 1.0e-11 ! min BG pressure in Pa
    real(rprec), parameter :: nmax = 1.0e7 ! dens in ple/m^3
    real(rprec), parameter :: nmin = 1.0e4 ! min dens in ple/m^3
    real(rprec), parameter :: potmax = 5.0e4 ! potential max
    real(rprec), parameter :: re = 6380.e3
    real(rprec) :: mhd_time_start
    real(rprec) :: mhd_time_end 
    real(rprec) :: mhd_dt 
    real(rprec) :: sc,sg
    real(rprec) :: colat_boundary
    real(rprec) :: rcm_boundary_s,rcm_boundary_e
    type(rcm_mhd_T) :: RM

    !Always start with fresh directory
    CALL SYSTEM("rm -rf RCMFiles > /dev/null 2>&1")

    !Get some XML stuff
    call getIDeckStr(XMLStr)
    inpXML = New_XML_Input(trim(XMLStr),"RCM",.true.)
    call inpXML%Set_Val(RunID,"sim/runid","rcmx")
    RM%rcm_runid = trim(RunID)

    call inpXML%Set_Val(mhd_time_start,"time/T0"  ,0.0)
    call inpXML%Set_Val(mhd_time_end  ,"time/tFin",36000.0)
    call inpXML%Set_Val(mhd_dt,        "time/dt"  ,500.0)

    call inpXML%Set_Val(doRestart,"restart/doRes",.false.)

    !Set planet and ionosphere radius for rid_torcm to use
    RM%planet_radius = re
    RM%iono_radius = 6.5e6
    if (doRestart) then
        call RCMRestartInfo(RM,inpXML,mhd_time_start,.true.)
        write(*,*) 'Restarting RCM @ t = ', mhd_time_start
        call rcm_mhd(mhd_time_start,mhd_dt,RM,RCMRESTART,iXML=inpXML)
        doColdstart = .false.
    else
        ! initialize
        call rcm_mhd(mhd_time_start,mhd_dt,RM,RCMINIT,iXML=inpXML)
    endif

    write(*,*) 'Start / End / dt = ', mhd_time_start,mhd_time_end,mhd_dt

    !Setup IO
    RM%rcm_nOut = 0
    call initRCMIO(RM,doRestart)

    !Set boundaries    
    rcm_boundary_s =35
    rcm_boundary_e =2
   
    mhdtime = mhd_time_start 
    ! now run 
!    do mhdtime=mhd_time_start,mhd_time_end-mhd_dt,mhd_dt
     do while (mhdtime <= mhd_time_end)
        IF(.not.doRestart)then
        rcmbndy = nint(rcm_boundary_s +&
            (rcm_boundary_e-rcm_boundary_s)*(mhdtime-mhd_time_start)/(mhd_time_end-mhd_time_start),iprec)
        write(*,'(a,g12.4,a,i5)')' At t =',mhdtime,' RCM boundary index =',rcmbndy
        colat_boundary = sin(RM%gcolat(rcmbndy))
        write(*,*)RM%nLat_ion,RM%nLon_ion
        ! compute flux tube volume and other items to pass to the RCM
        do i=1,RM%nLat_ion
            do j=1,RM%nLon_ion
                Lvalue = 1.0/sin(RM%gcolat(i))**2
                RM%Vol(i,j) =32./35.*Lvalue**4/mdipole
                RM%X_bmin(i,j,1) = Lvalue*cos(RM%glong(j))*re
                RM%x_bmin(i,j,2) = Lvalue*sin(RM%glong(j))*re
                RM%x_bmin(i,j,3) = 0.0
                RM%bmin(i,j) = mdipole/Lvalue**3
                RM%iopen(i,j) =-1  ! declare closed
                RM%beta_average(i,j) = 0.1
                RM%Pave(i,j) = pmax * exp(-(Lvalue-Lmax)**2) + pmin
                RM%Nave(i,j) = nmax * exp(-(Lvalue-Lmax)**2) + nmin
                ! add a potential that goes to zero near the inner boundary 5/20 frt
                sc = sin(colat_boundary)
                sg = sin(RM%gcolat(i))
                if(RM%gcolat(i) < colat_boundary)then
                    RM%pot(i,j) = -potmax/2.*sin(RM%glong(j))*sg/sc
                else
                    RM%pot(i,j) = -potmax/2.*sin(RM%glong(j))/sc/(1.-1./sc**2)*(sg-1./sg)
                end if

            end do
        end do

        ! set rcm boundary
        RM%Vol(1:rcmbndy,:) = -1.0
        RM%iopen(1:rcmbndy,:) = 1 ! declare open

        ELSE
        do i=1,RM%nLat_ion
            do j=1,RM%nLon_ion
               RM%Vol(i,j) = 1/abs(vm(i,j))**1.5 * sign(1.0d0,vm(i,j))*1.0e9
               RM%Bmin(i,j) = bmin(i,j)
               RM%X_bmin(i,j,1) = xmin(i,j)
               RM%X_bmin(i,j,2) = ymin(i,j)
               RM%X_bmin(i,j,3) = zmin(i,j)
               RM%iopen(i,j) = sign(vm(i,j),1.0d0)
               RM%pot(i,j) = v(i,j)
            end do
        end do
        END IF ! restart

        if (doColdstart)then
            write(*,'(2(a,g14.4))')' calling rcm_mhd at time: ',mhdtime,' delta t=',mhd_dt
            call rcm_mhd(mhdtime,mhd_dt,RM,RCMCOLDSTART)
            doColdstart = .false.
        else
            write(*,'(2(a,g14.4))')' calling rcm_mhd at time: ',mhdtime,' delta t=',mhd_dt
            call rcm_mhd(mhdtime,mhd_dt,RM,RCMADVANCE)
        end if
        
        !Commenting out old-style output for now
        !call write_2d(RM,mhdtime+mhd_dt) ! write out results

        call WriteRCM(RM,RM%rcm_nOut,mhdtime,mhdtime)
        write(*,*) 'Total pressure = ', sum(RM%Prcm)
        RM%rcm_nOut = RM%rcm_nOut+1

        mhdtime = mhdtime + mhd_dt
    end do 

    ! done now close out
    call rcm_mhd(mhdtime,mhd_dt,RM,RCMWRITETIMING)

    stop
end program rcmx

subroutine write_2d(RM,time)

    use rcm_mhd_interfaces
    use rcm_mhd_mod, ONLY: rcm_mhd
    USE Rcm_mod_subs, ONLY : iprec,rprec
    type(rcm_mhd_T),intent(in) :: RM
    real(rprec),intent(in) :: time
    character(len=15) :: fileout
    character(len=5) :: ctime
    integer(iprec) ::i,j

    write (ctime, '(i5.5)')int(time,iprec) 

    fileout = adjustr('tomhd') //ctime// '.dat'
    write(*,*)' writing file =',fileout

    open(unit=10,file=fileout,status='unknown')

    write(10,*)time
    write(10,*)RM%nLat_ion 
    write(10,*)RM%nLon_ion 

    do i=1,RM%nLat_ion 
        do j=1,RM%nLon_ion 
            write(10,'(4(g14.6,1x))')RM%X_bmin(i,j,1), RM%X_bmin(i,j,2), RM%Prcm(i,j), RM%Nrcm(i,j)
        end do
    end do

    close(10)
    return
end subroutine write_2d



