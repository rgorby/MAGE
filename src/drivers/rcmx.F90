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

    implicit none
    character(len=strLen) :: RunID
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
    real(rprec) :: colat_boundary
    real(rprec) :: rcm_boundary_s,rcm_boundary_e
    type(rcm_mhd_T) :: RM
    ! intialize transfer arrays

    !Always start with fresh directory
    CALL SYSTEM("rm -rf RCMFiles > /dev/null 2>&1")

    RunID = "rcmx"
    mhd_time_start = 0.0
    mhd_time_end   = 600.0
    mhd_dt = 5.0

    write(*,*) 'Start / End / dt = ', mhd_time_start,mhd_time_end,mhd_dt

    !write(*,'(a,$)')' input MHD time start, end and dt: '
    !read(5,*)mhd_time_start,mhd_time_end,mhd_dt

    ! initialize
    call rcm_mhd(mhd_time_start,mhd_dt,RM,RCMINIT)

    !Setup IO
    RM%rcm_runid = trim(RunID)
    RM%rcm_nOut = 0
    call initRCMIO(RM)

    !Set boundaries    
    rcm_boundary_s =35
    rcm_boundary_e =2
    
    ! now run 
    do mhdtime=mhd_time_start,mhd_time_end-mhd_dt,mhd_dt
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
                if(RM%gcolat(i) < colat_boundary)then
                    RM%pot(i,j) = -potmax/2.*sin(RM%glong(j))*sin(RM%gcolat(i))
                else
                    RM%pot(i,j) = -potmax/2.*sin(RM%glong(j))*sin(colat_boundary)/sin(RM%gcolat(i))
                end if
            end do
        end do

        ! set rcm boundary
        RM%Vol(1:rcmbndy,:) = -1.0
        RM%iopen(1:rcmbndy,:) = 1 ! declare open

        write(*,'(2(a,g14.4))')' calling rcm_mhd at time: ',mhdtime,' delta t=',mhd_dt
        call rcm_mhd(mhdtime,mhd_dt,RM,RCMADVANCE)
        call write_2d(RM,mhdtime+mhd_dt) ! write out results

        call WriteRCM(RM,RM%rcm_nOut,mhdtime,mhdtime)
        RM%rcm_nOut = RM%rcm_nOut+1
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



