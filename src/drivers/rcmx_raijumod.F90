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

    use planethelper

    implicit none

    character(len=strLen) :: RunID
    character(len=strLen) :: XMLStr
    type(XML_Input_T) :: inpXML
    logical :: doRestart
    integer(iprec) :: i,j,n, status
    real(rprec) :: Lvalue
    real(rprec) ::mdipole = 3.0e-5 ! dipole moment in T
    real(rprec) :: mhdtime
    integer(iprec) :: rcmbndy 

    real(rprec) :: Lpeak = 4.0 ! location of Pressure max
    real(rprec) :: dL = 0.625 ! width of Pressure max
    real(rprec) :: D0 = 10.0 ! Density at pressure max
    real(rprec) :: T0 = 30.0 ! Temperature at pressure max
    real(rprec) :: D, P  ! Used in loop to set at given location
    !real(rprec), parameter :: pmax = 5.0e-8 ! pressure max in Pa
    !real(rprec), parameter :: pmin = 1.0e-11 ! min BG pressure in Pa
    !real(rprec), parameter :: nmax = 1.0e7 ! dens in ple/m^3
    !real(rprec), parameter :: nmin = 1.0e4 ! min dens in ple/m^3
    !real(rprec), parameter :: potmax = 5.0e4 ! potential max
    real(rprec) :: mhd_time_start
    real(rprec) :: mhd_time_end 
    real(rprec) :: mhd_dt 
    real(rprec) :: sc,sg
    real(rprec) :: colat_boundary
    real(rprec) :: rcm_boundary_s,rcm_boundary_e
    type(rcm_mhd_T) :: RM

    type(planet_T) :: planet
    real(rp):: cpcp, gamma

    !Always start with fresh directory
    status = SYSTEM("rm -rf RCMFiles > /dev/null 2>&1")

    !Get some XML stuff
    call getIDeckStr(XMLStr)
    inpXML = New_XML_Input(trim(XMLStr),"Kaiju/RCM",.true.)
    call inpXML%Set_Val(RunID,"sim/runid","rcmx")
    RM%rcm_runid = trim(RunID)

    call inpXML%Set_Val(mhd_time_start,"time/T0"  ,0.0)
    call inpXML%Set_Val(mhd_time_end  ,"time/tFin",36000.0)
    call inpXML%Set_Val(mhd_dt,        "time/dt"  ,500.0)

    call inpXML%Set_Val(doRestart,"restart/doRes",.false.)

    ! Set planet and ionosphere radius for rid_torcm to use
    call getPlanetParams(planet, inpXML, doLoudO=.true.)
    RM%planet_radius = planet%rp_m  ! [m]
    RM%iono_radius   = planet%ri_m  ! [m]
    mdipole = planet%magMoment*G2T  ! [T]
    !RM%planet_radius = re  ! [m]
    !RM%iono_radius   = 6.5e6  ! [m]
    call inpXML%Set_Val(cpcp ,"SAprob/cpcp"   , 50.0)  ! cross polar cap potential in kV
    call inpXML%Set_Val(gamma,"SAprob/fShield", 1.0 )  ! Shielding factor. 1 = no shielding, 2-3 is more realistic
    call inpXML%Set_Val(Lpeak,"SAprob/L" , Lpeak)
    call inpXML%Set_Val(dL   ,"SAprob/dL", dL   )
    call inpXML%Set_Val(D0,"SAprob/D" , D0 )  ! [#/cc]
    call inpXML%Set_Val(T0,"SAprob/T" , T0 )  ! [keV]

    
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
   
    mhdtime = mhd_time_start 
    ! now run 
!    do mhdtime=mhd_time_start,mhd_time_end-mhd_dt,mhd_dt
    do while (mhdtime <= mhd_time_end)
        IF(.not.doRestart)then
            ! compute flux tube volume and other items to pass to the RCM
            do i=1,RM%nLat_ion
                do j=1,RM%nLon_ion
                    Lvalue = 1.0/sin(RM%gcolat(i))**2
                    RM%Vol(i,j) = 32./35.*Lvalue**4/mdipole
                    RM%X_bmin(i,j,1) = Lvalue*cos(RM%glong(j))*RM%planet_radius
                    RM%X_bmin(i,j,2) = Lvalue*sin(RM%glong(j))*RM%planet_radius
                    !RM%X_bmin(i,j,1) = Lvalue*cos(RM%glong(j))*re
                    !RM%X_bmin(i,j,2) = Lvalue*sin(RM%glong(j))*re
                    RM%X_bmin(i,j,3) = 0.0
                    RM%bmin(i,j) = mdipole/Lvalue**3
                    RM%iopen(i,j) =-1  ! declare closed
                    RM%beta_average(i,j) = 0.1
                    !RM%Pave(i,j) = pmax * exp(-(Lvalue-Lmax)**2) + pmin
                    !RM%Nave(i,j) = nmax * exp(-(Lvalue-Lmax)**2) + nmin
                    D = D0*exp(-abs(Lvalue-Lpeak)/dL)  ! [#/cc]
                    P = DkT2P(D, T0)  ! [nPa]
                    RM%Nave(i,j) = D*1e6  ! [#/m^3]
                    RM%Pave(i,j) = P*1e-9  ! [Pa]

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
            write(*,'(2(a,g14.4))')' Coldstarting rcm_mhd at time: ',mhdtime,' delta t=',mhd_dt
            call rcm_mhd(mhdtime,mhd_dt,RM,RCMCOLDSTART)
            doColdstart = .false.

            ! Baumjohann & Treumann potential
            call SetEspot_BT(RM, cpcp, gamma)
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

contains

    subroutine SetEspot_BT(RM, cpcp, gamma)
        !! Electrostatic potential described in "Basic Space Plasma Physics" - Baumjohann and Treumann
        !! (page 99, Eqs 5.14-5.15)
        type(rcm_mhd_T), intent(inout) :: RM
        real(rp), intent(in) :: cpcp, gamma

        integer :: i,j
        real(rp) :: Agamma, dy, L


        ! Get deltaY
        j = RM%nLon_ion/4
        dy = ( RM%X_bmin(1,j,2) - RM%X_bmin(RM%nLat_ion,j,2) ) / RM%planet_radius
        !write(*,*)RM%X_bmin(1,j,2)
        !write(*,*)RM%X_bmin(RM%nLat_ion,j,2)
        !write(*,*)dy
        Agamma = 0.5*cpcp*dy**(-1.0*gamma)*1.0D3  ! [V]
        
        do i=1,RM%nLat_ion
            do j=1,RM%nLon_ion
                L = 1.0/sin(RM%gcolat(i))**2
                RM%pot(i,j) =  -1*Agamma * L**gamma * sin(RM%glong(j))
            end do
        end do
        
    end subroutine SetEspot_BT


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