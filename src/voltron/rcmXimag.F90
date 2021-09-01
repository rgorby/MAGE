! A version of rcmimag which also calls sstimagLL and replaces RCM pressures with SST.
! Eventually, we want to merge the two pressures in some assimilative way.
! X stands for eXtended.
module rcmXimag
    use rcmimag
    use sstLLimag
    use mixdefs
    use mixgeom

    implicit none

    type, extends(innerMagBase_T) :: rcmXIMAG_T

        class(rcmIMAG_T), allocatable :: rcmApp
        class(empData_T), allocatable :: empApp
        type(mixGrid_T) :: rcmG
        type(Map_T) :: rcmMap  

        contains

        ! over-ride the base functions with RCM versions
        procedure :: doInit => initRCMX
        procedure :: doAdvance => advanceRCMX
        procedure :: doEval => evalRCMX
        procedure :: doIO => doRCMXIO
        procedure :: doRestart => doRCMXRestart

    end type

    contains 

    subroutine initRCMX(imag,iXML,isRestart,vApp)
        class(rcmXIMAG_T), intent(inout) :: imag
        type(XML_Input_T), intent(in) :: iXML
        logical, intent(in) :: isRestart
        type(voltApp_T), intent(inout) :: vApp

        allocate(rcmIMAG_T :: imag%rcmApp)
        allocate(empData_T :: imag%empApp)

        call imag%rcmApp%doInit(iXML,isRestart,vApp)
        call imag%empApp%doInit(iXML,isRestart,vApp)

        ! define rcm grid, store inside the rcmXIMAG class
        call rcmGrid(imag)
        ! set map (note, sstG is defined in empApp%doInit)
        call mix_set_map(imag%empApp%sstG,imag%rcmG,imag%rcmMap)

        contains

        ! note this is ugly, as it simply reuses the code from rcm_mix_interface
        ! FIXME: consider defining RCM mix-style grid in, e.g., initRCM and reusing here and in rcm_mix_interface
        subroutine rcmGrid(imagX) 
            class(rcmXIMAG_T), intent(inout) :: imagX
            integer :: i, j, Np, Nt
            real(rp), dimension(:,:), allocatable :: rcmp, rcmt ! remix-style 2-D arrays to hold the RCM grid

            Np = size(imagX%rcmApp%rcmCpl%glong)
            Nt = size(imagX%rcmApp%rcmCpl%gcolat)
               
        !Now do remix mapping
            if (.not.allocated(rcmp)) allocate(rcmp(Np,Nt))
            if (.not.allocated(rcmt)) allocate(rcmt(Np,Nt))

        ! construct the 2-D grid
            do j=1,Np
               rcmt(j,:) = imagX%rcmApp%rcmCpl%gcolat
            enddo

            do i=1,Nt
               rcmp(:,i) = imagX%rcmApp%rcmCpl%glong
            enddo

            ! call remix grid constructor
            call init_grid_fromTP(imagX%rcmG,rcmt,rcmp,isSolverGrid=.false.)

        end subroutine rcmGrid

    end subroutine initRCMX

    subroutine advanceRCMX(imag,vApp,tAdv)
        class(rcmXIMAG_T), intent(inout) :: imag
        type(voltApp_T), intent(inout) :: vApp
        real(rp), intent(in) :: tAdv
        real(rp), dimension(:,:), allocatable :: empPressureOnRCMGrid

        call imag%rcmApp%doAdvance(vApp,tAdv)
        call imag%empApp%doAdvance(vApp,tAdv)

        ! interpolate from emp to rcm here 
        call mix_map_grids(imag%rcmMap,imag%empApp%sstP,empPressureOnRCMGrid)

        ! replace RCM pressure for now but think about merging like this later
        !rcmX%Pressure = w1(x,y)*rcm%Pressure + w2(x,y)*sst%Pressure

        ! note, converting sst pressure (nPa) to rcm (Pa)
        ! doEval below converts back to nPa
        imag%rcmApp%rcmCpl%Prcm = 1.0e-9*transpose(empPressureOnRCMGrid)

        ! Manipulate "RCM's" density to be some combination of Nmhd and Npsph
        call setRCMXDensity(imag%rcmApp%rcmCpl, 2)


    end subroutine advanceRCMX

    subroutine evalRCMX(imag,x1,x2,t,imW,isEdible)
        class(rcmXIMAG_T), intent(inout) :: imag
        real(rp), intent(in) :: x1,x2,t
        real(rp), intent(out) :: imW(NVARIMAG)
        logical, intent(out) :: isEdible

        call imag%rcmApp%doEval(x1,x2,t,imW,isEdible)

    end subroutine evalRCMX

!IO wrappers -- just do RCM things
    subroutine doRCMXIO(imag,nOut,MJD,time)
        class(rcmXIMAG_T), intent(inout) :: imag
        integer, intent(in) :: nOut
        real(rp), intent(in) :: MJD,time

        call imag%rcmApp%doIO(nOut,MJD,time)
    end subroutine doRCMXIO

    subroutine doRCMXRestart(imag,nRes,MJD,time)
        class(rcmXIMAG_T), intent(inout) :: imag
        integer, intent(in) :: nRes
        real(rp), intent(in) :: MJD, time

        call imag%rcmApp%doRestart(nRes,MJD,time)
    end subroutine doRCMXRestart

    subroutine setRCMXDensity(rcmCpl,option)
        class(rcm_mhd_T), intent(inout) :: rcmCpl
        integer, intent(in) :: option

        integer :: i,j

        select case (option)
        case(1)
            ! Use only whatever's in MHD
            rcmCpl%Nrcm = rcmCpl%Nave
        case(2)
            ! Wherever Npsph>Nave, use Npsph
            ! Else, use Nave

            !$OMP PARALLEL DO default(shared) collapse(2) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do j=1,rcmCpl%nLon_ion
                do i=1,rcmCpl%nLat_ion
                    
                    if (rcmCpl%Npsph(i,j) > rcmCpl%Nave(i,j)) then
                        rcmCpl%Nrcm(i,j) = rcmCpl%Npsph(i,j)
                    else
                        rcmCpl%Nrcm(i,j) = rcmCpl%Nave(i,j)
                    endif
                enddo
            enddo

        case DEFAULT
            !Don't do anything, will use RCM's density
        end select

    end subroutine setRCMXDensity

end module rcmXimag
