! A version of rcmimag which also calls sstimagLL and replaces RCM pressures with SST.
! Eventually, we want to merge the two pressures in some assimilative way.
! X stands for eXtended.
module rcmXimag
    use rcmimag
    use sstLLimag
    use mixdefs
    use mixgeom
    USE Rcm_mod_subs, ONLY : isize, jsize
    use rcmdefs, only : RCMTOPOPEN,RCMTOPCLOSED

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
        real(rp), dimension(:,:), allocatable :: empBvolOnRCMGrid

        call imag%rcmApp%doAdvance(vApp,tAdv)
        call imag%empApp%doAdvance(vApp,tAdv)

        ! interpolate from emp to rcm here 
        !call mix_map_grids(imag%rcmMap,imag%empApp%sstP,empPressureOnRCMGrid)
        !call mix_map_grids(imag%rcmMap,imag%empApp%sstBvol,empBvolOnRCMGrid)
        

        ! replace RCM pressure for now but think about merging like this later
        !rcmX%Pressure = w1(x,y)*rcm%Pressure + w2(x,y)*sst%Pressure

        ! note, converting sst pressure (nPa) to rcm (Pa)
        ! doEval below converts back to nPa
        !imag%rcmApp%rcmCpl%Prcm = 1.0e-9*transpose(empPressureOnRCMGrid)
        
        ! Set RCM pressure via rcm(pV^gamma)=sst(pV^gamma)

        call setPressViaEntropy(imag)


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
        !Hijack mhdrcm file and include SST information
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

    subroutine setPressViaEntropy(imag)
        class(rcmXIMAG_T), intent(inout) :: imag

        real(rp), dimension(:,:), allocatable :: empIOpenOnRCMGrid, empTpIo
        real(rp), dimension(:,:), allocatable :: empPressureOnRCMGrid, empTpP
        real(rp), dimension(:,:), allocatable :: empBvolOnRCMGrid, empTpBvol
        real(rp) :: gamma = 5./3.
        integer :: i,j

        call mix_map_grids(imag%rcmMap,imag%empApp%Iopen,empIOpenOnRCMGrid)
        call mix_map_grids(imag%rcmMap,imag%empApp%sstP,empPressureOnRCMGrid)
        call mix_map_grids(imag%rcmMap,imag%empApp%sstBvol,empBvolOnRCMGrid)
        empTpIo = transpose(empIOpenOnRCMGrid)
        empTpP = transpose(empPressureOnRCMGrid)
        empTpBvol = transpose(empBvolOnRCMGrid)

        DO i = 1,isize
            DO j = 1,jsize

                if (empTpIo(i,j) > (-1+TINY)) then  ! If interpolated point is even somewhat influenced by an open line, kill it
                    ! Make sure mhd won't ingest this point
                    ! Probably only need to set one of these but idk where we are in RCMEval pipeline so set both to be safe
                    imag%rcmApp%rcmCpl%iopen(i,j) = 1
                    imag%rcmApp%rcmCpl%toMHD(i,j) = .false.
                else  ! Only other option here is that its a closed line (-1) according to SST
                    if (imag%rcmApp%rcmCpl%iopen(i,j) /= RCMTOPOPEN) then  ! Only overwrite closed and buffer region
                        imag%rcmApp%rcmCpl%Prcm(i,j) = 1.0e-9*empTpP(i,j)  &
                                    *empTpBvol(i,j)*1.0e9**gamma   &
                                    *imag%rcmApp%rcmCpl%Vol(i,j)**(-gamma)
                    end if
                end if
            END DO
        END DO

    end subroutine setPressViaEntropy

    subroutine doSSTIO(imag)
        class(rcmXIMAG_T), intent(in) :: imag

    end subroutine doSSTIO

end module rcmXimag
