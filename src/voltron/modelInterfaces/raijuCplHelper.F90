module raijuCplHelper

    use volttypes
    use raijutypes
    use remixReader
    use shellinterp
    use shellGridIO
    use shellUtils
    use ebtypes
    use raijugrids
    use ioh5
    use files
    
    use imagtubes
    use mixdefs
    

    implicit none

    integer, private, parameter :: MAXIOVARS = 50

    contains

    subroutine raijuCpl_init(raiCpl, inpXML)
        class(raijuCoupler_T), intent(inout) :: raiCpl
        type(XML_Input_T), intent(inout) :: inpXML

        type(XML_Input_T) :: iXML
        character(len=strLen) :: tmpStr
        integer, dimension(4) :: shGhosts
        integer :: i

        ! Make sure root is Kaiju/raiju
        call inpXML%GetFileStr(tmpStr)
        ! Create new XML reader w/ RAIJU as root
        iXML = New_XML_Input(trim(tmpStr),'Kaiju/RAIJU',.true.)

        ! Options
        call iXML%Set_Val(raiCpl%startup_blendTscl, "cpl/startupTscl", raiCpl%startup_blendTscl)
        call iXML%Set_Val(raiCpl%doColdstartCX,'prob/coldstartCX',raiCpl%doColdstartCX)

        ! Allocations
        associate(sh => raiCpl%raiApp%Grid%shGrid, nFluidIn => raiCpl%raiApp%Model%nFluidIn)


            ! Shell Grid inits
            shGhosts(NORTH) = sh%Ngn
            shGhosts(SOUTH) = sh%Ngs
            shGhosts(EAST)  = sh%Nge
            shGhosts(WEST)  = sh%Ngw
            call raijuGenGridFromShGrid(raiCpl%shGr, raiCpl%opt%voltGrid, iXML)
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%pot_total)
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%pot_corot)
            call initShellVar(raiCpl%shGr, SHGR_CC, raiCpl%bvol_cc)

            do i=1,NDIM
                call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%Bmin(i))
                call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%xyzMin(i))
                call initShellVar(raiCpl%shGr, SHGR_CC    , raiCpl%xyzMincc(i))
            enddo
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%thcon)
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%phcon)
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%bVol)
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%topo)
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%vaFrac)

            allocate(raiCpl%Pavg(0:nFluidIn))
            allocate(raiCpl%Davg(0:nFluidIn))
            allocate(raiCpl%Pstd(0:nFluidIn))
            allocate(raiCpl%Dstd(0:nFluidIn))
            do i=0,nFluidIn
                call initShellVar(raiCpl%shGr, SHGR_CC, raiCpl%Pavg(i))
                call initShellVar(raiCpl%shGr, SHGR_CC, raiCpl%Davg(i))
                call initShellVar(raiCpl%shGr, SHGR_CC, raiCpl%Pstd(i))
                call initShellVar(raiCpl%shGr, SHGR_CC, raiCpl%Dstd(i))
            enddo
            call initShellVar(raiCpl%shGr, SHGR_CC, raiCpl%tiote)
            call initShellVar(raiCpl%shGr, SHGR_CC, raiCpl%Tb)
        end associate
        
        ! Initial values
        raiCpl%tLastUpdate = -1.0*HUGE
        raiCpl%pot_total%data = 0.0
        raiCpl%pot_total%mask = .true.
        raiCpl%pot_corot%data = 0.0
        raiCpl%pot_corot%mask = .true.

    end subroutine raijuCpl_init


    subroutine tubeShell2RaiCpl(voltGrid, tubeShell, raiCpl)
        !! Takes voltron tubeShell data and stores what we need into our coupler
        !! Using the shell interp, we are mapping some quantities from corners to centers as expected by raijuApp
        type(ShellGrid_T), intent(in) :: voltGrid
        type(TubeShell_T), intent(in) :: tubeShell
        class(raijuCoupler_T), intent(inout) :: raiCpl

        type(ShellGridVar_T) :: tmpTopo
        logical, dimension(tubeShell%topo%isv:tubeShell%topo%iev,tubeShell%topo%jsv:tubeShell%topo%jev) :: topoSrcMask
        integer :: i,j,s

        where (tubeShell%topo%data == TUBE_CLOSED)
            topoSrcMask = .true.
        elsewhere
            topoSrcMask = .false.
        end where

        call initShellVar(raiCpl%shGr, SHGR_CORNER, tmpTopo)

        ! Corners
        do i=1,NDIM
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%X_bmin(i), raiCpl%shGr, raiCpl%xyzMin(i))
        enddo
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%bmin, raiCpl%shGr, raiCpl%Bmin(ZDIR))
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%latc, raiCpl%shGr, raiCpl%thcon)
        raiCpl%thcon%data = PI/2 - raiCpl%thcon%data
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%lonc, raiCpl%shGr, raiCpl%phcon)
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%bVol, raiCpl%shGr, raiCpl%bvol)
        !call InterpShellVar_ParentToChild(voltGrid, tubeShell%bVol, raiCpl%shGr, raiCpl%bvol)
        !call InterpShellVar_TSC_SG(voltGrid, tubeShell%bVol, raiCpl%shGr, raiCpl%bvol_cc)
        do j=raiCpl%shGr%jsg,raiCpl%shGr%jeg
            do i=raiCpl%shGr%isg,raiCpl%shGr%ieg
                raiCpl%bvol_cc%data(i,j) = toCenter2D(raiCpl%bvol%data(i:i+1,j:j+1))
            enddo
        enddo
        
        ! Get topo and then convert to RAIJU's definition
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%topo, raiCpl%shGr, tmpTopo)
        where (abs(tmpTopo%data - TUBE_CLOSED) < TINY)
            raiCpl%topo%data = RAIJUCLOSED
        elsewhere
            raiCpl%topo%data = RAIJUOPEN
        end where
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%wMAG, raiCpl%shGr, raiCpl%vaFrac)

        ! Centers
        do s=0,raiCpl%raiApp%Model%nFluidIn
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%avgP(s), raiCpl%shGr, raiCpl%Pavg(s), srcMaskO=topoSrcMask)
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%avgN(s), raiCpl%shGr, raiCpl%Davg(s), srcMaskO=topoSrcMask)
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%stdP(s), raiCpl%shGr, raiCpl%Pstd(s), srcMaskO=topoSrcMask)
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%stdN(s), raiCpl%shGr, raiCpl%Dstd(s), srcMaskO=topoSrcMask)
            !call InterpShellVar_TSC_SG(voltGrid, tubeShell%avgP(s), raiCpl%shGr, raiCpl%Pavg(s))
            !call InterpShellVar_TSC_SG(voltGrid, tubeShell%avgN(s), raiCpl%shGr, raiCpl%Davg(s))
            !call InterpShellVar_TSC_SG(voltGrid, tubeShell%stdP(s), raiCpl%shGr, raiCpl%Pstd(s))
            !call InterpShellVar_TSC_SG(voltGrid, tubeShell%stdN(s), raiCpl%shGr, raiCpl%Dstd(s))
        enddo
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%TioTe0, raiCpl%shGr, raiCpl%tiote)
        do i=1,NDIM
            call InterpShellVar_TSC_SG(raiCpl%shGr, raiCpl%xyzMin(i), raiCpl%shGr, raiCpl%xyzMincc(i), srcMaskO=topoSrcMask)
            !call InterpShellVar_TSC_SG(raiCpl%shGr, raiCpl%xyzMin(i), raiCpl%shGr, raiCpl%xyzMincc(i))
        enddo
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%Tb, raiCpl%shGr, raiCpl%Tb, srcMaskO=topoSrcMask)
        !call InterpShellVar_TSC_SG(voltGrid, tubeShell%Tb, raiCpl%shGr, raiCpl%Tb)
        

    end subroutine tubeShell2RaiCpl


    subroutine raiCpl2RAIJU(raiCpl)
        !! Take info at raijuCoupler level and put it into RAIJU proper
        !! raiCpl should have everything in the sizes we expect, so just need to do a bunch of copies
        class(raijuCoupler_T), intent(inout) :: raiCpl

        integer :: i,j,k
    
        associate(Model=>raiCpl%raiApp%Model, State=>raiCpl%raiApp%State, shGr=>raiCpl%shGr)

        ! Reset quantities we will only selectively write to
        State%Pavg = 0
        State%Davg = 0
        State%Pstd = 0
        State%Dstd = 0
        State%bvol_cc = 0
        State%Tb%data = 0

        !--- Ionosphere data ---!
        State%espot(:,:)     = raiCpl%pot_total%data(:,:) ! They live on the same grid so this is okay
        State%pot_corot(:,:) = raiCpl%pot_corot%data(:,:)


        !--- Mag data ---!

        ! Defaults
        State%Tb%data = HUGE

        State%topo = raiCpl%topo%data

        ! Copy no matter the topo value
        do i=1,NDIM
            State%xyzMin  (:,:,i) = raiCpl%xyzMin(i)%data
            State%xyzMincc(:,:,i) = raiCpl%xyzMincc(i)%data
        enddo
        State%thcon(:,:)     = raiCpl%thcon%data
        State%phcon(:,:)     = raiCpl%phcon%data
        State%Bmin(:,:,ZDIR) = raiCpl%bmin(ZDIR)%data
        State%bvol(:,:)      = raiCpl%bvol%data
        State%vaFrac(:,:)    = raiCpl%vaFrac%data

        ! Now only copy for good points
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j)
        do j=shGr%jsg,shGr%jeg
            do i=shGr%isg,shGr%ieg

                if (any(State%topo(i:i+1,j:j+1) .eq. RAIJUOPEN)) then
                    cycle
                endif

                do k=0,Model%nFluidIn
                    State%Pavg(i,j,k) = raiCpl%Pavg(k)%data(i,j)
                    State%Davg(i,j,k) = raiCpl%Davg(k)%data(i,j)
                    State%Pstd(i,j,k) = raiCpl%Pstd(k)%data(i,j) / max(State%Pavg(i,j,k), TINY)  ! Normalize
                    State%Dstd(i,j,k) = raiCpl%Dstd(k)%data(i,j) / max(State%Davg(i,j,k), TINY)
                enddo

                State%bvol_cc(i,j) = raiCpl%bvol_cc%data(i,j)
                State%tiote(i,j) = raiCpl%tiote%data(i,j)
                State%Tb%data(i,j) = raiCpl%Tb%data(i,j)
            enddo
        enddo

        end associate
    end subroutine

!------
! One-way driving from file helpers
!------

    !> This function takes updated model states and does the operations
    !> necessary to update the cplBase%fromV object
    subroutine packRaijuCoupler_OWD(raiCpl, vApp, rmReader)
        class(raijuCoupler_T), intent(inout) :: raiCpl
        type(voltApp_T), intent(inout) :: vApp
        type(rmReader_T) :: rmReader

        ! Update coupling time
        raiCpl%tLastUpdate = raiCpl%raiApp%State%t

        ! Using chimp, populate imagTubes
        call genImagTubes(raiCpl, vApp)
        ! Set potential
        call InterpShellVar_TSC_SG(rmReader%shGr, rmReader%nsPot(1), raiCpl%shGr, raiCpl%pot_total)

    end subroutine packRaijuCoupler_OWD


    subroutine genImagTubes(raiCpl, vApp)
        class(raijuCoupler_T), intent(inout) :: raiCpl
        type(voltApp_T), intent(in   ) :: vApp

        integer :: i,j
        real(rp) :: seedR, eqR
        type(magLine_T) :: magLine
        ! Get field line info and potential from voltron
        ! And put the data into RAIJU's fromV coupling object

        associate(sh=>raiCpl%raiApp%Grid%shGrid , &
            planet=>raiCpl%raiApp%Model%planet, &
            ebApp =>vApp%ebTrcApp)

            seedR =  planet%ri_m/planet%rp_m
            ! Do field line tracing, populate fromV%ijTubes
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,eqR)
            do i=sh%isg,sh%ieg+1
                do j=sh%jsg,sh%jeg+1
                    call CleanLine(raiCpl%magLines(i,j))

                    eqR = DipColat2L(raiCpl%raiApp%Grid%thRp(i))  ! Function assumes colat coming from 1 Rp, make sure we use the right theta value
                    if (eqR < raiCpl%opt%mhdRin) then
                        call DipoleTube(vApp, sh%th(i), sh%ph(j), raiCpl%ijTubes(i,j))
                    else
                        call MHDTube(ebApp, planet,   & !ebTrcApp, planet
                            sh%th(i), sh%ph(j), seedR, &  ! colat, lon, r
                            raiCpl%ijTubes(i,j), raiCpl%magLines(i,j), &  ! IMAGTube_T, magLine_T
                            doShiftO=.true.)
                    endif

                enddo
            enddo
        end associate

    end subroutine genImagTubes


!------
! Real-time coupling stuff
!------

    subroutine packRaijuCoupler_RT(raiCpl, vApp)
        !! Temporary imagTube generator for realtime coupling
        !! Eventually, someone else should probably be packing raiCpl objects for us
        class(raijuCoupler_T), intent(inout) :: raiCpl
        class(voltApp_T), intent(in) :: vApp

        raiCpl%tLastUpdate = vApp%time

        !call genImagTubes(raiCpl, vApp)
        call tubeShell2RaiCpl(vApp%shGrid, vApp%State%tubeShell, raiCpl)
        !call mixPot2Raiju_RT(raiCpl, vApp%remixApp)
        
        call InterpShellVar_ParentToChild(vApp%shGrid, vApp%State%potential_total, raiCpl%shGr, raiCpl%pot_total)
        call InterpShellVar_ParentToChild(vApp%shGrid, vApp%State%potential_corot, raiCpl%shGr, raiCpl%pot_corot)
        
    end subroutine


!------
! Coupler I/O
!------

    subroutine writeRaiCplRes(raiCpl, nRes)
        class(raijuCoupler_T), intent(in) :: raiCpl
        integer, intent(in) :: nRes

        character(len=strLen) :: ResF,lnResF !Name of restart file
        !logical :: fExist
        type(IOVAR_T), dimension(MAXIOVARS) :: IOVars

        write (ResF  , '(A,A,I0.5,A)') trim(raiCpl%raiApp%Model%RunID), ".raiCpl.Res.", nRes   , ".h5"
        write (lnResF, '(A,A,A,A)'   ) trim(raiCpl%raiApp%Model%RunID), ".raiCpl.Res.", "XXXXX", ".h5"
        call CheckAndKill(ResF)     
        
        call writeShellGrid(raiCpl%shGr, ResF,"/ShellGrid")

        call ClearIO(IOVars)
        call AddOutVar(IOVars, "tLastUpdate", raiCpl%tLastUpdate, uStr="s")
        call AddOutSGV(IOVars, "Pavg"       , raiCpl%Pavg       , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "Davg"       , raiCpl%Davg       , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "Pstd"       , raiCpl%Pstd       , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "Dstd"       , raiCpl%Dstd       , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "Bmin"       , raiCpl%Bmin       , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "xyzMin"     , raiCpl%xyzMin     , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "xyzMincc"   , raiCpl%xyzMincc   , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "topo"       , raiCpl%topo       , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "thcon"      , raiCpl%thcon      , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "phcon"      , raiCpl%phcon      , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "bvol"       , raiCpl%bvol       , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "bvol_cc"    , raiCpl%bvol_cc    , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "vaFrac"     , raiCpl%vaFrac     , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "Tb"         , raiCpl%Tb         , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "pot_total"  , raiCpl%pot_total  , doWriteMaskO=.true.)
        call AddOutSGV(IOVars, "pot_corot"  , raiCpl%pot_corot  , doWriteMaskO=.true.)
        call WriteVars(IOVars, .false., ResF)
        call MapSymLink(ResF,lnResF)
    end subroutine


    subroutine readRaiCplRes(raiCpl, resId, nRes)
        class(raijuCoupler_T), intent(inout) :: raiCpl
        character(len=*), intent(in) :: resId
        integer, intent(in) :: nRes

        character(len=strLen) :: ResF, nStr
        type(IOVAR_T), dimension(MAXIOVARS) :: IOVars

        if (nRes == -1) then
            nStr = "XXXXX"
        else
            write(nStr,'(I0.5)') nRes
        endif
        write(ResF,'(A,A,A,A)')trim(resId),".raiCpl.Res.",trim(nStr),".h5"

        call GenShellGridFromFile(raiCpl%shGr, "RAIJU", ResF,"/ShellGrid")

        call ClearIO(IOVars)
        call AddInVar(IOVars, "tLastUpdate",vTypeO=IOREAL)
        call ReadVars(IOVars, .false., ResF)
        raiCpl%tLastUpdate = GetIOReal(IOVars, "tLastUpdate")

        ! ShellGridVars
        call ReadInSGV(raiCpl%Pavg     ,ResF, "Pavg"     , doIOpO=.false.)
        call ReadInSGV(raiCpl%Davg     ,ResF, "Davg"     , doIOpO=.false.)
        call ReadInSGV(raiCpl%Pstd     ,ResF, "Pstd"     , doIOpO=.false.)
        call ReadInSGV(raiCpl%Dstd     ,ResF, "Dstd"     , doIOpO=.false.)
        call ReadInSGV(raiCpl%Bmin     ,ResF, "Bmin"     , doIOpO=.false.)
        call ReadInSGV(raiCpl%xyzMin   ,ResF, "xyzMin"   , doIOpO=.false.)
        call ReadInSGV(raiCpl%xyzMincc ,ResF, "xyzMincc" , doIOpO=.false.)
        call ReadInSGV(raiCpl%topo     ,ResF, "topo"     , doIOpO=.false.)
        call ReadInSGV(raiCpl%thcon    ,ResF, "thcon"    , doIOpO=.false.)
        call ReadInSGV(raiCpl%phcon    ,ResF, "phcon"    , doIOpO=.false.)
        call ReadInSGV(raiCpl%bvol     ,ResF, "bvol"     , doIOpO=.false.)
        call ReadInSGV(raiCpl%bvol_cc  ,ResF, "bvol_cc"  , doIOpO=.false.)
        call ReadInSGV(raiCpl%vaFrac   ,ResF, "vaFrac"   , doIOpO=.false.)
        call ReadInSGV(raiCpl%Tb       ,ResF, "Tb"       , doIOpO=.false.)
        call ReadInSGV(raiCpl%pot_total,ResF, "pot_total", doIOpO=.false.)
        call ReadInSGV(raiCpl%pot_corot,ResF, "pot_corot", doIOpO=.false.)

    end subroutine

end module raijuCplHelper