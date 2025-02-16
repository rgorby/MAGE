module raijuCplHelper

    use volttypes
    use raijutypes
    use remixReader
    use shellinterp
    use ebtypes
    use raijugrids
    
    use imagtubes
    use mixdefs
    

    implicit none

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

        associate(sh => raiCpl%raiApp%Grid%shGrid, nFluidIn => raiCpl%raiApp%Model%nFluidIn)

            ! Allocations
            !allocate(raiCpl%magLines (sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))
            !allocate(raiCpl%ijTubes( sh%isg:sh%ieg+1, sh%jsg:sh%jeg+1))

            ! Shell Grid inits
            shGhosts(NORTH) = sh%Ngn
            shGhosts(SOUTH) = sh%Ngs
            shGhosts(EAST) = sh%Nge
            shGhosts(WEST) = sh%Ngw
            !call GenChildShellGrid(sh, raiCpl%shGr, "raijuCpl", nGhosts=shGhosts)
            !call GenChildShellGrid(sh, raiCpl%opt%voltGrid, "raijuCpl", nGhosts=shGhosts)
            call raijuGenGridFromShGrid(raiCpl%shGr, raiCpl%opt%voltGrid, iXML)
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%pot_total)
            call initShellVar(raiCpl%shGr, SHGR_CORNER, raiCpl%pot_corot)
            call initShellVar(raiCpl%shGr, SHGR_CC, raiCpl%bvol_cc)
            associate(mask_corner=>raiCpl%pot_total%mask, mask_cc=>raiCpl%bvol_cc%mask)
                mask_corner = .true.
                mask_cc     = .true.

            do i=1,NDIM
                call initShellVar(sh, SHGR_CORNER, raiCpl%Bmin(i), maskO=mask_corner)
                call initShellVar(sh, SHGR_CORNER, raiCpl%xyzMin(i), maskO=mask_corner)
                call initShellVar(sh, SHGR_CC    , raiCpl%xyzMincc(i), maskO=mask_cc)
            enddo
            call initShellVar(sh, SHGR_CORNER, raiCpl%thcon, maskO=mask_corner)
            call initShellVar(sh, SHGR_CORNER, raiCpl%phcon, maskO=mask_corner)
            call initShellVar(sh, SHGR_CORNER, raiCpl%bVol, maskO=mask_corner)
            call initShellVar(sh, SHGR_CORNER, raiCpl%topo, maskO=mask_corner)
            call initShellVar(sh, SHGR_CORNER, raiCpl%vaFrac, maskO=mask_corner)

            allocate(raiCpl%Pavg(0:nFluidIn))
            allocate(raiCpl%Davg(0:nFluidIn))
            allocate(raiCpl%Pstd(0:nFluidIn))
            allocate(raiCpl%Dstd(0:nFluidIn))
            do i=0,nFluidIn
                call initShellVar(sh, SHGR_CC, raiCpl%Pavg(i), maskO=mask_cc)
                call initShellVar(sh, SHGR_CC, raiCpl%Davg(i), maskO=mask_cc)
                call initShellVar(sh, SHGR_CC, raiCpl%Pstd(i), maskO=mask_cc)
                call initShellVar(sh, SHGR_CC, raiCpl%Dstd(i), maskO=mask_cc)
            enddo
            call initShellVar(sh, SHGR_CC, raiCpl%Tb, maskO=mask_cc)
            end associate
        end associate
        
            ! Initial values
            raiCpl%tLastUpdate = -1.0*HUGE
            raiCpl%pot_total%data = 0.0
            raiCpl%pot_total%mask = .true.
            raiCpl%pot_corot%data = 0.0
            raiCpl%pot_corot%mask = .true.

    end subroutine raijuCpl_init


    subroutine imagTubes2RAIJU(Model, Grid, State, ijTubes)
        !! Map 2D array of IMAGTubes to RAIJU State
        type(raijuModel_T), intent(in) :: Model
        type(raijuGrid_T ), intent(in) :: Grid
        type(raijuState_T), intent(inout) :: State
        type(IMAGTube_T), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1,&
                                    Grid%shGrid%jsg:Grid%shGrid%jeg+1), intent(in) :: ijTubes
        

        integer :: i,j,s
        real(rp) :: P, D, Pstd, Dstd
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg+1) :: bVol_dip_corner
        real(rp), dimension(Grid%shGrid%isg:Grid%shGrid%ieg) :: bVol_dip_cc
        real(rp), dimension(2,2) :: dBVol
        !real(rp) :: VaMKS, Tiev, csMKS


        associate(sh=>Grid%shGrid)

            ! We'll use this later
            do i=sh%isg,sh%ieg+1
                bVol_dip_corner(i) = DipFTV_colat(Grid%thRp(i) , Model%planet%magMoment)
            enddo
            do i=sh%isg,sh%ieg
                bVol_dip_cc(i)     = DipFTV_colat(Grid%thcRp(i), Model%planet%magMoment)
            enddo

            ! Copy over all the tube info we want to have available to us

            ! Map ijTube's definition of topology to RAIJU's
            where (ijTubes%topo == 2)
                State%topo = RAIJUCLOSED
            elsewhere
                State%topo = RAIJUOPEN
            end where

            ! Assign corner quantities
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j)
            do i=sh%isg,sh%ieg+1
                do j=sh%jsg,sh%jeg+1

                    State%xyzMin(i,j,:)  = ijTubes(i,j)%X_bmin / Model%planet%rp_m  ! xyzMin in Rp
                    State%thcon(i,j)     = PI/2-ijTubes(i,j)%latc
                    State%phcon(i,j)     = ijTubes(i,j)%lonc
                    State%Bmin(i,j,ZDIR) = ijTubes(i,j)%bmin * 1.0e+9  ! Tesla -> nT
                    State%bvol(i,j)      = ijTubes(i,j)%Vol  * 1.0e-9  ! Rp/T -> Rp/nT
                    State%vaFrac(i,j)    = ijTubes(i,j)%wIMAG
                enddo
            enddo

            State%bvol_cc = 0.0
            ! Assign cell-centered quantities
            !$OMP PARALLEL DO default(shared) &
            !$OMP schedule(dynamic) &
            !$OMP private(i,j,s,P,D,Pstd,Dstd,dBVol)
            do i=sh%isg,sh%ieg
                do j=sh%jsg,sh%jeg
                    ! Note: we are doing this for all cells regardless of their goodness
                    State%xyzMincc(i,j,XDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,XDIR))  ! [Rp]
                    State%xyzMincc(i,j,YDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,YDIR))  ! [Rp]
                    State%xyzMincc(i,j,ZDIR) = toCenter2D(State%xyzMin(i:i+1,j:j+1,ZDIR))  ! [Rp]


                    ! Now we only calculate values for good cells
                    if (all(State%topo(i:i+1,j:j+1) .eq. RAIJUCLOSED)) then

                        do s=0,Model%nFluidIn
                            ! This means all 4 corners are good, can do cell centered stuff
                            P = 0.25*(ijTubes(i  ,j)%Pave(s) + ijTubes(i  ,j+1)%Pave(s) &
                                    + ijTubes(i+1,j)%Pave(s) + ijTubes(i+1,j+1)%Pave(s)) * 1.0D+9 ! [Pa -> nPa]
                            D = 0.25*(ijTubes(i  ,j)%Nave(s) + ijTubes(i  ,j+1)%Nave(s) &
                                    + ijTubes(i+1,j)%Nave(s) + ijTubes(i+1,j+1)%Nave(s)) * 1.0D-6 ! [#/m^3 --> #/cc]
                            Pstd = 0.25*(ijTubes(i  ,j)%Pstd(s) + ijTubes(i  ,j+1)%Pstd(s) &
                                       + ijTubes(i+1,j)%Pstd(s) + ijTubes(i+1,j+1)%Pstd(s)) * 1.0D+9 ! [Pa -> nPa]
                            Dstd = 0.25*(ijTubes(i  ,j)%Nstd(s) + ijTubes(i  ,j+1)%Nstd(s) &
                                       + ijTubes(i+1,j)%Nstd(s) + ijTubes(i+1,j+1)%Nstd(s)) * 1.0D+9 ! [Pa -> nPa]
                         
                            State%Pavg(i,j,s) = P
                            State%Davg(i,j,s) = D

                            State%Pstd(i,j,s) = Pstd / max(P, TINY)
                            State%Dstd(i,j,s) = Dstd / max(D, TINY)

                        enddo

                        State%Tb%data(i,j) = 0.25*(ijTubes(i  ,j)%Tb + ijTubes(i  ,j+1)%Tb &
                                                 + ijTubes(i+1,j)%Tb + ijTubes(i+1,j+1)%Tb)
                        State%bvol_cc(i,j) = toCenter2D(State%bvol(i:i+1,j:j+1))
                    endif
                enddo
            enddo

        end associate
        
    end subroutine imagTubes2RAIJU

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
        associate(Model=>raiCpl%raiApp%Model)

        ! Corners
        do i=1,NDIM
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%X_bmin(i), raiCpl%shGr, raiCpl%xyzMin(i))
        enddo
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%bmin, raiCpl%shGr, raiCpl%Bmin(ZDIR))
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%latc, raiCpl%shGr, raiCpl%thcon)
        raiCpl%thcon%data = PI/2 - raiCpl%thcon%data
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%lonc, raiCpl%shGr, raiCpl%phcon)
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%bVol, raiCpl%shGr, raiCpl%bvol)
        !call InterpShellVar_TSC_SG(voltGrid, tubeShell%bVol, raiCpl%shGr, raiCpl%bvol_cc)
        do j=raiCpl%shGr%jsg,raiCpl%shGr%jeg
            do i=raiCpl%shGr%isg,raiCpl%shGr%ieg
                raiCpl%bvol_cc%data(i,j) = toCenter2D(raiCpl%bvol%data(i:i+1,j:j+1))
            enddo
        enddo

        ! Get topo and then convert to RAIJU's definition
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%topo, raiCpl%shGr, tmpTopo)
        where (tmpTopo%data == TUBE_CLOSED)
            raiCpl%topo%data = RAIJUCLOSED
        elsewhere
            raiCpl%topo%data = RAIJUOPEN
        end where
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%wMAG, raiCpl%shGr, raiCpl%vaFrac)

        ! Centers
        
        do s=0,Model%nFluidIn
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%avgP(s), raiCpl%shGr, raiCpl%Pavg(s), srcMaskO=topoSrcMask)
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%avgN(s), raiCpl%shGr, raiCpl%Davg(s), srcMaskO=topoSrcMask)
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%stdP(s), raiCpl%shGr, raiCpl%Pstd(s), srcMaskO=topoSrcMask)
            call InterpShellVar_TSC_SG(voltGrid, tubeShell%stdN(s), raiCpl%shGr, raiCpl%Dstd(s), srcMaskO=topoSrcMask)
        enddo
        do i=1,NDIM
            call InterpShellVar_TSC_SG(raiCpl%shGr, raiCpl%xyzMin(i), raiCpl%shGr, raiCpl%xyzMincc(i), srcMaskO=topoSrcMask)
        enddo
        call InterpShellVar_TSC_SG(voltGrid, tubeShell%Tb, raiCpl%shGr, raiCpl%Tb, srcMaskO=topoSrcMask)
        

        end associate

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
        State%thcon     = raiCpl%thcon%data
        State%phcon     = raiCpl%phcon%data
        State%Bmin(:,:,ZDIR) = raiCpl%bmin(ZDIR)%data
        State%bvol      = raiCpl%bvol%data
        State%vaFrac    = raiCpl%vaFrac%data

        ! Now only copy for good points
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


!    subroutine mixPot2Raiju_RT(raiCpl, rmApp)
!        !! Take remix's potential, shove it into a remix ShellGrid, use InterpShellVar to get it onto raiju's ShellGrid
!        class(raijuCoupler_T), intent(inout) :: raiCpl
!        class(mixApp_T), intent(inout) :: rmApp
!
!        real(rp), dimension(rmApp%ion(NORTH)%shGr%Nt,rmApp%ion(NORTH)%shGr%Np) :: tmpPot
!
!        associate(rmHemi=>rmApp%ion(NORTH), Nt=>rmApp%ion(NORTH)%shGr%Nt, Np=>rmApp%ion(NORTH)%shGr%Np)       
!            rmHemi%St%pot_shGr%data(:,1:Np) = transpose(rmHemi%St%Vars(:,:,POT))
!            rmHemi%St%pot_shGr%data(:,Np+1) = rmHemi%St%pot_shGr%data(:,1)
!            rmHemi%St%pot_shGr%mask = .true.
!            call InterpShellVar_TSC_SG(rmHemi%shGr, rmHemi%St%pot_shGr, raiCpl%shGr, raiCpl%pot)
!        end associate
!
!    end subroutine mixPot2Raiju_RT

end module raijuCplHelper