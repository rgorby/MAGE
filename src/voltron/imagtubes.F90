
!Various routines to handle flux-tube calculations to feed inner magnetosphere models
module imagtubes
    use volttypes
    use voltcpltypes
    use streamline
    use chmpdbz, ONLY : DipoleB0
    use imaghelper
    use planethelper
    use shellGrid

    implicit none

    contains

    !subroutine init_TubeShell(sh, tubeShell, maskO, TioTeO)
    subroutine init_TubeShell(sh, tubeShell, TioTeO)
        type(ShellGrid_T), intent(in) :: sh
            !! Voltron shellgrid
        type(TubeShell_T), intent(inout) :: tubeShell
            !! IMAGTubeShell object we are initializing
        !logical, dimension(:,:), intent(in), optional :: maskO
        real(rp), intent(in), optional :: TioteO
            !! Default Ti/Te ratio

        integer :: i
        real(rp) :: tiote

        if (present(TioteO)) then
            tiote = TioteO
        else 
            tiote = 4.0
        endif


        do i=1,NDIM
            call initShellVar(sh, SHGR_CORNER, tubeShell%xyz0(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%X_bmin(i))
        enddo
        call initShellVar(sh, SHGR_CORNER, tubeShell%lat0   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%lon0   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%latc   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%lonc   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%invlat )
        call initShellVar(sh, SHGR_CORNER, tubeShell%topo   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%bmin   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%bVol   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%Lb     )
        call initShellVar(sh, SHGR_CORNER, tubeShell%Tb     )
        call initShellVar(sh, SHGR_CORNER, tubeShell%wMAG   )
        call initShellVar(sh, SHGR_CORNER, tubeShell%rCurv  )
        call initShellVar(sh, SHGR_CORNER, tubeShell%avgBeta)
        do i=0,MAXTUBEFLUIDS
            call initShellVar(sh, SHGR_CORNER, tubeShell%avgP(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%avgN(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%stdP(i))
            call initShellVar(sh, SHGR_CORNER, tubeShell%stdN(i))
        enddo
        call initShellVar(sh, SHGR_CORNER, tubeShell%losscone )
        call initShellVar(sh, SHGR_CORNER, tubeShell%lossconec)
        call initShellVar(sh, SHGR_CORNER, tubeShell%TioTe0   )

        tubeShell%TioTe0%data = tiote


    end subroutine init_TubeShell

    ! Dipole flux tube info
    subroutine DipoleTube(vApp,lat,lon,ijTube)
        type(voltApp_T), intent(in) :: vApp
        real(rp), intent(in) :: lat,lon
        type(IMAGTube_T), intent(out)   :: ijTube
        
        real(rp) :: L,colat
        real(rp) :: mdipole

        mdipole = ABS(vapp%planet%magMoment)*G2T ! dipole moment in T
        colat = PI/2 - lat
        L = 1.0/(sin(colat)**2.0)
        !ijTube%Vol = 32./35.*L**4.0/mdipole
        !Use full dipole FTV formula, convert Rx/nT => Rx/T
        ijTube%Vol = DipFTV_L(L,vapp%planet%magMoment)*1.0e+9
        
        ijTube%X_bmin(XDIR) = L*cos(lon)*vApp%planet%rp_m !Rp=>meters
        ijTube%X_bmin(YDIR) = L*sin(lon)*vApp%planet%rp_m !Rp=>meters
        ijTube%X_bmin(ZDIR) = 0.0
        ijTube%bmin = mdipole/L**3.0

        ijTube%topo = 2
        
        ijTube%beta_average = 0.0
        ijTube%Pave = 0.0
        ijTube%Nave = 0.0

        ijTube%latc = -lat
        ijTube%lonc = lon
        ijTube%Lb   = L !Just lazily using L shell
        ijTube%Tb   = 0.0
        ijTube%losscone = 0.0
        ijTube%rCurv = L/3.0

        ijTube%wIMAG = 1.0 !Much imag
        ijTube%TioTe0 = tiote_RCM
    end subroutine DipoleTube


    subroutine MHDTube(ebTrcApp,planet,colat,lon,r,ijTube,bTrc,nTrcO,doShiftO)
        type(ebTrcApp_T), intent(in) :: ebTrcApp
        type(planet_T), intent(in) :: planet
        real(rp), intent(in) :: colat, lon, r
        type(IMAGTube_T), intent(inout) :: ijTube
        type(magLine_T), intent(inout) :: bTrc
        integer, intent(in), optional :: nTrcO
        logical, intent(in), optional :: doShiftO

        integer :: s
        real(rp) :: t, bMin,bIon
        real(rp), dimension(NDIM) :: x0, bEq, xyzIon
        real(rp), dimension(NDIM) :: xyzC,xyzIonC
        integer :: OCb
        real(rp) :: bD,bP,dvB,bBeta,rCurv, bDRC, bPRC, bBetaRC, stdP, stdD
        real(rp) :: VaMKS,CsMKS,VebMKS !Speeds in km/s
        real(rp) :: TiEV !Temperature in ev

        logical :: doShift,doShue

        if (present(doShiftO)) then
            doShift = doShiftO
        else
            doShift = .false.
        endif

        ! Take seed point in spherical coordinates
        xyzIon(XDIR) = r*sin(colat)*cos(lon)
        xyzIon(YDIR) = r*sin(colat)*sin(lon)
        xyzIon(ZDIR) = r*cos(colat)
        if (doShift) then
            !! TODO: Decide if this is how we want to do it
            x0 = DipoleShift(xyzIon, planet%ri_m/planet%rp_m + TINY)
        else
            x0 = xyzIon
        endif
        
        bIon = norm2(DipoleB0(xyzIon))*oBScl*1.0e-9 !EB=>T, ionospheric field strength


        associate(ebModel=>ebTrcApp%ebModel,ebGr=>ebTrcApp%ebState%ebGr,ebState=>ebTrcApp%ebState)

    !Now do field line trace
        t = ebState%eb1%time !Time in CHIMP units
            !!TODO: What do we do when we want times in between steps? Somethign similar to what slice/chop do to output
        
        if (present(nTrcO)) then
            call genLine(ebModel,ebState,x0,t,bTrc,nTrcO,doShueO=.false.,doNHO=.true.)
        else
            call genLine(ebModel,ebState,x0,t,bTrc,      doShueO=.false.,doNHO=.true.)
        endif


        !Topology
        !OCB =  0 (solar wind), 1 (half-closed), 2 (both ends closed)
        OCb = FLTop(ebModel,ebGr,bTrc)
        ijTube%topo = OCb
        
        if (OCb /= 2) then
            ! If the field line hasn't closed after 15 minutes we are legally allowed to leave
            ijTube%X_bmin = 0.0
            ijTube%bmin = 0.0
            ijTube%Vol = -1.0
            ijTube%Pave = 0.0
            ijTube%Nave = 0.0
            ijTube%Pstd = 0.0
            ijTube%Nstd = 0.0
            ijTube%beta_average = 0.0
            ijTube%latc = 0.0
            ijTube%lonc = 0.0
            ijTube%Lb   = 0.0
            ijTube%Tb   = 0.0
            ijTube%losscone = 0.0
            ijTube%rCurv = 0.0
            ijTube%wIMAG = 0.0
            ijTube%Veb   = 0.0
            return
        endif

    !Get diagnostics from closed field line
        !Minimal surface (bEq in Rp, bMin in EB)
        call FLEq(ebModel,bTrc,bEq,bMin)
        bMin = bMin*oBScl*1.0e-9 !EB=>Tesla
        bEq = bEq*planet%rp_m ! Rp -> meters

        ! Plasma quantities
        do s=0,ebModel%nSpc
            !dvB = Flux-tube volume (Re/EB)
            !write(*,*)"FLThermo, s=",s
            call FLThermo(ebModel,ebGr,bTrc,bD,bP,dvB,bBeta,s)
            call FLStdev (ebModel,ebGr,bTrc,bD,bP,stdD,stdP,s)
            !call FLStdev (ebModel,ebGr,bTrc,stdD,stdP,s)
            bP = bP*1.0e-9 !nPa=>Pa
            bD = bD*1.0e+6 !#/cc => #/m3
            ijTube%Pave(s) = bP
            ijTube%Nave(s) = bD
            ijTube%Pstd(s) = stdP * 1.0e-9  ! nPa=>Pa
            ijTube%Nstd(s) = stdD * 1.0e+6  ! #/cc=>#/m3
            if (s .eq. RCFLUID) then
                bPRC = bP
                bDRC = bD
                bBetaRC = bBeta
            endif
        enddo

        !Converts Rp/EB => Rp/T
        dvB = dvB/(oBScl*1.0e-9)
        
        ijTube%X_bmin = bEq
        ijTube%bmin = bMin
        ijTube%Vol = dvB  ! Taken from last FLThermo call
        ijTube%beta_average = bBetaRC

        call FLConj(ebModel,ebGr,bTrc,xyzC)
        if (doShift) then
            xyzIonC = DipoleShift(xyzC, planet%ri_m)
        else
            xyzIonC = xyzC
        endif
        ijTube%latc = asin(xyzIonC(ZDIR)/norm2(xyzIonC))
        ijTube%lonc = modulo( atan2(xyzIonC(YDIR),xyzIonC(XDIR)),2*PI )
        ijTube%Lb = FLArc(ebModel,ebGr,bTrc)
        !NOTE: Bounce timescale may be altered to use IMAG hot density
        ijTube%Tb = FLAlfvenX(ebModel,ebGr,bTrc)
        !ijTube%Tb = FLFastX(ebModel,ebGr,bTrc)

        !write(*,*) 'Bounce compare: = ', FLFastX(ebModel,ebGr,bTrc)/FLAlfvenX(ebModel,ebGr,bTrc)
        
        ijTube%losscone = asin(sqrt(bMin/bIon))

        !Get curvature radius and ExB velocity [km/s]
        call FLCurvRadius(ebModel,ebGr,ebState,bTrc,rCurv,VebMKS)
        ijTube%rCurv = rCurv
        ijTube%Veb = VebMKS

    !Get confidence interval
        !VaMKS = flux tube arc length [km] / Alfven crossing time [s]
        VaMKS = (ijTube%Lb*planet%rp_m*1.0e-3)/ijTube%Tb 
        !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
        !TiEV = (1.0e+3)*DP2kT(bDRC*1.0e-6,bPRC*1.0e+9) !Temp in eV
        TiEV = (1.0e+3)*DP2kT(ijTube%Nave(BLK)*1.0D-6,ijTube%Pave(BLK)*1.0D9) !Temp in eV
        CsMKS = 9.79*sqrt((5.0/3)*TiEV)

        ijTube%wIMAG = VaMKS/( sqrt(VaMKS**2.0 + CsMKS**2.0) + VebMKS)

        end associate

    end subroutine MHDTube

end module imagtubes
