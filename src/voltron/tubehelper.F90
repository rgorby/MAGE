
!Helper routines to handle tube tracing and data extraction
module tubehelper
    use volttypes
    use voltcpltypes
    use streamline
    use chmpdbz
    use shellUtils

    implicit none

    contains

    subroutine FakeTube(P,xyz0,bTube)
        !! Given a seed xyz0 generate a fake tube that seems plausible
        type(planet_T), intent(in)  :: P
        real(rp)      , intent(in)  :: xyz0(NDIM)
        type(Tube_T)  , intent(out) :: bTube

        call FreshenTube(bTube)
        call DipoleTube(P,xyz0,bTube)
        !Add fake plasma stuff
        !TODO: Add calls to gallagher/TM03/etc

    end subroutine FakeTube

    
    subroutine DipoleTube(P,xyz0,bTube)
        !! Given a seed xyz0 generate an empty tube w/ dipole quantities
        type(planet_T), intent(in)  :: P
        real(rp)      , intent(in)  :: xyz0(NDIM)
        type(Tube_T)  , intent(out) :: bTube

        real(rp) :: L,bIon

        call FreshenTube(bTube)
        !Start w/ recording seed information and such
        bTube%xyz0(:) = xyz0
        call xyz2LL(bTube%xyz0,bTube%lat0,bTube%lon0)
        bTube%invlat = InvLatitude(bTube%xyz0)

        !Topology stuff
        bTube%topo = TUBE_CLOSED !Always closed
        bTube%latc = -bTube%lat0
        bTube%lonc =  bTube%lon0
        
        !Magnetic stuff
        L = DipoleL(xyz0)
        bTube%bmin = abs(P%magMoment)*G2T*1.0e+9/L**3.0 !Min B in nT
        bTube%X_bmin = L*[cos(bTube%lon0),sin(bTube%lon0),0.0_rp]
        bTube%bVol = DipFTV_L(L,P%magMoment) !Rx/nT
        bTube%Lb = L !Just lazily using L shell
        bTube%rCurv = L/3.0
        bTube%wMAG = 1.0
        
        !Loss cone stuff
        bIon = norm2(DipoleB0(xyz0))*oBScl !Ionospheric field strength [nT]
        bTube%losscone = asin(sqrt(bTube%bmin/bIon))
        bTube%lossconec = bTube%losscone !Assuming symmetry

        !Null plasma stuff
        bTube%Tb = 0.0
        bTube%avgBeta = 0.0
        bTube%avgP(:) = 0.0
        bTube%stdP(:) = 0.0
        bTube%avgN(:) = 0.0
        bTube%stdN(:) = 0.0

        bTube%TioTe0 = 1.0 !Meh
        
    end subroutine DipoleTube

    subroutine Line2Tube(ebTrcApp,P,bTrc,bTube)
        !! Given a traced field line, populate tube object
        type(ebTrcApp_T), intent(in) :: ebTrcApp
        type(planet_T)  , intent(in)  :: P
        type(magLine_T) , intent(in)  :: bTrc
        type(Tube_T)    , intent(out) :: bTube

        integer :: OCb, s
        real(rp) :: bD,bP,dvB,bBeta,stdD,stdP,bmin
        real(rp) :: bIon,bIonC,rCurv,VebMKS,VaMKS,TiEV,CsMKS

        real(rp), dimension(NDIM) :: bEq,xyzC,xyzIonC

        call FreshenTube(bTube)

        !Start w/ recording seed information and such
        bTube%xyz0(:) = bTrc%xyz(0,:)

        call xyz2LL(bTube%xyz0,bTube%lat0,bTube%lon0)

        bTube%invlat = InvLatitude(bTube%xyz0)
        
        associate(ebModel=>ebTrcApp%ebModel,ebGr=>ebTrcApp%ebState%ebGr,ebState=>ebTrcApp%ebState)
        !Figure out topology
        !OCB =  0 (solar wind), 1 (half-closed), 2 (both ends closed)
        OCb = FLTop(ebModel,ebGr,bTrc)

        if (OCb == 2) then
            bTube%topo = TUBE_CLOSED
        elseif (OCb == 1) then
            bTube%topo = TUBE_OPEN
        else
            !What the shit?
            bTube%topo = TUBE_DISASTER
            !Well, we've done everything we can let's just call it a day
            return

        endif

    !If we're still here, let's do some stuff
        !Start w/ stuff for both closed and open
        bTube%TioTe0 = 0.0 !For now

        !Right now doing species loop for open or closed
        !TODO: Rewrite open to get near-earth average instead of field line average

        do s=0,ebModel%nSpc
            call FLThermo(ebModel,ebGr,bTrc,bD,bP,dvB,bBeta,s)
            call FLStdev (ebModel,ebGr,bTrc,bD,bP,stdD,stdP,s)
            bTube%avgP(s) = bP
            bTube%stdP(s) = stdP
            bTube%avgN(s) = bD
            bTube%stdN(s) = stdD
            if (s == BLK) then
                bTube%avgBeta = bBeta
            endif
        enddo

        if (bTube%topo == TUBE_CLOSED) then
        !Do closed field line stuff
            bTube%bVol = dvB/oBScl !Rx/EB => Rx/nT

        !Minimal surface (bEq in Rp, bMin in EB)
            call FLEq(ebModel,bTrc,bEq,bmin)
            bTube%bmin = bmin*oBScl !EB => nT
            bTube%X_bmin = bEq

        !Line stuff
            bTube%Lb = FLArc(ebModel,ebGr,bTrc) !Rx
            !TODO: Check units of AlfvenX calculation
            bTube%Tb = FLAlfvenX(ebModel,ebGr,bTrc)
        !Conjugate stuff
            call FLConj(ebModel,ebGr,bTrc,xyzC)
            !Shift conjugate point to match radius of seed point
            xyzIonC = DipoleShift(xyzC, norm2(bTube%xyz0))

            call xyz2LL(xyzIonC,bTube%latc,bTube%lonc)
        !Loss cone stuff
            bIon  = norm2(DipoleB0(bTube%xyz0))*oBScl !Ionospheric field strength [nT]
            bIonC = norm2(DipoleB0(xyzIonC))*oBScl !Ionospheric field strength [nT]
            bTube%losscone   = asin(sqrt(bmin/bIon))
            bTube%lossconec  = asin(sqrt(bmin/bIonC)) !Keeping these separate for IGRF later

        !Get energy partition and curvature stuff
            !Get curvature radius and ExB velocity [km/s]
            call FLCurvRadius(ebModel,ebGr,ebState,bTrc,rCurv,VebMKS)
            bTube%rCurv = rCurv
            bD = bTube%avgN(BLK)
            bP = bTube%avgP(BLK)

            !VaMKS = flux tube arc length [km] / Alfven crossing time [s]
            VaMKS = (bTube%Lb*P%rp_m*1.0e-3)/bTube%Tb 
            !CsMKS = 9.79 x sqrt(5/3 * Ti) km/s, Ti eV
            !TiEV = (1.0e+3)*DP2kT(bDRC*1.0e-6,bPRC*1.0e+9) !Temp in eV
            TiEV = (1.0e+3)*DP2kT(bD,bP) !Temp in eV
            CsMKS = 9.79*sqrt((5.0/3)*TiEV)
            bTube%wMAG = VaMKS/( sqrt(VaMKS**2.0 + CsMKS**2.0) + VebMKS)
        endif

        end associate !ebTrc stuff

    end subroutine Line2Tube

    subroutine xyz2LL(xyz,lat,lon)
        real(rp), intent(in)  :: xyz(NDIM)
        real(rp), intent(out) :: lat,lon

        lon = katan2(xyz(YDIR),xyz(XDIR) )
        lat = asin  (xyz(ZDIR)/norm2(xyz))

    end subroutine xyz2LL

    !Just zero out values
    subroutine FreshenTube(bTube)
        type(Tube_T), intent(out) :: bTube
        bTube%xyz0    = 0.0
        bTube%lat0    = 0.0
        bTube%lon0    = 0.0
        bTube%invlat  = 0.0

        bTube%topo = TUBE_UNDEF
        bTube%bmin = 0.0
        bTube%X_bmin = 0.0
        bTube%bVol = 0.0
        bTube%Lb = 0.0
        bTube%Tb = 0.0

        bTube%wMAG = 0.0
        bTube%rCurv = 0.0
        bTube%avgBeta = 0.0

        bTube%avgP(:) = 0.0
        bTube%stdP(:) = 0.0
        bTube%avgN(:) = 0.0
        bTube%stdN(:) = 0.0

        bTube%losscone  = 0.0
        bTube%lossconec = 0.0

        bTube%TioTe0 = 0.0
    end subroutine FreshenTube


    subroutine tubes2Shell(shGr, ijTubes, tubeShell)
        type(ShellGrid_T), intent(in) :: shGr
            !! ShellGrid that ijTubes and tubeShell live on
        type(Tube_T), dimension(shGr%Nt,shGr%Np), intent(in) :: ijTubes
            !! Nt+1,Np+1 array of Tube_T's
        type(TubeShell_T), intent(inout) :: tubeShell
            !! ShellGridVar versions of 2D Tube_T data we are populating

        integer :: i,j,k

        ! Copy tubes to active tubeShell domain
        !$OMP PARALLEL DO default(shared) &
        !$OMP schedule(dynamic) &
        !$OMP private(i,j,k)
        do j=shGr%js,shGr%je+1
            do i=shGr%is,shGr%ie+1
                do k=1,NDIM
                    tubeShell%xyz0  (k)%data(i,j) = ijTubes(i,j)%xyz0(k)
                    tubeShell%X_bmin(k)%data(i,j) = ijTubes(i,j)%X_bmin(k)
                enddo
                tubeShell%lat0   %data(i,j) = ijTubes(i,j)%lat0
                tubeShell%lon0   %data(i,j) = ijTubes(i,j)%lon0
                tubeShell%invlat %data(i,j) = ijTubes(i,j)%invlat
                tubeShell%latc   %data(i,j) = ijTubes(i,j)%latc
                tubeShell%lonc   %data(i,j) = ijTubes(i,j)%lonc
                tubeShell%topo   %data(i,j) = ijTubes(i,j)%topo
                tubeShell%bmin   %data(i,j) = ijTubes(i,j)%bmin
                tubeShell%bVol   %data(i,j) = ijTubes(i,j)%bVol
                tubeShell%Lb     %data(i,j) = ijTubes(i,j)%Lb
                tubeShell%Tb     %data(i,j) = ijTubes(i,j)%Tb
                tubeShell%wMAG   %data(i,j) = ijTubes(i,j)%wMAG
                tubeShell%rCurv  %data(i,j) = ijTubes(i,j)%rCurv
                tubeShell%avgBeta%data(i,j) = ijTubes(i,j)%avgBeta
                do k=0,MAXTUBEFLUIDS
                    tubeShell%avgP(k)%data(i,j) = ijTubes(i,j)%avgP(k)
                    tubeShell%avgN(k)%data(i,j) = ijTubes(i,j)%avgN(k)
                    tubeShell%stdP(k)%data(i,j) = ijTubes(i,j)%stdP(k)
                    tubeShell%stdN(k)%data(i,j) = ijTubes(i,j)%stdN(k)
                enddo
                tubeShell%losscone %data(i,j) = ijTubes(i,j)%losscone
                tubeShell%lossconec%data(i,j) = ijTubes(i,j)%lossconec
                tubeShell%TioTe0   %data(i,j) = ijTubes(i,j)%TioTe0
            enddo
        enddo

        ! Now wrap everyone so ghosts are filled in
        do k=1,NDIM
            call wrapJ_SGV(shGr, tubeShell%xyz0(k))
            call wrapJ_SGV(shGr, tubeShell%X_bmin(k))
        enddo
        call wrapJ_SGV(shGr, tubeShell%lat0   )
        call wrapJ_SGV(shGr, tubeShell%lon0   )
        call wrapJ_SGV(shGr, tubeShell%invlat )
        call wrapJ_SGV(shGr, tubeShell%latc   )
        call wrapJ_SGV(shGr, tubeShell%lonc   )
        call wrapJ_SGV(shGr, tubeShell%topo   )
        call wrapJ_SGV(shGr, tubeShell%bmin   )
        call wrapJ_SGV(shGr, tubeShell%bVol   )
        call wrapJ_SGV(shGr, tubeShell%Lb     )
        call wrapJ_SGV(shGr, tubeShell%Tb     )
        call wrapJ_SGV(shGr, tubeShell%wMAG   )
        call wrapJ_SGV(shGr, tubeShell%rCurv  )
        call wrapJ_SGV(shGr, tubeShell%avgBeta)
        do k=0,MAXTUBEFLUIDS
            call wrapJ_SGV(shGr, tubeShell%avgP(k))
            call wrapJ_SGV(shGr, tubeShell%avgN(k))
            call wrapJ_SGV(shGr, tubeShell%stdP(k))
            call wrapJ_SGV(shGr, tubeShell%stdN(k))
        enddo
        call wrapJ_SGV(shGr, tubeShell%losscone )
        call wrapJ_SGV(shGr, tubeShell%lossconec)
        call wrapJ_SGV(shGr, tubeShell%TioTe0   )
    end subroutine tubes2Shell
end module tubehelper
