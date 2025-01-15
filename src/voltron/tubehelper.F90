
!Helper routines to handle tube tracing and data extraction
module tubehelper
    use volttypes
    use voltcpltypes
    use streamline
    use chmpdbz
    use shellUtils

    implicit none

    contains

    ! subroutine Line2Tube(bTrc,bTube)
    ! 	!! Given a traced field line, populate tube object

    ! 	type(magLine_T), intent(in)  :: bTrc
    !     type(Tube_T)   , intent(out) :: bTube

    !     !Start w/ recording seed information and such
    !     bTube%xyz0(:) = bTrc%xyz(0,:)

    !     bTube%lon0 = katan2(bTube%xyz0(YDIR),bTube%xyz0(XDIR))
    !     bTube%lat0 = asin( bTube%xyz0(ZDIR),norm2(bTube%xyz0) )

    !     bTube%invlat = InvLatitude(bTube%xyz0)
        
    ! end subroutine Line2Tube


    subroutine FakeTube(P,xyz0,bTube)
        !! Given a seed xyz0 generate a fake tube that seems plausible
        type(planet_T), intent(in)  :: P
        real(rp)      , intent(in)  :: xyz0(NDIM)
        type(Tube_T)  , intent(out) :: bTube

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

        !Start w/ recording seed information and such
        bTube%xyz0(:) = xyz0

        bTube%lon0 = katan2(bTube%xyz0(YDIR),bTube%xyz0(XDIR))
        bTube%lat0 = asin( bTube%xyz0(ZDIR)/norm2(bTube%xyz0) )
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
