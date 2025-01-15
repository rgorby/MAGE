
!Helper routines to handle tube tracing and data extraction
module tubehelper
    use volttypes
    use voltcpltypes
    use streamline
    use chmpdbz

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
end module tubehelper
