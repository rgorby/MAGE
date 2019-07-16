module fields
    use types
    use clocks
    use gamutils
    use recon
    use ringutils
    use gridutils

    implicit none

    logical, parameter :: doVa  = .true. !Use Alfven speed in diffusive velocity
    logical, parameter :: doRingRenorm = .false. !Do ring renormalization on faces/fluxes, pretty slow
    logical, parameter :: doVdA = .true. !Do area scaling for velocity->corner
    logical, parameter :: doBdA = .true. !Do area scaling for face flux->edge

    logical :: initField = .true. !Do we need to initialize module workspaces
    real(rp), dimension(:,:,:,:), allocatable, private :: Vf
    !Vf(i,j,k,XYZ-DIR), XYZ velocities pushed to dT1 faces
    !Vf = (Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM)

    contains

    !Use electric field to update face fluxes
    subroutine E2Flux(Model,Gr,magFlux,E,dtOpt)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(inout) :: magFlux(Gr%isg:Gr%ieg+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,NDIM)
        real(rp), intent(inout) :: E(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM)
        real(rp), intent(in), optional :: dtOpt

        integer :: i,j,k
        real(rp) :: dt,dfi,dfj,dfk

        !Use dt=1 if unspecified (ie for ring-avg)
        dt = 1.0
        if (present(dtOpt)) then
            dt = dtOpt
        endif

        !$OMP PARALLEL DO default (shared) collapse(2) &
        !$OMP private(dfi,dfj,dfk) 
        do k=Gr%ks,Gr%ke+1
            do j=Gr%js,Gr%je+1
                do i=Gr%is,Gr%ie+1
                    dfi = E(i,j+1,k,KDIR) - E(i,j,k,KDIR) - E(i,j,k+1,JDIR) + E(i,j,k,JDIR)
                    dfj = E(i,j,k+1,IDIR) - E(i,j,k,IDIR) - E(i+1,j,k,KDIR) + E(i,j,k,KDIR)
                    dfk = E(i+1,j,k,JDIR) - E(i,j,k,JDIR) - E(i,j+1,k,IDIR) + E(i,j,k,IDIR)
                    magFlux(i,j,k,IDIR) = magFlux(i,j,k,IDIR) - dt*dfi
                    magFlux(i,j,k,JDIR) = magFlux(i,j,k,JDIR) - dt*dfj
                    magFlux(i,j,k,KDIR) = magFlux(i,j,k,KDIR) - dt*dfk

                enddo
            enddo
        enddo   
        
    end subroutine E2Flux

    !Calculates electric field using variables in State (generally predictor state)
    !General structure of computation
    !Calculate all cell-centered velocities (XYZ)
    !Loop over E-field directions (eD)
    !--Set triad, (eD,dT1,dT2)
    !--Push cell-center Vxyz velocities to faces along dT1
    !-----Vary loop order for cache locality
    !-----NEED IMPLICIT OMP BARRIERS at end of loops 
    !--Loop over active i,j,k edges (block inner i dimension)
    !----Get corner Vxyz's from face-centers by pushing from face->edge along dT2
    !----Turn corner Vxyz's into V1/V2 using edge transforms
    !----Now have full v1/v2
    !----Get bT1,bT2.  Map to 1/2, calculate dB1,dB2
    !----Add B0 corrections
    !----Calculate diffusive velocity and finish EMF calculation, scale w/ edge length
    !Finally, do any other E field relevant calculations, ie resistivity

    subroutine CalcElecField(Model,Gr,State,E)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real(rp), dimension(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM), intent(inout) :: E

        !Vector buffers
        real(rp), dimension(vecLen) :: v1,v2,b1,b2,Jd,Dc,vDiff
        real(rp) :: VelB(vecLen,NDIM)
        real(rp) :: vA
        integer :: i,iB,ieB,j,k,iG,iMax
        integer :: ie,je,ke,ksg,keg
        integer :: eD,eD0,dT1,dT2

        !DIR$ ASSUME_ALIGNED E: ALIGN
        !DIR$ ATTRIBUTES align : ALIGN :: v1,v2,b1,b2,Jd,Dc,vDiff,VelB

        if (initField) then
            !Initialize workspaces
            !Vf(i,j,k,XYZ-DIR), XYZ velocities pushed to dT1 faces
            allocate(Vf(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM))
            Vf = 0.0
            initField = .false.
        endif

        !Prep bounds for this timestep
        eD0 = 1 !Starting direction for EMF

        ksg = Gr%ksg
        keg = Gr%keg

        !Open big parallel region
        !$OMP PARALLEL default(shared) &
        !$OMP private(v1,v2,b1,b2,Jd,Dc,vDiff,VelB) &
        !$OMP private(i,iB,ieB,j,k,iG,iMax,ie,je,ke,eD,dT1,dT2,vA)
        
        !Initialize thread-private blocks
        v1 = 0.0
        v2 = 0.0
        b1 = 0.0
        b2 = 0.0
        Jd = 0.0
        Dc = 0.0
        VelB = 0.0
        

        do eD=eD0,NDIM
            !$OMP SINGLE
            call Tic("Mom2Face")
            !$OMP END SINGLE NOWAIT

            !Set last active edges in IJK
            !TODO: Modify per direction
            ie = Gr%ie+1
            je = Gr%je+1
            ke = Gr%ke+1
            
            !Use edge normal direction to calculate local triad
            !Push cell-centered velocities along dT1 to face
            !TODO: Vary loop order for locality
            !TODO: Combine I/J center->face
            select case(eD)
                !Ei fields
                case(IDIR)
                    dT1 = JDIR; dT2 = KDIR
                    !$OMP DO collapse(2)
                    do k=ksg,keg
                        do iB=Gr%isg,Gr%ieg,vecLen !Block loop
                            do j=Gr%js,Gr%je+1 !J reconstruction
                            
                                iMax = min(vecLen,Gr%ieg-iB+1)
                                ieB = iB+iMax-1

                                call Mom2Face(Model,Gr,State%Gas(:,:,:,:,BLK),VelB,iB,j,k,iMax,dT1)

                                Vf(iB:ieB,j,k,:) = VelB(1:iMax,:)
                            enddo
                        enddo
                    enddo

                case(JDIR)
                    !Ej fields
                    dT1 = KDIR; dT2 = IDIR
                    !$OMP DO collapse(2)
                    do j=Gr%jsg,Gr%jeg
                        do iB=Gr%isg,Gr%ieg,vecLen !Block loop
                            do k=Gr%ks,Gr%ke+1 !K reconstruction
                            
                                iMax = min(vecLen,Gr%ieg-iB+1)
                                ieB = iB+iMax-1

                                call Mom2Face(Model,Gr,State%Gas(:,:,:,:,BLK),VelB,iB,j,k,iMax,dT1)

                                Vf(iB:ieB,j,k,:) = VelB(1:iMax,:)
                            enddo
                        enddo
                    enddo

                case(KDIR)
                    !Ek fields
                    dT1 = IDIR; dT2 = JDIR
                    !$OMP DO collapse(2)
                    do k=ksg,keg
                        do j=Gr%jsg,Gr%jeg
                            do iB=Gr%is,Gr%ie+1,vecLen !Block loop/I reconstruction
                            
                                iMax = min(vecLen,Gr%ie+1-iB+1)
                                ieB = iB+iMax-1

                                call Mom2Face(Model,Gr,State%Gas(:,:,:,:,BLK),VelB,iB,j,k,iMax,dT1)

                                Vf(iB:ieB,j,k,:) = VelB(1:iMax,:)
                            enddo
                        enddo
                    enddo
                !NOTE: Each OMP DO has an implicit barrier at the end
            end select

            !$OMP SINGLE
            call Toc("Mom2Face")
            call Tic("VxB")
            !$OMP END SINGLE NOWAIT

            !Loop over active edges
            !$OMP DO collapse(2)
            do k=Gr%ks, ke
                do j=Gr%js, je
                    do iB=Gr%is,ie,vecLen
                        !Get size of this vector brick
                        iMax = min(vecLen,ie-iB+1)

                        vDiff = 0.0 !Zero out by default
                        !Push from face to edge (v1/v2) and get diffusive flow speed (vDiff)
                        call GetCornerV(Model,Gr,Vf,v1,v2,vDiff,iB,j,k,iMax,eD,dT1,dT2)

                        call GetCornerD(Model,Gr,State%Gas(:,:,:,DEN,BLK),Dc,iB,j,k,iMax,eD)
                        call GetCornerB(Model,Gr,State%magFlux,b1,b2,Jd,iB,j,k,iMax,eD,dT1,dT2)

                        !Now we have everything, calculate diffusive speed and do field
                        do i=1,iMax
                            iG =iB+i-1
                            
                            !vDiff = 0.0
                            if (doVa) then
                                !Add Alfven speed to diffusive speed
                                vA = sqrt(b1(i)**2.0 + b2(i)**2.0)/sqrt(Dc(i))
                                !Boris correct
                                if (Model%doBoris) then
                                    vA = Model%Ca*vA/sqrt(Model%Ca**2.0 + vA**2.0)
                                endif
                                vDiff(i) = vDiff(i) + vA
                                vDiff(i) = min(vDiff(i),Model%CFL*Gr%edge(iG,j,k,eD)/Model%dt)
                            endif

                            !Final field
                            E(iG,j,k,eD) = -( v1(i)*b2(i) - v2(i)*b1(i) ) + Model%Vd0*vDiff(i)*Jd(i)
                            !Scale by edge length
                            E(iG,j,k,eD) = Gr%edge(iG,j,k,eD)*E(iG,j,k,eD)
                        enddo
                    enddo !iB loop
                enddo

            enddo !k loop

            !$OMP SINGLE
            call Toc("VxB")
            !$OMP END SINGLE !IMPLICIT BARRIER 
            
        enddo !eD loop, EMF direction

        !$OMP END PARALLEL

        if(Model%useResistivity) call resistivity(Model,Gr,State,E)
        
    end subroutine CalcElecField


    subroutine Mom2Face(Model,Gr,W,VfB,iB,j,k,iMax,dT)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(in) :: W(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NVAR)
        integer, intent(in) :: iB,j,k,iMax,dT

        real(rp), intent(out) :: VfB(vecLen,NDIM)

        integer :: i,n,d
        real(rp), dimension(vecLen,recLen,NDIM) :: MomB,VelB
        real(rp), dimension(vecLen,recLen) :: DenB,VolB
        real(rp), dimension(vecLen) :: dV

        !DIR$ ASSUME_ALIGNED W: ALIGN
        !DIR$ ASSUME_ALIGNED VfB: ALIGN
        !DIR$ attributes align : ALIGN :: MomB,VelB,DenB,VolB,dV

        MomB = 0.0
        DenB = 0.0
        VfB  = 0.0
        VolB = 0.0

        !Get stencils
        call LoadBlock(Model,Gr,VolB          ,Gr%volume    ,iB,j,k,iMax,dT)
        call LoadBlock(Model,Gr,DenB          ,W(:,:,:,DEN ),iB,j,k,iMax,dT)
        call LoadBlock(Model,Gr,MomB(:,:,XDIR),W(:,:,:,MOMX),iB,j,k,iMax,dT)
        call LoadBlock(Model,Gr,MomB(:,:,YDIR),W(:,:,:,MOMY),iB,j,k,iMax,dT)
        call LoadBlock(Model,Gr,MomB(:,:,ZDIR),W(:,:,:,MOMZ),iB,j,k,iMax,dT)

        !Get stencil for VdV
        do n=1,recLen
            do i=1,iMax
                VelB(i,n,XDIR) = VolB(i,n)*MomB(i,n,XDIR)/max( DenB(i,n),dFloor )
                VelB(i,n,YDIR) = VolB(i,n)*MomB(i,n,YDIR)/max( DenB(i,n),dFloor )
                VelB(i,n,ZDIR) = VolB(i,n)*MomB(i,n,ZDIR)/max( DenB(i,n),dFloor )
            enddo
        enddo

        !Interpolate dV to corner
        do i=1,iMax
            dV(i) = dot_product(interpWgt,VolB(i,:))
        enddo

        !Now interpolate, <VdV>/<dV>
        do d=1,NDIM
            do i=1,iMax
                !VfB(i,d) = dot_product(interpWgt,MomB(i,:,d)/DenB(i,:))
                Vfb(i,d) = dot_product(interpWgt,VelB(i,:,d))/dV(i)
            enddo
        enddo

    end subroutine Mom2Face

    !Get Edge LRs for field/face areas in triad dN,dT1,dT2
    subroutine EdgeLRs(Model,Gr,bFlux,bL,bR,fA,iB,j,k,iMax,dN,dT1,dT2)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(in) :: bFlux(Gr%isg:Gr%ieg+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,NDIM)
        integer, intent(in) :: iB,j,k,iMax,dN,dT1,dT2
        real(rp), dimension(vecLen), intent(out) :: bL,bR,fA

        real(rp), dimension(vecLen,recLen) :: AreaB,VolB
        real(rp), dimension(vecLen,recLen,1) :: MagB
        real(rp), dimension(vecLen,1) :: MagLB,MagRB
        integer :: i,n
        !Note, bFlux is *NOT* aligned
        !DIR$ ASSUME_ALIGNED bL: ALIGN
        !DIR$ ASSUME_ALIGNED bR: ALIGN
        !DIR$ ASSUME_ALIGNED fA: ALIGN
        !DIR$ attributes align : ALIGN :: AreaB,VolB,MagB,MagLB,MagRB

        fA   = 0.0
        bL   = 0.0
        bR   = 0.0
        VolB = 0.0

        !Get stencils for this sweep
        !dT1 faces in dT2 direction
        call LoadBlockI(Model,Gr,AreaB      ,Gr%Face(:,:,:,dT1),iB,j,k,iMax,dT2)
        call LoadBlockI(Model,Gr,MagB(:,:,1),  bFlux(:,:,:,dT1),iB,j,k,iMax,dT2)

        if (Model%doRing .and. doRingRenorm) then
            !Note flipped dT2/dT1 for different component/stencil ordering
            call RingRenorm(Model,Gr,AreaB      ,iB,j,k,dT2,dT1)
            call RingRenorm(Model,Gr,MagB(:,:,1),iB,j,k,dT2,dT1)
        endif
        
        !Split into L/Rs
        if (doBdA) then
            !Do area weighting
            do i=1,iMax
                VolB(i,:) = AreaB(i,:)
                MagB(i,:,1) = MagB(i,:,1)/max(AreaB(i,:),TINY)
            enddo
        else
            !Fake volume-weighting for BlockLR routine
            VolB = 1.0 
        endif !doBdA

        call BlockLRs(VolB,MagB,MagLB,MagRB,1)
        bL = MagLB(:,1)
        bR = MagRB(:,1)

        !Interpolate (no need to split) face areas
        if (doBdA) then
            do i=1,iMax
                fA(i) = dot_product(interpWgt,AreaB(i,:))
                bL(i) = bL(i)*fA(i)
                bR(i) = bR(i)*fA(i)
            enddo
        else
            do i=1,iMax
                fA(i) = dot_product(interpWgt,AreaB(i,:))
            enddo
        endif !doBdA            

    end subroutine EdgeLRs

    !Get corner fields and diffusive current
    subroutine GetCornerB(Model,Gr,bFlux,b1,b2,Jd,iB,j,k,iMax,dN,dT1,dT2)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(in) :: bFlux(Gr%isg:Gr%ieg+1,Gr%jsg:Gr%jeg+1,Gr%ksg:Gr%keg+1,NDIM)
        real(rp), dimension(vecLen), intent(out) :: b1,b2,Jd
        integer, intent(in) :: iB,j,k,iMax,dN,dT1,dT2

        real(rp), dimension(vecLen) :: bT1L,bT1R,fA1,bT1
        real(rp), dimension(vecLen) :: bT2L,bT2R,fA2,bT2
        real(rp) :: db1,db2,xni,yni,xnj,ynj,detT
        real(rp) :: b1L,b2L,b1R,b2R
        integer :: i,iG

        !DIR$ ASSUME_ALIGNED b1: ALIGN
        !DIR$ ASSUME_ALIGNED b2: ALIGN
        !DIR$ ASSUME_ALIGNED Jd: ALIGN
        !DIR$ attributes align : ALIGN :: bT1L,bT1R,fA1,bT1,bT2L,bT2R,fA2,bT2

        !Start by getting bT1LRs and bT2LRs
        call EdgeLRs(Model,Gr,bFlux,bT1L,bT1R,fA1,iB,j,k,iMax,dN,dT1,dT2)
        call EdgeLRs(Model,Gr,bFlux,bT2L,bT2R,fA2,iB,j,k,iMax,dN,dT2,dT1)

        
        do i=1,iMax
            iG =iB+i-1

            !Turn bTx LR's into bTx's and diffusive currents

            !Turn bTx's into b1,b2 w/ transform
            xni = Gr%Teb(iG,j,k,XNQI,dN)
            yni = Gr%Teb(iG,j,k,YNQI,dN)
            xnj = Gr%Teb(iG,j,k,XNQJ,dN)
            ynj = Gr%Teb(iG,j,k,YNQJ,dN)
            detT = 1.0/(xni*ynj - xnj*yni)

            !Pre-transform version
            ! b1L = detT*( ynj*bT1L(i)/fA1(i) - yni*bT2L(i)/fA2(i) )
            ! b1R = detT*( ynj*bT1R(i)/fA1(i) - yni*bT2R(i)/fA2(i) )

            ! b2L = detT*(-xnj*bT1L(i)/fA1(i) + xni*bT2L(i)/fA2(i) )
            ! b2R = detT*(-xnj*bT1R(i)/fA1(i) + xni*bT2R(i)/fA2(i) )

            ! b1(i) = 0.5*(b1L + b1R)/fA1(i)
            ! b2(i) = 0.5*(b2L + b2R)/fA2(i)

            ! db1 = b1R-b1L
            ! db2 = b2R-b2L
            ! Jd(i) = db2 - db1

            !TODO: Consider transforming to b1/b2 before calculating diffusive current

            !Post-transform version
            bT1(i) = 0.5*( bT1L(i) + bT1R(i) )/fA1(i)
            bT2(i) = 0.5*( bT2L(i) + bT2R(i) )/fA2(i)

            !Note sign difference between 1st and 2nd term for curl
            db1 = (bT1R(i) - bT1L(i))/fA1(i)
            db2 = (bT2R(i) - bT2L(i))/fA2(i)
            Jd(i) = db2-db1

            b1(i) = detT*( ynj*bT1(i) - yni*bT2(i) )
            b2(i) = detT*(-xnj*bT1(i) + xni*bT2(i) )

            !Incorporate background field if necessary
            !Use edgB0(i,j,k,1/2,dN), already have mapped XYZ->1/2 system
            if (Model%doBackground) then
                b1(i) = b1(i) + Gr%edgB0(iG,j,k,1,dN)
                b2(i) = b2(i) + Gr%edgB0(iG,j,k,2,dN)
            endif


        enddo

    end subroutine GetCornerB

    !Get corner density at dN edge
    subroutine GetCornerD(Model,Gr,ccD,Dc,iB,j,k,iMax,dN)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(in) :: ccD(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg)
        real(rp), intent(out) :: Dc(vecLen)
        integer, intent(in) :: iB,j,k,iMax,dN

        integer, dimension(NDIM) :: e1,e2
        integer :: i,iG
        integer :: i1,j1,k1,i2,j2,k2,i12,j12,k12
        !DIR$ ASSUME_ALIGNED ccD: ALIGN
        !DIR$ ASSUME_ALIGNED Dc: ALIGN
        !Get normal plane
        call getNormalPlane(dN,e1,e2)

        do i=1,iMax
            iG = iB+i-1 !Global index
            !4-point average coordinates
            i1  = iG-e1(IDIR); j1 = j-e1(JDIR); k1 = k-e1(KDIR)
            i2  = iG-e2(IDIR); j2 = j-e2(JDIR); k2 = k-e2(KDIR)
            i12 = iG-e1(IDIR)-e2(IDIR)
            j12 = j -e1(JDIR)-e2(JDIR)
            k12 = k -e1(KDIR)-e2(KDIR)
            Dc(i) = 0.25*( ccD(iG,j,k) + ccD(i1,j1,k1) + ccD(i2,j2,k2) + ccD(i12,j12,k12) )
        enddo
    end subroutine GetCornerD


    subroutine GetCornerV(Model,Gr,Vf,v1,v2,vDiff,iB,j,k,iMax,dN,dT1,dT2)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        real(rp), intent(in) :: Vf(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM)
        real(rp), intent(inout) :: vDiff(vecLen)
        real(rp), intent(out) :: v1(vecLen), v2(vecLen)
        integer, intent(in) :: iB,j,k,iMax,dN,dT1,dT2

        integer :: i,n,d,iG
        real(rp) :: VelB(vecLen,recLen,NDIM) !Recon stencils
        real(rp) :: VxyzC(vecLen,NDIM) !Corner xyz velocities
        real(rp), dimension(vecLen,recLen) :: AreaB !Area scaling
        real(rp), dimension(vecLen) :: dA !Interpolated area
        !DIR$ ASSUME_ALIGNED Vf: ALIGN
        !DIR$ ASSUME_ALIGNED v1: ALIGN
        !DIR$ ASSUME_ALIGNED v2: ALIGN
        !DIR$ ASSUME_ALIGNED vDiff: ALIGN
        !DIR$ attributes align : ALIGN :: VelB,VxyzC,AreaB,dA
        !Initialize
        v1 = 0.0
        v2 = 0.0
        VxyzC = 0.0
        AreaB = 0.0
        dA = 0.0

        if (doVdA) then
            !Get face areas for scaling, dT1 faces in dT2 direction
            call LoadBlockI(Model,Gr,AreaB      ,Gr%Face(:,:,:,dT1),iB,j,k,iMax,dT2)
            !Interpolate first
            do i=1,iMax
                dA(i) = dot_product(interpWgt,AreaB(i,:))
            enddo
            !Loop over face Vxyz
            do d=1,NDIM
                call LoadBlock(Model,Gr,VelB(:,:,d),Vf(:,:,:,d),iB,j,k,iMax,dT2)
                do i=1,iMax
                    VxyzC(i,d) = dot_product(interpWgt,AreaB(i,:)*VelB(i,:,d))/dA(i)
                enddo
            enddo            
        else
            !Loop over face Vxyz
            do d=1,NDIM
                call LoadBlock(Model,Gr,VelB(:,:,d),Vf(:,:,:,d),iB,j,k,iMax,dT2)
                do i=1,iMax
                    VxyzC(i,d) = dot_product(interpWgt,VelB(i,:,d))
                enddo
            enddo   
        endif !do VdA scaling

        do i=1,iMax
            iG = iB+i-1
            v1(i) = VxyzC(i,XDIR)*Gr%Te(iG,j,k,TAN1X,dN) + &
                    VxyzC(i,YDIR)*Gr%Te(iG,j,k,TAN1Y,dN) + &
                    VxyzC(i,ZDIR)*Gr%Te(iG,j,k,TAN1Z,dN)
            v2(i) = VxyzC(i,XDIR)*Gr%Te(iG,j,k,TAN2X,dN) + &
                    VxyzC(i,YDIR)*Gr%Te(iG,j,k,TAN2Y,dN) + &
                    VxyzC(i,ZDIR)*Gr%Te(iG,j,k,TAN2Z,dN)
            vDiff(i) = sqrt( v1(i)**2.0 + v2(i)**2.0 )
        enddo

    end subroutine GetCornerV


    subroutine resistivity(Model,Gr,State,E)
        type(Model_T), intent(in) :: Model
        type(Grid_T), intent(in) :: Gr
        type(State_T), intent(in) :: State
        real (rp), dimension(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,NDIM), intent(inout) :: E

        integer :: i,j,k
        ! Cell centerd values
        real(rp), dimension(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM) :: Btot, Jcc
        real(rp), dimension(NDIM)  :: Jedge, e1,e2

        
        call TIC("EtaBtot")
        !$OMP PARALLEL DO default (shared) collapse(2) &
        !$OMP private(i,j,k)
        do k=Gr%ksg, Gr%keg
            do j=Gr%jsg, Gr%jeg
                do i=Gr%isg, Gr%ieg
                    Btot(i,j,k,:) = Gr%B0(i,j,k,:) + State%Bxyz(i,j,k,:)
                end do
            end do
        end do

        call TOC("EtaBtot")
        
        call TIC("Jcell")
        call bFld2Jxyz(Model,Gr,Btot, Jcc)
        call TOC("Jcell")
        

        call TIC("Eeta")
        !Loop over active edges and calculate edge current
        !Use resistivity and edge length to correct electric field
        !NOTE: cell-centered current is Jxyz, but need IJK E fields for CT update
        !$OMP PARALLEL DO default (shared) collapse(2) &
        !$OMP private(i,j,k,e1,e2,Jedge) 
        do k=Gr%ks-1, Gr%ke+1
            do j=Gr%js-1, Gr%je+1
                do i=Gr%is-1,Gr%ie +1
                    !For each direction, use 4-point cell-centered average for edge Jxyz

                    !---- Ei calculation
                    Jedge(XDIR:ZDIR) = 0.25*( Jcc(i  ,j-1,k-1,:) + Jcc(i  ,j  ,k-1,:) + &
                                              Jcc(i  ,j  ,k  ,:) + Jcc(i  ,j-1,k  ,:) )
                    call edgeCoords(Model,Gr,i,j,k,IDIR,e1,e2)
                    E(i,j,k,IDIR) = E(i,j,k,IDIR) + State%Deta(i,j,k,IDIR)*dot_product(Jedge,e2-e1)

                    !---- Ej calculation
                    Jedge(XDIR:ZDIR) = 0.25*( Jcc(i-1,j  ,k-1,:) + Jcc(i  ,j  ,k-1,:) + &
                                              Jcc(i  ,j  ,k  ,:) + Jcc(i-1,j  ,k  ,:) )
                    call edgeCoords(Model,Gr,i,j,k,JDIR,e1,e2)
                    E(i,j,k,JDIR) = E(i,j,k,JDIR) + State%Deta(i,j,k,JDIR)*dot_product(Jedge,e2-e1)

                    !---- Ek calculation
                    Jedge(XDIR:ZDIR) = 0.25*( Jcc(i-1,j-1,k  ,:) + Jcc(i  ,j-1,k  ,:) + &
                                              Jcc(i  ,j  ,k  ,:) + Jcc(i-1,j  ,k  ,:) )
                    call edgeCoords(Model,Gr,i,j,k,KDIR,e1,e2)
                    E(i,j,k,KDIR) = E(i,j,k,KDIR) + State%Deta(i,j,k,KDIR)*dot_product(Jedge,e2-e1)

                 
                enddo
            enddo
        enddo

        call TOC("Eeta")
        
    end subroutine resistivity
end module fields
