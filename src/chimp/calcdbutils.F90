module calcdbutils
	use kdefs
	use chmpdefs
	use ebtypes
    use calcdbtypes

	implicit none

	real(rp) :: dzGG = 30.0 !Default height spacing [km]
    logical, private, parameter :: doHall = .true.
    logical, private, parameter :: doPed  = .true.
    logical, private, parameter :: doAmm  = .true.

    integer, private, parameter :: TDIR=1,PDIR=2
    integer, private, parameter :: Ngm = 4 !Stencil for gradient

	contains

	subroutine facGridInit(Model,ebState,rmState,ionGrid,facGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(ionGrid_T)  , intent(in) :: ionGrid
        type(facGrid_T)  , intent(inout) :: facGrid

        integer :: Np,Nth,n,i,j,l
        real(rp) :: rMHD,rIon,dr,dV
        real(rp), dimension(:), allocatable :: rC
        real(rp) :: phi,theta,ionlat,c2lat,lat
        real(rp) :: latC,thGap,thIon,dl
        real(rp) :: latP,latM,dlat
        Np = rmState%Np
        Nth = rmState%Nth

        !Allocate and zero out facGrid
        allocate(facGrid%XYZcc(Np,Nth,rSegs,2,NDIM))
        allocate(facGrid%Jxyz (Np,Nth,rSegs,2,NDIM))
        allocate(facGrid%dV   (Np,Nth,rSegs,2))
        allocate(facGrid%pcc  (Np,Nth,2)) !Cell-centered phi
        allocate(facGrid%tcc  (Np,Nth,2)) !Cell-centered theta

        facGrid%XYZcc = 0.0
        facGrid%Jxyz  = 0.0
        facGrid%dV    = 0.0
        facGrid%Np = Np
        facGrid%Nth = Nth
        facGrid%rSegs = rSegs
        facGrid%pcc = ionGrid%pcc
        facGrid%tcc = ionGrid%tcc
        facGrid%dp = ionGrid%dp
        facGrid%dt = ionGrid%dt

        !Get radial spacing
        rMHD = norm2(ebState%ebGr%xyz(1,1,1,XDIR:ZDIR))
        !Using the current closure in the 2nd cell center
        !rMHD = norm2(ebState%ebGr%xyzcc(2,1,1,XDIR:ZDIR))

        rIon = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880

        dr = (rMHD-rIon)/rSegs
        allocate(rC(rSegs))

        do n=1,rSegs
            rC(n) = rIon + 0.5*dr + (n-1)*dr
        enddo

        !Get cell centers of segments
        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,n,l,phi,theta,ionlat,c2lat,lat)
        do l=1,2
            do n=1,rSegs
                do j=1,Nth
                    do i=1,Np
                        phi   = facGrid%pcc(i,j,l) !N/S hemisphere
                        theta = facGrid%tcc(i,j,l) !N/S hemisphere
                        ionlat = PI/2 - theta

                        c2lat = (rIon/rC(n)) * (cos(ionlat)**2.0)
                        lat = acos(sqrt(c2lat))
                        if ( (l==SOUTH) .and. (lat>0) ) then
                            lat = -lat
                        endif
                        theta = PI/2-lat
                  
                        facGrid%XYZcc(i,j,n,l,XDIR) = rC(n)*sin(theta)*cos(phi)
                        facGrid%XYZcc(i,j,n,l,YDIR) = rC(n)*sin(theta)*sin(phi)
                        facGrid%XYZcc(i,j,n,l,ZDIR) = rC(n)*cos(theta)
                        
                    enddo !i

                enddo !j
            enddo !n
        enddo !l

        !J dA dl = J_ion dA_ion dl = J11_ion dA_ion cos(dip) dl
        !dV = dA_ion cos(dip) dl
        !Get dV, do northern hemisphere then copy

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(j,n,latC,thGap,thIon,dl,dV)
        do n=1,rSegs
            do j=1,Nth
                latC = xyz2lat(facGrid%XYZcc(1,j,n,NORTH,:))
                thGap = PI/2 - latC !theta of this cell-center in gap region
                thIon = facGrid%tcc(1,j,NORTH) !Ionospheric theta

                !Arc length along dipole at this radius/theta
                dl = dr * sqrt(1.0 + 3.0*(cos(thGap)**2.0))/2/cos(thGap)
                dV = dl*ionGrid%dS(1,j,NORTH)*abs(CosDip(thIon))

                facGrid%dV(:,j,n,:) = dV
            enddo
        enddo

	end subroutine facGridInit

    subroutine ionGridInit(Model,ebState,rmState,ionGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(ionGrid_T)  , intent(inout) :: ionGrid

        integer :: Np,Nth
        integer :: i,j
        real(rp), dimension(:,:), allocatable :: Z
        real(rp) :: R0,theta,phi,thesp,phisp
        real(rp) :: phi0,dp,th0,dth,th1,th2

        Np = rmState%Np
        Nth = rmState%Nth

        allocate(ionGrid%XYZcc(Np,Nth,2,NDIM))
        allocate(ionGrid%Jxyz (Np,Nth,2,NDIM))
        allocate(ionGrid%Etp  (Np,Nth,2,TDIR:PDIR))
        allocate(ionGrid%hJ   (Np,Nth,2,TDIR:PDIR))
        allocate(ionGrid%pJ   (Np,Nth,2,TDIR:PDIR))
        allocate(ionGrid%dS   (Np,Nth,2)) !Surface area per patch, Re^2
        allocate(ionGrid%pcc  (Np,Nth,2)) !Cell-centered phi
        allocate(ionGrid%tcc  (Np,Nth,2)) !Cell-centered theta

        ionGrid%XYZcc = 0.0
        ionGrid%Jxyz  = 0.0
        ionGrid%dS    = 0.0
        ionGrid%Np = Np
        ionGrid%Nth = Nth

        allocate(Z(Np,Nth))
        Z = sqrt(1.0 - rmState%XY(:,:,XDIR)**2.0 - rmState%XY(:,:,YDIR)**2.0)
        !Note: Z is in Rion, not Re

        !Get angular spacing & starting point
        phi0 = atan2(rmState%XY(1,1,YDIR),rmState%XY(1,1,XDIR))
        dp = 2*PI/Np

        th0 = acos(Z(1,1))
        dth = acos(Z(1,2)) - acos(Z(1,1))
        !Doing some dumb angle rounding to force sensible values from janky remix grid
        th0 = RoundAng(th0)
        dth = RoundAng(dth)

        write(*,*) 'ReMIX Low-latitude = ', ang2deg(acos(Z(1,Nth))),ang2deg(RoundAng(acos(Z(1,Nth))))

        ionGrid%dt = dth
        ionGrid%dp = dp

        R0 = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880

        !Get coordinates and dS
        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,phi,theta,thesp,phisp,th1,th2)
        do j=1,ionGrid%Nth !Theta

            !Get NH theta (cc)
            theta = th0 + 0.5*dth + dth*(j-1)
            th1 = theta-0.5*dth
            th2 = theta+0.5*dth

            !Fix hole at pole
            if (j == 1) then
                th1 = 0.0 !Make sure highest cell goes up to pole
            endif

            do i=1,ionGrid%Np !phi
                !cell-centered phi
                phi = phi0 + 0.5*dp + dp*(i-1)

                !Store northern hemisphere quantities
                ionGrid%XYZcc(i,j,NORTH,XDIR) = R0*sin(theta)*cos(phi)
                ionGrid%XYZcc(i,j,NORTH,YDIR) = R0*sin(theta)*sin(phi)
                ionGrid%XYZcc(i,j,NORTH,ZDIR) = R0*cos(theta)

                !Get southern hemisphere angles
                thesp = PI - theta !Southern hemisphere theta
                phisp = 2*PI-phi !Southern hemisphere MLT

                !Store southern hemisphere quantities
                ionGrid%XYZcc(i,j,SOUTH,XDIR) = R0*sin(thesp)*cos(phisp)
                ionGrid%XYZcc(i,j,SOUTH,YDIR) = R0*sin(thesp)*sin(phisp)
                ionGrid%XYZcc(i,j,SOUTH,ZDIR) = R0*cos(thesp)

                !Store surface area elements
                ionGrid%dS(i,j,NORTH) = (R0**2.0)*sin(theta)*(th2-th1)*dp
                ionGrid%dS(i,j,SOUTH) = (R0**2.0)*sin(thesp)*(th2-th1)*dp

                !Store cell-centered angles
                ionGrid%pcc(i,j,NORTH) = phi
                ionGrid%pcc(i,j,SOUTH) = phisp
                ionGrid%tcc(i,j,NORTH) = theta
                ionGrid%tcc(i,j,SOUTH) = thesp

            enddo
        enddo

        contains 
            function RoundAng(inAng)
                real(rp), intent(in) :: inAng
                real(rp) :: RoundAng
                real(rp) :: inAngD
                inAngD = inAng*180.0/PI
                RoundAng = (PI/180.0)*nint(inAngD*8.0)/8.0
            end function RoundAng

    end subroutine ionGridInit

    !Initialize holders for Bios-Savart contributions
    subroutine BSGridInit(Model,ebState,rmState,magBS,ionBS,facBS)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(BSGrid_T), intent(inout) :: magBS,ionBS,facBS

        integer :: NMag,NIon,NFac
        real(rp) :: B0
    !Mag BS grid
        B0 = bScale()
        !Number of cells
        NMag = (ebState%ebGr%Nip)*(ebState%ebGr%Njp)*(ebState%ebGr%Nkp)
        call BSSubInit(magBS,NMag)
        magBS%jScl = B0/(4.0*PI) !Scaling factor

    !Ion BS grid
        !Number of cells
        NIon = (rmState%Np)*(rmState%Nth)*(2) !Include N/S hemispheres
        call BSSubInit(ionBS,NIon)
        ionBS%jScl = (1.0e+9)*Mu0/(4.0*PI) !Ensure final db is nT

    !FAC BS grid
        !Number of cells
        NFac = NIon*rSegs
        call BSSubInit(facBS,NFac)
        facBS%jScl = (1.0e+9)*REarth*Mu0/(4.0*PI) !Ensure final db is nT & needs extra factor of Re

    end subroutine BSGridInit

    !Initialize BS grid w/ N points
    subroutine BSSubInit(xBS,N)
        type(BSGrid_T), intent(inout) :: xBS
        integer, intent(in) :: N
        xBS%NumP = N
        allocate(xBS%XYZcc(N,NDIM))
        allocate(xBS%Jxyz (N,NDIM))
        allocate(xBS%dV(N))
        xBS%XYZcc = 0.0
        xBS%Jxyz = 0.0
        xBS%dV = 0.0
        
    end subroutine BSSubInit

	!Using a rmState (remix data), fill facGrid Jxyz
	subroutine facGridUpdate(Model,ebState,rmState,facGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(facGrid_T)  , intent(inout) :: facGrid

        integer :: i,j,n,l
        real(rp), dimension(NDIM) :: xC,bhat
        real(rp) :: rIon,J11i,J11m,Biom

        facGrid%Jxyz = 0.0

        rIon = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880

        !$OMP PARALLEL DO default(shared) collapse (3) &
        !$OMP private(i,j,n,l,xC,bhat) &
        !$OMP private(J11i,J11m,Biom)
        do l=1,2
            do n=1,facGrid%rSegs
                do j=1,facGrid%Nth
                    do i=1,facGrid%Np
                        xC = facGrid%XYZcc(i,j,n,l,XDIR:ZDIR)
                        
                        Biom = BionoBm(xC) !Ratio of ionosphere B field and this point
                        !Get J11,ion stored in muA/m2
                        if (l == NORTH) then
                            J11i = rmState%nFac(i,j)
                        else
                            J11i = rmState%sFac(i,j)
                        endif

                        !Calculate current in A/m2
                        J11i = (1.0e-6)*J11i
                        !Use same gap region J11 as ion, scaling has been included in dV
                        J11m = J11i

                        !Now get local bhat
                        bhat = bhatXYZ(xC)

                        !J11,bhat,amplification factor
                        facGrid%Jxyz(i,j,n,l,:) = J11m*bhat

                    enddo !i
                enddo !j
            enddo !k
        enddo !l

        contains

            function bhatXYZ(xyz)
                real(rp), intent(in) :: xyz(NDIM)
                real(rp) :: bhatXYZ(NDIM)

                real(rp) :: x,y,z,r,r2,A

                x = xyz(XDIR)
                y = xyz(YDIR)
                z = xyz(ZDIR)
                r = norm2(xyz)
                r2 = r**2.0
                A = sqrt( 1 + 3.0*(z/r)**2.0 )

                bhatXYZ(XDIR) = -(1/A)*(3*x*z)/r2
                bhatXYZ(YDIR) = -(1/A)*(3*y*z)/r2
                bhatXYZ(ZDIR) =  (1/A)*( -2 + 3*(x**2.0+y**2.0)/r2 )

            end function bhatXYZ

            !Bion/Bm ratio, see mhd2mix code for explanation
            function BionoBm(xyz)
                real(rp), intent(in) :: xyz(NDIM)
                real(rp) :: BionoBm

                real(rp) :: x,y,z,rMag,Rp,zor2,mZ,pZ

                x = xyz(XDIR)
                y = xyz(YDIR)
                z = xyz(ZDIR)
                rMag = norm2(xyz)
                Rp = rMag/rIon
                
                zor2 = (z/rMag)**2.0
                mZ = (1.0-zor2)/Rp
                pZ = 1.0+3.0*zor2

                BionoBm = (Rp**3.0)*sqrt( 1.0 + 3.0*(1.0-mZ) )/sqrt(pZ)

            end function BionoBm

    end subroutine facGridUpdate

    !Using a rmState (remix data), fill ionGrid Jxyz
    subroutine ionGridUpdate(Model,ebState,rmState,ionGrid)
        type(chmpModel_T), intent(in) :: Model
        type(ebState_T)  , intent(in) :: ebState
        type(rmState_T)  , intent(in) :: rmState
        type(ionGrid_T)  , intent(inout) :: ionGrid

        real(rp) :: nEp,nEt,sEp,sEt
        integer :: i,j
        real(rp) :: R0,dth,dph,thnh,thsh,phnh,phsh
        real(rp) :: nhcdip,shcdip,nJt,sJt,nJp,sJp

        real(rp), dimension(NDIM) :: nJ,sJ

        !Need to calculate E field from potential in both N/S hemispheres
        !Then get J from E and SigP/SigH
        ionGrid%Jxyz = 0.0

        R0 = (RionE*1.0e+6)/REarth !Ionospheric radius in units of Re, ~1.01880
        dth = ionGrid%dt
        dph = ionGrid%dp

        !$OMP PARALLEL DO default(shared) &
        !$OMP private(i,j,nJ,sJ) &
        !$OMP private(nEp,nEt,sEp,sEt,thnh,thsh,phnh,phsh) &
        !$OMP private(nhcdip,shcdip,nJt,sJt,nJp,sJp)
        do j=1,ionGrid%Nth !Theta
            do i=1,ionGrid%Np !phi
            !Start w/ E field via gradient    
            !E = -grad(pot), kV/Re

                !Theta derivatives
                nEt = -GradTheta(rmState%nPot(i,:),j,ionGrid%Nth)/(R0*dth)
                sEt = -GradTheta(rmState%sPot(i,:),j,ionGrid%Nth)/(R0*dth)
                sEt = -sEt !Flip direction

                !Phi derivatives
                thnh = ionGrid%tcc(i,j,NORTH)
                thsh = ionGrid%tcc(i,j,SOUTH)

                phnh = ionGrid%pcc(i,j,NORTH)
                phsh = ionGrid%pcc(i,j,SOUTH)

                nEp = -GradPhi(rmState%nPot(:,j),i,ionGrid%Np)/(R0*sin(thnh)*dph)
                sEp = -GradPhi(rmState%sPot(:,j),i,ionGrid%Np)/(R0*sin(thsh)*dph)
                sEp = -sEp !Flip direction

                !Convert E [kV/Re] to [V/m]
                nEp = (1.0e+3)*nEp/REarth
                nEt = (1.0e+3)*nEt/REarth
                sEp = (1.0e+3)*sEp/REarth
                sEt = (1.0e+3)*sEt/REarth

                !Store E fields
                ionGrid%Etp(i,j,NORTH,TDIR) = nEt
                ionGrid%Etp(i,j,NORTH,PDIR) = nEp
                ionGrid%Etp(i,j,SOUTH,TDIR) = sEt
                ionGrid%Etp(i,j,SOUTH,PDIR) = sEp

            !Now have E fields [V/m], calculate currents
                !Get cos of dip angles for both
                !CALCDB-TODO: Check dip angle formula in SH?
                if (doAmm) then
                    nhcdip = -2*cos(thnh)/sqrt(1.0 + 3*cos(thnh)*cos(thnh))
                    shcdip = -2*cos(thsh)/sqrt(1.0 + 3*cos(thsh)*cos(thsh))
                else
                    nhcdip = -1.0
                    shcdip = +1.0
                endif

                !Get theta/phi currents from Hall/Pederson
                !Conductance units are S = A/V, so currents are J = A/m
                ionGrid%hJ(i,j,NORTH,TDIR) = -rmState%nSigH(i,j)*nEp/nhcdip
                ionGrid%hJ(i,j,NORTH,PDIR) =  rmState%nSigH(i,j)*nEt/nhcdip
                ionGrid%pJ(i,j,NORTH,TDIR) =  rmState%nSigP(i,j)*nEt/nhcdip**2.0
                ionGrid%pJ(i,j,NORTH,PDIR) =  rmState%nSigP(i,j)*nEp

                ionGrid%hJ(i,j,SOUTH,TDIR) = -rmState%sSigH(i,j)*sEp/shcdip
                ionGrid%hJ(i,j,SOUTH,PDIR) =  rmState%sSigH(i,j)*sEt/shcdip
                ionGrid%pJ(i,j,SOUTH,TDIR) =  rmState%sSigP(i,j)*sEt/shcdip**2.0
                ionGrid%pJ(i,j,SOUTH,PDIR) =  rmState%sSigP(i,j)*sEp


                nJt = 0.0 ; nJp = 0.0 ; sJt = 0.0 ; sJp = 0.0
                if (doHall) then
                    nJt = nJt + ionGrid%hJ(i,j,NORTH,TDIR)
                    nJp = nJp + ionGrid%hJ(i,j,NORTH,PDIR)
                    sJt = sJt + ionGrid%hJ(i,j,SOUTH,TDIR)
                    sJp = sJp + ionGrid%hJ(i,j,SOUTH,PDIR)
                endif
                if (doPed) then
                    nJt = nJt + ionGrid%pJ(i,j,NORTH,TDIR)
                    nJp = nJp + ionGrid%pJ(i,j,NORTH,PDIR)
                    sJt = sJt + ionGrid%pJ(i,j,SOUTH,TDIR)
                    sJp = sJp + ionGrid%pJ(i,j,SOUTH,PDIR)
                endif


            !Now have currents (theta/phi) [A/m], convert to XYZ and store
                nJ = tp2xyz(phnh,thnh,nJt,nJp)
                sJ = tp2xyz(phsh,thsh,sJt,sJp)

                ionGrid%Jxyz(i,j,NORTH,:) = nJ
                ionGrid%Jxyz(i,j,SOUTH,:) = sJ
            enddo

        enddo !j,Nth

        ! write(*,*) "North hJ-Theta"
        ! write(*,*) ionGrid%hJ(:,ionGrid%Nth,NORTH,TDIR)

        write(*,*) "North pJ-Theta"
        write(*,*) minval(ionGrid%pJ(:,ionGrid%Nth,NORTH,TDIR)),maxval(ionGrid%pJ(:,ionGrid%Nth,NORTH,TDIR))

        ! write(*,*) "South hJ-Theta"
        ! write(*,*) ionGrid%hJ(:,ionGrid%Nth,SOUTH,TDIR)

        write(*,*) "South pJ-Theta"
        write(*,*) minval(ionGrid%pJ(:,ionGrid%Nth,SOUTH,TDIR)),maxval(ionGrid%pJ(:,ionGrid%Nth,SOUTH,TDIR))

        ! write(*,*) 'Hall '
        ! write(*,*) '   Theta: ',minval(ionGrid%hJ(:,:,NORTH,TDIR)),maxval(ionGrid%hJ(:,:,NORTH,TDIR))
        ! write(*,*) '   Phi  : ',minval(ionGrid%hJ(:,:,NORTH,PDIR)),maxval(ionGrid%hJ(:,:,NORTH,PDIR))

        ! write(*,*) 'Pede '
        ! write(*,*) '   Theta: ',minval(ionGrid%pJ(:,:,NORTH,TDIR)),maxval(ionGrid%pJ(:,:,NORTH,TDIR))
        ! write(*,*) '   Phi  : ',minval(ionGrid%pJ(:,:,NORTH,PDIR)),maxval(ionGrid%pJ(:,:,NORTH,PDIR))

        contains
            function GradTheta(Q,i0,Ni) result(Qp)
                real(rp), intent(in) :: Q(Ni)
                integer , intent(in) :: i0,Ni
                real(rp) :: Qp

                real(rp) :: Qblk(Ngm),c(Ngm)
                integer :: is,ie
                is = 1
                ie = Ni
                if (i0 == is) then
                    !Forward
                    Qblk = [Q(is),Q(is+1),Q(is+2),Q(is+3)]
                    c = [-11.0,18.0,-9.0,2.0]/6.0
                else if (i0 == is+1) then
                    !1 back
                    Qblk = [Q(is),Q(is+1),Q(is+2),Q(is+3)]
                    c = [-2.0,-3.0,6.0,-1.0]/6.0
                else if (i0 == ie) then
                    Qblk = [Q(ie-3),Q(ie-2),Q(ie-1),Q(ie)]
                    c = [-2.0,9.0,-18.0,11.0]/6.0
                else if (i0 == ie-1) then
                    Qblk = [Q(ie-3),Q(ie-2),Q(ie-1),Q(ie)]
                    c = [1.0,-6.0,3.0,2.0]/6.0
                else
                    !Centered
                    Qblk = [Q(i0-2),Q(i0-1),Q(i0+1),Q(i0+2)]
                    c = [1.0,-8.0,8.0,-1.0]/12.0
                endif
                Qp = dot_product(Qblk,c)
            end function GradTheta

            function GradPhi(Q,i0,Ni) result(Qp)
                real(rp), intent(in) :: Q(Ni)
                integer , intent(in) :: i0,Ni
                real(rp) :: Qp

                real(rp) :: Qblk(Ngm),c(Ngm)
                integer :: ip,ipp,im,imm
                integer :: is,ie
                is = 1
                ie = Ni

                !Always centered b/c of periodicity
                c = [1.0,-8.0,8.0,-1.0]/12.0
                if (i0 == is) then
                    ip  = i0+1
                    ipp = i0+2
                    im  = ie
                    imm = ie-1
                    
                else if (i0 == is+1) then
                    ip  = i0+1
                    ipp = i0+2
                    im  = is
                    imm = ie
                else if (i0 == ie) then
                    ip  = is
                    ipp = is+1
                    im  = i0-1
                    imm = i0-2

                else if (i0 == ie-1) then
                    ip  = ie
                    ipp = is
                    im  = i0-1
                    imm = i0-2
                else
                    !Centered
                    ip  = i0+1
                    ipp = i0+2
                    im  = i0-1
                    imm = i0-2

                    Qblk = [Q(i0-2),Q(i0-1),Q(i0+1),Q(i0+2)]
                    c = [1.0,-8.0,8.0,-1.0]/12.0
                endif
                Qblk = [Q(imm),Q(im),Q(ip),Q(ipp)]
                Qp = dot_product(Qblk,c)
            end function GradPhi

    end subroutine ionGridUpdate

    !Convert theta-phi vector to Jxyz
    function tp2xyz(phi,theta,Jt,Jp) result(Jxyz)
        real(rp), intent(in) :: phi,theta,Jt,Jp

        real(rp), dimension(NDIM) :: Jxyz
        real(rp), dimension(NDIM) :: phat,that
        
        phat = [-sin(phi),cos(phi),0.0_rp]
        that = [cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)]

        Jxyz = that*Jt + phat*Jp

    end function tp2xyz

    function CosDip(theta) result(cosd)
        real(rp), intent(in) :: theta
        real(rp) :: cosd
        cosd = -2*cos(theta)/sqrt(1.0 + 3*cos(theta)*cos(theta))
    end function CosDip

    !Set rmState given properly set 4 hemispheres and temporal weights
    subroutine hemi2rm(rmState,w1,w2)
        type(rmState_T)  , intent(inout) :: rmState
        real(rp), intent(in) :: w1,w2

        rmState%nFac  = w1*rmState%rmN1%xFac  + w2*rmState%rmN2%xFac 
        rmState%nPot  = w1*rmState%rmN1%xPot  + w2*rmState%rmN2%xPot 
        rmState%nSigP = w1*rmState%rmN1%xSigP + w2*rmState%rmN2%xSigP
        rmState%nSigH = w1*rmState%rmN1%xSigH + w2*rmState%rmN2%xSigH
        rmState%sFac  = w1*rmState%rmS1%xFac  + w2*rmState%rmS2%xFac 
        rmState%sPot  = w1*rmState%rmS1%xPot  + w2*rmState%rmS2%xPot 
        rmState%sSigP = w1*rmState%rmS1%xSigP + w2*rmState%rmS2%xSigP
        rmState%sSigH = w1*rmState%rmS1%xSigH + w2*rmState%rmS2%xSigH

    end subroutine hemi2rm

    !Lazy function to return scaling, should be replaced
    function bScale() result(B0)
        real(rp) :: B0
        real(rp) :: mu0,Mp,x0,u0,t0,d0,p0
        mu0 = 4*PI*1e-7
        Mp = 1.67e-27 ![kg]
        x0 = 1*6.38e6 ![m]   - RE
        u0 = 100e3    ![m/s] - 100 km/s
        t0 = x0/u0 ![s]   -
        d0 = Mp*1e6 ! [kg/m^3] - 1 particle/cc
        p0 = d0*u0*u0 ![N/m^2]
        B0 = sqrt(mu0*d0*u0*u0)*1e9 ! [nT]
    end function bScale

    function xyz2lat(xyz) result(lat)
        real(rp), dimension(NDIM), intent(in) :: xyz

        real(rp) :: lat
        real(rp) :: z,r
        r = norm2(xyz)
        z = xyz(ZDIR)
        lat = abs(asin(z/r))

    end function xyz2lat

end module calcdbutils
