!Routines for dst-type calculations

module dstutils
    use gamtypes
    use earthhelper
    use volttypes
    use gridutils
    use gamutils
    
    implicit none

    integer, parameter, private :: i0 = 1 !What shell to start at (i0+1) for BS-Dst    
    real(rp), private :: B0 = 1.0 !Scaling factor for Jxyz=>Dst
    real(rp), parameter, private :: BSMaxR = 40.0 !Sphere to do quick BS integral over

    !Some lazy hard-coded stations (SM lat/lon [deg])
    !NOTE/TODO: Right now just calculating dBz (instead of dBn), so stick to equatorial stations
    integer , parameter, private :: NumStat = 4 !Number of stations
    real(rp), parameter, private :: SR0 = 1.0 !Station radius in Re
    real(rp), dimension(NumStat) :: SLats = [  0.0,  0.0,  0.0,  0.0]
    real(rp), dimension(NumStat) :: SLons = [  0.0, 90.0,180.0,270.0]

    ! integer, parameter, private :: NumStat = 11 !Number of stations
    ! real(rp), dimension(NumStat) :: SLats = [+28.04,+48.14,+48.24,+39.73,+21.71,+35.63,+34.34,+10.37,-46.22,-34.08,+49.75]
    ! real(rp), dimension(NumStat) :: SLons = [  6.54,353.93,321.28,316.74,270.27,211.74,162.53,146.55,144.93, 84.63, 85.80]

    contains

    !Use Gamera data to estimate DST
    !(Move this to msphutils?)
    !BSSMRs holds 4 quadrant estimate: 12,18,00,06
    subroutine EstDST(Model,Gr,State,BSDst0,AvgBSDst,BSSMRs,DPSDst)
        type(Model_T), intent(in)  :: Model
        type(Grid_T) , intent(in)  :: Gr
        type(State_T), intent(in)  :: State
        real(rp)     , intent(out) :: BSDst0,DPSDst,AvgBSDst,BSSMRs(4)

        integer :: i,j,k,n
        real (rp), dimension(:,:,:,:), allocatable :: dB,Jxyz !Full-sized arrays
        real(rp), dimension(NDIM) :: xyz0
        real(rp) :: d0,u0,phi,lat
        real(rp), dimension(NumStat) :: StatDSTs

        BSDst0   = 0.0
        DPSDst   = 0.0
        AvgBSDst = 0.0

        !Very lazy scaling (TODO: Pull from Model%xxx)
        d0 = (1.67e-27)*1.0e+6
        u0 = 1.0e+5
        B0 = sqrt(Mu0*d0*u0*u0)*1.0e+9 !nT

        call allocGridVec(Model,Gr,dB  )
        call allocGridVec(Model,Gr,Jxyz)

    !Calculate current
        !Subtract dipole before calculating current
        !$OMP PARALLEL DO default(shared) collapse(2)
        do k=Gr%ksg,Gr%keg
            do j=Gr%jsg,Gr%jeg
                do i=Gr%isg,Gr%ieg
                    dB(i,j,k,:) = State%Bxyz(i,j,k,:) - MagsphereDipole(Gr%xyzcc(i,j,k,:),Model%MagM0)
                    if (Model%doBackground) then
                        dB(i,j,k,:) = dB(i,j,k,:) + Gr%B0(i,j,k,:)
                    endif
                enddo
            enddo
        enddo
        call bFld2Jxyz(Model,Gr,dB,Jxyz)

    !Get Dst's
        !DPS Dst from ring current energy density
        DPSDst = CalcDPSDst(Model,Gr)
        !Quick Dst from center of earth
        xyz0 = 0.0 !Measure at center of Earth
        BSDst0  = BSDstAt(Model,Gr,Jxyz,xyz0)

        !Get simple station-averaged Dst (using Bn ~ zhat at equator)         
        do n=1,NumStat
            phi = (PI/180.0)*SLons(n)
            lat = (PI/180.0)*SLats(n)
            xyz0(XDIR) = SR0*cos(lat)*cos(phi)
            xyz0(YDIR) = SR0*cos(lat)*sin(phi)
            xyz0(ZDIR) = SR0*sin(lat)
            StatDSTs(n) = BSDstAt(Model,Gr,Jxyz,xyz0)
        enddo
        AvgBSDst = sum(StatDSTs)/NumStat
        
        if (NumStat /= 4) then
            write(*,*) "Fix EstDST ..."
            stop
        else
            !Assuming ordering, 12/18/00/06
            BSSMRs = StatDSTs
        endif
    end subroutine EstDST

    !Calculate Bios-Savart Dst at given point
    function BSDstAt(Model,Gr,Jxyz,xyz0) result(BSDst)
        type(Model_T), intent(in)  :: Model
        type(Grid_T) , intent(in)  :: Gr
        real(rp)     , intent(in) :: Jxyz(Gr%isg:Gr%ieg,Gr%jsg:Gr%jeg,Gr%ksg:Gr%keg,1:NDIM)
        real(rp)     , intent(in) :: xyz0(NDIM)
        real(rp) :: BSDst

        integer :: i,j,k
        real(rp) :: dV,r,bs1,bs2,bScl,dBz
        real(rp), dimension(NDIM) :: xyz

        BSDst  = 0.0 !Bios-Savart

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,xyz,dV,r,bs1,bs2,bScl,dBz) &
        !$OMP reduction(+:BSDst)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is+i0,Gr%ie
                    xyz = Gr%xyzcc (i,j,k,:)
                    dV  = Gr%volume(i,j,k)
                    r = norm2(xyz-xyz0)
                    if (r <= BSMaxR) then
                        bs1 = Jxyz(i,j,k,XDIR)*(xyz(YDIR)-xyz0(YDIR))
                        bs2 = Jxyz(i,j,k,YDIR)*(xyz(XDIR)-xyz0(XDIR))
                        bScl = B0*dV/(4*PI)

                        dBz = -(bs1 - bs2)/(r**3.0)
                    else
                        dBz = 0.0
                    endif

                    !Bios-Savart Dst
                    BSDst = BSDst + bScl*dBz
                enddo ! i loop
            enddo
        enddo !k loop

    end function BSDstAt

    !Calculate DPS Dst
    function CalcDPSDst(Model,Gr) result(DPSDst)
        type(Model_T), intent(in)  :: Model
        type(Grid_T) , intent(in)  :: Gr
        real(rp) :: DPSDst

        integer :: i,j,k
        real(rp) :: dV,KTot

        DPSDst = 0.0 !DPS
        KTot = 0.0 !Total energy content

        if (.not. Model%doSource) return

        !$OMP PARALLEL DO default(shared) collapse(2) &
        !$OMP private(i,j,k,dV) &
        !$OMP reduction(+:KTot)
        do k=Gr%ks,Gr%ke
            do j=Gr%js,Gr%je
                do i=Gr%is+i0,Gr%ie
                    dV  = Gr%volume(i,j,k)
                    if (Gr%Gas0(i,j,k,IMPR ,BLK)>TINY) then
                        KTot = KTot + dV*Gr%Gas0(i,j,k,IMPR ,BLK) !Code units
                    endif

                enddo ! i loop
            enddo
        enddo !k loop

        !Scale code units to energy
        !KTot = x0^3 x P0
        KTot = (Model%Units%gx0**3.0)*(Model%Units%gP0)*KTot
        !KTot = nPa*m3, now go to keV
        KTot = ((1.0e-9)/kev2J)*KTot
        DPSDst = -4.2*(1.0e-30)*KTot
    end function CalcDPSDst

end module dstutils