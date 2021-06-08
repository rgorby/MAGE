!Utilities for loss calculations

MODULE lossutils
    USE kdefs, ONLY : TINY,PI,Mp_cgs,kev2J
    USE rcm_precision
    USE rcmdefs
    use math, ONLY : SmoothOpTSC,SmoothOperator33,ClampValue,LinRampUp

    implicit none

    contains

    !Loss functions
    FUNCTION CXKaiju(isp,enrg,rloc) result(cxrate)
        IMPLICIT NONE

        integer(iprec), intent(in) :: isp
        real(rprec), intent(in) :: enrg,rloc

        real(rprec) :: cxrate
        real(rprec) :: K,Ngeo,Sig,M,Kj,V,Tau,tScl

        K = enrg*1.0e-3 !Energy in kev
    !Geocoronal density [#/cc]
        Ngeo = OstgaardGeocorona(rloc)
      
    !K in keV, Sig in cm2
        !Using Lindsay & Stebbings 2005
        if (isp == RCMPROTON) then
            Sig = CXSigma(K,isp)
            M = (Mp_cgs)*(1.0e-3) !Proton mass
        else
            write(*,*) 'Unknown charge exchange species, bailing!'
            stop
        endif

    !Get velocity [cm/s] from energy [keV]
        Kj = kev2J*K !Joules
        V = sqrt(2*Kj/M)*100.0 !m/s->cm/s

    !Timescale
        tScl = cos(45*PI/180.0)**3.5 !Using Smith & Bewtra 1976 scaling
        Tau = tScl*1.0/(Ngeo*V*Sig)

        cxrate = 1.0/Tau
    END FUNCTION CXKaiju

    !Charge exchange cross-section for K [keV] and species ispc
    !Sig in cm2
    !Using Lindsay & Stebbings 2005
    FUNCTION CXSigma(K,ispc) result(Sig)
        real(rprec), intent(in) :: K
        integer(iprec), intent(in) :: ispc
        real(rprec) :: Sig

        real(rprec) :: Sig0, KSig,a1,a2,a3,B1,B2

        Sig0 = 1.0e-16

        select case(ispc)
            case(RCMPROTON)
            !Charge exchange cross-section for H+/H
                !Cap for validity of CX cross-section
                KSig = K
                call ClampValue(KSig,0.005_rprec,250.0_rprec)

                a1 = 4.15
                a2 = 0.531
                a3 = 67.3

                B1 = (a1-a2*log(KSig))**2.0
                B2 = 1.0-exp(-a3/KSig) 
                Sig =  Sig0*B1*(B2**(4.5))
            case(RCMOXYGEN)
            !Charge exchange cross-section for O+/H
                !Cap for validity of CX cross-section
                KSig = K
                call ClampValue(KSig,0.025_rprec,600.0_rprec)
                a1 = 3.13
                a2 = 0.170
                a3 = 87.5

                B1 = (a1-a2*log(KSig))**2.0
                B2 = 1.0-exp(-a3/KSig) 
                Sig =  Sig0*B1*(B2**(0.8))
            case default
                Sig = 0.0
        end select

    END FUNCTION CXSigma

    !Geocoronal density afa L [#/cc], Taken from Ostgaard 2003 
    FUNCTION OstgaardGeocorona(L) result(Ngeo)
        real(rprec), intent(in) :: L
        real(rprec) :: Ngeo

        Ngeo = 10000.0*exp(-L/1.02) + 70.0*exp(-L/8.2)

    END FUNCTION OstgaardGeocorona

    !Simple mock-up for FLC losses
    FUNCTION FLCRat(ie,alam,vm,beq,rcurv,lossc) result(lossFLC)
        use constants, only : radius_earth_m
        use kdefs, only : TINY
        use math, only : RampUp
        IMPLICIT NONE
        integer(iprec), intent(in) :: ie
        real(rprec), intent(in) :: alam,vm,beq,rcurv,lossc
        real(rprec) :: lossFLC
        real(rprec) :: Np,bfp,ftv,K,V,TauSS,Rgyro,eps,xSS,earg

        bfp = beq/(sin(lossc)**2.0) !Foot point field strength, nT
        ftv = (1.0/vm)**(3.0/2.0) !flux-tube volume Re/nT
        K = alam*vm*1.0e-3 !Energy [keV]

        if (ie == RCMPROTON) then
            Np = 1 !Number of nucleons
        else
            lossFLC = 0.0
            return
        endif

        V = sqrt(2*K/Np)*sqrt(kev2J/(Mp_cgs*1.0e-3)) !V in m/s

        !Convert V from m/s to Re/s
        V = V/radius_earth_m

        TauSS = 3*2*ftv*bfp/V !Strong scattering lifetime [s], assuming ion w/ gamma=1

        Rgyro = (4.6e+3)*sqrt(K)/beq !Gyroradius of proton [km], assuming K in keV and beq in nT
        Rgyro = Rgyro/(radius_earth_m*1.0e-3) !In terms of Re
        eps = Rgyro/rcurv

        !Chen+ 2019, w/ correction
        if (eps>TINY) then
            !1/TauSS = Strong scattering loss rate
            earg = eps**5.0
            lossFLC = min(1.0,100.0*earg)*(1/TauSS) !Rate, 1/s
        else
            lossFLC = 0.0
        endif

        ! !Gilson criteria, use Chen Tau_SS
        ! !Ramp up to LossSS between kappa <= sqrt(8) and 1
        ! xC = 1.0/8.0
        ! if (eps >= xC) then
        !     xSS = LinRampUp(eps,xC,1.0-xC)
        !     lossFLC = xSS*(1.0/TauSS) !Rate, 1/s
        ! else
        !     lossFLC = 0.0 !None
        ! endif
        ! eps = Rgyro/rcurv

        ! !Chen+ 2019

        ! !K: Mockup between Chen/Gilson, transition between eps^-5 dep. and strong scattering at kappa = sqrt(8)
        ! !xSS = max( (8.0*eps)**(-5.0), 1.0 )
        ! earg = eps**(-5.0)
        ! xSS = max(100.0*earg,1.0)

        ! TauFLC = xSS*TauSS
        ! lossFLC = 1.0/TauFLC !Rate, 1/s

    END FUNCTION FLCRat

    FUNCTION RatefnC_tau_s(alam,vm,beq,lossc) result(TauSS)
        use constants, only : radius_earth_m
        use kdefs, only : TINY
        use math, only : RampUp
        IMPLICIT NONE
        real(rprec), intent(in) :: alam,vm,beq,lossc
        real(rprec) :: TauSS
        real(rprec) :: bfp,ftv,K,V
        bfp = beq/(sin(lossc)**2.0) !Foot point field strength, nT
        ftv = (1.0/vm)**(3.0/2.0) !flux-tube volume Re/nT
        K = alam*vm*1.0e-3 !Energy [keV]
        V = (3.1e+2)*sqrt(abs(K)) !km/s
        !Convert V from km/s to Re/s
        V = V/(radius_earth_m*1.0e-3)
        TauSS = 3.D0*2.D0*ftv*bfp/V !Strong scattering lifetime [s], assuming gamma=1 and eta=2/3. Need to include relativisitc factor.
!        write(*,"(a,1x,5e25.10)") "In RatefnC_tau_s",bfp,ftv,K,V,TauSS
    END FUNCTION RatefnC_tau_s

    FUNCTION RatefnC_lambda(mltx,engx,Lshx) result(lambda)
    ! default scattering rate lambda based on Chen et al. 2005.
    ! tau = 1/(1+a1*sin(phi+phi0)+a2*cos(2*(phi+phi0)))/lambda0, 
    ! where lambda0 = min(0.08*E[MeV]^(-1.32), 0.4*10^(2L-6+0.4*log_2(E)))/day, a1=1.2, a2=-0.25*a1, phi0=pi/6.
    ! lambda max = 2.6*lambda0 at 04MLT, and min = 0.6 lambda0 at 22 MLT.
        IMPLICIT NONE
        REAL (rprec), INTENT (IN) :: mltx,engx,Lshx
        REAL (rprec) :: lambda, lambda0, a1, a2, phi0, phi, tau

        lambda0 = min(0.08*engx**(-1.32), 0.4*10**(2*Lshx-6.D0+0.4*dlog(engx)/dlog(2.D0)))/86400.D0 ! 1/s.
        a1 = 1.2D0
        a2 = -0.25D0*a1
        phi0 = pi/6.D0
        phi = (mltx-0.D0)/12.D0*pi+phi0
        lambda = (1.D0+a1*sin(phi)+a2*cos(2.D0*phi))*lambda0
    END FUNCTION RatefnC_lambda
    
    FUNCTION RatefnC_lambda_w(mltx,engx,Lshx,kpx) result(lambda)
    ! empirical scattering rate lambda based on Orlova and Shprits, 2014.
        IMPLICIT NONE
        REAL (rprec), INTENT (IN) :: mltx,engx,kpx,Lshx ! engx in MeV.
        REAL (rprec) :: lambda, R, K, E, RK, KE, RE, RKE, R2, K2, E2
        REAL (rprec), DIMENSION(33) :: a1_33, rke_pol

        lambda = 0.0
        if(mltx<=21.0 .and. mltx>15.0) then
            return
        endif
        a1_33 = 0.D0
        if((mltx>21.0 .or. mltx<=3.0)) then ! Night side
            if(engx<=10e-3) then ! Night side Ek <= 10 keV
                a1_33 = [ & 
               -6.29800e+00,  2.47510e+00, -7.26390e-01,  8.69730e+01, -1.57180e+03, -1.53970e-02, -1.25300e+00,  1.16620e-03, -3.19940e+02, -5.98660e+00,  & 
                0.00000e+00, -3.27510e-01,  9.05410e-02, -4.77420e+00,  1.56930e-03,  1.18960e+01,  1.07320e+03,  1.50100e-02, -3.96060e-03,  0.00000e+00,  & 
                2.11370e+00, -4.24700e+02,  2.64410e+04, -9.88010e+05, -8.82250e-02,  8.66970e+00,  4.19100e+02,  4.70090e-03, -7.94090e-01,  2.02540e+03,  & 
               -1.47880e+05,  5.62920e+06,  0.00000e+00]
            else if(engx>10e-3 .and. engx<0.5) then ! Night side 10 keV < Ek < 0.5 MeV
                if(kpx<=3.0) then ! Night side 10 keV < Ek < 0.5 MeV and Kp<=3 
                    a1_33 = [ & 
                    9.24560e+00, -3.76850e+00,  1.25590e+00,  4.02960e-02,  1.18330e+00, -1.08430e-01,  0.00000e+00,  1.36980e-02,  1.60750e+01, -1.98660e+01,  & 
                    8.13400e+00,  5.44060e-01, -1.66160e-01, -9.18410e-02,  0.00000e+00, -1.69710e+00,  9.91020e-01, -2.81450e-02,  9.26060e-03,  5.88920e-02,  & 
                   -3.56740e+00, -9.05730e+00, -1.76860e+01,  1.00290e+01,  6.32200e-01,  7.36800e+00,  0.00000e+00, -8.55330e-02, -1.00050e+00, -3.13980e+01,  & 
                    8.88950e+01, -5.97370e+01,  0.00000e+00]
                else !% Night side 10 keV < Ek < 0.5 MeV and Kp>3 
                    a1_33 = [ & 
                   -3.80640e+00, -1.06160e+00, -5.02820e-01, -9.34780e-01,  0.00000e+00,  1.59690e-01, -8.14210e-02, -8.10140e-03,  2.72620e+01, -1.05540e+01,  & 
                    2.35650e+00,  4.42100e-01, -5.90540e-02,  1.55420e-01, -2.73460e-03, -4.05000e+00,  9.05600e-01, -2.84240e-02,  4.96020e-03,  1.48180e-01,  & 
                    4.10080e+00, -5.58120e+01,  1.03920e+01,  0.00000e+00, -8.51250e-01,  9.90870e+00, -8.34090e-01,  4.40160e-02, -5.34040e-01,  8.23600e+01,  & 
                   -7.64720e+01,  1.58860e+02, -1.48750e+02]
                endif
            else !% Night side Ek >= 0.5 MeV
                if(kpx<=3.0) then !% Night side Ek >= 0.5 MeV and Kp<=3 
                    a1_33 = [ & 
                    9.49410e+00, -4.60800e+00,  1.98140e+00, -3.54580e-01,  1.07230e-01,  1.54310e-01, -1.54080e-02, -1.13010e-02,  3.12460e+00, -6.73170e-02,  & 
                   -8.89080e-02,  9.11370e-01, -3.74740e-01,  6.32320e-03, -6.35990e-03, -3.68530e-01,  0.00000e+00, -6.11110e-02,  2.10430e-02,  1.80110e-02,  & 
                   -7.95890e+00, -2.00000e+00,  9.35410e-01, -4.25710e-01,  2.03700e+00,  1.06340e+00,  0.00000e+00, -3.05470e-01, -1.32250e-01,  9.44100e+00,  & 
                   -1.22310e+01,  5.00080e+00, -5.05700e-01]
                else !% Night side Ek >= 0.5 MeV and Kp>3 
                    a1_33 = [ & 
                   -2.26310e+01,  1.16930e+01, -4.69090e+00, -2.47010e-01,  3.34580e-02,  9.81050e-01,  0.00000e+00, -5.58680e-02,  3.83790e+00, -8.52390e-01,  & 
                    0.00000e+00, -9.06580e-02, -1.58690e-01,  9.00520e-03, -7.96750e-03, -3.59650e-01,  6.80280e-02, -3.44430e-02,  1.87950e-02,  5.11320e-03,  & 
                    1.06470e+01, -1.13320e+01,  1.39530e-01, -6.84300e-02, -2.16830e+00,  2.09080e+00,  0.00000e+00,  1.28550e-01, -1.21910e-01,  1.99560e+01,  & 
                    3.90780e-01, -2.71180e-01,  1.92860e-01]
                endif
            endif
        else if((mltx>3.0 .and. mltx<=9.0)) then !% Dawn side
            if(engx<7e-3) then !% Dawn side Ek < 7 keV 
                a1_33 = [ & 
                1.30230e+01, -5.13970e+00,  3.66290e-01, -1.14330e+02,  5.23120e+03, -4.32600e-03,  1.79870e+00,  0.00000e+00,  2.36470e+03, -2.37100e+05,  & 
                8.16800e+06,  6.61540e-01, -3.53860e-02,  4.47500e+00,  0.00000e+00, -2.23630e+02,  9.03540e+03, -2.98050e-02,  1.35540e-03,  8.04240e+00,  & 
               -1.78800e+00,  5.68060e+02, -6.35360e+04,  2.30770e+06,  7.67790e-02, -1.68380e+01,  1.33480e+03,  0.00000e+00,  0.00000e+00, -6.51880e+03,  & 
                9.23510e+05, -4.27880e+07,  0.00000e+00]
            else if(engx>=7e-3 .and. engx<90e-3) then !% Dawn side 7 keV <= Ek < 90 keV 
                a1_33 = [ & 
                3.08720e+00,  4.40920e-01, -3.48240e-01,  2.83500e+00, -1.08380e+01,  1.37250e-02, -6.88500e-02,  0.00000e+00,  2.32830e+01, -4.26740e+02,  & 
                1.21980e+03, -1.41190e-01,  3.97480e-02, -1.01110e-01, -7.46140e-04, -1.04010e+00,  2.47740e+01,  1.07430e-02, -1.61820e-03, -9.16570e-02,  & 
                1.56840e-02, -1.22210e+01,  2.26340e+01,  3.96600e+02,  3.73680e-03,  7.93640e-01, -3.80690e+00,  1.16140e-03,  0.00000e+00, -1.70600e+02,  & 
                4.94760e+03, -4.94390e+04,  1.81330e+05]
            else !% Dawn side Ek >= 90 keV 
                a1_33 = [ & 
                4.01120e+00, -3.90940e-01, -3.25600e-02, -3.79590e-02, -6.26790e-03,  5.22710e-03,  1.26300e-03, -1.76010e-04,  1.00320e-01,  7.42100e-02,  & 
                1.13470e-02,  4.03350e-02,  5.57050e-03,  3.34900e-03, -3.85910e-04, -2.60900e-02, -6.55570e-03, -2.42570e-03, -2.71600e-04,  1.59150e-03,  & 
               -9.36750e-01,  5.75970e-02,  5.61290e-02, -8.76050e-03,  5.90570e-02, -6.86030e-03,  1.47720e-03,  1.07560e-03, -3.19070e-04,  3.41930e+00,  & 
               -2.95490e+00,  1.25700e+00, -2.41230e-01]
            endif
        else if((mltx>9.0 .and. mltx<=12.0)) then !% Prenoon side
            if(engx<7e-3) then !% Prenoon side Ek < 7 keV 
                a1_33 = [ & 
                9.37910e+00, -3.38010e+00,  2.29500e-01, -2.10440e+01,  3.23730e+03, -1.49830e-02, -1.37450e+00,  0.00000e+00,  1.95630e+03, -2.60190e+05,  & 
                1.25510e+07,  4.05490e-01, -1.64040e-02,  0.00000e+00,  1.40030e-03, -1.80880e+02,  6.90820e+03, -1.82190e-02,  0.00000e+00,  7.95340e+00,  & 
               -1.44820e+00,  2.80490e+02, -3.29820e+04,  0.00000e+00,  1.09630e-01, -1.94360e+01,  2.54150e+03,  0.00000e+00,  0.00000e+00, -4.98720e+03,  & 
                8.43540e+05, -4.44550e+07,  0.00000e+00]
            else if(engx>=7e-3 .and. engx<0.1) then !% Prenoon side 7 keV <= Ek < 0.1 MeV 
                a1_33 = [ & 
                6.88650e+00, -1.66320e+00, -1.87690e-02, -4.04800e-01, -4.01940e+00, -4.79730e-03,  0.00000e+00,  4.97590e-04,  3.74930e+01, -1.84030e+02,  & 
                0.00000e+00,  1.48470e-01,  9.90130e-03,  8.86630e-02,  0.00000e+00, -3.66650e+00,  1.61190e+01, -2.52420e-03, -9.83090e-04,  4.54560e-02,  & 
               -3.00720e-01,  3.96140e+00, -4.95480e+01,  3.76000e+02,  1.89350e-02, -7.31370e-02,  0.00000e+00, -1.26880e-03,  0.00000e+00, -2.17920e+02,  & 
                3.72960e+03, -3.41250e+04,  1.27070e+05]
            else !% Prenoon side Ek >= 0.1 MeV 
                a1_33 = [ & 
                3.77720e-01,  3.19090e-01, -9.72410e-02,  5.28290e-02, -5.49410e-03,  1.25200e-02, -4.55740e-03, -1.44610e-04, -1.16590e-01,  9.03070e-02,  & 
                4.56800e-03, -5.65820e-02,  7.94530e-03, -1.95700e-03, -9.02940e-04, -1.63340e-02, -4.90140e-03,  1.45370e-03,  0.00000e+00,  3.07330e-03,  & 
                5.18830e-02, -2.32240e-01,  8.30020e-02, -1.89770e-02, -3.38290e-02,  1.88780e-02, -1.07900e-03,  4.20820e-04,  0.00000e+00,  4.70710e+00,  & 
               -4.06830e+00,  1.84240e+00, -3.35050e-01]
            endif
        else if((mltx>12.0 .and. mltx<=15.0)) then !% Postnoon side
            if(engx<6e-3) then !% Postnoon side Ek < 6 keV 
                a1_33 = [ & 
                5.49270e+00, -1.22390e+00,  1.67600e-01, -3.24120e+01,  4.52400e+03, -1.39930e-02, -1.13610e+00,  0.00000e+00,  1.18180e+03, -2.11750e+05,  & 
                1.39230e+07,  7.95320e-02, -7.59730e-03,  0.00000e+00,  9.95900e-04, -8.14630e+01,  1.78080e+03, -3.80660e-03,  0.00000e+00,  4.75400e+00,  & 
               -1.40150e+00,  4.78460e+02, -8.44970e+04,  3.56970e+06,  1.40880e-01, -2.63470e+01,  3.59790e+03, -2.24530e-03,  0.00000e+00, -3.68170e+03,  & 
                8.58490e+05, -6.09680e+07,  0.00000e+00]
            else if(engx>=6e-3 .and. engx<0.1) then !% Postnoon side 6 keV <= Ek < 0.1 MeV 
                a1_33 = [ & 
                5.25830e+00, -9.93980e-01, -1.03430e-01,  2.15500e-01, -5.89310e+00, -6.02980e-03,  1.83200e-02,  5.13720e-04,  3.45720e+01, -1.67450e+02,  & 
               -2.08860e+02,  6.29100e-02,  2.26720e-02,  4.07410e-02,  0.00000e+00, -3.92410e+00,  1.98090e+01,  1.48700e-03, -1.54070e-03,  5.62200e-02,  & 
               -5.01280e-02, -1.31250e-01, -1.27710e+01,  2.93380e+02,  2.57500e-02, -1.60430e-02, -1.25290e+00, -1.76190e-03,  0.00000e+00, -1.74200e+02,  & 
                3.26740e+03, -3.24120e+04,  1.30380e+05]
            else !% Postnoon side Ek >= 0.1 MeV 
                a1_33 = [ & 
                8.56720e-01,  4.13830e-01, -1.14870e-01,  5.24810e-02, -5.71760e-03,  1.33180e-02, -5.20970e-03, -1.26850e-04, -6.66380e-02,  1.57840e-01,  & 
               -4.22140e-03, -7.55980e-02,  1.08770e-02, -1.55610e-03, -9.77010e-04, -3.50350e-02, -8.20200e-03,  2.57510e-03, -1.58690e-04,  4.45860e-03,  & 
                1.40840e-01, -2.57110e-01,  9.20020e-02, -1.76770e-02, -3.68530e-02,  2.54710e-02, -3.38180e-03,  3.24390e-04,  0.00000e+00,  4.83550e+00,  & 
               -4.46000e+00,  1.98870e+00, -3.56050e-01]
            endif
        endif
        R  = Lshx
        K  = kpx+1.0
        E  = engx ! [MeV]
        RK = R*K
        KE = K*E
        RE = R*E
        RKE= RK*E
        R2 = R*R
        K2 = K*K
        E2 = E*E
        rke_pol = (/1.D0,R,RK,RKE,RKE*E,RK*K,RK*KE,RK*K2,RE,RE*E,RE*E2,R2,R*RK,R*RKE,RK*RK,R*RE,RE*RE,R*R2,R2*RK,R2*RE,K,KE,KE*E,KE*E2,K2,K*KE,KE*KE,K*K2,K2*KE,E,E2,E*E2,E2*E2/)
        ! log(tau) = a1+a2*R+...
        lambda = 10.0**(-dot_product(a1_33,rke_pol))
    END FUNCTION RatefnC_lambda_w

    FUNCTION RatefnC_lambda_h(mltx,engx,Lshx,kpx) result(lambda)
    ! empirical scattering rate lambda based on Orlova and Shprits, 2014.
    ! Inside plasmasphere, lambda_h=Dhaa, Kp and MLT parameterized PA diffusion coef 
    ! against Pspheric hiss for energies of 1 keV to 10 MeV from R0 of 3 to 6.
        IMPLICIT NONE
        REAL (rprec), INTENT (IN) :: mltx,engx,kpx,Lshx ! engx in MeV.
        REAL (rprec) :: lambda, tau, gEL, yKp, fL
        REAL (rprec) :: L, E, K, L2, L3, L4, E2, E3, E4, E5, K2, LE, LK, EK, LEK
        REAL (rprec), DIMENSION(28) :: a1_28, le_pol
        REAL (rprec), DIMENSION(59) :: a1_59, lek_pol

        lambda = 0.0
        L = Lshx ! L=3-6
        E = engx
        K = kpx
        L2 = L*L
        L3 = L2*L
        L4 = L3*L
        fL = -0.2573*L4 + 4.2781*L3 - 25.9348*L2 + 66.8113*L - 66.1182
        if(L>6.0 .or. L<3.0 .or. K>6.0 .or. E>10.0 .or. E<1.0e-3 .or. E<fL) then
            return
        endif
        E2 = E*E
        E3 = E2*E
        E4 = E3*E
        E5 = E4*E
        K2 = K*K
        LE = L*E
        LK = L*K
        EK = E*K
        LEK= LE*K
        if(mltx<=21 .and. mltx>=6) then ! Dayside, Ek = 1 keV - 10 MeV and > f(L), Kp<=6, Coefficients for g(E,L) parameterization on the dayside (Eq. 7 in the main text).
            if(E>fL .and. kpx<=6 .and. E<=10.0 .and. E>=1.0e-3) then
                a1_28 = [ &
                4.32340e+01, -3.70570e+01, -1.09340e+02,  1.16160e+01,  8.29830e+01,  3.57260e+01, -1.58030e+00, -2.26630e+01, -1.60060e+01,  3.51250e+01,  & 
                7.99240e-02,  2.67740e+00,  2.13960e+00, -2.41070e+01, -3.50920e+00, -1.15310e-01, -7.43760e-02,  5.16590e+00, -2.04880e+00, -4.52050e+00,  & 
               -3.52900e-01,  1.25340e+00,  2.32240e+00,  7.00390e-01, -1.36660e-01, -2.71340e-01, -1.41940e-01, -2.09270e-02 ]
                le_pol = (/1.D0,L,E,L2,LE,E2,L3,L2*E,L*E2,E3,L4,L3*E,L2*E2,L*E3,E4,L3*LE,L3*E2,L2*E3,L*E4,E5,L3*E3,L2*E4,L*E5,E3*E3,L3*E4,L2*E5,LE*E5,E2*E5/)
                gEL = dot_product(a1_28,le_pol) ! only valid for E>f(L)
                yKp = 0.015465*kpx*kpx - 0.26074*kpx + 1.0077
                tau = 10**(gEL+yKp)
            endif
        else ! Nightside, 21<MLT<6, Kp<=6, f(L) < Ek < 10 MeV, Eq. S1 in the supporting material.
            a1_59 = [ &
            5.35120e+01, -4.62820e+01, -1.30630e+02, -1.69970e+00,  1.51270e+01,  1.08050e+02,  5.42210e+01,  1.22940e+00, -1.66680e+00, -4.14470e-02,  & 
           -2.16170e+00, -3.27550e+01, -2.86350e+01,  3.67450e+00, -2.36980e-01,  1.75970e+00,  4.50400e+00,  1.14840e-01,  4.38180e+00,  4.99730e+00,  & 
           -4.06960e+00, -1.31850e+01,  1.49090e-02, -4.94890e-01, -3.07880e+00,  6.97390e-01, -2.18870e-01, -2.88410e-01,  9.19740e-01,  4.05090e+00,  & 
           -3.50100e+00,  4.15720e-02,  6.55190e-01, -1.16650e+00, -6.87900e-02, -5.47790e-02, -1.45610e-01,  1.89180e+00,  2.82220e+00, -4.41970e-02,  & 
            3.92090e-01,  1.03970e-01,  4.47600e-01, -2.25270e-02, -2.16760e-01, -5.81220e-01,  1.02370e+00, -3.64160e-02, -1.80690e-02, -1.86230e-01,  & 
            9.45330e-04, -2.49520e-01,  1.06160e-01,  1.82980e-02,  1.79470e-03,  1.89430e-02, -3.98110e-02, -1.36410e-02, -2.73280e-03 ]
            lek_pol = (/1.D0,L,E,K,L2,LE,E2,LK,EK,K2,L3,L2*E,L*E2,E3,L2*K,LE*K,E2*K,L4,L3*E,L2*E2,L*E3,E4,L3*K,L2*EK,LE*EK,E3*K,L4*E,L3*E2,L2*E3,L*E4,E5,L3*EK,L2*E2*K,LK*E3,E4*K &
                ,L3*E3,L2*E4,L*E5,E3*E3,L3*E2*K,L2*E3*K,LK*E4,E5*K,L3*E4,L2*E5,LE*E5,E3*E4,L3*E3*K,L2*E4*K,LEK*E4,E5*EK,L*E3*E4,E4*E4,L2*E5*K,LEK*E5,E3*E3*EK,LE*E3*E4,E4*E5,LEK*E3*E3/)
            tau = 10**(dot_product(a1_59,lek_pol))
        endif
        lambda = 1.D0/tau
    END FUNCTION RatefnC_lambda_h
END MODULE lossutils
