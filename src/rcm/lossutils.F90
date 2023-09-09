!Utilities for loss calculations

MODULE lossutils
    USE kdefs, ONLY : TINY,PI,Mp_cgs,kev2J,Me_cgs,vc_cgs
    USE rcm_precision
    USE rcmdefs
    USE math, ONLY : SmoothOpTSC,SmoothOperator33,ClampValue,LinRampUp
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

    END FUNCTION FLCRat

    FUNCTION RatefnC_tau_s(alam,vm,beq,lossc) result(TauSS)
    ! Strong diffusion lifetime based on Schulz, 1974, 1998.
    ! tau_s ~ [2*Psi*Bh/(1-eta)](gamma*m/p), where Psi is flux tube volume, Bh is |B| at field line foot point.
    ! eta is backscatter rate at alitude h, here eta=2/3.
    ! gamma = m/m0 is relativisitc factor, p is particle momentum.
    !       = mc2/m0c2 = (m0c2+K)/m0c2 = 1+K/mec2 ! mec2=0.511 is me*c^2 in MeV
    ! m = m0/sqrt(1-v^2/c^2)
    ! V = c*1/sqrt(1-1/gammar2)
        use kdefs, only : TINY,vc_cgs,Re_cgs,mec2
        use math, only : RampUp
        IMPLICIT NONE
        real(rprec), intent(in) :: alam,vm,beq,lossc
        real(rprec) :: TauSS
        real(rprec) :: bfp,ftv,K,V,gammar
        bfp = beq/(sin(lossc)**2.0)               ! Foot point field strength, nT
        ftv = (1.0/vm)**(3.0/2.0)                 ! flux-tube volume Re/nT
        K = abs(alam)*vm*1.0e-6                   ! Energy [MeV]
        gammar = 1.0+K/mec2
        V = vc_cgs*sqrt(1.0-1.0/gammar**2)/Re_cgs ! Re/s
        TauSS = 3.D0*2.D0*ftv*bfp/V*gammar        ! Strong scattering lifetime [s], assuming eta=2/3.

    END FUNCTION RatefnC_tau_s

    FUNCTION RatefnDW_tau_c(Kpx,mltx,Lx,Ekx) result(tau)
    ! linearly interpolate tau from EWMTauInput to current MLT,L,Kp,Ek value
        USE rice_housekeeping_module, ONLY: EWMTauInput
        IMPLICIT NONE
        REAL (rprec), INTENT (IN) :: Kpx, mltx,Lx,Ekx
        REAL(rprec) :: tau
        REAL(rprec) :: tauKMLE(2,2,2,2),tauMLE(2,2,2),tauLE(2,2),tauE(2)! tauKMLE(1,2,2,2) means tauKlMuLuEu, l:lower bound, u: upper bound in the NN methond
        REAL(rprec) :: dK,wK,dM,wM,dL,wL,dE,wE
        INTEGER :: iK,kL,kU,mL,mU,lL,lU,eL,eU
        LOGICAL :: Lflag = .false.


        associate(Nm=>EWMTauInput%ChorusTauInput%Nm,Nl=>EWMTauInput%ChorusTauInput%Nl,Nk=>EWMTauInput%ChorusTauInput%Nk,Ne=>EWMTauInput%ChorusTauInput%Ne,&
                  Kpi=>EWMTauInput%ChorusTauInput%Kpi,MLTi=>EWMTauInput%ChorusTauInput%MLTi,Li=>EWMTauInput%ChorusTauInput%Li,Eki=>EWMTauInput%ChorusTauInput%Eki,&
                  taui=>EWMTauInput%ChorusTauInput%taui)

        ! Find the nearest neighbors in Kp
        if (Kpx >= maxval(Kpi)) then
            kL = Nk !use Kp maximum 
            kU = Nk
        else if (Kpx <= minval(Kpi)) then
            kL = 1  !use Kp minimum
            kU = 1
        else
            kL = maxloc(Kpi,dim=1,mask=(Kpi<Kpx))
            kU = kL+1
        endif
        
        ! Find the nearest neighbours in MLT
        if ((mltx >= maxval(MLTi)) .or. (mltx <= minval(MLTi)))  then ! maxval of MLT is 24, minval of MLT is 0
            mL = 1 !use MLT = 0
            mU = 1
        else
            mL = maxloc(MLTi,dim=1,mask=(MLTi<mltx))
            mU = mL+1
        endif

        ! Find the nearest neighbours in L
        if (Lx >= maxval(Li)) then
            lL = Nl !use L maximum
            lU = Nl
            Lflag = .true.
        else if (Lx <= minval(Li)) then
            lL = 1 ! use L minimum
            lU = 1
        else
            lL = maxloc(Li,dim=1,mask=(Li<Lx))
            lU = lL+1
        endif

         ! Find the nearest neighbours in Ek
        if (Ekx < minval(Eki)) then
            tau = 1.D10 ! For low energies, assign a huge lifetime is 10^10s ~ 10^3 years.
        else if (Ekx >= maxval(Eki)) then
            eL = Ne !use Ek maximum
            eU = Ne
        else
            eL = maxloc(Eki,dim=1,mask=(Eki<Ekx))
            eU = eL + 1
        endif

        !linear interpolation in Kp
        if (kL == kU) then
            tauMLE(1,1,1) = log10(taui(kL,mL,lL,eL))!Interpolation in log10(taui) space
            tauMLE(1,1,2) = log10(taui(kL,mL,lL,eU))
            tauMLE(1,2,1) = log10(taui(kL,mL,lU,eL))
            tauMLE(1,2,2) = log10(taui(kL,mL,lU,eU))
            tauMLE(2,1,1) = log10(taui(kL,mU,lL,eL))
            tauMLE(2,1,2) = log10(taui(kL,mU,lL,eU))
            tauMLE(2,2,1) = log10(taui(kL,mU,lU,eL))
            tauMLE(2,2,2) = log10(taui(kL,mU,lU,eU))
        else
            dK = Kpi(kU)-Kpi(kL)
            wK = (Kpx-Kpi(kL))/dK
            tauKMLE(1,1,1,1) = log10(taui(kL,mL,lL,eL))
            tauKMLE(2,1,1,1) = log10(taui(kU,mL,lL,eL))
            tauMLE(1,1,1) = tauKMLE(1,1,1,1) + wK*(tauKMLE(2,1,1,1)-tauKMLE(1,1,1,1))
            tauKMLE(1,1,1,2) = log10(taui(kL,mL,lL,eU))
            tauKMLE(2,1,1,2) = log10(taui(kU,mL,lL,eU))
            tauMLE(1,1,2) = tauKMLE(1,1,1,2) + wK*(tauKMLE(2,1,1,2)-tauKMLE(1,1,1,2))
            tauKMLE(1,1,2,1) = log10(taui(kL,mL,lU,eL))
            tauKMLE(2,1,2,1) = log10(taui(kU,mL,lU,eL))
            tauMLE(1,2,1) = tauKMLE(1,1,2,1) + wK*(tauKMLE(2,1,2,1)-tauKMLE(1,1,2,1))
            tauKMLE(1,1,2,2) = log10(taui(kL,mL,lU,eU))
            tauKMLE(2,1,2,2) = log10(taui(kU,mL,lU,eU))
            tauMLE(1,2,2) = tauKMLE(1,1,2,2) + wK*(tauKMLE(2,1,2,2)-tauKMLE(1,1,2,2))
            tauKMLE(1,2,1,1) = log10(taui(kL,mU,lL,eL))
            tauKMLE(2,2,1,1) = log10(taui(kU,mU,lL,eL))
            tauMLE(2,1,1) = tauKMLE(1,2,1,1) + wK*(tauKMLE(2,2,1,1)-tauKMLE(1,2,1,1))
            tauKMLE(1,2,1,2) = log10(taui(kL,mU,lL,eU))
            tauKMLE(2,2,1,2) = log10(taui(kU,mU,lL,eU))
            tauMLE(2,1,2) = tauKMLE(1,2,1,2) + wK*(tauKMLE(2,2,1,2)-tauKMLE(1,2,1,2))
            tauKMLE(1,2,2,1) = log10(taui(kL,mU,lU,eL))
            tauKMLE(2,2,2,1) = log10(taui(kU,mU,lU,eL))
            tauMLE(2,2,1) = tauKMLE(1,2,2,1) + wK*(tauKMLE(2,2,2,1)-tauKMLE(1,2,2,1))
            tauKMLE(1,2,2,2) = log10(taui(kL,mU,lU,eU))
            tauKMLE(2,2,2,2) = log10(taui(kU,mU,lU,eU))
            tauMLE(2,2,2) = tauKMLE(1,2,2,2) + wK*(tauKMLE(2,2,2,2)-tauKMLE(1,2,2,2))
        end if

        ! linear interpolation in mlt
        if (mL == mU) then
            tauLE(1,1) = tauMLE(2,1,1)
            tauLE(1,2) = tauMLE(2,1,2)
            tauLE(2,1) = tauMLE(2,2,1)
            tauLE(2,2) = tauMLE(2,2,2)
        else
            dM = MLTi(mU)-MLTi(mL)
            wM = (mltx-MLTi(mL))/dM
            tauLE(1,1) = tauMLE(1,1,1) + wM*(tauMLE(2,1,1)-tauMLE(1,1,1))
            tauLE(1,2) = tauMLE(1,1,2) + wM*(tauMLE(2,1,2)-tauMLE(1,1,2))
            tauLE(2,1) = tauMLE(1,2,1) + wM*(tauMLE(2,2,1)-tauMLE(1,2,1))
            tauLE(2,2) = tauMLE(1,2,2) + wM*(tauMLE(2,2,2)-tauMLE(1,2,2))
        end if

        ! linear interpolation in L
        if (lL == lU) then
            tauE(1) = tauLE(2,1)
            tauE(2) = tauLE(2,2)
            if (Lflag) then ! use gaussian decay for L > maxval(Li) (7Re)
               tauE(1)=tauE(1)/exp(-(Lx-maxval(Li))**2)
               tauE(2)=tauE(2)/exp(-(Lx-maxval(Li))**2)
            endif
        else
            dL = Li(lU)-Li(lL)
            wL = (Lx-Li(lL))/dL
            tauE(1) = tauLE(1,1)+ wL*(tauLE(2,1)-tauLE(1,1))
            tauE(2) = tauLE(1,2)+ wL*(tauLE(2,2)-tauLE(1,2))
        end if

        ! linear interpolation in Ek
        if (eL == eU) then
            tau = tauE(1)
        else
            dE = Eki(eU)-Eki(eL)
            wE = (Ekx-Eki(eL))/dE
            tau = tauE(1) + wE*(tauE(2)-tauE(1))
        end if

        tau = 10.0**tau !convert back to tau in seconds

        end associate

    END FUNCTION RatefnDW_tau_c


    FUNCTION RatefnC_tau_h16(mltx,engx,Lshx,kpx) result(tau)
    ! Empirical lifetime against plasmaspheric hiss pitch angle diffusion, based on Orlova et al. 2015JA021878.
    ! Improvements relative to 2014GL060100: 1. Hiss wave intensity distribution model is based on new data 
    ! (O14 was based on single-component E field in CRRES data. O16 used Spasojevic+2015 model based on EMFISIS B data on VAP); 
    ! 2. Wave spectrum is assumed differently (O14 assume Gaussian spectrum based on CRRES data).
    ! Electron lifetime tau(L,E,MLT,Kp) = tau_av(L,E)/g(MLT)/h(Kp), 
    ! where 1.5<L<5.5, E=log10(Ek[MeV]) for 1 KeV < Ek < 10 MeV.
    ! log10(tau_av(L,E)) = a1+a2*L+a3*E+...+a20*E^3, when E >= f(L).
    ! f(L) = 0.1328*L^2-2.1463*L+3.7857.
    ! g(MLT) = 10^g0(MLT)/G0
    ! h(Kp) = 10^h0(Kp)/H0
    !   G0 = int_0^24(10^g0(MLT))dMLT / 24 = 782.3.
    !   g0(MLT) = b2*MLT^2 + b1*MLT + b0
    !   H0 = 1315.
    !   h0(Kp) = c2*Kp^2 + c1*Kp + c0

        IMPLICIT NONE
        REAL (rprec), INTENT (IN) :: mltx,engx,kpx,Lshx ! engx in MeV.
        REAL (rprec) :: lambda, tau, tau_av
        REAL (rprec) :: MLT, L, E, K, L2, L3, L4, fL, E2, E3, E4, E5, LE
        REAL (rprec) :: b0, b1, b2, G0, g0_MLT, g_MLT, c0, c1, c2, H0, h0_Kp, h_Kp
        REAL (rprec), DIMENSION(20) :: le_pol
        REAL(rprec), dimension(20), parameter :: a1_20 = [77.323, -92.641, -55.754, 44.497, 48.981, 8.9067, -10.704, &
                                                         -15.711, -3.3326, 1.5189, 1.294, 2.2546, 0.31889, -0.85916, & 
                                                         -0.22182, 0.034318, 0.097248, -0.12192, -0.062765, 0.0063218]
     
        lambda = 0.D0
        tau = 1.D10
        MLT = mltx
        L = Lshx ! L=3-6
        E = log10(engx) ! engx is Ek in MeV
        L2 = L*L
        fL = 0.1328*L2 - 2.1463*L + 3.7857
        !if(L>5.5 .or. L<1.5 .or. E>1.0 .or. E<-3.0 .or. E<fL) then 
        !! Both sectors are only valid for log10(Ek)>=f(L), 1keV<Ek<10MeV, 1.5<=L<=5.5.
        !    return
        !endif

        call ClampValue(L,1.5_rprec,5.5_rprec)
      
        L2 = L*L !renew L2
        
        call ClampValue(E,max(-3.0_rprec,fL),1.0_rprec) 

        b0 = 2.080
        b1 = 0.1773
        b2 = -0.007338
        G0 = 782.3
        g0_MLT = b2*MLT*MLT + b1*MLT + b0
        g_MLT = 10**g0_MLT/G0
        c0 = 2.598
        c1 = 0.2321
        c2 = -0.01414
        H0 = 1315.0
        K = min(kpx,5.0) ! 0<Kp<5, 1.2% data in Kp bin of 4.3-7.6
        h0_Kp = c2*K*K + c1*K + c0
        h_Kp = 10**h0_Kp/H0

        E2 = E*E
        E3 = E2*E
        E4 = E3*E
        E5 = E4*E
        L3 = L2*L
        L4 = L3*L
        LE = L*E
        le_pol = (/1.D0,L,E,L2,LE,E2,L3,L2*E,L*E2,E3,L4,L3*E,L2*E2,L*E3,E4,L*E4,L2*E3,L4*E,L2*L3,E5/)
        tau_av = 10.0**(dot_product(a1_20,le_pol))*86400.D0 ! seconds
        tau = tau_av/g_MLT/h_KP
        lambda = 1.D0/tau ! 1/s
    END FUNCTION RatefnC_tau_h16
END MODULE lossutils
