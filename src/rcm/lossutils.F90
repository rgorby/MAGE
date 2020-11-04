!Utilities for loss calculations

MODULE lossutils
  USE kdefs, ONLY : TINY
  USE rcm_precision
  
  implicit none

  contains

    !Loss functions
    FUNCTION CXKaiju(isp,enrg,rloc) result(cxrate)
      IMPLICIT NONE

      integer(iprec), intent(in) :: isp
      real(rprec), intent(in) :: enrg,rloc

      real(rprec) :: cxrate
      real(rprec) :: K,L,Ngeo,KSig,Sig0,a1,a2,a3,B1,B2,Sig,M,Kj,V,Tau,tScl

      K = enrg*1.0e-3 !Energy in kev
    !Geocoronal density afa L [#/cc], Taken from Ostgaard 2003 
      L = rloc
      Ngeo = 10000.0*exp(-L/1.02) + 70.0*exp(-L/8.2)

    !Charge exchange cross-section for H+/H
      !K in keV, Sig in cm2
      !Using Lindsay & Stebbings 2005
      KSig = min(K,250.0) !Cap for validity of CX cross-section
      
      Sig0 = 1.0e-16
      a1 = 4.15
      a2 = 0.531
      a3 = 67.3

      B1 = (a1-a2*log(KSig))**2.0
      B2 = 1.0-exp(-a3/KSig) 
      Sig =  Sig0*B1*(B2**(4.5))
    !Get velocity [cm/s] from energy [keV]
      M = 1.67*1.0e-27 !Proton mass
      Kj = K*1000.0*1.6*1.0e-19 !Joules
      V = sqrt(2*Kj/M)*100.0 !m/s->cm/s

    !Timescale
      tScl = cos(45*PI/180.0)**3.5 !Using Smith & Bewtra 1976 scaling
      Tau = tScl*1.0/(Ngeo*V*Sig)

      cxrate = 1.0/Tau
    END FUNCTION CXKaiju


    !Really quick test of simple FLCRat
    FUNCTION FLCRat(ie,alam,vm,beq,rcurv,lossc) result(lossFLC)
      use constants, only : radius_earth_m
      use kdefs, only : TINY
      use math, only : RampUp
      IMPLICIT NONE
      integer(iprec), intent(in) :: ie
      real(rprec), intent(in) :: alam,vm,beq,rcurv,lossc
      real(rprec) :: lossFLC
      
      real(rprec) :: bfp,ftv,K,V,TauSS,Rgyro,eps,xSS,TauFLC,earg

      bfp = beq/(sin(lossc)**2.0) !Foot point field strength, nT
      ftv = (1.0/vm)**(3.0/2.0) !flux-tube volume Re/nT
      K = alam*vm*1.0e-3 !Energy [keV]

      if (ie == RCMPROTON) then
        V = (3.1e+2)*sqrt(K) !km/s
      else
        lossFLC = 0.0
        return
      endif

      !Convert V from km/s to Re/s
      V = V/(radius_earth_m*1.0e-3)

      TauSS = 3*2*ftv*bfp/V !Strong scattering lifetime [s], assuming ion w/ gamma=1

      Rgyro = (4.6e+3)*sqrt(K)/beq !Gyroradius of proton [km], assuming K in keV and beq in nT
      Rgyro = Rgyro/(radius_earth_m*1.0e-3) !In terms of Re

      eps = Rgyro/rcurv

      !Chen+ 2019
      
      !K: Mockup between Chen/Gibson, transition between eps^-5 dep. and strong scattering at kappa = sqrt(8)
      !xSS = max( (8.0*eps)**(-5.0), 1.0 )
      earg = eps**(-5.0)
      xSS = max(100.0*earg,1.0)
      !xSS = max(10.0*earg,1.0)

      TauFLC = xSS*TauSS
      lossFLC = 1.0/TauFLC !Rate, 1/s

    END FUNCTION FLCRat

END MODULE lossutils