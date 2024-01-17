"""
A collection of functions to help debug and validate RAIJU
"""

import os
import h5py as h5
import numpy as np

import kaipy.kdefs as kd


#------
# Analytic functions
#------

def calc_bVol_ana(colat: float, b_surf: float) -> float:
    """
    Analytic calculation of flux tube volume [Rp/nT]
    colat: colatitude [radians]
    b_surf: Surface magnetic field [nT]
    """
    cSum =  (35.0      *np.cos(1.0*colat) -      7.0 *np.cos(3.0*colat) \
            +(7.0/5.0) *np.cos(5.0*colat) - (1.0/7.0)*np.cos(7.0*colat))/64.0
    S8 = np.sin(colat)**8.0
    bVol = 2*cSum/S8/b_surf  # [Rp/nT]
    return bVol


def calc_grad_bVol_ana(colat: float, b_surf:float , Ri_m: float) -> float:
    """
    Analytic calculation of the gradient of flux tube volume [1/nT w.r.t. arc length]
    colat: Colatitude [radians]
    b_surf: Surface magnetic field [nT]
    Ri_m: Ionospheric radius [meters]
    """
    cSum =  (35.0      *np.cos(1.0*colat) -      7.0 *np.cos(3.0*colat) \
            +(7.0/5.0) *np.cos(5.0*colat) - (1.0/7.0)*np.cos(7.0*colat))/64.0
    S8 = np.sin(colat)**8.0

    dSum = (-35.0*np.sin(1.0*colat) + 21.0*np.sin(3.0*colat) \
            - 7.0*np.sin(5.0*colat) +  1.0*np.sin(7.0*colat) )/64.0
    grad_bVol = (2.0/S8/b_surf) * (-8.0/np.tan(colat)*cSum + dSum) / Ri_m  # 1/Ri_m  w.r.t. arc length
    return grad_bVol


def vCorot(Rp_m: float, Ri_m: float, cPot: float, b_surf: float, r: float, doIono=False) -> float:
    """
    Calculate corotation velocity [m/s] (default to equatorial velocity)
    Rp_m: Planetary radius [m]
    Ri_m: Ionospheric radius [m]
    cPot: Corotation potential [kV]
    b_surf: surface B-field [nT]
    r: eval equatorial radius [Rp]
    doIono: If true, returns velocity in the ionosphere instead

    vCorot [rad/s] = cPot / b_surf / Rp^2
    vCorot [m/s] = cPot / b_surf * r [Rp] / Rp^2 = cPos / b_surf * r / Rp
    """
    v_eq = (cPot*1e3) / (b_surf * 1e-9) * (r*Rp_m) / Rp_m**2  # [m/s]

    if not doIono:
        return v_eq
    else:
        sinColat = np.sqrt(1/r)
        r_iono = Ri_m * sinColat

        v_iono = v_eq * r_iono / (r*Rp_m)
        return v_iono
    

def vGradEq(Rp_m: float, Ri_m: float, signQ: float, b_surf: float, E: float, r: float, doIono=False) -> float:
    """
    Calculate gradient-drift velocity of equatorial particle [m/s] (default to equatorial velocity)
    Rp_m: Planetary radius [m]
    Ri_m: Ionospheric radius [m]
    signQ: charge sign
    b_surf: surface B-field [nT]
    E: Particle equatorial energy [keV]
    r: radial distance [Rp]
    doIono: If true, returns velocity in the ionosphere instead

    vd = 1/2* v_perp * r_g * (B x grad-B)/B^2

    1/2 * v_perp * r_g = E/qB
    grad-B/B = -3/r

    vd = -3*E*r^2 / (q * B-surf)
    
    If E in eV, then = -3*E*r^2/B-surf

    """
    sinColat = np.sqrt(1/r)
    #cosColat = np.cos(np.arcsin(sinColat))
    r_iono = Ri_m * sinColat

    #b_iono = b_surf * np.sqrt(1 + 3*sinColat**2)
    v_eq = signQ * -3 * (E * 1e3) * r**2 /(b_surf*1e-9 * rp_m)  # [m/s]
    
    if not doIono:
        return v_eq
    else:
        v_iono = v_eq * r_iono / (r*Rp_m)
        return v_iono
    

def vGCiso(Rp_m: float, Ri_m: float, b_surf: float, alamc: float, r: float, doIono=False) -> float:
    """
    Calculate gradient-curvature drift of an isotropic distribution in a dipole field
    eq. 8.1.40 from chapter 8 of Wolf's unpublished book :)

    Rp_m: [m]
    Ri_m: [m]
    b_surf: [nT]
    alamc: eV * (Rp/nT)^(2/3)
    r: Radial distance [Rp]
    doIono: If true, returns velocity in the ionosphere instead

    v_d = B x grad-Wk / qB^2  --> grad-Wk / |B|

    lambda in eV so can drop the q

    grad-Wk = grad(lambda * bVol^-2/3) = lambda*(-2/3*bVol^(-5/3)) * grad(bVol)
    """

    sinColat = np.sqrt(1/r)
    r_iono = Ri_m * sinColat
    colat = np.arcsin(sinColat)

    bVol      = calc_bVol_ana     (colat, b_surf)  # [Rp/nT]
    grad_bVol = calc_grad_bVol_ana(colat, b_surf, Ri_m)  # [Rp/nT/m]
    

    grad_Wk = alamc*(-2./3.*bVol**(-5./3.)) * grad_bVol  # [eV/m]

    cosdip = 2.0*np.cos(colat)/np.sqrt(1.0 + 3.0*np.cos(colat)**2.0) 
    Bmag = (b_surf*1e-9) / (Ri_m/Rp_m)**3 * np.sqrt(1.0+3.0*np.cos(colat)**2.0)

    v_drift = -1*grad_Wk / cosdip / Bmag  # [m/s]

    # Already ionospheric speed. Map to equator if needed
    if doIono:
        return v_drift
    else:
        v_eq = v_drift * (r*Rp_m) / r_iono
        return v_eq