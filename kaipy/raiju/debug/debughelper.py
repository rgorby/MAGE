"""
A collection of functions to help debug and validate RAIJU
"""

import os
import h5py as h5
import numpy as np

import kaipy.kdefs as kd
import kaipy.kaiTools as kt


#------
# File/run info
#------

def getAttrs(f5name: str) -> dict:

    with h5.File(f5name) as f5:
        attrs = {}
        for k in f5.attrs.keys():
            attrs[k] = f5.attrs[k]
    return attrs

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
    

#------
# Plotting
#------
    
def drawCell(Ax, x, y, var, norm=None, cmap=None, active=None, iono=True):
    # Draws values of faces on grid lines
    # Assumes:
    # x.shape and y.shape: (Ni+1, Nj+1)
    # var.shape: (Ni, Nj)
    """
    col = cm.jet((c-np.min(c))/(np.max(c)-np.min(c)))
    ax = plt.gca()
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]], c=col[i])
    """

    xcc = kt.to_center2D(x)
    ycc = kt.to_center2D(y)

    # Face lines
    Ax.pcolormesh(x, y, var, norm=norm, cmap=cmap)

    # Draw contours for active domain if we have it and its actually interesting
    if active is not None and np.min(active) != np.max(active):
        lvls = [-1, 0, 1.0]
        Ax.contour(xcc, ycc, active, levels=lvls)


    if iono:
        ylim = Ax.get_ylim()
        Ax.set_ylim(np.max(ylim), np.min(ylim))
        Ax.set_ylabel('colat [deg]')
        Ax.set_xlabel('lon [deg]')


def drawFaces(Ax, x, y, var, direction=None, norm=None, cmap=None, active=None, iono=True):
    # Draws values of faces on grid lines
    # Assumes:
    # x.shape and y.shape: (Ni+1, Nj+1)
    # var.shape: (Ni+1, Nj+1, 2)
    import matplotlib.cm as cm
    from matplotlib.collections import LineCollection

    Ni, Nj = x.shape
    Ni -= 1; Nj -=1

    eCol = 'white'
    alpha = 1.0
    eLW = 4
    Ax.plot(x  , y  , c=eCol, alpha=alpha, linewidth=eLW)
    Ax.plot(x.T, y.T, c=eCol, alpha=alpha, linewidth=eLW)

    if cmap is None:
        cmap = 'viridis'
    pltcmap = getattr(cm, cmap)
    clrs = pltcmap( (var - np.min(var))/(np.max(var)-np.min(var))  )

    xcc = kt.to_center2D(x)
    ycc = kt.to_center2D(y)

    # Face lines
    iLW = 2
    
    # Theta dir
    if direction is None or direction=="theta":
        segments = np.zeros((2,2))[np.newaxis,:,:]  # dead first point just to initialize
        for i in range(Ni+1):
            tmp = [ [[x[i,j],y[i,j]], [x[i,j+1],y[i,j+1]]] for j in range(Nj)   ]
            segments = np.concatenate((segments, tmp), axis=0)
        segments = segments[1:,:,:]
        lc = LineCollection(segments, cmap=cmap, norm=norm)
        # Set the values used for colormapping
        #lc.set_array(var[:-1,:,0].reshape((Nj+1)*Ni))
        lc.set_array(var[:,:Nj,0].flatten())
        lc.set_linewidth(iLW)
        line = Ax.add_collection(lc)

    # Phi dir
    if direction is None or direction=="phi":
        segments = np.zeros((2,2))[np.newaxis,:,:]  # dead first point just to initialize
        for j in range(Nj+1):
            tmp = [ [[x[i,j],y[i,j]], [x[i+1,j],y[i+1,j]]] for i in range(Ni)   ]
            segments = np.concatenate((segments, tmp), axis=0)
        segments = segments[1:,:,:]
        lc = LineCollection(segments, cmap=cmap, norm=norm)
        # Set the values used for colormapping
        #lc.set_array(var[:-1,:,0].reshape((Nj+1)*Ni))
        lc.set_array(var[:Ni,:,1].T.flatten())
        lc.set_linewidth(iLW)
        line = Ax.add_collection(lc)
    

    # Draw contours for active domain if we have it and its actually interesting
    if active is not None and np.min(active) != np.max(active):
        lvls = [-1, 0, 1.0]
        Ax.contour(xcc, ycc, active, levels=lvls)


    if iono:
        ylim = Ax.get_ylim()
        Ax.set_ylim(np.max(ylim), np.min(ylim))
        Ax.set_ylabel('colat [deg]')
        Ax.set_xlabel('lon [deg]')

    return line