import numpy as np

def L_to_bVol(L):  # L shell [Re] to V [Re/nT]
    bsurf_nT = 3.11E4
    colat = np.arcsin(np.sqrt(1.0/L))

    cSum = 35*np.cos(colat) - 7*np.cos(3*colat) +(7./5.)*np.cos(5*colat) - (1./7.)*np.cos(7*colat)
    cSum /= 64.
    s8 = np.sin(colat)**8
    V = 2*cSum/s8/bsurf_nT
    return V
