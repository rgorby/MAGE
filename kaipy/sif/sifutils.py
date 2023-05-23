import numpy as np
import matplotlib.pyplot as plt


import kaipy.kaiTools as kt

domain = {"INACTIVE" : -1,
          "BUFFER" : 0,
          "ACTIVE" : 1}


def getMask(s5, dom="ACTIVE"):
    # s5 = step#X group object

    mask = s5['active'][:] != domain[dom]
        # Only include domain specified by caller
    return mask

# TODO:
#   Calc vcorot, vgc, veffective

def plotLonColat(Ax, lon, colat, var, unit="deg", norm=None, cmap='viridis'):
    """ Simple plot function in longitude - colatitude projection
        lon and colat inputs should be in radians
        lon and colat should be size Ni+1,Nj+1
    """

    if unit=="deg":
        uStr = "[deg]"
        colat = colat*180/np.pi
        lon   = lon  *180/np.pi
    else:
        uStr = "[rad]"

    cMin = np.min(colat)
    cMax = np.max(colat)

    Ax.pcolormesh(lon, colat, var, norm=norm, cmap=cmap)
    Ax.set_xlabel("Longitude " + uStr)
    Ax.set_ylabel("Colatitude " + uStr)
    Ax.set_ylim([cMin, cMax])

def plotIono(Ax, lon, colat, var, norm=None, cmap='viridis', Riono=1):
    """ Simple plot function to show ionosphere projection
        lon and colat inputs should be in radians
        lon and colat should be size Ni+1,Nj+1
    """

    # First get coordinates in polar projection
    r, theta = kt.rtp2rt(Riono, colat, lon)
    Ax.pcolormesh(theta, r, var, norm=norm, cmap=cmap)
    Ax.axis([0,2*np.pi,0,r.max()])
    Ax.grid(True)

    # Re-do Colat labels
    colatStride = 10*np.pi/180  # rad
    colatMax = np.max(colat) # rad
    colatTickLocs_rad = np.arange(10*np.pi/180, colatMax, colatStride)
    rTickLocs = Riono*np.sin(colatTickLocs_rad)
    rTickLabels = ["{:2n}".format(c*180/np.pi) for c in colatTickLocs_rad]
    Ax.set_yticks(rTickLocs)
    Ax.set_yticklabels(rTickLabels)

def plotXYMin(Ax, xmin, ymin, var, norm=None, cmap='viridis'):
    """ Simple plot function to show xy bmin projection
        If xBnds or yBnds is not none, it will use x or y bnds and rescale the other so that aspect ratio is 1
    """

    Ax.pcolormesh(xmin, ymin, var, norm=norm, cmap=cmap)
    Ax.set_xlabel('X [R$_p$]')
    Ax.set_ylabel('Y [R$_p$]')