import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from dataclasses_json import dataclass_json
from typing import List
from bidict import bidict

import kaipy.kaiTools as kt
import kaipy.kaiH5 as kh5

import kaipy.sif.lambdautils.AlamParams as aP

topo = {"OPEN" : 0,
        "CLOSED" : 1}

domain = {"INACTIVE" : -1,
          "BUFFER" : 0,
          "ACTIVE" : 1}

flavs_s = {"PSPH" : 0,  # Flav dict, lookup by string name
           "HOTE" : 1,
           "HOTP" : 2}
flavs_n = bidict(flavs_s).inv  # Flav dict, lookup by index


@dataclass_json
@dataclass
class SpeciesInfo:
    N: int
    flav: int
    kStart: int
    kEnd: int
    numNuc_p: int
    numNuc_n: int
    amu: float
    q: int
    alami: List[float]
    name: str = None

@dataclass_json
@dataclass
class SIFInfo(kh5.H5Info):
    species: List[aP.SpecParams] = None

    def getInfo(h5fname, noSubSec=True):
        # Base 
        fi = kh5.H5Info.getInfo(h5fname, noSubSec)
        specs = []
        with h5.File(h5fname) as f5:
            for spc_key in f5['Species'].keys():
                s5 = f5['Species'][spc_key]
                att = s5.attrs
                # Determine name if it exists
                flav = att['flav']
                if flav in flavs_n.keys():
                    name = flavs_n[flav]
                else:
                    name = None
                spc = SpeciesInfo(att['N'], att['flav'],
                                  att['kStart'], att['kEnd']+1,
                                  att['numNuc_p'], att['numNuc_p'],
                                  att['amu'], att['q'],
                                  s5['alami'][:], name)
                specs.append(spc)
        # Now make our final object
        return SIFInfo(fi.fname, fi.Nt, fi.steps, fi.stepStrs, fi.times, fi.MJDs,
                        fi.UTs, specs)
    
def spcIdx(spcList: List[SpeciesInfo], flav: int) -> int:
    # Get index of a certain species based on its flavor
    # spcList: list of SpeciesInfo
    for idx, s in enumerate(spcList):
        if s.flav == flav:
            return idx
    # If here, we didn't find index. Complain
    print("Warning: spcIdx didn't find flav '{}' in spcList".format(flav))
    return -1

def getSpcFromNkArr():
    pass

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
    #Ax.set_ylim([cMin, cMax])

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