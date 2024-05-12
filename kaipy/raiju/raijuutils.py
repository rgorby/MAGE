import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from dataclasses_json import dataclass_json
from typing import List
import scipy.special as sp

#from bidict import bidict

import kaipy.kdefs as kd
import kaipy.kaiTools as kt
import kaipy.kaiH5 as kh5

import kaipy.raiju.lambdautils.AlamParams as aP

dim = {"THETA": 0,
       "PHI": 1}

topo = {"OPEN" : 0,
        "CLOSED" : 1}

domain = {"INACTIVE" : -1,
          "BUFFER" : 0,
          "ACTIVE" : 1}

flavs_s = {"PSPH" : 0,  # Flav dict, lookup by string name
           "HOTE" : 1,
           "HOTP" : 2}
#flavs_n = bidict(flavs_s).inv  # Flav dict, lookup by index
flavs_n = {0 : "PSPH",
           1 : "HOTE",
           2 : "HOTP"}

spcs_s = {"IDK": 0,
          "ELE": 1,
          "H+" : 2,
          "O+" : 3}
spcs_n = {0: "IDK",
          1: "ELE",
          2: "H+" ,
          3: "O+" }

#------
# Containers
#------

@dataclass_json
@dataclass
class SpeciesInfo:
    N: int
    flav: int
    spcType: int
    kStart: int
    kEnd: int
    numNuc_p: int
    numNuc_n: int
    amu: float
    q: int
    alami: List[float]
    alamc: List[float]
    name: str = None

@dataclass_json
@dataclass
class RAIJUInfo(kh5.H5Info):
    nSpc: int = None
    species: List[SpeciesInfo] = None
    Nk: int = 0
    planet: dict = None
    extra: dict = None # Extra info we can add later
    

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
                alami = s5['alami'][:]
                alamc = 0.5*(alami[1:]+alami[:-1])
                spc = SpeciesInfo(att['N'], att['flav'], att['spcType'],
                                  att['kStart'], att['kEnd']+1,
                                  att['numNuc_p'], att['numNuc_p'],
                                  att['amu'], att['q'],
                                  alami, alamc, 
                                  name)
                specs.append(spc)
            
            planetInfo = {}
            for k in f5['Planet'].attrs.keys():
                v = f5['Planet'].attrs[k]
                if isinstance(v, bytes):
                    v = v.decode('utf-8')
                planetInfo[k] = v

        Nk = 0
        for s in specs: Nk += s.N
        # Now make our final object
        return RAIJUInfo(fi.fname, fi.Nt, fi.steps, fi.stepStrs, fi.times, fi.MJDs, fi.UTs,
                       len(specs), specs, Nk, planetInfo, {})


#------
# Data handlers
#------

def getVar(grp: h5.Group, varName: str, mask:bool=None, broadcast_dims=None) -> np.ndarray:
    """ Get vars from h5 file this way so that we're sure everyone agrees on type and shape
    """
    try:
        var = grp[varName][:].T  # .T to go from Fortran to python indexing order
        if mask is not None:
            if broadcast_dims is not None:
                for d in broadcast_dims:
                    mask = np.expand_dims(mask, axis=d)
                mask = np.broadcast_to(mask, var.shape)
            var = np.ma.masked_where(mask, var)
        return var
    except KeyError:
        print("Error: {} not in keys".format(varName))
        return None

#------
# Species helpers
#------

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

    mask = getVar(s5, 'active') != domain[dom]
        # Only include domain specified by caller
    return mask


#------
# Some analytic stuff
#------

def intensity_maxwell(n, mass, E, kT):
    """ Intensity for energy E in an analytic Maxwellian profile
        n: density in #/cc
        mass: mass in kg
        E: energy in keV
        kT: temp in keV
        j: Intensity [1/(s*sr*keV*cm^2)]
    """
    f = n * (mass/(2*np.pi*kT))**(3/2) * np.exp(-E/kT)
    j = 2*E/mass**2 * f  * 1e2 * np.sqrt(kd.kev2J)
    return j

def intensity_kappa(n, mass, E, kT, kappa=6):
    """ Intensity for energy E in an analytic Maxwellian profile
        n: density in #/cc
        mass: mass in kg
        E: energy in keV
        kT: temp in keV
        j: Intensity [1/(s*sr*keV*cm^2)]
    """
    gamfac = sp.gamma(kappa+1)/sp.gamma(kappa-0.5)

    #f = n * (mass/(2*np.pi*kappa*kT))**(3/2) * gamfac * (1+(E/kappa/kT))**(-kappa-1)  # From Baumjohann & Treumann
    kap15 = kappa-1.5
    E0 = kT*kap15/kappa
    kArg = 1 + (E/E0)/kap15
    f = n * (mass/(2*np.pi*E0*kap15))**(3/2) * gamfac * kArg**(-kappa-1)

    j = 2*E/mass**2 * f  * 1e2 * np.sqrt(kd.kev2J)
    return j

#------
# Conversions
#------

def etak2Press(etak, alamc, bVol):
    """
    etak [#/cc * Rx/T]
    alamc: [eV*(Rx/nT)^(2/3)]
    bVol [Rx/nT]

    Returns: pressure [nPa]
    """
    return 2./3.*etak*alamc*bVol**(-5./3.) * kd.ev2J * 1.e6

def etak2Den(etak, bVol):
    """
    etak [#/cc * Rx/T]
    bVol [Rx/nT]

    Returns: density [#/cc]
    """

    return (etak*1.0E-9)/bVol  # [#/cc]

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

def plotXYMin(Ax, xmin, ymin, var, norm=None, cmap='viridis', shading='auto', lblsize=None):
    """ Simple plot function to show xy bmin projection
        If xBnds or yBnds is not none, it will use x or y bnds and rescale the other so that aspect ratio is 1
    """
    if (type(xmin)==np.ma.core.MaskedArray):
        Ax.pcolor(xmin, ymin, var, norm=norm, cmap=cmap, shading=shading, rasterized=True)
    else:
        Ax.pcolormesh(xmin, ymin, var, norm=norm, cmap=cmap, shading=shading, rasterized=True)

    if lblsize is not None:
        Ax.set_xlabel('X [R$_p$]', fontsize=lblsize)
        Ax.set_ylabel('Y [R$_p$]', fontsize=lblsize)
