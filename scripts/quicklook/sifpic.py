import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
from argparse import RawTextHelpFormatter

import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv
import kaipy.kaiTools as kt
import kaipy.kaijson as kj

import kaipy.sif.sifutils as su



if __name__=="__main__":

    fdir = "."
    outdir = "."
    jdir = "jstore"
    ftag = "msphere"

    fcIDs = ['info']  # IDs for flrceCalc checking

    # Arg parse
    MainS = """Creates visualization of ground dB
    NOTE: Assumes ground dB has been calculated using calcdb.x on simulation data.
    """

    parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
    parser.add_argument('-o',type=str,metavar="directory",default=outdir,help="Subdirectory to write to (default: %(default)s)")
    parser.add_argument('-jdir',type=str,metavar="directory",default=jdir,help="Json directory (default: %(default)s)")
    parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
    parser.add_argument('-forceCalc',type=str,metavar=fcIDs,default="",help="Things to force recalc of (default: %(default)s)")
    
    #Finalize parsing
    args = parser.parse_args()
    fdir = args.d
    outdir = args.o
    jdir = args.jdir
    ftag = args.id
    fcStr = args.forceCalc

    # Build forceCalc array
    if fcStr == "":
        fcIDs = []
    elif "all" in fcStr:
        fcIDs = fcIDs
    else:
        fcIDs = fcStr.split(',')

    # Init file stuff
    sif5fname = ftag+".sif.h5"
    kh5.CheckOrDie(sif5fname)
    kh5.CheckDirOrMake(outdir)
    kh5.CheckDirOrMake(jdir)

    # Scrape info
    jfname = os.path.join(jdir, "sifInfo.json")
    if os.path.exists(jfname) and 'info' not in fcIDs:
        sifInfo = kh5.H5Info.from_dict(kj.load(jfname))
    else:
        sifInfo = kh5.H5Info.getInfo(sif5fname)
        kj.dump(jfname, sifInfo.to_dict())

    sifInfo.printStepInfo()

    # Start with just the first step
    sif5 = h5.File(sif5fname)
    s5 = sif5[sifInfo.stepStrs[-1]]
    colat = sif5['X'][:]
    lon = sif5['Y'][:]

    mask = su.getMask(s5)

    espot_m = np.ma.masked_where(mask, s5['espot'][:])

    bVol_m = np.ma.masked_where(mask, s5['bVol'][:])
    bVolNorm = kv.genNorm(np.min(bVol_m), np.max(bVol_m),doLog=True)

    Pmhd_e = np.ma.masked_where(mask, s5['Pavg_in'][:][1])
    pmhd_eNorm = kv.genNorm(np.min(Pmhd_e), np.max(Pmhd_e))

    Pmhd_p = np.ma.masked_where(mask, s5['Pavg_in'][:][2])
    pmhd_pNorm = kv.genNorm(np.min(Pmhd_p), np.max(Pmhd_p))

    # Plot setup
    xBnd = [-20,10]

    fig = plt.figure(1,figsize=(14,10))
    gsMain = gridspec.GridSpec(3,3, wspace=1, hspace=0.3)

    AxLC = fig.add_subplot(gsMain[0,0])
    AxPolar = fig.add_subplot(gsMain[1,0], projection='polar')

    # Electron pressure
    gsxymin1 = gsMain[2,0].subgridspec(1,8)
    AxXYMin1   = fig.add_subplot(gsxymin1[0,:-1])
    AxXYMin1CB = fig.add_subplot(gsxymin1[0,-1])

    # Proton pressure
    gsxymin2 = gsMain[2,1].subgridspec(1,8)
    AxXYMin2   = fig.add_subplot(gsxymin2[0,:-1])
    AxXYMin2CB = fig.add_subplot(gsxymin2[0,-1])


    su.plotLonColat(AxLC, lon, colat, s5['OCBDist'][:])
    su.plotIono(AxPolar, lon, colat, espot_m)

    su.plotXYMin(AxXYMin1, s5['xmin'][:], s5['ymin'][:], Pmhd_e, norm=pmhd_eNorm)
    kv.setBndsByAspect(AxXYMin1, xBnd)
    kv.genCB(AxXYMin1CB,pmhd_eNorm,"Electron\npressure [nPa]",doVert=True)

    su.plotXYMin(AxXYMin2, s5['xmin'][:], s5['ymin'][:], Pmhd_p, norm=pmhd_pNorm)
    kv.setBndsByAspect(AxXYMin2, xBnd)
    kv.genCB(AxXYMin2CB,pmhd_pNorm,"Proton\npressure [nPa]",doVert=True)
    
    
    plt.show()
