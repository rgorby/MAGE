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
    s5 = sif5[sifInfo.stepStrs[1]]
    colat = sif5['X'][:]
    lon = sif5['Y'][:]

    espot = s5['espot'][:]
    mask = su.getMask(s5)
    espot_m = np.ma.masked_where(mask, espot)

    # Plot setup
    fig = plt.figure(1,figsize=(14,10))
    
    gs = gridspec.GridSpec(3,3, wspace=1, hspace=0.3)

    AxLC = fig.add_subplot(gs[0,0])
    AxPolar = fig.add_subplot(gs[1,0], projection='polar')
    AxXYMin = fig.add_subplot(gs[2,0])

    su.plotLonColat(AxLC, lon, colat, s5['OCBDist'][:])
    su.plotIono(AxPolar, lon, colat, s5['espot'][:])
    su.plotXYMin(AxXYMin, s5['xmin'][:], s5['ymin'][:], espot_m)
    kv.setBndsByAspect(AxXYMin, [-20, 20])
    #plt.colorbar()
    plt.show()

    """
    colat = sif5['X'][:]
    lon = sif5['Y'][:]
    colat_deg = colat*180/np.pi
    lon_deg = lon*180/np.pi
    plt.pcolormesh(lon_deg,colat_deg,s5['OCBDist'][:])
    plt.ylim([np.max(colat_deg), np.min(colat_deg)])
    plt.colorbar()
    plt.show()
    """