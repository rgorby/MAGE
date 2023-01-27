#!/usr/bin/env python


"""Make a quick figure of a Gamera magnetosphere run.

Make a quick figure of a Gamera magnetosphere run.

Author
------
Kareem Sorathia (kareem.sorathia@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import os

# Import supplemental modules.
import matplotlib as mpl
mpl.use('Agg')  # Create figures in memory.
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

# Import project-specific modules.
import kaipy.gamera.gampp as gampp
import kaipy.gamera.rcmpp as rcmpp
import kaipy.kaiH5 as kh5
import kaipy.kaiTools as ktools
import kaipy.kaiViz as kv


# Program constants and defaults

# Program description.
description = """Creates simple multi-panel figure for RCM magnetosphere run
    Top Panel - XXX
    Bottom Panel - XXX
"""

# Default identifier for results to read.
default_runid = "msphere"

# Plot the last step by default.
default_step = -1

# Name of plot output file.
fOut = "qkrcmpic.png"


def create_command_line_parser():
    """Create the command-line argument parser.
    
    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.
    """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-beta', action='store_true', default=False,
        help="Show beta instead of FTE (default: %(default)s)"
    )
    parser.add_argument(
        "-big", action="store_true", default=False,
        help="Show entire RCM grid (default: %(default)s)."
    )
    parser.add_argument(
        '-bmin', action='store_true', default=False,
        help="Show B-min (default: %(default)s)"
    )
    parser.add_argument(
        "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-d", type=str, metavar="directory", default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        '-elec', action='store_true', default=False,
        help="Show electron pressure (default: %(default)s)"
    )
    parser.add_argument(
        '-fac', action='store_true', default=False,
        help="Show FAC (default: %(default)s)"
    )
    parser.add_argument(
        "-id", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        '-kt' , action='store_true', default=False,
        help="Show temperature instead of FTE (default: %(default)s)"
    )
    parser.add_argument(
        "-n", type=int, metavar="step", default=default_step,
        help="Time slice to plot (default: %(default)s)"
    )
    parser.add_argument(
        '-tbnc', action='store_true', default=False,
        help="Show Tb instead of FTE (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        '-vol', action='store_true', default=False,
        help="Show FTV instead of FTE (default: %(default)s)"
    )
    parser.add_argument(
        '-wgt', action='store_true', default=False,
        help="Show wRCM instead of FTE (default: %(default)s)"
    )
    return parser


if __name__ == "__main__":
    """Make a quick figure of a Gamera magnetosphere run."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    doBeta = args.beta
    doBig = args.big
    doBMin = args.bmin
    fdir = args.d
    debug = args.debug
    doElec = args.elec
    doFAC = args.fac
    ftag = args.id + ".mhdrcm"
    doT   = args.kt
    nStp = args.n
    doTb   = args.tbnc
    verbose = args.verbose
    doVol = args.vol
    doWgt = args.wgt
    if debug:
        print("args = %s" % args)

    # Defaults from external modules.
    MHDCol = rcmpp.MHDCol
    MHDLW = rcmpp.MHDLW
    rcmpp.doEll = not doBig

    # Figure parameters
    xTail = -20.0
    xSun = 10.0
    yMax = 15.0
    xyBds = [xTail, xSun, -yMax, yMax]
    figSz = (12, 6)
    eCol = "slategrey"
    eLW = 0.15
    cLW = 0.5
    vP = kv.genNorm(1.0e-1, 1.0e+2, doLog=True)
    vS = kv.genNorm(0.0, 0.25)
    vW = kv.genNorm(0, 1)
    vV = kv.genNorm(1.0e-2, 1.0, doLog=True)
    vT = kv.genNorm(0, 50)
    vB = kv.genNorm(1.0e-2, 1.0e+2, doLog=True)
    vI = kv.genNorm(0, 180)
    vBM = kv.genNorm(0, 100)
    vFAC = kv.genNorm(2.0)
    Nc = 10
    nMin = 1.0
    nMax = 1.0e+3
    vD = kv.genNorm(nMin, nMax, doLog=True)
    cVals = np.logspace(1.0, 3.0, Nc)
    pCMap = "viridis"
    sCMap = "terrain"
    dCMap = "cool"
    wCMap = "bwr_r"
    vCMap = "gnuplot2"

    # Read the RCM results.
    rcmdata = gampp.GameraPipe(fdir, ftag)
    if debug:
        print("rcmdata = %s" % rcmdata)
    if nStp < 0:
        nStp = rcmdata.sFin
        print("Using Step %d"%(nStp))
    if debug:
        print("nStp = %s" % nStp)

    # Create the figure.
    fig = plt.figure(figsize=figSz)

    # Create the grid for laying out the subplots.
    gs = gridspec.GridSpec(3, 3, height_ratios=[20, 1.0, 1.0], hspace=0.025)

    # Create the Axes objects for the individual plots.
    AxL = fig.add_subplot(gs[0, 0])
    AxM = fig.add_subplot(gs[0, 1])
    AxR = fig.add_subplot(gs[0, -1])

    # Create the Colorbar Axes.
    AxC1 = fig.add_subplot(gs[-1, 0])
    AxC2 = fig.add_subplot(gs[-1, 1])
    AxC3 = fig.add_subplot(gs[-1, -1])

    # Create the Colorbars``.
    kv.genCB(AxC1, vP, "Pressure [nPa]", cM=pCMap)
    kv.genCB(AxC2, vD, "Density [#/cc]", cM=dCMap)
    if doWgt:
        kv.genCB(AxC3, vW, r"wRCM", cM=wCMap)
    elif doVol:
        kv.genCB(AxC3, vV, r"Flux-Tube Volume [Re/nT]", cM=vCMap)
    elif doT:
        kv.genCB(AxC3, vT, r"Temperature [keV]", cM=vCMap)
    elif doBeta:
        kv.genCB(AxC3, vB, r"Beta", cM=wCMap)
    elif doTb:
        kv.genCB(AxC3, vI, r"Tb", cM=sCMap)
    elif doBMin:
        kv.genCB(AxC3, vBM, r"B-Minimum [nT]", cM=sCMap)
    elif doFAC:
        kv.genCB(AxC3, vFAC, r"FAC [uA/m2]", cM=wCMap)
    else:
        kv.genCB(AxC3, vS, r"Flux-Tube Entropy [nPa (R$_{E}$/nT)$^{\gamma}$]", cM=sCMap)

    # Clear the subplots (why is this needed?)
    AxL.clear()
    AxM.clear()
    AxR.clear()

    # Fetch the coordinates to plot.
    bmX, bmY = rcmpp.RCMEq(rcmdata, nStp, doMask=True)
    I = rcmpp.GetMask(rcmdata, nStp)
    Ni = (~I).sum()
    if debug:
        print("bmX = %s" % bmX)
        print("bmY = %s" % bmY)
        print("I = %s" % I)
        print("Ni = %s" % Ni)

    # Abort if a closed-field region does not exist.
    if Ni == 0:
        print("No closed field region in RCM, exiting ...")
        exit()

    # Fetch the data to plot.
    if doElec:
        Prcm = rcmpp.GetVarMask(rcmdata, nStp, "Pe", I)
    else:
        Prcm = rcmpp.GetVarMask(rcmdata, nStp, "P", I)
    Nrcm = rcmpp.GetVarMask(rcmdata, nStp, "N", I)
    Pmhd = rcmpp.GetVarMask(rcmdata, nStp, "Pmhd", I)
    Nmhd = rcmpp.GetVarMask(rcmdata, nStp, "Nmhd", I)
    S = rcmpp.GetVarMask(rcmdata,nStp, "S", I)
    toMHD = rcmpp.GetVarMask(rcmdata, nStp, "toMHD", I)
    pot, pVals = rcmpp.GetPotential(rcmdata, nStp, I)
    wRCM = None
    if doWgt:
        wRCM = rcmpp.GetVarMask(rcmdata, nStp, "wIMAG", I)
    bVol = None
    if doVol:
        bVol = rcmpp.GetVarMask(rcmdata, nStp, "bVol", I)
    beta = None
    if doBeta:
        beta = rcmpp.GetVarMask(rcmdata, nStp, "beta", I)
    Tb = None
    if doTb:
        Tb = rcmpp.GetVarMask(rcmdata,nStp,"Tb", I)
    Bmin = None
    if doBMin:
        Bmin = rcmpp.GetVarMask(rcmdata, nStp, "bMin", I)
    toRCM = None
    if doBig:
        toRCM = rcmpp.GetVarMask(rcmdata, nStp, "IOpen", I)
    jBirk = None
    if doFAC:
        jBirk = rcmpp.GetVarMask(rcmdata, nStp, "birk", I)
    if debug:
        print("Prcm = %s" % Prcm)
        print("Nrcm = %s" % Nrcm)
        print("Pmhd = %s" % Pmhd)
        print("Nmhd = %s" % Nmhd)
        print("S = %s" % S)
        print("toMHD = %s" % toMHD)
        print("pot = %s" % pot)
        print("pVals = %s" % pVals)
        print("wRCM = %s" % wRCM)
        print("bVol = %s" % bVol)
        print("beta = %s" % beta)
        print("Tb = %s" % Tb)
        print("Bmin = %s" % Bmin)
        print("toRCM = %s" % toRCM)
        print("jBirk = %s" % jBirk)

    # Read the dates the data file.
    fStr = os.path.join(fdir, ftag + '.h5')
    if debug:
        print("fStr = %s" % fStr)
    MJD = kh5.tStep(fStr, nStp, aID="MJD")
    if debug:
        print("MJD = %s" % MJD)
    utS = ktools.MJD2UT([MJD])
    if debug:
        print("utS = %s" % utS)
    utDT= utS[0]
    if debug:
        print("utDT = %s" % utDT)

    # Assemble the left-hand plot.
    AxL.set_title("RCM Pressure")
    AxL.pcolor(bmX, bmY, Prcm, norm=vP, cmap=pCMap)
    AxL.contour(bmX, bmY, pot, pVals, colors='grey', linewidths=cLW)
    kv.addEarth2D(ax=AxL)
    kv.SetAx(xyBds, AxL)

    # Assemble the middle plot.
    AxM.set_title("MHD Pressure")
    AxM.pcolor(bmX, bmY, Pmhd, norm=vP, cmap=pCMap)
    AxM.contour(bmX, bmY, Nmhd, cVals, norm=vD, cmap=dCMap, linewidths=cLW)
    kv.addEarth2D(ax=AxM)
    kv.SetAx(xyBds, AxM)

    # Gather all Axes.
    Axs = [AxL, AxM, AxR]

    # Plot the full domain if needed.
    if nStp > 0 and doBig:
        for Ax in Axs:
            CS1 = Ax.contour(bmX, bmY, toMHD, [0.5], colors=MHDCol, linewidths=MHDLW)
            manloc = [(0.0, 8.0)]
            fmt = {}
            fmt[0.5] = 'MHD'  # Key on a float is a bad idea.
            Ax.clabel(
                CS1, CS1.levels[::2], inline=True, fmt=fmt, fontsize=5,
                inline_spacing=25, manual=manloc
            )
            CS2 = Ax.contour(
                bmX, bmY, toRCM, [-0.5], colors=rcmpp.rcmCol,
                linewidths=MHDLW, linestyles='solid'
            )

    # Assemble the right-hand plot.
    if doWgt:
        AxR.set_title("RCM Weight")
        AxR.pcolor(bmX, bmY, wRCM, norm=vW, cmap=wCMap)
    elif doVol:
        AxR.set_title("Flux-tube Volume")
        AxR.pcolor(bmX, bmY, bVol, norm=vV, cmap=vCMap)
    elif doT:
        kT = 6.25*Prcm/Nrcm
        AxR.set_title("RCM Temperature")
        AxR.pcolor(bmX, bmY, kT, norm=vT, cmap=vCMap)
    elif doBeta:
        AxR.set_title("Average Beta")
        AxR.pcolor(bmX, bmY, beta, norm=vB, cmap=wCMap)
    elif doTb:
        AxR.set_title("Ingestion timescale")
        AxR.pcolor(bmX, bmY, Tb, norm=vI, cmap=sCMap)
    elif doBMin:
        AxR.set_title("B Minimum")
        AxR.pcolor(bmX, bmY, 1.0e+9*Bmin, norm=vBM, cmap=sCMap)
    elif doFAC:
        AxR.set_title("Vasyliunas FAC")
        AxR.pcolor(bmX, bmY, jBirk, norm=vFAC, cmap=wCMap)
    else:
        AxR.set_title("Flux-Tube Entropy")
        AxR.pcolor(bmX, bmY, S, norm=vS, cmap=sCMap)    
    AxR.plot(bmX, bmY, color=eCol, linewidth=eLW)
    AxR.plot(bmX.T, bmY.T, color=eCol, linewidth=eLW)
    kv.addEarth2D(ax=AxR)
    kv.SetAx(xyBds, AxR)

    # Create the title for the complete figure.
    tStr = "\n\n\n" + utDT.strftime("%m/%d/%Y, %H:%M:%S")
    plt.suptitle(tStr, fontsize="x-large")

    # Save the figure to a file.
    kv.savePic(fOut)
