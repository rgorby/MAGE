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
import sys
from multiprocessing import Pool
from psutil import cpu_count
import warnings

# Import project-specific modules.
import kaipy.cdaweb_utils as cdaweb_utils
import kaipy.gamera.gampp as gampp
import kaipy.gamera.rcmpp as rcmpp
import kaipy.kaiH5 as kh5
import kaipy.kaiTools as ktools
import kaipy.kaiViz as kv
import kaipy.kdefs as kdefs
import kaipy.kdefs as kd


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
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Names of spacecraft to plot trajectories, separated by commas (default: %(default)s)"
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
    parser.add_argument(
        '-vid', action='store_true', default=False,
        help="Make a video and store in mixVid directory (default: %(default)s)"
    )
    parser.add_argument(
        '-overwrite', action='store_true', default=False,
        help="Overwrite existing vid files (default: %(default)s)"
    )
    parser.add_argument(
        '--ncpus', type=int, metavar="ncpus", default=1,
        help="Number of threads to use with --vid (default: %(default)s)"
    )
    parser.add_argument(
        '-nohash', action='store_true', default=False,
        help="Don't display branch/hash info (default: %(default)s)"
    )
    return parser

def makePlot(i,rcmdata,nStp):

    if not debug:
        # Suppress the warning
        warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
        warnings.filterwarnings("ignore", message="The input coordinates to pcolor are interpreted as cell centers.*")
        warnings.filterwarnings("ignore", message="Log scale: values of z <= 0 have been masked.*")
        warnings.filterwarnings("ignore", message="No contour levels were found within the data range.*")

    # Name of plot output file.
    if do_vid:
        fOut = "{}.{:0>{n}d}.png".format("rcmpic", i, n=n_pad)
        outPath = os.path.join(outDir, fOut)
    else:
        outPath = "qkrcmpic.png"

    # Skip this file if it already exists and we're not supposed to overwrite
    if not do_overwrite and os.path.exists(outPath) and do_vid:
        return

    plt.clf()

    # Create the grid for laying out the subplots.
    gs = gridspec.GridSpec(2, 3, height_ratios=[20, 1.0], hspace=0.025)

    # Create the Axes objects for the individual plots.
    AxL = fig.add_subplot(gs[0, 0])
    AxM = fig.add_subplot(gs[0, 1])
    AxR = fig.add_subplot(gs[0, -1])

    # Create the Colorbar Axes.
    AxC1 = fig.add_subplot(gs[-1, 0])
    AxC2 = fig.add_subplot(gs[-1, 1])
    AxC3 = fig.add_subplot(gs[-1, -1])

    # Adjust the positions of the individual subplots.
    AxL.set_position([0.05, 0.1, 0.25, 1.0])
    AxM.set_position([0.35, 0.1, 0.25, 1.0])
    AxR.set_position([0.65, 0.1, 0.25, 1.0])

    # Update the subplot parameters for visibility.
    AxL.tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, labelbottom=True, labelleft=True)
    AxM.tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, labelbottom=True, labelleft=True)
    AxR.tick_params(axis='both', which='both', bottom=True, top=False, left=True, right=False, labelbottom=True, labelleft=True)

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
    try:
        bmX, bmY = rcmpp.RCMEq(rcmdata, nStp, doMask=True)
    except:
        print(f"Step #{nStp} does not exist!")
        sys.exit(1)
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

    # This is not working yet. A blank bar still shows up
    doVerb = debug
    doVerb = True

    # Fetch the data to plot.
    if doElec:
        Prcm = rcmpp.GetVarMask(rcmdata, nStp, "Pe", I, doVerb=doVerb)
    else:
        Prcm = rcmpp.GetVarMask(rcmdata, nStp, "P", I, doVerb=doVerb)
    Nrcm = rcmpp.GetVarMask(rcmdata, nStp, "N", I, doVerb=doVerb)
    Pmhd = rcmpp.GetVarMask(rcmdata, nStp, "Pmhd", I, doVerb=doVerb)
    Nmhd = rcmpp.GetVarMask(rcmdata, nStp, "Nmhd", I, doVerb=doVerb)
    S = rcmpp.GetVarMask(rcmdata,nStp, "S", I, doVerb=doVerb)
    toMHD = rcmpp.GetVarMask(rcmdata, nStp, "toMHD", I, doVerb=doVerb)
    pot, pVals = rcmpp.GetPotential(rcmdata, nStp, I, doVerb=doVerb)
    wRCM = None
    if doWgt:
        wRCM = rcmpp.GetVarMask(rcmdata, nStp, "wIMAG", I, doVerb=doVerb)
    bVol = None
    if doVol:
        bVol = rcmpp.GetVarMask(rcmdata, nStp, "bVol", I, doVerb=doVerb)
    beta = None
    if doBeta:
        beta = rcmpp.GetVarMask(rcmdata, nStp, "beta", I, doVerb=doVerb)
    Tb = None
    if doTb:
        Tb = rcmpp.GetVarMask(rcmdata,nStp,"Tb", I, doVerb=doVerb)
    Bmin = None
    if doBMin:
        Bmin = rcmpp.GetVarMask(rcmdata, nStp, "bMin", I, doVerb=doVerb)
    toRCM = None
    if doBig:
        toRCM = rcmpp.GetVarMask(rcmdata, nStp, "IOpen", I, doVerb=doVerb)
    jBirk = None
    if doFAC:
        jBirk = rcmpp.GetVarMask(rcmdata, nStp, "birk", I, doVerb=doVerb)
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
    if debug or do_vid:
        print("fStr = %s" % fStr)
    MJD = kh5.tStep(fStr, nStp, aID="MJD")
    if debug:
        print("MJD = %s" % MJD)
    utS = ktools.MJD2UT([MJD])
    if debug:
        print("utS = %s" % utS)
    utDT = utS[0]
    if debug:
        print("utDT = %s" % utDT)

    # Fetch the date of the earliest available step.
    n_steps, step_numbers = kh5.cntSteps(fStr)
    MJD = kh5.tStep(fStr, step_numbers[0], aID="MJD")
    if debug:
        print("MJD = %s" % MJD)
    utS = ktools.MJD2UT([MJD])
    if debug:
        print("utS = %s" % utS)
    ut0 = utS[0]
    if debug:
        print("ut0 = %s" % ut0)

    # If needed, fetch the trajectory of each spacecraft from CDAWeb.
    sc_X = []
    sc_Y = []
    sc_Z = []
    if spacecraft:
        spacecrafts = spacecraft.split(',')
        for (i_sc, sc) in enumerate(spacecrafts):
            if debug:
                print("i_sc, sc = %s, %s" % (i_sc, sc))

            # Fetch the spacecraft trajectory in Solar Magnetic (SM)
            # Cartesian coordinates between the start and end times.
            sc_x, sc_y, sc_z = cdaweb_utils.fetch_spacecraft_SM_trajectory(
                sc, ut0, utDT
            )
            if debug:
                print("sc_x, sc_y, sc_z = %s, %s, %s" % (sc_x, sc_y, sc_z))
            sc_X.append(sc_x)
            sc_Y.append(sc_y)
            sc_Z.append(sc_z)

            # Skip if no trajectory found.
            if sc_x is None:
                print("No trajectory found for spacecraft %s." % sc)
                continue

            # Convert coordinates to units of Earth radius.
            CM_TO_KM = 1e-5  # Centimeters to kilometers
            Re_km = kdefs.Re_cgs*CM_TO_KM  # Earth radius in kilometers
            sc_X[-1] = sc_x/Re_km
            sc_Y[-1] = sc_y/Re_km
            sc_Z[-1] = sc_z/Re_km
            if debug:
                print("sc_X, sc_Y, sc_Z = %s, %s, %s" % (sc_X, sc_Y, sc_Z))

    # Assemble the left-hand plot.
    AxL.set_title("RCM Pressure")
    AxL.pcolor(bmX, bmY, Prcm, norm=vP, cmap=pCMap, shading='auto')
    AxL.contour(bmX, bmY, pot, pVals, colors='grey', linewidths=cLW)
    kv.addEarth2D(ax=AxL)
    kv.SetAx(xyBds, AxL)

    # Assemble the middle plot.
    AxM.set_title("MHD Pressure")
    AxM.pcolor(bmX, bmY, Pmhd, norm=vP, cmap=pCMap, shading='auto')
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
        AxR.pcolor(bmX, bmY, wRCM, norm=vW, cmap=wCMap, shading='auto')
    elif doVol:
        AxR.set_title("Flux-tube Volume")
        AxR.pcolor(bmX, bmY, bVol, norm=vV, cmap=vCMap, shading='auto')
    elif doT:
        kT = 6.25*Prcm/Nrcm
        AxR.set_title("RCM Temperature")
        AxR.pcolor(bmX, bmY, kT, norm=vT, cmap=vCMap, shading='auto')
    elif doBeta:
        AxR.set_title("Average Beta")
        AxR.pcolor(bmX, bmY, beta, norm=vB, cmap=wCMap, shading='auto')
    elif doTb:
        AxR.set_title("Ingestion timescale")
        AxR.pcolor(bmX, bmY, Tb, norm=vI, cmap=sCMap, shading='auto')
    elif doBMin:
        AxR.set_title("B Minimum")
        AxR.pcolor(bmX, bmY, 1.0e+9*Bmin, norm=vBM, cmap=sCMap, shading='auto')
    elif doFAC:
        AxR.set_title("Vasyliunas FAC")
        AxR.pcolor(bmX, bmY, jBirk, norm=vFAC, cmap=wCMap, shading='auto')
    else:
        AxR.set_title("Flux-Tube Entropy")
        AxR.pcolor(bmX, bmY, S, norm=vS, cmap=sCMap, shading='auto')
    AxR.plot(bmX, bmY, color=eCol, linewidth=eLW)
    AxR.plot(bmX.T, bmY.T, color=eCol, linewidth=eLW)
    kv.addEarth2D(ax=AxR)
    kv.SetAx(xyBds, AxR)

    # Overlay spacecraft trajectories if requested.
    if spacecraft:
        spacecrafts = spacecraft.split(',')
        for (i_sc, sc) in enumerate(spacecrafts):
            if debug:
                print("i_sc, sc = %s, %s" % (i_sc, sc))

            # Skip this spacecraft if no trajectory is available.
            if sc_X[i_sc] is None:
                continue

            # Plot a labelled trajectory of the spacecraft. Also plot a larger
            # dot at the last point in the trajectory.

            # NOTE: Need to add a filter to not plot spacecraft positions
            # outside of the plot limits. Many X values are > +10, which puts
            # them off the right side of the plots.

            # Left plot
            SPACECRAFT_COLORS = list(mpl.colors.TABLEAU_COLORS.keys())
            color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
            AxL.plot(sc_X[i_sc], sc_Y[i_sc], marker=None, linewidth=1, c=color)
            AxL.plot(sc_X[i_sc][-1], sc_Y[i_sc][-1], 'o', c=color)
            x_nudge = 1.0
            y_nudge = 1.0
            AxL.text(sc_X[i_sc][-1] + x_nudge, sc_Y[i_sc][-1] + y_nudge, sc, c=color)

            # Middle plot
            AxM.plot(sc_X[i_sc], sc_Y[i_sc], marker=None, linewidth=1, c=color)
            AxM.plot(sc_X[i_sc][-1], sc_Y[i_sc][-1], 'o', c=color)
            x_nudge = 1.0
            y_nudge = 1.0
            AxM.text(sc_X[i_sc][-1] + x_nudge, sc_Y[i_sc][-1] + y_nudge, sc, c=color)

            # Right plot
            AxR.plot(sc_X[i_sc], sc_Y[i_sc], marker=None, linewidth=1, c=color)
            AxR.plot(sc_X[i_sc][-1], sc_Y[i_sc][-1], 'o', c=color)
            x_nudge = 1.0
            y_nudge = 1.0
            AxR.text(sc_X[i_sc][-1] + x_nudge, sc_Y[i_sc][-1] + y_nudge, sc, c=color)

    # Set left labels for subplots.
    ylabel = 'X [R$_E$]'
    AxL.set_ylabel(ylabel)

    # Set bottom labels for subplots.
    xlabel = 'Y [R$_E$]'
    AxL.set_xlabel(xlabel)
    AxM.set_xlabel(xlabel)
    AxR.set_xlabel(xlabel)

    # Create the title for the complete figure.
    tStr = "\n\n\n" + utDT.strftime("%m/%d/%Y, %H:%M:%S")
    plt.suptitle(tStr, fontsize="x-large")

    # Add Branch and Hash info
    if do_hash:
        fig.text(0.1,0.85,f"branch/commit: {branch}/{githash}", fontsize=4)

    # Adjust layout to reduce white space
    plt.subplots_adjust(top=0.9, bottom=0.1, hspace=0.025)

    # Save the figure to a file.
    kv.savePic(outPath, dpiQ=300)

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
    spacecraft = args.spacecraft
    doTb   = args.tbnc
    verbose = args.verbose
    doVol = args.vol
    doWgt = args.wgt
    do_vid = args.vid
    do_overwrite = args.overwrite
    do_hash = not args.nohash
    ncpus = args.ncpus
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
    #sCMap = "terrain"
    sCMap = "turbo"
    dCMap = "cool"
    wCMap = "bwr_r"
    vCMap = "gnuplot2"

    # Read the RCM results.
    rcmdata = gampp.GameraPipe(fdir, ftag)
    fnrcm = os.path.join(fdir, f'{ftag}.h5')

    # Get branch/hash info
    if do_hash:
        branch = kh5.GetBranch(fnrcm)
        githash = kh5.GetHash(fnrcm)
        if debug:
            print(f'branch/commit: {branch}/{githash}')

    if debug:
        print("rcmdata = %s" % rcmdata)

    # Set up Figure Parameters

    # Set global plot font options.
    mpl.rc('mathtext', fontset='stixsans', default='regular')
    mpl.rc('font', size=10)

    # Init figure.
    fig = plt.figure(figsize=figSz)

    if not do_vid: # Then we are making a single image, keep original functionality
        if nStp < 0:  # ANY negative index gets the last step.
            nStp = rcmdata.sFin
            print("Using Step %d"%(nStp))
        if debug:
            print("nStp = %s" % nStp)
        makePlot(nStp,rcmdata,nStp)

    else: # Then we make a video, i.e. series of images saved to rcmVid

        # Get video loop parameters
        s0 = max(rcmdata.s0,1) # Skip Step#0
        sFin = rcmdata.sFin
        nsteps = sFin - s0
        sIds = np.array(range(s0,sFin))
        outDir = 'rcmVid'
        kh5.CheckDirOrMake(outDir)

        # How many 0's do we need for filenames?
        n_pad = int(np.log10(nsteps)) + 1

        if ncpus == 1:
            for i, nStp in enumerate(sIds):
                makePlot(i,rcmdata, nStp)
        else:
            # Make list of parallel arguments
            ag = ((i,rcmdata,nStp) for i, nStp in enumerate(sIds) )
            # Check we're not exceeding cpu_count on computer
            ncpus = min(int(ncpus),cpu_count(logical=False))
            print('Doing multithreading on ',ncpus,' threads')
            # Do parallel job
            with Pool(processes=ncpus) as pl:
                pl.starmap(makePlot,ag)
            print("Done making all the images. Go to mixVid folder")

