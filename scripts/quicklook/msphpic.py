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
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import warnings
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from multiprocessing import Pool
from psutil import cpu_count

# Import project-specific modules.
from kaipy import cdaweb_utils
import kaipy.gamera.gampp as gampp
import kaipy.gamera.magsphere as msph
import kaipy.gamera.msphViz as mviz
import kaipy.gamera.rcmpp as rcmpp
import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv
import kaipy.kaiTools as ktools
import kaipy.kdefs as kdefs
import kaipy.remix.remix as remix


# Program constants and defaults

# Program description.
description = """Creates simple multi-panel figure for Gamera magnetosphere run
Top Panel - Residual vertical magnetic field
Bottom Panel - Pressure (or density) and hemispherical insets
NOTE: There is an optional -size argument for domain bounds options
(default: std), which is passed to kaiViz functions.
"""

# Default identifier for results to read.
default_runid = "msphere"

# Plot the last step by default.
default_step = -1

# Color to use for spacecraft position symbols.
SPACECRAFT_COLOR = 'red'


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
        "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-d", type=str, metavar="directory", default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "-id", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-n", type=int, metavar="step", default=default_step,
        help="Time slice to plot (default: %(default)s)"
    )
    parser.add_argument(
        "-bz", action="store_true", default=False,
        help="Show Bz instead of dBz (default: %(default)s)."
    )
    parser.add_argument(
        "-den", action="store_true", default=False,
        help="Show density instead of pressure (default: %(default)s)."
    )
    parser.add_argument(
        "-jy", action="store_true", default=False,
        help="Show Jy instead of pressure (default: %(default)s)."
    )
    parser.add_argument(
        "-ephi", action="store_true", default=False,
        help="Show Ephi instead of pressure (default: %(default)s)."
    )
    parser.add_argument(
        "-noion", action="store_true", default=False,
        help="Don't show ReMIX data (default: %(default)s)."
    )
    parser.add_argument(
        "-nompi", action="store_true", default=False,
        help="Don't show MPI boundaries (default: %(default)s)."
    )
    parser.add_argument(
        "-norcm", action="store_true", default=False,
        help="Don't show RCM data (default: %(default)s)."
    )
    parser.add_argument(
        "-bigrcm", action="store_true", default=False,
        help="Show entire RCM domain (default: %(default)s)."
    )
    parser.add_argument(
        "-src", action="store_true", default=False,
        help="Show source term (default: %(default)s)."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Names of spacecraft to plot positions, separated by commas (default: %(default)s)"
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
    # Add an option for plot domain size.
    mviz.AddSizeArgs(parser)
    return parser

def makePlot(i,spacecraft,nStp):

    # Disable some warning spam if not debug
    if not debug:
        warnings.filterwarnings("ignore", message="The input coordinates to pcolor are interpreted as cell centers.*")

    # Name of plot output file.
    if do_vid:
        fOut = "{}.{:0>{n}d}.png".format("msphpic", i, n=n_pad)
        outPath = os.path.join(outDir, fOut)
    else:
        # Name of plot output file.
        fOut = "qkmsphpic.png"
        outPath = fOut

    # Skip this file if it already exists and we're not supposed to overwrite
    if not do_overwrite and os.path.exists(outPath) and do_vid:
        return

    # Open remix data if available.
    if doMIX:
        print("Found ReMIX data")
        ion = remix.remix(rmxChk, nStp)
        if debug:
            print("ion = %s" % ion)

    gs = gridspec.GridSpec(3, 6, height_ratios=[20, 1, 1], hspace=0.025)
    if debug:
        print("fig = %s" % fig)
        print("gs = %s" % gs)

    # Create the plotting Axes objects.
    AxL = fig.add_subplot(gs[0, 0:3])
    AxR = fig.add_subplot(gs[0, 3:])
    AxC1 = fig.add_subplot(gs[-1, 0:2])
    AxC2 = fig.add_subplot(gs[-1, 2:4])
    AxC3 = fig.add_subplot(gs[-1, 4:6])
    if debug:
        print("AxL = %s" % AxL)
        print("AxR = %s" % AxR)
        print("AxC1 = %s" % AxC1)
        print("AxC2 = %s" % AxC2)
        print("AxC3 = %s" % AxC3)

    # Create the field-aligned current colorbar on Axes #2.
    cbM = kv.genCB(
        AxC2, kv.genNorm(remix.facMax), "FAC", cM=remix.facCM, Ntk=4
    )
    AxC2.xaxis.set_ticks_position('top')
    if debug:
        print("cbM = %s" % cbM)

    # On the left, plot the z-component of the residual magnetic field.
    Bz = mviz.PlotEqB(gsph, nStp, xyBds, AxL, AxC1, doBz=doBz)

    # Make any requested optional plots, or just pressure.
    if doJy:
        mviz.PlotJyXZ(gsph, nStp, xyBds, AxR, AxC3)
    elif doEphi:
        mviz.PlotEqEphi(gsph, nStp, xyBds, AxR, AxC3)
    else:
        mviz.PlotMerid(gsph, nStp, xyBds, AxR, doDen, doRCM, AxC3, doSrc=doSrc)

    # Add the date and time for the plot.
    gsph.AddTime(nStp, AxL, xy=[0.025, 0.89], fs="x-large")

    # Add the solar wind description text.
    gsph.AddSW(nStp, AxL, xy=[0.625, 0.025], fs="small")

    # If available, add the inset RCM plot.
    if not noRCM:
        AxRCM = inset_axes(AxL, width="30%", height="30%", loc=3)
        rcmpp.RCMInset(AxRCM, rcmdata, nStp, mviz.vP)
        # Add dBz contours.
        AxRCM.contour(
            kv.reWrap(gsph.xxc), kv.reWrap(gsph.yyc), kv.reWrap(Bz), [0.0],
            colors=mviz.bz0Col, linewidths=mviz.cLW
        )
        # Show the RCM region as a box.
        rcmpp.AddRCMBox(AxL)

    # Plot the REMIX data, if requested.
    if doIon:
        gsph.AddCPCP(nStp, AxR, xy=[0.610, 0.925])
        mviz.AddIonBoxes(gs[0, 3:], ion)

    # Show the MPI decomposition, if requested.
    if doMPI:
        mviz.PlotMPI(gsph, AxL)
        mviz.PlotMPI(gsph, AxR)

    # If requested, overlay the spacecraft locations.
    if spacecraft:
        print("Overplotting spacecraft trajectories of %s." % spacecraft)

        # Split the list into individual spacecraft names.
        spacecraft = spacecraft.split(',')
        if debug:
            print("spacecraft = %s" % spacecraft)

        # Fetch the MJD start and end time of the model results.
        fname = gsph.f0
        if debug:
            print("fname = %s" % fname)
        MJD_start = kh5.tStep(fname, gsph.s0, aID="MJD")
        if debug:
            print("MJD_start = %s" % MJD_start)
        MJD_end = kh5.tStep(fname, gsph.sFin, aID="MJD")
        if debug:
            print("MJD_end = %s" % MJD_end)

        # Convert the statrt and stop MJD to a datetime object in UT.
        ut_start = ktools.MJD2UT(MJD_start)
        if debug:
            print("ut_start = %s" % ut_start)
        ut_end = ktools.MJD2UT(MJD_end)
        if debug:
            print("ut_end = %s" % ut_end)

        # Fetch and plot the trajectory of each spacecraft from CDAWeb.
        for (i_sc, sc) in enumerate(spacecraft):

            # Fetch the spacecraft trajectory in Solar Magnetic (SM)
            # Cartesian coordinates between the start and end times.
            sc_x, sc_y, sc_z = cdaweb_utils.fetch_spacecraft_SM_trajectory(
                sc, ut_start, ut_end
            )
            if debug:
                print("sc_x, sc_y, sc_z = %s, %s, %s" % (sc_x, sc_y, sc_z))

            # Skip if no trajectory found.
            if sc_x is None:
                print("No trajectory found for spacecraft %s." % sc)
                continue

            # Convert coordinates to units of Earth radius.
            CM_TO_KM = 1e-5  # Centimeters to kilometers
            Re_km = kdefs.Re_cgs*CM_TO_KM  # Earth radius in kilometers
            sc_x_Re = sc_x/Re_km
            sc_y_Re = sc_y/Re_km
            sc_z_Re = sc_z/Re_km
            if debug:
                print("sc_x_Re, sc_y_Re, sc_z_Re = %s, %s, %s" %
                (sc_x_Re, sc_y_Re, sc_z_Re))

            # Plot a labelled trajectory of the spacecraft. Also plot a larger
            # dot at the last point in the trajectory.
            # Left plot
            SPACECRAFT_COLORS = list(mpl.colors.TABLEAU_COLORS.keys())
            color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
            AxL.plot(sc_x_Re, sc_y_Re, marker=None, linewidth=1, c=color)
            AxL.plot(sc_x_Re[-1], sc_y_Re[-1], 'o', c=color)
            x_nudge = 1.0
            y_nudge = 1.0
            AxL.text(sc_x_Re[-1] + x_nudge, sc_y_Re[-1] + y_nudge, sc, c=color)
            # Right plot
            AxR.plot(sc_x_Re, sc_z_Re, marker=None, linewidth=1, c=color)
            AxR.plot(sc_x_Re[-1], sc_z_Re[-1], 'o', c=color)
            x_nudge = 1.0
            z_nudge = 1.0
            AxR.text(sc_x_Re[-1] + x_nudge, sc_z_Re[-1] + z_nudge, sc, c=color)

    # Add Branch and Hash info
    if do_hash:
        fig.text(0.1,0.87,f"branch/commit: {branch}/{githash}", fontsize=6)

    # Save the plot to a file.
    kv.savePic(outPath, bLenX=45)


if __name__ == "__main__":
    """Make a quick figure of a Gamera magnetosphere run."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    verbose = args.verbose
    fdir = args.d
    ftag = args.id
    nStp = args.n
    doDen = args.den
    noIon = args.noion
    noMPI = args.nompi
    doMPI = not noMPI
    doJy = args.jy
    doEphi = args.ephi
    doSrc = args.src
    doBz = args.bz
    noRCM = args.norcm
    doBigRCM = args.bigrcm
    do_vid = args.vid
    do_overwrite = args.overwrite
    do_hash = not args.nohash
    ncpus = args.ncpus
    spacecraft = args.spacecraft
    if debug:
        print("args = %s" % args)

    # Get the domain size in Re.
    xyBds = mviz.GetSizeBds(args)
    if debug:
        print("xyBds = %s" % xyBds)

    # Set figure parameters.
    doFast = False
    doIon = not noIon
    figSz = (12, 7.5)

    # Open the gamera results pipe.
    gsph = msph.GamsphPipe(fdir, ftag, doFast=doFast)

    # Check for the presence of RCM results.
    rcmChk = os.path.join(fdir, "%s.mhdrcm.h5" % ftag)
    doRCM = os.path.exists(rcmChk)
    if debug:
        print("rcmChk = %s" % rcmChk)
        print("doRCM = %s" % doRCM)

    # Check for the presence of remix results.
    rmxChk = os.path.join(fdir, "%s.mix.h5" % ftag)
    doMIX = os.path.exists(rmxChk)
    if debug:
        print("rmxChk = %s" % rmxChk)
        print("doMIX = %s" % doMIX)

    # Get branch/hash info
    if doMIX:
        branch = kh5.GetBranch(rmxChk)
        githash = kh5.GetHash(rmxChk)
        if debug:
            print(f'branch/commit: {branch}/{githash}')

    # Open RCM data if available, and initialize visualization.
    if doRCM:
        print("Found RCM data")
        rcmdata = gampp.GameraPipe(fdir, ftag + ".mhdrcm")
        mviz.vP = kv.genNorm(1.0e-2, 100.0, doLog=True)
        rcmpp.doEll = not doBigRCM
        if debug:
            print("rcmdata = %s" % rcmdata)
    else:
        rcmdata = None

    # Setup the figure.
    mpl.use('Agg')  # Plot in memory buffer.

    # Set global plot font options.
    mpl.rc('mathtext', fontset='stixsans', default='regular')
    mpl.rc('font', size=10)

    # Init figure
    fig = plt.figure(figsize=figSz)

    if not do_vid: # If we are making a single image, keep original functionality
        # If needed, fetch the number of the last step.
        if nStp < 0:
            nStp = gsph.sFin
            print("Using Step %d" % nStp)
        makePlot(nStp,spacecraft,nStp)

    else: # then we make a video, i.e. series of images saved to msphVid

        # Get video loop parameters
        s0 = max(gsph.s0,1) # Skip Step#0
        sFin = gsph.sFin
        nsteps = sFin - s0
        sIds = np.array(range(s0,sFin))
        outDir = 'msphVid'
        kh5.CheckDirOrMake(outDir)

        # How many 0's do we need for filenames?
        n_pad = int(np.log10(nsteps)) + 1

        if ncpus == 1:
            for i, nStp in enumerate(sIds):
                makePlot(i, spacecraft, nStp)
        else:
            # Make list of parallel arguments
            ag = ((i,spacecraft,nStp) for i, nStp in enumerate(sIds) )
            # Check we're not exceeding cpu_count on computer
            ncpus = min(int(ncpus),cpu_count(logical=False))
            print('Doing multithreading on ',ncpus,' threads')
            # Do parallel job
            with Pool(processes=ncpus) as pl:
                pl.starmap(makePlot,ag)
            print("Done making all the images. Go to mixVid folder")
