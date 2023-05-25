#!/usr/bin/env python


"""Make a quick figure of a Gamera heliosphere run.

Make a quick figure of a Gamera heliosphere run.

Five different sets of plots are supported, and are distinguished by the
value of the "pic" argument.

pic1 (default): A 4-panel display showing radial speed, number
density*(r/r0)**2, temperature*(r/r0), and radial magnetic field*(r/r0)**2.
These plots are done in the XY plane of the gamhelio frame (which is a
Heliographic Stonyhurst (HGS) frame, modified with the +x reversed from
the usual HGS definition)

pic2:

pic3:

pic4:

pic5:

Author
------
Elena Provornikova (elena.provornikova@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import os
import time

# Import supplemental modules.
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np

# Import project-specific modules.
from kaipy import cdaweb_utils
import kaipy.gamhelio.helioViz as hviz
import kaipy.gamhelio.heliosphere as hsph
import kaipy.kaiH5 as kh5
import kaipy.kaiTools as ktools
import kaipy.kaiViz as kv
import kaipy.kdefs as kdefs
from kaipy.satcomp import scutils


# Program constants and defaults

# Program description.
description = """Creates multi-panel figure for Gamera heliosphere run
Upper left - Solar wind speed
Upper right - Solar wind number density
Lower left - Solar wind temperature
Lower right - Solar wind radial magnetic field
"""

# Default identifier for results to read.
default_runid = "wsa"

# List of steps
default_steps = "-1"

# Default slices
default_slice = "1:2:1"

# Code for default picture type.
default_pictype = "pic1"

# Color to use for spacecraft position symbols.
SPACECRAFT_COLOR = "red"

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
        "--nlist", type=lambda n: [int(item) for item in n.split(',')], metavar="list of steps", default=default_steps,
        help="List of time slice(s) to plot (default: %(default)s)"
    )
    parser.add_argument(
        "--nslice", type=lambda n: [int(item) for item in n.split(':')], metavar="step slice", default=default_slice,
        help="Slice for range of time slice(s) to plot (default: %(default)s)"
    )
    parser.add_argument(
        "-nompi", action="store_true", default=False,
        help="Don't show MPI boundaries (default: %(default)s)."
    )
    parser.add_argument(
        "-p", type=str, metavar="pictype", default=default_pictype,
        help="Code for plot type (default: %(default)s)"
    )
    parser.add_argument(
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Names of spacecraft to plot positions, separated by commas (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        "-o", "--outfile", type=str, metavar="outFile", default=fOut,
        help="Output file name (default: %(default)s)."
    )
    return parser

def initFig(pic):
    # Determine figure size (width, height) (inches) based on the picture type.
    if pic == "pic1" or pic == "pic2":
        figSz = (10, 12.5)
    elif pic == "pic3":
        figSz = (10, 6.5)
    elif pic == "pic4":
        figSz = (10, 6)
    elif pic == "pic5":
        figSz = (12, 12)
    elif  pic == "pic6":
        figSz = (12.5, 12.5)

    # Create the figure.
    fig = plt.figure(figsize=figSz)
    return fig

# Name of plot output file.
def fOut(id, pic, nStp):
    return "qkpic_{}_{}_n{}.png".format(id, pic, nStp)

if __name__ == "__main__":
    from alive_progress import alive_bar

    """Make a quick figure of a Gamera heliosphere run."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    fdir = args.d
    ftag = args.id
    noMPI = args.nompi
    steps = args.nlist
    slices = args.nslice
    pic = args.p
    spacecraft = args.spacecraft
    verbose = args.verbose
    fOut = args.outfile
    if debug:
        print("args = %s" % args)
    if slices:
        print("Slice selected {}".format(slice(slices[0],slices[1],slices[2])))
    # Invert the MPI flag for convenience.
    doMPI = not noMPI

    tic = time.perf_counter()
    # Fetch the plot domain based on the picture type.
    xyBds = hviz.GetSizeBds(pic)
    toc = time.perf_counter()
    print(xyBds)
    print(f"Get bounds took {toc-tic} s")
    # Do work?ß
    doFast = False

    # Create figures in a memory buffer.
    mpl.use("Agg")

    # Open a pipe to the results data.
    tic = time.perf_counter()
    gsph = hsph.GamsphPipe(fdir, ftag, doFast=doFast)
    toc = time.perf_counter()

    print(f"Open pipe took {toc-tic} s")

    if(slices and steps[0] == -1):
        steps = range(gsph.sFin)[slice(slices[0],slices[1],slices[2])]

    nsteps = len(steps)

    for nStp in steps: 
        tic = time.perf_counter()
        print(f"Generating {pic} for time step {nStp}")
        fig = initFig(pic)

        # Lay out the subplots.
        if pic == "pic1" or pic == "pic2" or pic == "pic3" or pic == "pic6":
            gs = gridspec.GridSpec(4, 6, height_ratios=[20, 1, 20, 1])
            # Axes for plots.
            AxL0 = fig.add_subplot(gs[0, 0:3])
            AxR0 = fig.add_subplot(gs[0, 3:])
            AxL1 = fig.add_subplot(gs[2, 0:3])
            AxR1 = fig.add_subplot(gs[2, 3:])
            # Axes for colorbars.
            AxC1_0 = fig.add_subplot(gs[1, 0:3])
            AxC2_0 = fig.add_subplot(gs[1, 3:])
            AxC1_1 = fig.add_subplot(gs[3, 0:3])
            AxC2_1 = fig.add_subplot(gs[3, 3:])
        elif pic == "pic4":
            gs = gridspec.GridSpec(2, 1, height_ratios=[20, 1])
            Ax = fig.add_subplot(gs[0, 0])
            AxC = fig.add_subplot(gs[1, 0])
        elif pic == "pic5":
            gs = gridspec.GridSpec(2, 2)
            Ax = fig.add_subplot(gs[0, 0])
            AxC = fig.add_subplot(gs[0, 1])
            AxC1 = fig.add_subplot(gs[1, 0])

        if nStp < 0:
            nStp = gsph.sFin
            print("Using Step %d" % nStp)

        # Now create the actual plots.
        if pic == "pic1":
            # These are all equatorial plots in the XY plane of the HGS frame
            # used by gamhelio.
            hviz.PlotEqMagV(gsph, nStp, xyBds, AxL0, AxC1_0)
            hviz.PlotEqD(gsph, nStp, xyBds, AxR0, AxC2_0)
            hviz.PlotEqTemp(gsph, nStp, xyBds, AxL1, AxC1_1)
            hviz.PlotEqBr(gsph, nStp, xyBds, AxR1, AxC2_1)
        elif pic == "pic2":
            hviz.PlotMerMagV(gsph ,nStp, xyBds, AxL0, AxC1_0,indx=(None,2.95))
            hviz.PlotMerDNorm(gsph, nStp, xyBds, AxR0, AxC2_0,indx=(None,2.95))
            hviz.PlotMerTemp(gsph, nStp, xyBds, AxL1, AxC1_1,indx=(None,2.95))
            hviz.PlotMerBrNorm(gsph, nStp, xyBds, AxR1, AxC2_1,indx=(None,2.95))
        elif pic == "pic3":
            hviz.PlotiSlMagV(gsph, nStp, xyBds, AxL0, AxC1_0,idx=0)
            hviz.PlotiSlD(gsph, nStp, xyBds, AxR0, AxC2_0,idx=0)
            hviz.PlotiSlTemp(gsph, nStp, xyBds, AxL1, AxC1_1,idx=0)
            hviz.PlotiSlBr(gsph, nStp, xyBds, AxR1, AxC2_1,idx=0)
        elif pic == "pic4":
            hviz.PlotiSlBrRotatingFrame(gsph, nStp, xyBds, Ax, AxC)
        elif pic == "pic5":
            hviz.PlotDensityProf(gsph, nStp, xyBds, Ax)
            hviz.PlotSpeedProf(gsph, nStp, xyBds, AxC)
            hviz.PlotFluxProf(gsph, nStp, xyBds, AxC1)
        elif pic == "pic6":
            hviz.PlotEqBr(gsph, nStp, xyBds, AxL0, AxC1_0)
            hviz.PlotEqBx(gsph, nStp, xyBds, AxR0, AxC2_0)
            hviz.PlotEqBy(gsph, nStp, xyBds, AxL1, AxC1_1)
            hviz.PlotEqBz(gsph, nStp, xyBds, AxR1, AxC2_1)
        else:
            print ("Pic is empty. Choose pic1 or pic2 or pic3")

        # Add time in the upper left.
        if pic == "pic1" or pic == "pic2" or pic == "pic6":
            gsph.AddTime(nStp, AxL0, xy=[0.025, 0.875], fs="x-large")
        elif pic == "pic3":
            gsph.AddTime(nStp, AxL0, xy=[0.015, 0.82], fs="small")
        elif pic == "pic4" or pic == "pic5":
            gsph.AddTime(nStp, Ax, xy=[0.015, 0.92], fs="small")
        else:
            print ("Pic is empty. Choose pic1 or pic2 or pic3")

        # Overlay the spacecraft trajectory, if needed.
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
            MJD_start = kh5.tStep(fname, 0, aID="MJD")
            if debug:
                print("MJD_start = %s" % MJD_start)
            MJD_end = kh5.tStep(fname, gsph.sFin, aID="MJD")
            if debug:
                print("MJD_end = %s" % MJD_end)

            # Convert the start and stop MJD to a datetime object in UT.
            ut_start = ktools.MJD2UT(MJD_start)
            if debug:
                print("ut_start = %s" % ut_start)
            ut_end = ktools.MJD2UT(MJD_end)
            if debug:
                print("ut_end = %s" % ut_end)

            # Get the MJDc value for use in computing the gamhelio frame.
            MJDc = scutils.read_MJDc(fname)
            if debug:
                print("mjdc = %s" % MJDc)

            # Fetch and plot the trajectory of each spacecraft from CDAWeb.
            for (i_sc, sc_id) in enumerate(spacecraft):
                if verbose:
                    print("Fetching trajectory for %s." % sc_id)

                # Fetch the spacecraft trajectory in whatever frame is available
                # from CDAWeb.
                sc_data = cdaweb_utils.fetch_helio_spacecraft_trajectory(
                    sc_id, ut_start, ut_end
                )
                if sc_data is None:
                    print("No trajectory found for %s." % sc_id)
                    continue

                # Ingest the trajectory by converting it to the GH(MJDc) frame.
                if verbose:
                    print("Converting ephemeris for %s into gamhelio format." % sc_id)
                x, y, z = cdaweb_utils.ingest_helio_spacecraft_trajectory(sc_id, sc_data, MJDc)
                if debug:
                    print("x, y, z = %s, %s, %s" % (x, y, z))

                # If needed, compute heliocentric spherical coordinates.
                if pic == "pic3" or pic == "pic4":
                    r = np.sqrt(x**2 + y**2 + z**2)
                    lon = np.degrees(np.arccos(x/(x**2 + y**2)))
                    lat = np.degrees(-np.arccos(z/r) + np.pi/2)

                # Plot a labelled trajectory of the spacecraft. Also plot a larger
                # dot at the last point in the trajectory.
                # Left plot
                SPACECRAFT_COLORS = list(mpl.colors.TABLEAU_COLORS.keys())
                color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
                x_nudge = 5.0
                y_nudge = 5.0
                if pic == "pic1":
                    for ax in (AxL0, AxR0, AxL1, AxR1):
                        ax.plot(x, y, marker=None, linewidth=1, c=color)
                        ax.plot(x[-1], y[-1], 'o', c=color)
                        ax.text(x[-1] + x_nudge, y[-1] + y_nudge, sc_id, c=color)
                elif pic == "pic2":
                    for ax in (AxL0, AxR0, AxL1, AxR1):
                        ax.plot(x, z, marker=None, linewidth=1, c=color)
                        ax.plot(x[-1], z[-1], 'o', c=color)
                        ax.text(x[-1] + x_nudge, z[-1] + y_nudge, sc_id, c=color)
                elif pic == "pic3":
                    for ax in (AxL0, AxR0, AxL1, AxR1):
                        ax.plot(lon, lat, marker=None, linewidth=1, c=color)
                        ax.plot(lon[-1], lat[-1], 'o', c=color)
                        ax.text(lon[-1] + x_nudge, lat[-1] + y_nudge, sc_id, c=color)
                elif pic == "pic4":
                    ax = Ax
                    ax.plot(lon, lat, marker=None, linewidth=1, c=color)
                    ax.plot(lon[-1], lat[-1], 'o', c=color)
                    ax.text(lon[-1] + x_nudge, lat[-1] + y_nudge, sc_id, c=color)
                elif pic == "pic5":
                    pass

        # Save the figure to a file.
        path = os.path.join(fdir, fOut(ftag, pic, nStp))
        kv.savePic(path, bLenX=40)
        fig.clear()
        toc = time.perf_counter()
        print(f"Step {nStp} took {toc-tic} s")
