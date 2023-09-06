#!/usr/bin/env python


"""Make a quick-look figure of a Gamera heliosphere run.

Make a quick-look figure of a Gamera heliosphere run.

Five different sets of plots are supported, and are distinguished by the
value of the "pic" argument.

pic1 (default): A 4-panel display showing pcolormesh plots in the z = 0
(equatorial) plane of the gamhelio frame used in the simulation. The plots
are:

    Upper left: Solar wind speed (km/s)
    Upper right: Solar wind number density scaled by (r/r0)**2 (cm**-3)
    Lower left: Solar wind temperature scaled by r/r0 (MK)
    Lower right: Solar wind radial magnetic field scaled by r/r0 (nT)

pic2: A 4-panel display showing pcolormesh plots in the y = 0 (meridional,
containing Earth) plane of the gamhelio frame used in the simulation. The
plots are:

    Upper left: Solar wind speed (km/s)
    Upper right: Solar wind number density scaled by (r/r0)**2 (cm**-3)
    Lower left: Solar wind temperature scaled by r/r0 (MK)
    Lower right: Solar wind radial magnetic field scaled bryr/r0 (nT)

pic3: A 4-panel display showing pcolormesh plots in the r = 1 AU slice of the
gamhelio frame used in the simulation. The plots are:

    Upper left: Solar wind speed (km/s)
    Upper right: Solar wind number density (cm**-3)
    Lower left: Solar wind temperature (MK)
    Lower right: Solar wind radial magnetic field (nT)

pic4: A pcolormesh plot in the innermost radial slice (r = 22 Rsun) of the
gamhelio frame used in the simulation. The plot shows the radial magnetic
field in nT, in a coordinate frame rotating with the Sun.

pic5: A 3-panel display showing line as a function of radius,
22 Rsun <= r <= 220 Rsun. The plots are:

    Upper left: Solar wind number density (cm**-3)
    Upper right: Solar wind speed (km/s)
    Lower left: Solar wind radial momentum flux (km**2/s**2/cm**3)

All plots can optionally display the contemporary location of relevant
spacecraft.

Author
------
Elena Provornikova (elena.provornikova@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import os

# Import supplemental modules.
import astropy.time
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
from kaipy.satcomp import scutils


# Program constants and defaults

# Program description.
description = "Create quick-look plots for a Gamera heliosphere run"

# Default identifier for results to read.
default_runid = "wsa"

# Plot the last step by default.
default_step = -1

# Code for default picture type.
default_pictype = "pic1"

# Name of plot output file.
fOut = "qkpic.png"


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
        "-p", type=str, metavar="pictype", default=default_pictype,
        help="Code for plot type (default: %(default)s)"
    )
    parser.add_argument(
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Names of spacecraft to plot positions, separated by commas"
        " (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


if __name__ == "__main__":
    """Make a quick figure of a Gamera heliosphere run."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    fdir = args.d
    ftag = args.id
    nStp = args.n
    pic = args.p
    spacecraft = args.spacecraft
    verbose = args.verbose
    if debug:
        print("args = %s" % args)

    # Fetch the plot domain based on the picture type.
    xyBds = hviz.GetSizeBds(pic)
    print(xyBds)

    # Do work?
    doFast = False

    # Create figures in a memory buffer.
    mpl.use("Agg")

    # Read the CDAWeb spacecraft database.
    sc_metadata_path = os.path.join(
        os.environ["KAIJUHOME"], "kaipy", "satcomp", "sc_helio.json"
    )
    sc_metadata = scutils.getScIds(spacecraft_data_file=sc_metadata_path)

    # Determine figure size (width, height) (inches) based on the picture type.
    if pic == "pic1" or pic == "pic2":
        figSz = (10, 12.5)
    elif pic == "pic3":
        figSz = (10, 6.5)
    elif pic == "pic4":
        figSz = (10, 6)
    elif pic == "pic5":
        figSz = (12, 12)

    # Create the figure.
    fig = plt.figure(figsize=figSz)

    # Lay out the subplots.
    if pic == "pic1" or pic == "pic2" or pic == "pic3":
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

    # Open a pipe to the results data.
    gsph = hsph.GamsphPipe(fdir, ftag, doFast=doFast)
    if nStp < 0:
        nStp = gsph.sFin
        print("Using Step %d" % nStp)

    # Extract the date/time of the plot.
    mjd = gsph.MJDs[nStp]
    if debug:
        print(f"mjd = {mjd}")

    # Now create the actual plots.
    if pic == "pic1":
        # These are all equatorial plots in the XY plane of the HGS frame
        # used by gamhelio.
        hviz.PlotEqMagV(gsph, nStp, xyBds, AxL0, AxC1_0)
        hviz.PlotEqD(gsph, nStp, xyBds, AxR0, AxC2_0)
        hviz.PlotEqTemp(gsph, nStp, xyBds, AxL1, AxC1_1)
        hviz.PlotEqBr(gsph, nStp, xyBds, AxR1, AxC2_1)
    elif pic == "pic2":
        hviz.PlotMerMagV(gsph, nStp, xyBds, AxL0, AxC1_0)
        hviz.PlotMerDNorm(gsph, nStp, xyBds, AxR0, AxC2_0)
        hviz.PlotMerTemp(gsph, nStp, xyBds, AxL1, AxC1_1)
        hviz.PlotMerBrNorm(gsph, nStp, xyBds, AxR1, AxC2_1)
    elif pic == "pic3":
        hviz.PlotiSlMagV(gsph, nStp, xyBds, AxL0, AxC1_0)
        hviz.PlotiSlD(gsph, nStp, xyBds, AxR0, AxC2_0)
        hviz.PlotiSlTemp(gsph, nStp, xyBds, AxL1, AxC1_1)
        hviz.PlotiSlBr(gsph, nStp, xyBds, AxR1, AxC2_1)
    elif pic == "pic4":
        hviz.PlotiSlBrRotatingFrame(gsph, nStp, xyBds, Ax, AxC)
    elif pic == "pic5":
        hviz.PlotDensityProf(gsph, nStp, xyBds, Ax)
        hviz.PlotSpeedProf(gsph, nStp, xyBds, AxC)
        hviz.PlotFluxProf(gsph, nStp, xyBds, AxC1)
    else:
        print("Pic is empty. Choose pic1 or pic2 or pic3")

    # Add time in the upper left.
    if pic == "pic1" or pic == "pic2":
        gsph.AddTime(nStp, AxL0, xy=[0.025, 0.875], fs="x-large")
    elif pic == "pic3":
        gsph.AddTime(nStp, AxL0, xy=[0.015, 0.82], fs="small")
    elif pic == "pic4" or pic == "pic5":
        gsph.AddTime(nStp, Ax, xy=[0.015, 0.92], fs="small")
    else:
        print("Pic is empty. Choose pic1 or pic2 or pic3")

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
                print("Converting ephemeris for %s into gamhelio format." %
                      sc_id)
            x, y, z = cdaweb_utils.ingest_helio_spacecraft_trajectory(
                sc_id, sc_data, MJDc
            )
            if debug:
                print("x, y, z = %s, %s, %s" % (x, y, z))

            # Convert the datetime objects from the trajectory to MJD.
            t_strings = np.array([str(t) for t in sc_data["Epoch"]])
            t = astropy.time.Time(t_strings, scale='utc').mjd

            # Interpolate the spacecraft position at the time for the plot.
            t_sc = mjd
            x_sc = np.interp(t_sc, t, x)
            y_sc = np.interp(t_sc, t, y)
            z_sc = np.interp(t_sc, t, z)

            # If needed, compute heliocentric spherical coordinates.
            if pic == "pic3" or pic == "pic4":
                rxy = np.sqrt(x**2 + y**2)
                theta = np.arctan2(rxy, z)
                phi = np.arctan2(y, x)
                lat = np.degrees(np.pi/2 - theta)
                lon = np.degrees(phi)
                lat_sc = np.interp(t_sc, t, lat)
                lon_sc = np.interp(t_sc, t, lon)

            # Plot a labelled trajectory of the spacecraft. Also plot a larger
            # dot at the last point in the trajectory.
            # Left plot
            SPACECRAFT_COLORS = list(mpl.colors.TABLEAU_COLORS.keys())
            color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
            x_nudge = 0.0
            y_nudge = 8.0
            sc_label = sc_metadata[sc_id]["label"]
            if pic == "pic1":
                for ax in (AxL0, AxR0, AxL1, AxR1):
                    ax.plot(x_sc, y_sc, 'o', c=color)
                    ax.plot(x_sc, y_sc, 'o', c="black", fillstyle="none")
                    ax.text(x_sc + x_nudge, y_sc + y_nudge, sc_label,
                            c="black", horizontalalignment="center")
            elif pic == "pic2":
                for ax in (AxL0, AxR0, AxL1, AxR1):
                    ax.plot(x_sc, z_sc, 'o', c=color)
                    ax.plot(x_sc, z_sc, 'o', c="black", fillstyle="none")
                    ax.text(x_sc + x_nudge, z_sc + y_nudge, sc_label,
                            c="black", horizontalalignment="center")
            elif pic == "pic3":
                for ax in (AxL0, AxR0, AxL1, AxR1):
                    ax.plot(lon_sc, lat_sc, 'o', c=color)
                    ax.plot(lon_sc, lat_sc, 'o', c="black", fillstyle="none")
                    ax.text(lon_sc + x_nudge, lat_sc + y_nudge, sc_label,
                            c="black", horizontalalignment="center")
            elif pic == "pic4":
                ax = Ax
                ax.plot(lon_sc, lat_sc, 'o', c=color)
                ax.plot(lon_sc, lat_sc, 'o', c="black", fillstyle="none")
                ax.text(lon_sc + x_nudge, lat_sc + y_nudge, sc_label,
                        c="black", horizontalalignment="center")
            elif pic == "pic5":
                pass

    # Save the figure to a file.
    path = os.path.join(fdir, fOut)
    kv.savePic(path, bLenX=45)
