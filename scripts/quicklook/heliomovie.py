#!/usr/bin/env python


"""Make a movie from a gamhelio run.

Make a movie from a gamhelio run. The frames in the movie are based on the
individual plots made by heliopic.py.

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
import subprocess

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
description = "Make a movie from a gamhelio run"

# Default identifier for results to read.
default_runid = "wsa"

# Plot all steps by default.
default_first_step = 1
default_last_step = -1

# Code for default picture type.
default_pictype = "pic1"

# Valid plot type strings.
valid_pictypes = ("pic1", "pic2", "pic3", "pic4", "pic5")

# Path to spacecraft metadata file.
sc_metadata_path = os.path.join(
    os.environ["KAIJUHOME"], "kaipy", "satcomp", "sc_helio.json"
)

# Figure sizes by plot type, (width, height) in inches.
figure_sizes = {
    "pic1": (10, 12.5),
    "pic2": (10, 12.5),
    "pic3": (10, 6.5),
    "pic4": (10, 6),
    "pic5": (12, 12),
}

# List of colors to use for spacecraft position dots.
SPACECRAFT_COLORS = list(mpl.colors.TABLEAU_COLORS.keys())


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

    Raises
    ------
    None
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
        "-d", "--directory", type=str, metavar="directory",
        default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "-id", "--runid", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-n0", "--first_step", type=int, metavar="n0",
        default=default_first_step,
        help="First time step to plot (default: %(default)s)"
    )
    parser.add_argument(
        "-n1", "--last_step", type=int, metavar="n1",
        default=default_last_step,
        help="Last time step to plot (default: %(default)s)"
    )
    parser.add_argument(
        "-p", "--pictype", type=str, metavar="pictype",
        default=default_pictype,
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


def fetch_spacecraft_trajectories(spacecraft_names, gsph):
    """Fetch spacecraft trajectories from CDAWeb.

    Fetch spacecraft trajectories from CDAWeb.

    Parameters
    ----------
    spacecraft_names : list of str
        List of spacecraft names.
    gsph : GamsphPipe
        Pipe object for reading simulation results.

    Returns
    -------
    sc_t, sc_x, sc_y, sc_z : dict of np.ndarray
        Dictionaries of (t, x, y, z), keyed by spacecraft names.

    Raises
    ------
    None
    """
    # Initialize the coordinate dictionaries.
    sc_t = {}
    sc_x = {}
    sc_y = {}
    sc_z = {}

    # Fetch the MJD start and end time of the model results.
    fname = gsph.f0
    MJD_start = kh5.tStep(fname, 0, aID="MJD")
    MJD_end = kh5.tStep(fname, gsph.sFin, aID="MJD")

    # Convert the start and stop MJD to a datetime object in UT.
    ut_start = ktools.MJD2UT(MJD_start)
    ut_end = ktools.MJD2UT(MJD_end)

    # Get the MJDc value for use in computing the gamhelio frame.
    MJDc = scutils.read_MJDc(fname)

    # Fetch the trajectory of each spacecraft from CDAWeb.
    for (i_sc, sc_id) in enumerate(spacecraft_names):

        # Fetch the spacecraft trajectory in whatever frame is available
        # from CDAWeb.
        sc_data = cdaweb_utils.fetch_helio_spacecraft_trajectory(
            sc_id, ut_start, ut_end
        )
        if sc_data is None:
            print(f"No trajectory found for {sc_id}.")
            continue

        # Ingest the trajectory by converting it to the GH(MJDc) frame.
        x, y, z = cdaweb_utils.ingest_helio_spacecraft_trajectory(
            sc_id, sc_data, MJDc
        )

        # Convert the datetime objects from the trajectory to MJD.
        t_strings = np.array([str(t) for t in sc_data["Epoch"]])
        t = astropy.time.Time(t_strings, scale='utc').mjd

        # Save the trajectory for this spacecraft.
        sc_t[sc_id] = t
        sc_x[sc_id] = x
        sc_y[sc_id] = y
        sc_z[sc_id] = z

    # Return the trajectory dictionaries.
    return sc_t, sc_x, sc_y, sc_z


def create_pic1_movie(args):
    """Create a pic1-style gamhelio movie.

    Create a pic1-style gamhelio movie.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options.

    Returns
    -------
    movie_file : str
        Path to movie file.

    Raises
    ------
    None
    """
    # Extract arguments.
    debug = args.debug
    pictype = args.pictype
    spacecraft = args.spacecraft
    verbose = args.verbose

    # Fetch the plot limits based on the picture type.
    plot_limits = hviz.GetSizeBds(pictype)
    if debug:
        print(f"plot_limits = {plot_limits}")

    # Create all plot images in a memory buffer.
    mpl.use("Agg")

    # Fetch the figure size.
    figsize = figure_sizes[pictype]
    if debug:
        print(f"figsize = {figsize}")

    # Create figures in a memory buffer.
    mpl.use("Agg")

    # Create the figure.
    fig = plt.figure(figsize=figsize)
    if debug:
        print(f"fig = {fig}")

    # Lay out the subplots for this figure. The grid contains separate axes
    # to use for the color bars (the rows with relative heights of 1).
    nrows = 4
    ncols = 6
    gs = mpl.gridspec.GridSpec(nrows, ncols, height_ratios=[20, 1, 20, 1])
    if debug:
        print(f"gs = {gs}")

    # Create the Axes objects for the individual subplots.
    # Each subplot is 1 row x 3 columns in the grid.
    ax_v = fig.add_subplot(gs[0, :3])   # Upper left
    ax_n = fig.add_subplot(gs[0, 3:])   # Upper right
    ax_T = fig.add_subplot(gs[2, :3])   # Lower left
    ax_Br = fig.add_subplot(gs[2, 3:])  # Lower right
    if debug:
        print(f"ax_v = {ax_v}")
        print(f"ax_n = {ax_n}")
        print(f"ax_T = {ax_T}")
        print(f"ax_Br = {ax_Br}")

    # Create the Axes objects for the individual color bars.
    # Each color bar is 1 (thin) row x 3 columns in the grid.
    ax_cb_v = fig.add_subplot(gs[1, :3])   # Upper left
    ax_cb_n = fig.add_subplot(gs[1, 3:])   # Upper right
    ax_cb_T = fig.add_subplot(gs[3, :3])   # Lower left
    ax_cb_Br = fig.add_subplot(gs[3, 3:])  # Lower right
    if debug:
        print(f"ax_cb_v = {ax_cb_v}")
        print(f"ax_cb_n = {ax_cb_n}")
        print(f"ax_cb_T = {ax_cb_T}")
        print(f"ax_cb_Br = {ax_cb_Br}")

    # Open a "pipe" to the data for this run.
    fdir = args.directory
    ftag = args.runid
    gsph = hsph.GamsphPipe(fdir, ftag)

    # If spacecraft positions will be plotted, read the spacecraft metadata
    # and fetch trajectories.
    if spacecraft:
        sc_metadata = scutils.getScIds(spacecraft_data_file=sc_metadata_path)
        if debug:
            print(f"scPmetadata = {sc_metadata}")

        # Split the list into individual spacecraft names.
        spacecraft_names = spacecraft.split(',')
        if debug:
            print(f"spacecraft_names = {spacecraft_names}")

        # Fetch the spacecraft trajectories.
        sc_t, sc_x, sc_y, sc_z = fetch_spacecraft_trajectories(
            spacecraft_names, gsph
        )

    # Create and save frame images for each step.
    first_step = args.first_step
    last_step = args.last_step
    if last_step == -1:
        last_step = gsph.Nt
    else:
        last_step += 1
    if debug:
        print(f"first_step, last_step = {first_step, last_step}")
    frame_files = []
    for i_step in range(first_step, last_step):
        if verbose:
            print(f"Creating {pictype} frame for step {i_step}.")

        # Extract the MJD for the frame.
        mjd = gsph.MJDs[i_step]
        if debug:
            print(f"mjd = {mjd}")

        # Create the individual plots for this frame.
        hviz.PlotEqMagV(gsph, i_step, plot_limits, ax_v, ax_cb_v)
        hviz.PlotEqD(gsph, i_step, plot_limits, ax_n, ax_cb_n)
        hviz.PlotEqTemp(gsph, i_step, plot_limits, ax_T, ax_cb_T)
        hviz.PlotEqBr(gsph, i_step, plot_limits, ax_Br, ax_cb_Br)

        # Add time in the upper left.
        gsph.AddTime(i_step, ax_v, xy=[0.025, 0.875], fs="x-large")

        # Overlay spacecraft positions (optional).
        if spacecraft:
            for (i_sc, sc_id) in enumerate(spacecraft_names):

                # Skip this spacecraft if no trajectory available.
                if sc_id not in sc_t:
                    continue

                # Interpolate the spacecraft position at the time for the plot.
                t_sc = mjd
                x_sc = np.interp(t_sc, sc_t[sc_id], sc_x[sc_id])
                y_sc = np.interp(t_sc, sc_t[sc_id], sc_y[sc_id])
                z_sc = np.interp(t_sc, sc_t[sc_id], sc_z[sc_id])

                # Plot the spacecraft position as a colored circle with black
                # outline and a label.
                x_nudge = 0.0
                y_nudge = 8.0
                sc_label = sc_metadata[sc_id]["label"]
                color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
                for ax in (ax_v, ax_n, ax_T, ax_Br):
                    ax.plot(x_sc, y_sc, 'o', c=color)
                    ax.plot(x_sc, y_sc, 'o', c="black", fillstyle="none")
                    ax.text(x_sc + x_nudge, y_sc + y_nudge, sc_label,
                            c="black", horizontalalignment="center")

        # Save the figure to a file.
        path = os.path.join(fdir, f"{pictype}-{i_step}.png")
        if debug:
            print(f"path = {path}")
        kv.savePic(path, bLenX=45)
        frame_files.append(path)

    if debug:
        print(f"frame_files = {frame_files}")

    # Assemble the frames into a movie.
    cmd = ["convert",  "-delay", "10", "-loop", "0"]
    cmd += frame_files
    movie_file = os.path.join(fdir, f"{pictype}.gif")
    cmd.append(movie_file)
    if debug:
        print(f"cmd = {cmd}")
    if verbose:
        print(f"Assembling frames into {movie_file}.")
    subprocess.run(cmd, check=True)

    # Return the path to the movie file.
    return movie_file


def create_pic2_movie(args):
    """Create a pic2-style gamhelio movie.

    Create a pic2-style gamhelio movie.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options.

    Returns
    -------
    movie_file : str
        Path to movie file.

    Raises
    ------
    None
    """
    # Extract arguments.
    debug = args.debug
    pictype = args.pictype
    spacecraft = args.spacecraft
    verbose = args.verbose

    # Fetch the plot limits based on the picture type.
    plot_limits = hviz.GetSizeBds(pictype)
    if debug:
        print(f"plot_limits = {plot_limits}")

    # Create all plot images in a memory buffer.
    mpl.use("Agg")

    # Fetch the figure size.
    figsize = figure_sizes[pictype]
    if debug:
        print(f"figsize = {figsize}")

    # Create figures in a memory buffer.
    mpl.use("Agg")

    # Create the figure.
    fig = plt.figure(figsize=figsize)
    if debug:
        print(f"fig = {fig}")

    # Lay out the subplots for this figure. The grid contains separate axes
    # to use for the color bars (the rows with relative heights of 1).
    nrows = 4
    ncols = 6
    gs = mpl.gridspec.GridSpec(nrows, ncols, height_ratios=[20, 1, 20, 1])
    if debug:
        print(f"gs = {gs}")

    # Create the Axes objects for the individual subplots.
    # Each subplot is 1 row x 3 columns in the grid.
    ax_v = fig.add_subplot(gs[0, :3])   # Upper left
    ax_n = fig.add_subplot(gs[0, 3:])   # Upper right
    ax_T = fig.add_subplot(gs[2, :3])   # Lower left
    ax_Br = fig.add_subplot(gs[2, 3:])  # Lower right
    if debug:
        print(f"ax_v = {ax_v}")
        print(f"ax_n = {ax_n}")
        print(f"ax_T = {ax_T}")
        print(f"ax_Br = {ax_Br}")

    # Create the Axes objects for the individual color bars.
    # Each color bar is 1 (thin) row x 3 columns in the grid.
    ax_cb_v = fig.add_subplot(gs[1, :3])   # Upper left
    ax_cb_n = fig.add_subplot(gs[1, 3:])   # Upper right
    ax_cb_T = fig.add_subplot(gs[3, :3])   # Lower left
    ax_cb_Br = fig.add_subplot(gs[3, 3:])  # Lower right
    if debug:
        print(f"ax_cb_v = {ax_cb_v}")
        print(f"ax_cb_n = {ax_cb_n}")
        print(f"ax_cb_T = {ax_cb_T}")
        print(f"ax_cb_Br = {ax_cb_Br}")

    # Open a "pipe" to the data for this run.
    fdir = args.directory
    ftag = args.runid
    gsph = hsph.GamsphPipe(fdir, ftag)

    # If spacecraft positions will be plotted, read the spacecraft metadata
    # and fetch trajectories.
    if spacecraft:
        sc_metadata = scutils.getScIds(spacecraft_data_file=sc_metadata_path)
        if debug:
            print(f"sc_metadata = {sc_metadata}")

        # Split the list into individual spacecraft names.
        spacecraft_names = spacecraft.split(',')
        if debug:
            print(f"spacecraft_names = {spacecraft_names}")

        # Fetch the spacecraft trajectories.
        sc_t, sc_x, sc_y, sc_z = fetch_spacecraft_trajectories(
            spacecraft_names, gsph
        )

    # Create and save frame images for each step.
    first_step = args.first_step
    last_step = args.last_step
    if last_step == -1:
        last_step = gsph.Nt
    else:
        last_step += 1
    if debug:
        print(f"first_step, last_step = {first_step, last_step}")
    frame_files = []
    for i_step in range(first_step, last_step):
        if verbose:
            print(f"Creating {pictype} frame for step {i_step}.")

        # Extract the MJD for the frame.
        mjd = gsph.MJDs[i_step]
        if debug:
            print(f"mjd = {mjd}")

        # Create the individual plots for this frame.
        hviz.PlotMerMagV(gsph, i_step, plot_limits, ax_v, ax_cb_v)
        hviz.PlotMerDNorm(gsph, i_step, plot_limits, ax_n, ax_cb_n)
        hviz.PlotMerTemp(gsph, i_step, plot_limits, ax_T, ax_cb_T)
        hviz.PlotMerBrNorm(gsph, i_step, plot_limits, ax_Br, ax_cb_Br)

        # Add time in the upper left.
        gsph.AddTime(i_step, ax_v, xy=[0.025, 0.875], fs="x-large")

        # Overlay spacecraft positions (optional).
        if spacecraft:
            for (i_sc, sc_id) in enumerate(spacecraft_names):

                # Skip this spacecraft if no trajectory available.
                if sc_id not in sc_t:
                    continue

                # Interpolate the spacecraft position at the time for the plot.
                t_sc = mjd
                x_sc = np.interp(t_sc, sc_t[sc_id], sc_x[sc_id])
                y_sc = np.interp(t_sc, sc_t[sc_id], sc_y[sc_id])
                z_sc = np.interp(t_sc, sc_t[sc_id], sc_z[sc_id])

                # Plot the spacecraft position as a colored circle with black
                # outline and a label.
                x_nudge = 0.0
                y_nudge = 8.0
                sc_label = sc_metadata[sc_id]["label"]
                color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
                for ax in (ax_v, ax_n, ax_T, ax_Br):
                    ax.plot(x_sc, z_sc, 'o', c=color)
                    ax.plot(x_sc, z_sc, 'o', c="black", fillstyle="none")
                    ax.text(x_sc + x_nudge, z_sc + y_nudge, sc_label,
                            c="black", horizontalalignment="center")

        # Save the figure to a file.
        path = os.path.join(fdir, f"{pictype}-{i_step}.png")
        if debug:
            print(f"path = {path}")
        kv.savePic(path, bLenX=45)
        frame_files.append(path)

    if debug:
        print(f"frame_files = {frame_files}")

    # Assemble the frames into a movie.
    cmd = ["convert",  "-delay", "10", "-loop", "0"]
    cmd += frame_files
    movie_file = os.path.join(fdir, f"{pictype}.gif")
    cmd.append(movie_file)
    if debug:
        print(f"cmd = {cmd}")
    if verbose:
        print(f"Assembling frames into {movie_file}.")
    subprocess.run(cmd, check=True)

    # Return the path to the movie file.
    return movie_file


def create_pic3_movie(args):
    """Create a pic3-style gamhelio movie.

    Create a pic3-style gamhelio movie.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options.

    Returns
    -------
    movie_file : str
        Path to movie file.

    Raises
    ------
    None
    """
    # Extract arguments.
    debug = args.debug
    pictype = args.pictype
    spacecraft = args.spacecraft
    verbose = args.verbose

    # Fetch the plot limits based on the picture type.
    plot_limits = hviz.GetSizeBds(pictype)
    if debug:
        print(f"plot_limits = {plot_limits}")

    # Create all plot images in a memory buffer.
    mpl.use("Agg")

    # Fetch the figure size.
    figsize = figure_sizes[pictype]
    if debug:
        print(f"figsize = {figsize}")

    # Create figures in a memory buffer.
    mpl.use("Agg")

    # Create the figure.
    fig = plt.figure(figsize=figsize)
    if debug:
        print(f"fig = {fig}")

    # Lay out the subplots for this figure. The grid contains separate axes
    # to use for the color bars (the rows with relative heights of 1).
    nrows = 4
    ncols = 6
    gs = mpl.gridspec.GridSpec(nrows, ncols, height_ratios=[20, 1, 20, 1])
    if debug:
        print(f"gs = {gs}")

    # Create the Axes objects for the individual subplots.
    # Each subplot is 1 row x 3 columns in the grid.
    ax_v = fig.add_subplot(gs[0, :3])   # Upper left
    ax_n = fig.add_subplot(gs[0, 3:])   # Upper right
    ax_T = fig.add_subplot(gs[2, :3])   # Lower left
    ax_Br = fig.add_subplot(gs[2, 3:])  # Lower right
    if debug:
        print(f"ax_v = {ax_v}")
        print(f"ax_n = {ax_n}")
        print(f"ax_T = {ax_T}")
        print(f"ax_Br = {ax_Br}")

    # Create the Axes objects for the individual color bars.
    # Each color bar is 1 (thin) row x 3 columns in the grid.
    ax_cb_v = fig.add_subplot(gs[1, :3])   # Upper left
    ax_cb_n = fig.add_subplot(gs[1, 3:])   # Upper right
    ax_cb_T = fig.add_subplot(gs[3, :3])   # Lower left
    ax_cb_Br = fig.add_subplot(gs[3, 3:])  # Lower right
    if debug:
        print(f"ax_cb_v = {ax_cb_v}")
        print(f"ax_cb_n = {ax_cb_n}")
        print(f"ax_cb_T = {ax_cb_T}")
        print(f"ax_cb_Br = {ax_cb_Br}")

    # Open a "pipe" to the data for this run.
    fdir = args.directory
    ftag = args.runid
    gsph = hsph.GamsphPipe(fdir, ftag)

    # If spacecraft positions will be plotted, read the spacecraft metadata
    # and fetch trajectories.
    if spacecraft:
        sc_metadata = scutils.getScIds(spacecraft_data_file=sc_metadata_path)
        if debug:
            print(f"sc_metadata = {sc_metadata}")

        # Split the list into individual spacecraft names.
        spacecraft_names = spacecraft.split(',')
        if debug:
            print(f"spacecraft_names = {spacecraft_names}")

        # Fetch the spacecraft trajectories.
        sc_t, sc_x, sc_y, sc_z = fetch_spacecraft_trajectories(
            spacecraft_names, gsph
        )

    # Create and save frame images for each step.
    first_step = args.first_step
    last_step = args.last_step
    if last_step == -1:
        last_step = gsph.Nt
    else:
        last_step += 1
    if debug:
        print(f"first_step, last_step = {first_step, last_step}")
    frame_files = []
    for i_step in range(first_step, last_step):
        if verbose:
            print(f"Creating {pictype} frame for step {i_step}.")

        # Extract the MJD for the frame.
        mjd = gsph.MJDs[i_step]
        if debug:
            print(f"mjd = {mjd}")

        # Create the individual plots for this frame.
        hviz.PlotiSlMagV(gsph, i_step, plot_limits, ax_v, ax_cb_v)
        hviz.PlotiSlD(gsph, i_step, plot_limits, ax_n, ax_cb_n)
        hviz.PlotiSlTemp(gsph, i_step, plot_limits, ax_T, ax_cb_T)
        hviz.PlotiSlBr(gsph, i_step, plot_limits, ax_Br, ax_cb_Br)

        # Add time in the upper left.
        gsph.AddTime(i_step, ax_v, xy=[0.025, 0.875], fs="x-large")

        # Overlay spacecraft positions (optional).
        if spacecraft:
            for (i_sc, sc_id) in enumerate(spacecraft_names):

                # Skip this spacecraft if no trajectory available.
                if sc_id not in sc_t:
                    continue

                # Interpolate the spacecraft position at the time for the plot.
                t_sc = mjd
                # x_sc = np.interp(t_sc, sc_t[sc_id], sc_x[sc_id])
                # y_sc = np.interp(t_sc, sc_t[sc_id], sc_y[sc_id])
                # z_sc = np.interp(t_sc, sc_t[sc_id], sc_z[sc_id])

                # Convert Cartesian location to heliocentric lon/lat.
                rxy = np.sqrt(sc_x[sc_id]**2 + sc_y[sc_id]**2)
                theta = np.arctan2(rxy, sc_z[sc_id])
                phi = np.arctan2(sc_y[sc_id], sc_x[sc_id])
                lat = np.degrees(np.pi/2 - theta)
                lon = np.degrees(phi)
                lat_sc = np.interp(t_sc, sc_t[sc_id], lat)
                lon_sc = np.interp(t_sc, sc_t[sc_id], lon)

                # Plot the spacecraft position as a colored circle with black
                # outline and a label.
                x_nudge = 0.0
                y_nudge = 8.0
                sc_label = sc_metadata[sc_id]["label"]
                color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
                for ax in (ax_v, ax_n, ax_T, ax_Br):
                    ax.plot(lon_sc, lat_sc, 'o', c=color)
                    ax.plot(lon_sc, lat_sc, 'o', c="black", fillstyle="none")
                    ax.text(lon_sc + x_nudge, lat_sc + y_nudge, sc_label,
                            c="black", horizontalalignment="center")

        # Save the figure to a file.
        path = os.path.join(fdir, f"{pictype}-{i_step}.png")
        if debug:
            print(f"path = {path}")
        kv.savePic(path, bLenX=45)
        frame_files.append(path)

    if debug:
        print(f"frame_files = {frame_files}")

    # Assemble the frames into a movie.
    cmd = ["convert",  "-delay", "10", "-loop", "0"]
    cmd += frame_files
    movie_file = os.path.join(fdir, f"{pictype}.gif")
    cmd.append(movie_file)
    if debug:
        print(f"cmd = {cmd}")
    if verbose:
        print(f"Assembling frames into {movie_file}.")
    subprocess.run(cmd, check=True)

    # Return the path to the movie file.
    return movie_file


def create_gamhelio_movie(args):
    """Create a gamhelio movie.

    Create a gamhelio movie.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options.

    Returns
    -------
    None

    Raises
    ------
    TypeError : If an invalid type code is provided.
    """
    # Extract arguments.
    debug = args.debug
    pictype = args.pictype

    # Check that a valid plot code was provided.
    if pictype not in valid_pictypes:
        raise TypeError(f"Invalid plot type ({pictype})!")

    # Make the movie for the selected plot type.
    if pictype == "pic1":
        movie_file = create_pic1_movie(args)
    elif pictype == "pic2":
        movie_file = create_pic2_movie(args)
    elif pictype == "pic3":
        movie_file = create_pic3_movie(args)
    else:
        raise TypeError(f"Invalid plot type ({pictype})!")
    if debug:
        print(f"movie_file = {movie_file}")

#     # Fetch the plot domain based on the picture type.
#     xyBds = hviz.GetSizeBds(pic)
#     print(xyBds)

#     # Do work?
#     doFast = False

#     # Create figures in a memory buffer.
#     mpl.use("Agg")

#     # Read the CDAWeb spacecraft database.
#     sc_metadata_path = os.path.join(
#         os.environ["KAIJUHOME"], "kaipy", "satcomp", "sc_helio.json"
#     )
#     sc_metadata = scutils.getScIds(spacecraft_data_file=sc_metadata_path)

#     # Determine figure size (width, height) (inches) based on the picture type.
#     if pic == "pic1" or pic == "pic2":
#         figSz = (10, 12.5)
#     elif pic == "pic3":
#         figSz = (10, 6.5)
#     elif pic == "pic4":
#         figSz = (10, 6)
#     elif pic == "pic5":
#         figSz = (12, 12)

#     # Create the figure.
#     fig = plt.figure(figsize=figSz)

#     # Lay out the subplots.
#     if pic == "pic1" or pic == "pic2" or pic == "pic3":
#         gs = gridspec.GridSpec(4, 6, height_ratios=[20, 1, 20, 1])
#         # Axes for plots.
#         AxL0 = fig.add_subplot(gs[0, 0:3])
#         AxR0 = fig.add_subplot(gs[0, 3:])
#         AxL1 = fig.add_subplot(gs[2, 0:3])
#         AxR1 = fig.add_subplot(gs[2, 3:])
#         # Axes for colorbars.
#         AxC1_0 = fig.add_subplot(gs[1, 0:3])
#         AxC2_0 = fig.add_subplot(gs[1, 3:])
#         AxC1_1 = fig.add_subplot(gs[3, 0:3])
#         AxC2_1 = fig.add_subplot(gs[3, 3:])
#     elif pic == "pic4":
#         gs = gridspec.GridSpec(2, 1, height_ratios=[20, 1])
#         Ax = fig.add_subplot(gs[0, 0])
#         AxC = fig.add_subplot(gs[1, 0])
#     elif pic == "pic5":
#         gs = gridspec.GridSpec(2, 2)
#         Ax = fig.add_subplot(gs[0, 0])
#         AxC = fig.add_subplot(gs[0, 1])
#         AxC1 = fig.add_subplot(gs[1, 0])

#     # Open a pipe to the results data.
#     gsph = hsph.GamsphPipe(fdir, ftag, doFast=doFast)
#     if nStp < 0:
#         nStp = gsph.sFin
#         print("Using Step %d" % nStp)

#     # Extract the date/time of the plot.
#     mjd = gsph.MJDs[nStp]
#     if debug:
#         print(f"mjd = {mjd}")

#     # Now create the actual plots.
#     if pic == "pic1":
#         # These are all equatorial plots in the XY plane of the HGS frame
#         # used by gamhelio.
#         hviz.PlotEqMagV(gsph, nStp, xyBds, AxL0, AxC1_0)
#         hviz.PlotEqD(gsph, nStp, xyBds, AxR0, AxC2_0)
#         hviz.PlotEqTemp(gsph, nStp, xyBds, AxL1, AxC1_1)
#         hviz.PlotEqBr(gsph, nStp, xyBds, AxR1, AxC2_1)
#     elif pic == "pic2":
#         hviz.PlotMerMagV(gsph, nStp, xyBds, AxL0, AxC1_0)
#         hviz.PlotMerDNorm(gsph, nStp, xyBds, AxR0, AxC2_0)
#         hviz.PlotMerTemp(gsph, nStp, xyBds, AxL1, AxC1_1)
#         hviz.PlotMerBrNorm(gsph, nStp, xyBds, AxR1, AxC2_1)
#     elif pic == "pic3":
#         hviz.PlotiSlMagV(gsph, nStp, xyBds, AxL0, AxC1_0)
#         hviz.PlotiSlD(gsph, nStp, xyBds, AxR0, AxC2_0)
#         hviz.PlotiSlTemp(gsph, nStp, xyBds, AxL1, AxC1_1)
#         hviz.PlotiSlBr(gsph, nStp, xyBds, AxR1, AxC2_1)
#     elif pic == "pic4":
#         hviz.PlotiSlBrRotatingFrame(gsph, nStp, xyBds, Ax, AxC)
#     elif pic == "pic5":
#         hviz.PlotDensityProf(gsph, nStp, xyBds, Ax)
#         hviz.PlotSpeedProf(gsph, nStp, xyBds, AxC)
#         hviz.PlotFluxProf(gsph, nStp, xyBds, AxC1)
#     else:
#         print("Pic is empty. Choose pic1 or pic2 or pic3")

#     # Add time in the upper left.
#     if pic == "pic1" or pic == "pic2":
#         gsph.AddTime(nStp, AxL0, xy=[0.025, 0.875], fs="x-large")
#     elif pic == "pic3":
#         gsph.AddTime(nStp, AxL0, xy=[0.015, 0.82], fs="small")
#     elif pic == "pic4" or pic == "pic5":
#         gsph.AddTime(nStp, Ax, xy=[0.015, 0.92], fs="small")
#     else:
#         print("Pic is empty. Choose pic1 or pic2 or pic3")

#     # Overlay the spacecraft trajectory, if needed.
#     if spacecraft:
#         print("Overplotting spacecraft trajectories of %s." % spacecraft)

#         # Split the list into individual spacecraft names.
#         spacecraft = spacecraft.split(',')
#         if debug:
#             print("spacecraft = %s" % spacecraft)

#         # Fetch the MJD start and end time of the model results.
#         fname = gsph.f0
#         if debug:
#             print("fname = %s" % fname)
#         MJD_start = kh5.tStep(fname, 0, aID="MJD")
#         if debug:
#             print("MJD_start = %s" % MJD_start)
#         MJD_end = kh5.tStep(fname, gsph.sFin, aID="MJD")
#         if debug:
#             print("MJD_end = %s" % MJD_end)

#         # Convert the start and stop MJD to a datetime object in UT.
#         ut_start = ktools.MJD2UT(MJD_start)
#         if debug:
#             print("ut_start = %s" % ut_start)
#         ut_end = ktools.MJD2UT(MJD_end)
#         if debug:
#             print("ut_end = %s" % ut_end)

#         # Get the MJDc value for use in computing the gamhelio frame.
#         MJDc = scutils.read_MJDc(fname)
#         if debug:
#             print("mjdc = %s" % MJDc)

#         # Fetch and plot the trajectory of each spacecraft from CDAWeb.
#         for (i_sc, sc_id) in enumerate(spacecraft):
#             if verbose:
#                 print("Fetching trajectory for %s." % sc_id)

#             # Fetch the spacecraft trajectory in whatever frame is available
#             # from CDAWeb.
#             sc_data = cdaweb_utils.fetch_helio_spacecraft_trajectory(
#                 sc_id, ut_start, ut_end
#             )
#             if sc_data is None:
#                 print("No trajectory found for %s." % sc_id)
#                 continue

#             # Ingest the trajectory by converting it to the GH(MJDc) frame.
#             if verbose:
#                 print("Converting ephemeris for %s into gamhelio format." %
#                       sc_id)
#             x, y, z = cdaweb_utils.ingest_helio_spacecraft_trajectory(
#                 sc_id, sc_data, MJDc
#             )
#             if debug:
#                 print("x, y, z = %s, %s, %s" % (x, y, z))

#             # Convert the datetime objects from the trajectory to MJD.
#             t_strings = np.array([str(t) for t in sc_data["Epoch"]])
#             t = astropy.time.Time(t_strings, scale='utc').mjd

#             # Interpolate the spacecraft position at the time for the plot.
#             t_sc = mjd
#             x_sc = np.interp(t_sc, t, x)
#             y_sc = np.interp(t_sc, t, y)
#             z_sc = np.interp(t_sc, t, z)

#             # If needed, compute heliocentric spherical coordinates.
#             if pic == "pic3" or pic == "pic4":
#                 rxy = np.sqrt(x**2 + y**2)
#                 theta = np.arctan2(rxy, z)
#                 phi = np.arctan2(y, x)
#                 lat = np.degrees(np.pi/2 - theta)
#                 lon = np.degrees(phi)
#                 lat_sc = np.interp(t_sc, t, lat)
#                 lon_sc = np.interp(t_sc, t, lon)

#             # Plot a labelled trajectory of the spacecraft. Also plot a larger
#             # dot at the last point in the trajectory.
#             # Left plot
#             SPACECRAFT_COLORS = list(mpl.colors.TABLEAU_COLORS.keys())
#             color = SPACECRAFT_COLORS[i_sc % len(SPACECRAFT_COLORS)]
#             x_nudge = 0.0
#             y_nudge = 8.0
#             sc_label = sc_metadata[sc_id]["label"]
#             if pic == "pic1":
#                 for ax in (AxL0, AxR0, AxL1, AxR1):
#                     ax.plot(x_sc, y_sc, 'o', c=color)
#                     ax.plot(x_sc, y_sc, 'o', c="black", fillstyle="none")
#                     ax.text(x_sc + x_nudge, y_sc + y_nudge, sc_label,
#                             c="black", horizontalalignment="center")
#             elif pic == "pic2":
#                 for ax in (AxL0, AxR0, AxL1, AxR1):
#                     ax.plot(x_sc, z_sc, 'o', c=color)
#                     ax.plot(x_sc, z_sc, 'o', c="black", fillstyle="none")
#                     ax.text(x_sc + x_nudge, z_sc + y_nudge, sc_label,
#                             c="black", horizontalalignment="center")
#             elif pic == "pic3":
#                 for ax in (AxL0, AxR0, AxL1, AxR1):
#                     ax.plot(lon_sc, lat_sc, 'o', c=color)
#                     ax.plot(lon_sc, lat_sc, 'o', c="black", fillstyle="none")
#                     ax.text(lon_sc + x_nudge, lat_sc + y_nudge, sc_label,
#                             c="black", horizontalalignment="center")
#             elif pic == "pic4":
#                 ax = Ax
#                 ax.plot(lon_sc, lat_sc, 'o', c=color)
#                 ax.plot(lon_sc, lat_sc, 'o', c="black", fillstyle="none")
#                 ax.text(lon_sc + x_nudge, lat_sc + y_nudge, sc_label,
#                         c="black", horizontalalignment="center")
#             elif pic == "pic5":
#                 pass

#     # Save the figure to a file.
#     path = os.path.join(fdir, fOut)
#     kv.savePic(path, bLenX=45)


def main():
    """Make a movie from a gamhelio run."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    if debug:
        print(f"args = {args}")

    # Create the movie based on the selected picture type.
    create_gamhelio_movie(args)


if __name__ == "__main__":
    """Begin main program."""
    main()
