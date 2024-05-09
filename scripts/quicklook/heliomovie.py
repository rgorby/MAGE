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
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import spacepy.datamodel as dm
from sunpy.coordinates import frames

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
valid_pictypes = ("pic1", "pic2", "pic3", "pic4", "pic5", "pic6", "pic7")

# Default movie format.
default_movie_format = "mp4"

# Valid movie format strings.
valid_movie_formats = ("gif", "mp4")

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
    "pic6": (12.5, 12.5),
    "pic7": (10, 12.5),
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
        "--clobber", action="store_true",
        help="Overwrite existing frame and movie files (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--directory", "-d", type=str, metavar="directory",
        default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "--first_step", "-n0", type=int, metavar="n0",
        default=default_first_step,
        help="First time step to plot (default: %(default)s)"
    )
    parser.add_argument(
        "--hgsplot", action="store_true",
        help="Plot in the Heliographic Stonyhurst frame corresponding to the "
             "date of the plot (default: %(default)s)."
    )
    parser.add_argument(
        "--jslice", "-jslice", type=int, metavar="jSlice", default=None,
        help="Index of j-slice for pic7 (default: Nj/2-1)"
    )
    parser.add_argument(
        "--last_step", "-n1", type=int, metavar="n1",
        default=default_last_step,
        help="Last time step to plot (default: %(default)s)"
    )
    parser.add_argument(
        "--movie_format", type=str, metavar="movie_format",
        default=default_movie_format,
        help="Output movie format (default: %(default)s)"
    )
    parser.add_argument(
        "--pictype", "-p", type=str, metavar="pictype",
        default=default_pictype,
        help="Code for plot type (default: %(default)s)"
    )
    parser.add_argument(
        "--runid", "-id", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Names of spacecraft to plot positions, separated by commas"
        " (default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", default=False,
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
    for sc_id in spacecraft_names:

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


def assemble_frames_into_gif(frame_files, args):
    """Assemble individual frames into an animated GIF movie.

    Assemble individual frames into an animated GIF movie.

    Parameters
    ----------
    frame_files : list of str
        List of paths to frame images, in frame order.
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
    # Fetch arguments.
    movie_directory = args.directory
    pictype = args.pictype

    # Create the movie.
    cmd = ["convert",  "-delay", "10", "-loop", "0"]
    cmd += frame_files
    if args.hgsplot:
        movie_file = os.path.join(movie_directory, f"{pictype}-HGS.gif")
    else:
        movie_file = os.path.join(movie_directory, f"{pictype}.gif")
    cmd.append(movie_file)
    # NOTE: movie_file is overwritten by default.
    subprocess.run(cmd, check=True)

    # Return the path to the movie file.
    return movie_file


def assemble_frames_into_mp4(frame_files, args):
    """Assemble individual frame images into a MP4 movie.

    Assemble individual frame images into a MP4 movie.

    Parameters
    ----------
    frame_files : list of str (ignored)
        List of paths to frame images, in frame order.
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
    # Fetch arguments.
    movie_directory = args.directory
    pictype = args.pictype

    # Create the movie.
    frame_directory = os.path.split(frame_files[0])[0]
    if args.hgsplot:
        movie_file = os.path.join(movie_directory, f"{pictype}-HGS.mp4")
        frame_pattern = os.path.join(frame_directory, f"{pictype}-HGS-%06d.png")
    else:
        movie_file = os.path.join(movie_directory, f"{pictype}.mp4")
        frame_pattern = os.path.join(frame_directory, f"{pictype}-%06d.png")
    cmd = [
        "ffmpeg", "-r", "4", "-s", "1920x1080",
        "-i", frame_pattern,
        "-vcodec", "libx264", "-crf", "25", "-pix_fmt", "yuv420p",
        movie_file
    ]
    if args.clobber:
        cmd += ["-y"]
    subprocess.run(cmd, check=True)

    # Return the path to the movie file.
    return movie_file


def assemble_frames_into_movie(frame_files, args):
    """Assemble individual frame images into a movie.

    Assemble individual frame images into a movie.

    Parameters
    ----------
    frame_files : list of str
        List of paths to frame images, in frame order.
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
    # Extract command-line options.
    movie_format = args.movie_format

    # Assemble the movie.
    if movie_format == "gif":
        movie_file = assemble_frames_into_gif(frame_files, args)
    elif movie_format == "mp4":
        movie_file = assemble_frames_into_mp4(frame_files, args)
    else:
        raise TypeError(f"Invalid movie format: {movie_format}")

    # Return the path to the movie file.
    return movie_file


def GHtoHGS(mjd_gh, x_gh, y_gh, z_gh, mjd_hgs):
    """Convert Cartesin GH coordinates to HGS.

    Convert Cartesian coordinates in the gamhelio frame at time mjdc to
    the Heliographic Sonyhurst frame at time mjd.

    NOTE: The gamhelio frame at time t is related to the Heliographic
    Stonyhurst frame at time t by the reflection of the x- and y-axes:

    x_gh(t) = -x_hgs(t)
    y_gh(t) = -y_hgs(t)
    z_gh(t) = z_hgs(t)

    Since HGS is a time-dependent frame, a time must be provided for each set
    of coordinates.

    Parameters
    ----------
    mjd_gh : float
        MJD of source gamhelio frame
    x_gh, y_gh, z_gh : np.array of float (any shape) or scalar float
        Cartesian coordinates in GH(mjdc) frame. All three arrays x, y, z must
        have identical shapes.
    mjd_hgs : float
        MJD of target HGS frame

    Returns
    -------
    x_hgs, y_hgs, z_hgs : np.array of float (same shape as x_gh, y_gh, z_gh)
        Cartesian coordinates converted to HGS(mjd) frame.

    Raises
    ------
    None
    """
    # Load the source coordinates (originially in the GH(mjd_gh) frame) into
    #  the equivalent HGS(mjd_gh) frame.
    c_gh = SkyCoord(
        -x_gh*u.Rsun, -y_gh*u.Rsun, z_gh*u.Rsun,
        frame=frames.HeliographicStonyhurst,
        obstime=ktools.MJD2UT(mjd_gh),
        representation_type="cartesian"
    )

    # Create the target Heliographic Stonyhurst frame.
    hgs_frame = frames.HeliographicStonyhurst(
        obstime=ktools.MJD2UT(mjd_hgs)
    )

    # Convert the coordinates from GH(mjd_gh) to HGS(mjd_hgs).
    c_hgs = c_gh.transform_to(hgs_frame)

    # Extract and return the converted coordinates.
    x_hgs = dm.dmarray(c_hgs.cartesian.x)
    y_hgs = dm.dmarray(c_hgs.cartesian.y)
    z_hgs = dm.dmarray(c_hgs.cartesian.z)
    return x_hgs, y_hgs, z_hgs


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
    hgsplot=args.hgsplot
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

    # Compute the path to the frame directory.
    if hgsplot:
        frame_directory = os.path.join(fdir, f"frames-{pictype}-HGS")
    else:
        frame_directory = os.path.join(fdir, f"frames-{pictype}")
    try:
        os.mkdir(frame_directory)
    except FileExistsError as e:
        if args.clobber:
            pass
        else:
            raise e

    # Get the MJDc value for use in computing the gamhelio frame.
    fname = gsph.f0
    MJDc = scutils.read_MJDc(fname)

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
        hviz.PlotEqMagV(gsph, i_step, plot_limits, ax_v, ax_cb_v,
                        hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotEqD(gsph, i_step, plot_limits, ax_n, ax_cb_n,
                     hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotEqTemp(gsph, i_step, plot_limits, ax_T, ax_cb_T,
                        hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotEqBr(gsph, i_step, plot_limits, ax_Br, ax_cb_Br,
                      hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        if hgsplot:
            fig.suptitle("Heliographic Stonyhurst frame for "
                            f"{ktools.MJD2UT(mjd)}")
        else:
            fig.suptitle(f"GAMERA-Helio frame for {ktools.MJD2UT(mjd)}")

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

                # If needed, convert the position to HGS(mjd).
                if hgsplot:
                    x_sc, y_sc, z_sc = GHtoHGS(MJDc, x_sc, y_sc, z_sc, mjd)

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
        if hgsplot:
            path = os.path.join(frame_directory, f"{pictype}-HGS-{i_step:06d}.png")
        else:
            path = os.path.join(frame_directory, f"{pictype}-{i_step:06d}.png")
        if debug:
            print(f"path = {path}")
        kv.savePic(path, bLenX=45)
        frame_files.append(path)

    if debug:
        print(f"frame_files = {frame_files}")

    # Assemble the frames into a movie.
    if verbose:
        print("Assembling frames into movie.")
    movie_file = assemble_frames_into_movie(frame_files, args)
    if verbose:
        print(f"Movie file {movie_file} created.")

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
    hgsplot=args.hgsplot
    pictype = args.pictype
    spacecraft = args.spacecraft
    verbose = args.verbose

    # Fetch the plot limits based on the picture type.
    plot_limits = hviz.GetSizeBds(pictype)
    if debug:
        print(f"plot_limits = {plot_limits}")

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

    # Compute the path to the frame directory.
    if hgsplot:
        frame_directory = os.path.join(fdir, f"frames-HGS-{pictype}-HGS")
    else:
        frame_directory = os.path.join(fdir, f"frames-{pictype}")
    try:
        os.mkdir(frame_directory)
    except FileExistsError as e:
        if args.clobber:
            pass
        else:
            raise e

    # Get the MJDc value for use in computing the gamhelio frame.
    fname = gsph.f0
    MJDc = scutils.read_MJDc(fname)

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
        hviz.PlotMerMagV(gsph, i_step, plot_limits, ax_v, ax_cb_v,
                         hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotMerDNorm(gsph, i_step, plot_limits, ax_n, ax_cb_n,
                          hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotMerTemp(gsph, i_step, plot_limits, ax_T, ax_cb_T,
                         hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotMerBrNorm(gsph, i_step, plot_limits, ax_Br, ax_cb_Br,
                           hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        if hgsplot:
            fig.suptitle("Heliographic Stonyhurst frame for "
                            f"{ktools.MJD2UT(mjd)}")
        else:
            fig.suptitle(f"GAMERA-Helio frame for {ktools.MJD2UT(mjd)}")

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

                # If needed, convert the position to HGS(mjd).
                if hgsplot:
                    x_sc, y_sc, z_sc = GHtoHGS(MJDc, x_sc, y_sc, z_sc, mjd)

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
        if hgsplot:
            path = os.path.join(frame_directory, f"{pictype}-HGS-{i_step:06d}.png")
        else:
            path = os.path.join(frame_directory, f"{pictype}-{i_step:06d}.png")
        if debug:
            print(f"path = {path}")
        kv.savePic(path, bLenX=45)
        frame_files.append(path)

    if debug:
        print(f"frame_files = {frame_files}")

    # Assemble the frames into a movie.
    if verbose:
        print("Assembling frames into movie.")
    movie_file = assemble_frames_into_movie(frame_files, args)
    if verbose:
        print(f"Movie file {movie_file} created.")

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
    hgsplot=args.hgsplot
    pictype = args.pictype
    spacecraft = args.spacecraft
    verbose = args.verbose

    # Fetch the plot limits based on the picture type.
    plot_limits = hviz.GetSizeBds(pictype)
    if debug:
        print(f"plot_limits = {plot_limits}")

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

    # Compute the path to the frame directory.
    if hgsplot:
        frame_directory = os.path.join(fdir, f"frames-HGS-{pictype}-HGS")
    else:
        frame_directory = os.path.join(fdir, f"frames-{pictype}")
    try:
        os.mkdir(frame_directory)
    except FileExistsError as e:
        if args.clobber:
            pass
        else:
            raise e

    # Get the MJDc value for use in computing the gamhelio frame.
    fname = gsph.f0
    MJDc = scutils.read_MJDc(fname)

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
        AU_RSUN = 215.0
        radius = AU_RSUN
        hviz.PlotiSlMagV(gsph, i_step, plot_limits, ax_v, ax_cb_v, idx=radius,
                         idx_is_radius=True,
                         hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotiSlD(gsph, i_step, plot_limits, ax_n, ax_cb_n, idx=radius,
                      idx_is_radius=True,
                      hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotiSlTemp(gsph, i_step, plot_limits, ax_T, ax_cb_T, idx=radius,
                         idx_is_radius=True,
                         hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        hviz.PlotiSlBr(gsph, i_step, plot_limits, ax_Br, ax_cb_Br, idx=radius,
                       idx_is_radius=True,
                       hgsplot=hgsplot, MJDc=MJDc, MJD_plot=mjd)
        if hgsplot:
            fig.suptitle("Heliographic Stonyhurst frame at 1 AU for "
                            f"{ktools.MJD2UT(mjd)}")
        else:
            fig.suptitle("GAMERA-Helio frame at 1 AU for "
                            f"{ktools.MJD2UT(mjd)}")

        # Overlay spacecraft positions (optional).
        if spacecraft:
            for (i_sc, sc_id) in enumerate(spacecraft_names):

                # Skip this spacecraft if no trajectory available.
                if sc_id not in sc_t:
                    continue

                # Interpolate the spacecraft position at the time for the plot.
                t_sc = mjd

                # If needed, convert the position to HGS(mjd).
                if hgsplot:
                    sc_x[sc_id], sc_y[sc_id], sc_z[sc_id] = (
                        GHtoHGS(MJDc, sc_x[sc_id], sc_y[sc_id], sc_z[sc_id],
                                mjd)
                    )

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
        if hgsplot:
            path = os.path.join(frame_directory, f"{pictype}-HGS-{i_step:06d}.png")
        else:
            path = os.path.join(frame_directory, f"{pictype}-{i_step:06d}.png")
        if debug:
            print(f"path = {path}")
        kv.savePic(path, bLenX=45)
        frame_files.append(path)

    if debug:
        print(f"frame_files = {frame_files}")

    # Assemble the frames into a movie.
    if verbose:
        print("Assembling frames into movie.")
    movie_file = assemble_frames_into_movie(frame_files, args)
    if verbose:
        print(f"Movie file {movie_file} created.")

    # Return the path to the movie file.
    return movie_file


def create_pic4_movie(args):
    """Create a pic4-style gamhelio movie.

    Create a pic4-style gamhelio movie.

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
    # to use for the color bar (the row with relative heights of 1).
    nrows = 2
    ncols = 1
    gs = mpl.gridspec.GridSpec(nrows, ncols, height_ratios=[20, 1])
    if debug:
        print(f"gs = {gs}")

    # Create the Axes objects for the individual subplots.
    # Each subplot is 1 row x 3 columns in the grid.
    ax_Br = fig.add_subplot(gs[0, 0])
    if debug:
        print(f"ax_Br = {ax_Br}")

    # Create the Axes objects for the individual color bars.
    # Each color bar is 1 (thin) row x 3 columns in the grid.
    ax_cb_Br = fig.add_subplot(gs[1, 0])
    if debug:
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

    # Compute the path to the frame directory.
    frame_directory = os.path.join(fdir, f"frames-{pictype}")
    try:
        os.mkdir(frame_directory)
    except FileExistsError as e:
        if args.clobber:
            pass
        else:
            raise e

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
        hviz.PlotiSlBrRotatingFrame(gsph, i_step, plot_limits, ax_Br,
                                    ax_cb_Br)

        # Add time in the upper left.
        gsph.AddTime(i_step, ax_Br, xy=[0.015, 0.92], fs="small")

        # Overlay spacecraft positions (optional).
        if spacecraft:
            for (i_sc, sc_id) in enumerate(spacecraft_names):

                # Skip this spacecraft if no trajectory available.
                if sc_id not in sc_t:
                    continue

                # Interpolate the spacecraft position at the time for the plot.
                t_sc = mjd

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
                for ax in (ax_Br,):
                    ax.plot(lon_sc, lat_sc, 'o', c=color)
                    ax.plot(lon_sc, lat_sc, 'o', c="black", fillstyle="none")
                    ax.text(lon_sc + x_nudge, lat_sc + y_nudge, sc_label,
                            c="black", horizontalalignment="center")

        # Save the figure to a file.
        path = os.path.join(frame_directory, f"{pictype}-{i_step:06d}.png")
        if debug:
            print(f"path = {path}")
        kv.savePic(path, bLenX=45)
        frame_files.append(path)

    if debug:
        print(f"frame_files = {frame_files}")

    # Assemble the frames into a movie.
    if verbose:
        print("Assembling frames into movie.")
    movie_file = assemble_frames_into_movie(frame_files, args)
    if verbose:
        print(f"Movie file {movie_file} created.")

    # Return the path to the movie file.
    return movie_file


def create_pic5_movie(args):
    """Create a pic5-style gamhelio movie.

    Create a pic5-style gamhelio movie.

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
    # to use for the color bar (the row with relative heights of 1).
    nrows = 2
    ncols = 2
    gs = mpl.gridspec.GridSpec(nrows, ncols)
    if debug:
        print(f"gs = {gs}")

    # Create the Axes objects for the individual subplots.
    # Each subplot is 1 row x 3 columns in the grid.
    ax_n = fig.add_subplot(gs[0, 0])
    ax_v = fig.add_subplot(gs[0, 1])
    ax_mf = fig.add_subplot(gs[1, 0])
    if debug:
        print(f"ax_n = {ax_n}")
        print(f"ax_v = {ax_v}")
        print(f"ax_mf = {ax_mf}")

    # Open a "pipe" to the data for this run.
    fdir = args.directory
    ftag = args.runid
    gsph = hsph.GamsphPipe(fdir, ftag)

    # Compute the path to the frame directory.
    frame_directory = os.path.join(fdir, f"frames-{pictype}")
    try:
        os.mkdir(frame_directory)
    except FileExistsError as e:
        if args.clobber:
            pass
        else:
            raise e

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
        hviz.PlotDensityProf(gsph, i_step, plot_limits, ax_n)
        hviz.PlotSpeedProf(gsph, i_step, plot_limits, ax_v)
        hviz.PlotFluxProf(gsph, i_step, plot_limits, ax_mf)

        # Add time in the upper left.
        gsph.AddTime(i_step, ax_n, xy=[0.015, 0.92], fs="small")

        # Save the figure to a file.
        path = os.path.join(frame_directory, f"{pictype}-{i_step:06d}.png")
        if debug:
            print(f"path = {path}")
        kv.savePic(path, bLenX=45)
        frame_files.append(path)

    if debug:
        print(f"frame_files = {frame_files}")

    # Assemble the frames into a movie.
    if verbose:
        print("Assembling frames into movie.")
    movie_file = assemble_frames_into_movie(frame_files, args)
    if verbose:
        print(f"Movie file {movie_file} created.")

    # Return the path to the movie file.
    return movie_file


def create_pic6_movie(args):
    """Create a pic6-style gamhelio movie.

    Create a pic6-style gamhelio movie.

    Parameters
    ----------
    args : dict
        Dictionary of command-line arguments.

    Returns
    -------
    movie_file : str
        Path to movie file.

    Raises
    ------
    None
    """
    # Return the path to the movie file.
    movie_file = None
    return movie_file


def create_pic7_movie(args):
    """Create a pic7-style gamhelio movie.

    Create a pic7-style gamhelio movie.

    Parameters
    ----------
    args : dict
        Dictionary of command-line arguments.

    Returns
    -------
    movie_file : str
        Path to movie file.

    Raises
    ------
    None
    """
    # Return the path to the movie file.
    movie_file = None
    return movie_file


def create_gamhelio_movie(args):
    """Create a gamhelio movie.

    Create a gamhelio movie.

    Parameters
    ----------
    args : dict
        Dictionary of command-line arguments.

    Returns
    -------
    movie_file : str
        Path to movie file.

    Raises
    ------
    TypeError : If an invalid type code is provided.
    """
    # Extract arguments.
    debug = args.debug
    pictype = args.pictype

    # Check that a valid picture type code was provided.
    if pictype not in valid_pictypes:
        raise TypeError(f"Invalid plot type ({pictype})!")

    # Make the movie for the selected plot type.
    if pictype == "pic1":
        movie_file = create_pic1_movie(args)
    elif pictype == "pic2":
        movie_file = create_pic2_movie(args)
    elif pictype == "pic3":
        movie_file = create_pic3_movie(args)
    elif pictype == "pic4":
        movie_file = create_pic4_movie(args)
    elif pictype == "pic5":
        movie_file = create_pic5_movie(args)
    elif pictype == "pic6":
        movie_file = create_pic6_movie(args)
    elif pictype == "pic7":
        movie_file = create_pic7_movie(args)
    else:
        raise TypeError(f"Invalid plot type ({pictype})!")
    if debug:
        print(f"movie_file = {movie_file}")

    # Return the path to the movie file.
    return movie_file


def main():
    """Make a movie from a gamhelio run."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")

    # Create the movie based on the selected picture type.
    movie_file = create_gamhelio_movie(args)
    if args.verbose:
        print(f"The movie is available in {movie_file}.")


if __name__ == "__main__":
    main()
