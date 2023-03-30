#!/usr/bin/env python


"""Make a quick-look plot for the bw3d example run.

The quick-look plot displays the magnetic pressure:

Pb = 0.5*(Bx**2 + By**2 + Bz**2)

from the first and nth steps in the HDF file, with an overlay of the
magnetic field lines.
"""


# Import standard modules.
import argparse
import os
from pathlib import Path
import re

# Import 3rd-party modules.
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Import project-specific modules.
import kaipy.gamera.gampp as gampp


# Program constants and defaults

# Default identifier for model to run,
default_runid = "bw3d"

# Program description.
description = "Create a quick-look plot for the %s example." % default_runid

# Some RGB colors
light_grey = (0.5, 0.5, 0.5)
transparent_red = (1.0, 0, 0, 0.5)


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the command-line argument parser.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-d", "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--directory", type=str, metavar="directory", default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "--runid", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-s", "--step", type=int, metavar="step", default=7,
        help="Step to use for second plot (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def determine_mpi_grid_shape(directory, runid):
    """Determine the shape of the MPI grid.

    Determine the numbers of planes, rows, and columns of the MPI grid.

    Parameters
    ----------
    directory : str
        Path to the directory containing the results.
    runid : str
        ID string for the run to examine.

    Returns
    -------
    nx, ny, nz : int
        Length of each dimension of the MPI grid.
    """

    # Get the path object for the data directory.
    path = Path(directory)

    # Start with a null shape.
    nx = None
    ny = None
    nz = None

    # Get a list of the files for this runID.
    pattern =  r"%s_(\d{4})_(\d{4})_(\d{4})_\d{4}_\d{4}_\d{4}.*.h5" % runid
    regexp = re.compile(pattern)
    for f in path.iterdir():
        filename = os.path.split(f)[-1]
        m = regexp.search(filename)
        if m:
            nx = int(m.group(1))
            ny = int(m.group(2))
            nz = int(m.group(3))
            break

    return (nx, ny, nz)


def create_quicklook_plot(directory, runid):
    """Create the quicklook plot for the run.

    Create the quicklook plot for the run.

    Parameters
    ----------
    directory : str
        Path to directory containing results.
    runid : str
        ID string for results to examine.

    Returns
    -------
    figure_file_name : str
        Path to quicklook plot file.
    """
    # Determine the shape of the MPI grid.
    (n_x, n_y, n_z) = determine_mpi_grid_shape(directory, runid)

    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(directory, runid, doVerbose=False)

    # Find the layer closest to Z=0.
    Z = data_pipe.Z[...]
    kz_0 = np.argmin(np.abs(Z[0, 0, :]))
    z_0 = Z[0, 0, kz_0]

    # Read the grid coordinates from the Z=0 layer.
    X = data_pipe.X[..., kz_0]
    Y = data_pipe.Y[..., kz_0]

    # Determine the X and Y limits.
    X_min = X.min()
    X_max = X.max()
    Y_min = Y.min()
    Y_max = Y.max()

    # Compute a modified grid with the coordinates of the cell centers.
    XT = X.T
    (n_rows, n_cols) = XT.shape
    n_rows -= 1
    n_cols -= 1
    XTC = np.empty((n_rows, n_cols))
    for row in range(n_rows):
        for col in range(n_cols):
            XTC[row, col] = (XT[row, col] + XT[row, col + 1])/2
    YT = Y.T
    (n_rows, n_cols) = YT.shape
    n_rows -= 1
    n_cols -= 1
    YTC = np.empty((n_rows, n_cols))
    for row in range(n_rows):
        for col in range(n_cols):
            YTC[row, col] = (YT[row, col] + YT[row + 1, col])/2

    # Compute the MPI tile boundaries assuming they are equally-sized along
    # each dimension.
    mpi_tiles_x = np.linspace(X_min, X_max, n_x + 1)
    mpi_tiles_y = np.linspace(Y_min, Y_max, n_y + 1)

    # Read the pressure and B-field components at the first and nth steps.
    P_first = data_pipe.GetVar("P", data_pipe.s0, doVerb=False)[..., kz_0]
    Bx_first = data_pipe.GetVar("Bx", data_pipe.s0, doVerb=False)[..., kz_0]
    By_first = data_pipe.GetVar("By", data_pipe.s0, doVerb=False)[..., kz_0]
    P_last = data_pipe.GetVar("P", step, doVerb=False)[..., kz_0]
    Bx_last = data_pipe.GetVar("Bx", step, doVerb=False)[..., kz_0]
    By_last = data_pipe.GetVar("By", step, doVerb=False)[..., kz_0]

    # Plot parameters
    units = "code units"
    vmin = 0
    vmax = 1

    # Create the figure in a memory buffer.
    matplotlib.use("Agg")
    fig, axes = plt.subplots(2, 1)

    # Plot the pressure and XY-plane B-field from the first step,
    # with MPI tiling.
    axes[0].set_aspect("equal")
    axes[0].set_ylabel("Y")
    values = axes[0].pcolormesh(X, Y, P_first, cmap="viridis", vmin=vmin,
                                vmax=vmax)
    axes[0].streamplot(XTC, YTC, Bx_first.T, By_first.T, linewidth=0.5,
                       color="tab:orange", density=[0.5, 1])
    for x in mpi_tiles_x[1:-1]:
        axes[0].axvline(x, linestyle="--", linewidth=0.5, color=light_grey)
    for y in mpi_tiles_y[1:-1]:
        axes[0].axhline(y, linestyle="--", linewidth=0.5, color=light_grey)
    axes[0].text(0.5, 0.4, "Step 0", color="white")

    # Plot the pressure and XY-plane B-field from the nth step,
    # with MPI tiling.
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    values = axes[1].pcolormesh(X, Y, P_last, cmap="viridis", vmin=vmin,
                                vmax=vmax)
    axes[1].streamplot(XTC, YTC, Bx_last.T, By_last.T, linewidth=0.5,
                       color="tab:orange", density=[0.5, 1])
    for x in mpi_tiles_x[1:-1]:
        axes[1].axvline(x, linestyle="--", linewidth=0.5, color=light_grey)
    for y in mpi_tiles_y[1:-1]:
        axes[1].axhline(y, linestyle="--", linewidth=0.5, color=light_grey)
    axes[1].text(0.5, 0.4, "Step %s" % step, color="white")

    # Show the shared colorbar.
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(values, cax=cbar_ax, label="%s [%s]" % ("Pressure", units))

    # Set the plot title.
    plt.suptitle("Pressure at start and step %s for %s" % (step, runid))

    # Save the quicklook plot.
    figure_file_name = "%s_quicklook.png" % (runid)
    plt.savefig(figure_file_name)

    return figure_file_name


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    directory = args.directory
    runid = args.runid
    step = args.step
    verbose = args.verbose

    if verbose:
        print("Creating quicklook plot.")
    quicklook_file = create_quicklook_plot(directory, runid)
    if verbose:
        print("The quicklook plot is in %s." % quicklook_file)
