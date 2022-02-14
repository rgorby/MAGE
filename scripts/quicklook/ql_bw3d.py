#!/usr/bin/env python

"""Make a quick-look plot for the bw3d example run.

The quick-look plot displays the magnetic pressure:

Pb = 0.5*(Bx**2 + By**2 + Bz**2)

from the first and nth steps in the HDF file.
"""


import argparse
import glob
import os
from pathlib import Path
import re

import matplotlib
from matplotlib import style
import matplotlib.pyplot as plt
import numpy as np

import kaipy.kaiH5 as kaih5
import kaipy.gamera.gampp as gampp


# Some RGB colors
light_grey = (0.5, 0.5, 0.5)
transparent_red = (1.0, 0, 0, 0.5)

def determine_mpi_grid_shape(directory, runid):
    """Determine the shape of the MPI grid."""

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


if __name__ == "__main__":
    """Make a quick-look plot for the bw3d example run."""

    # Defaults for argument parser

    # Directory to read data from.
    directory = os.getcwd()

    # Data tag
    runid = "blast3D"

    # Step number to show in lower plot.
    step = 7

    # Program description.
    description = "Create a quick-look plot (Pb at start and end) for the bw3d test case."

    # Create the argument parser.
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-d", "--directory", type=str, metavar="directory",
        default=directory,
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "-id", "--runid", type=str, metavar="runid",
        default=runid,
        help="RunID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-s", "--step", type=int, metavar="step",
        default=step,
        help="RunID of data (default: %(default)s)"
    )

    # Parse the command-line arguments.
    args = parser.parse_args()
    directory = args.directory
    runid = args.runid
    step = args.step

    # Determine the shape of the MPI grid.
    (n_x, n_y, n_z) = determine_mpi_grid_shape(directory, runid)

    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(directory, runid)

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
    if step is None:
        step = data_pipe.sFin
    P_first = data_pipe.GetVar("P", data_pipe.s0, doVerb=False)[..., kz_0]
    Bx_first = data_pipe.GetVar("Bx", data_pipe.s0, doVerb=False)[..., kz_0]
    By_first = data_pipe.GetVar("By", data_pipe.s0, doVerb=False)[..., kz_0]
    P_last = data_pipe.GetVar("P", step, doVerb=False)[..., kz_0]
    Bx_last = data_pipe.GetVar("Bx", step, doVerb=False)[..., kz_0]
    By_last = data_pipe.GetVar("By", step, doVerb=False)[..., kz_0]

    # Plot parameters
    name = "Pressure at z = %s" % z_0
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
    values = axes[0].pcolormesh(X, Y, P_first, cmap="viridis", vmin=vmin, vmax=vmax)
    axes[0].streamplot(XTC, YTC, Bx_first.T, By_first.T, linewidth=0.5, color=transparent_red, density=[0.5, 1])
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
    values = axes[1].pcolormesh(X, Y, P_last, cmap="viridis", vmin=vmin, vmax=vmax)
    axes[1].streamplot(XTC, YTC, Bx_last.T, By_last.T, linewidth=0.5, color=transparent_red, density=[0.5, 1])
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
