#!/usr/bin/env python

"""Make a quick-look plot for the bw3d example run.

The quick-look plot displays the magnetic pressure:

Pb = 0.5*(Bx**2 + By**2 + Bz**2)

from the first and last steps in the HDF file.
"""


import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import kaipy.kaiH5 as kaih5
import kaipy.gamera.gampp as gampp


if __name__ == "__main__":
    """Make a quick-look plot for the bw3d example run."""

    # Defaults for argument parser

    # Directory to read data from.
    directory = os.getcwd()

    # Data tag
    runid = "blast3D"

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

    # Parse the command-line arguments.
    args = parser.parse_args()
    directory = args.directory
    runid = args.runid

    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(directory, runid)

    # Find the layer closest to Z=0.
    Z = data_pipe.Z[...]
    kz_0 = np.argmin(np.abs(Z[0, 0, :]))
    z_0 = Z[0, 0, kz_0]

    # Read the grid coordinates.
    X = data_pipe.X[..., kz_0]
    Y = data_pipe.Y[..., kz_0]

    # Read the magnetic field components for the first step, and compute
    # the corresponding magnetic pressure.
    Bx = data_pipe.GetVar("Bx", data_pipe.s0, doVerb=False)[..., kz_0]
    By = data_pipe.GetVar("By", data_pipe.s0, doVerb=False)[..., kz_0]
    Bz = data_pipe.GetVar("Bz", data_pipe.s0, doVerb=False)[..., kz_0]
    Pb_first = (Bx**2 + By**2 + Bz**2)/2

    # Read the magnetic field components for the last step, and compute
    # the corresponding magnetic pressure.
    Bx = data_pipe.GetVar("Bx", data_pipe.sFin, doVerb=False)[..., kz_0]
    By = data_pipe.GetVar("By", data_pipe.sFin, doVerb=False)[..., kz_0]
    Bz = data_pipe.GetVar("Bz", data_pipe.sFin, doVerb=False)[..., kz_0]
    Pb_last = (Bx**2 + By**2 + Bz**2)/2

    # Plot parameters
    name = "Magnetic pressure at z = %s" % z_0
    units = "code units"
    vmin = 0
    vmax = 1

    # Create the figure in a memory buffer.
    matplotlib.use("Agg")
    fig, axes = plt.subplots(2, 1)

    # Plot the magnetic pressure from the first step.
    axes[0].set_aspect("equal")
    axes[0].set_ylabel("Y")
    values = axes[0].pcolormesh(X, Y, Pb_first, cmap="viridis", vmin=vmin, vmax=vmax)
    axes[0].text(0.5, 0.4, "Step 0", color="white")

    # Plot the magnetic pressure from the last step.
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    values = axes[1].pcolormesh(X, Y, Pb_last, cmap="viridis", vmin=vmin, vmax=vmax)
    axes[1].text(0.5, 0.4, "Step %s" % data_pipe.sFin, color="white")

    # Show the shared colorbar.
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(values, cax=cbar_ax, label="%s [%s]" % ("Magnetic pressure", units))

    # Set the plot title.
    plt.suptitle("Magnetic pressure at start and end for %s" % (runid))

    # Save the quicklook plot.
    figure_file_name = "%s_quicklook.png" % (runid)
    plt.savefig(figure_file_name)
