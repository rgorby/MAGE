#!/usr/bin/env python

"""Make a quick-look plot for the loop2d example run.

The quick-look plot displays the magnetic pressure:

Pb = 0.5*(Bx**2 + By**2 + Bz**2)

from the first and last steps in the HDF file.
"""


import argparse
from curses import meta
import os

import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    """Make a quick-look plot for the loop2d example run."""

    # Defaults for argument parser

    # Directory to read data from.
    directory = os.getcwd()

    # Data tag
    runid = "loop2d"

    # Name of variable to plot.
    variable_name = "Pb"

    # Program description.
    description = "Create a quick-look plot (Pb at start and end) for the loop2d test case."

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

    # Open the data file.
    file_name = runid + ".gam.h5"
    path = os.path.join(directory, file_name)
    f = h5py.File(path, "r")

    # Load and the top-level datasets.
    X = f["X"]
    Y = f["Y"]

    # Fetch the coordinate limits.
    X_min = X[:].min()
    X_max = X[:].max()
    Y_min = Y[:].min()
    Y_max = Y[:].max()

    # Determine the number of time steps in the file.
    groups = f.keys()
    n_steps = 0
    for g in groups:
        if g.startswith("Step#"):
            n_steps += 1

    # Compute the name of the first and last steps.
    step_first = "Step#0"
    step_last = "Step#%s" % (n_steps - 1)

    # Compute the magnetic pressures.
    group = f[step_first]
    Bx = group["Bx"][:]
    By = group["By"][:]
    Bz = group["Bz"][:]
    Pb_first = 0.5*(Bx**2 + By**2 + Bz**2)
    group = f[step_last]
    Bx = group["Bx"][:]
    By = group["By"][:]
    Bz = group["Bz"][:]
    Pb_last = 0.5*(Bx**2 + By**2 + Bz**2)

    # Plot parameters
    name = "Magnetic pressure"
    units = "code units"
    vmin = 0
    vmax = 6e-7

    # Create the figure in a memory buffer.
    matplotlib.use("Agg")
    fig, axes = plt.subplots(2, 1)
    fig.subplots_adjust(hspace=0.3)

    # Pb at start
    values = axes[0].imshow(
        Pb_first, extent=(X_min, X_max, Y_min, Y_max),
        cmap="viridis", vmin=vmin, vmax=vmax)
    axes[0].set_ylabel("Y")
    axes[0].text(0.5, 0.4, "Step 0", color="white")

    # Pb at end
    values = axes[1].imshow(
        Pb_last, extent=(X_min, X_max, Y_min, Y_max),
        cmap="viridis", vmin=vmin, vmax=vmax)
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    axes[1].text(0.5, 0.4, "Step %s" % (n_steps - 1), color="white")

    plt.suptitle("Magnetic pressure at start and end for %s" % (runid))

    # Show the shared colorbar.
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(values, cax=cbar_ax, label="%s [%s]" % ("Magnetic pressure", units))

    # Save the quicklook plot.
    figure_file_name = "%s_quicklook.png" % (runid)
    plt.savefig(figure_file_name)

    # Close the data file.
    f.close()
