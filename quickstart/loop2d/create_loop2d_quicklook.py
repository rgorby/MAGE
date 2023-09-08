#!/usr/bin/env python


"""Make a quick-look plot for the loop2d example run.

The quick-look plot displays the magnetic pressure:

Pb = 0.5*(Bx**2 + By**2 + Bz**2)

from the first and last steps in the HDF file.
"""


# Import standard modules.
import argparse
import os

# Import 3rd-party modules.
import matplotlib
import matplotlib.pyplot as plt

# Import project-specific modules.
import kaipy.gamera.gampp as gampp


# Program constants and defaults

# Default identifier for model to run,
default_runid = "loop2d"

# Program description.
description = "Create a quick-look plot for the loop2d quickstart case."


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

    Raises
    ------
    None
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--debug", "-d", action="store_true", default=False,
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
        "--verbose", "-v", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def create_quicklook_plot(directory, runid):
    """Create the quicklook plot for the loop2d run.

    Create the quicklook plot for the loop2d run.

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

    Raises
    ------
    None
    """
    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(directory, runid, doVerbose=False)

    # Read the grid coordinates.
    X = data_pipe.X[...]
    Y = data_pipe.Y[...]

    # Read the magnetic field components for the first step, and
    # compute the corresponding magnetic pressure.
    Bx = data_pipe.GetVar("Bx", data_pipe.s0, doVerb=False)[...]
    By = data_pipe.GetVar("By", data_pipe.s0, doVerb=False)[...]
    Bz = data_pipe.GetVar("Bz", data_pipe.s0, doVerb=False)[...]
    Pb_first = (Bx**2 + By**2 + Bz**2)/2

    # Read the magnetic field components for the last step, and
    # compute the corresponding magnetic pressure.
    Bx = data_pipe.GetVar("Bx", data_pipe.sFin, doVerb=False)[...]
    By = data_pipe.GetVar("By", data_pipe.sFin, doVerb=False)[...]
    Bz = data_pipe.GetVar("Bz", data_pipe.sFin, doVerb=False)[...]
    Pb_last = (Bx**2 + By**2 + Bz**2)/2

    # Plot parameters
    units = "code units"
    vmin = 0
    vmax = 6e-7

    # Create the figure in a memory buffer.
    matplotlib.use("Agg")
    fig, axes = plt.subplots(2, 1)

    # Plot the magnetic pressure from the first step.
    axes[0].set_aspect("equal")
    axes[0].set_ylabel("Y")
    values = axes[0].pcolormesh(X, Y, Pb_first, cmap="viridis",
                                vmin=vmin, vmax=vmax)
    axes[0].text(0.5, 0.4, "Step 0", color="white")

    # Plot the magnetic pressure from the last step.
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    values = axes[1].pcolormesh(X, Y, Pb_last, cmap="viridis",
                                vmin=vmin, vmax=vmax)
    axes[1].text(0.5, 0.4, f"Step {data_pipe.sFin}", color="white")

    # Show the shared colorbar.
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(values, cax=cbar_ax, label=f"Magnetic pressure [{units}]")

    # Set the plot title.
    plt.suptitle(f"Magnetic pressure at start and end for {runid}")

    # Save the quicklook plot.
    figure_file_name = f"{runid}_quicklook.png"
    fig.savefig(figure_file_name)

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
    verbose = args.verbose
    if debug:
        print(f"args = {args}")

    if verbose:
        print("Creating quicklook plot.")
    quicklook_file = create_quicklook_plot(directory, runid)
    if verbose:
        print(f"The quicklook plot is in {quicklook_file}")
