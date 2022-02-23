#!/usr/bin/env python


"""Perform quality checks on the results of running the gamera loop2d case..

Examine the results of running the loop2d example through gamera, and
apply consistency checks.
"""


# Import standard modules.
import argparse
import subprocess

# Import 3rd-party modules.
import numpy as np

# Import project-specific modules.
import kaipy.gamera.gampp as gampp


# Program constants and defaults

# Program description.
description = "Perform consistency checks on the gamera loop2d test case."

# Location of quick-look plots file.
quicklook_file = "loop2d_quicklook.png"

    
def create_command_line_parser():
    """Create the command-line argument parser.
    
    Ceate the parser for command-line arguments.

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
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def compute_volume_integrated_magnetic_pressure():
    """Compute the volume-integrated magnetic pressure at start and end.
    
    Compute the volume-integrated magnetic pressre for the first and last
    steps:

    Pb = (Bx**2 + By**2 + Bz**2)/2
    Pb_integrated = SUM(Pb*dV)

    Parameters
    ----------
    None

    Returns
    -------
    Pb_integrated_first, Pb_integrated_last : float
        Volume-integrated magnetic pressure (in code units) for first and
        last steps.
    """

    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(".", "loop2d", doVerbose=False)

    # Load the grid cell volumes.
    dV = data_pipe.GetVar("dV", None, doVerb=False)[...]

    # Load the magnetic field components from the first step, and use
    # to compute the magnetic pressure.
    Bx = data_pipe.GetVar("Bx", data_pipe.s0, doVerb=False)[...]
    By = data_pipe.GetVar("By", data_pipe.s0, doVerb=False)[...]
    Bz = data_pipe.GetVar("Bz", data_pipe.s0, doVerb=False)[...]
    Pb_first = (Bx**2 + By**2 + Bz**2)/2
    Pb_integrated_first = np.sum(Pb_first*dV)

    # Load the magnetic field components from the last step, and use
    # to compute the magnetic pressure.
    Bx = data_pipe.GetVar("Bx", data_pipe.sFin, doVerb=False)[...]
    By = data_pipe.GetVar("By", data_pipe.sFin, doVerb=False)[...]
    Bz = data_pipe.GetVar("Bz", data_pipe.sFin, doVerb=False)[...]
    Pb_last = (Bx**2 + By**2 + Bz**2)/2
    Pb_integrated_last = np.sum(Pb_last*dV)

    return Pb_integrated_first, Pb_integrated_last


def create_quicklook_plot():
    """Create the quick-look plots in a file.
    
    Create the quicklook plots summarizing the loop2d run.

    Parameters
    ----------
    None

    Returns
    -------
    quicklook_file : str
        Path to the file containing the quick-look plots.
    """
    cmd = "ql_loop2d.py"
    subprocess.run([cmd])
    return quicklook_file


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    verbose = args.verbose

    # Verify the volume-integrated magnetic pressure.
    if verbose:
        print("Computing volume-integrated magnetic pressure.")
    Pb_integrated_first, Pb_integrated_last = compute_volume_integrated_magnetic_pressure()
    print("Volume-integrated magnetic pressure:")
    print("At start: %s (code units)" % Pb_integrated_first)
    print("At end: %s (code units)" % Pb_integrated_last)

    # Generate the quick-look plot.
    if verbose:
        print("Creating quick-look plot.")
    quicklook_file = create_quicklook_plot()
    print("The quicklook plot is in %s." % quicklook_file)