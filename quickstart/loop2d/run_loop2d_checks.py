#!/usr/bin/env python


"""Perform sample computations on the results of the loop2d example.

Perform sample computations on the results of the loop2d example.
"""


# Import standard modules.
import argparse
import os

# Import 3rd-party modules.
import numpy as np

# Import project-specific modules.
import kaipy.gamera.gampp as gampp


# Program constants and defaults

# Default identifier for model to run,
default_runid = "loop2d"

# Program description.
description = (
    "Perform sample computations on the gamera %s test case." % default_runid
)


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
        help="Directory to contain files generated for the run (default: %(default)s)"
    )
    parser.add_argument(
        "--runid", type=str, metavar="runid", default=default_runid,
        help="ID string of the run (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def compute_volume_integrated_magnetic_pressure(directory, runid):
    """Compute the volume-integrated magnetic pressure at start and endxs.

    Compute the volume-integrated magnetic pressure for the first and
    last steps.

    Parameters
    ----------
    directory : str
        Path to directory containing model results.
    runid : str
        ID string for model to run.

    Returns
    -------
    Pb_integrated_first, Pb_integrated_last : float
        Volume-integrated magnetic pressure for first and last steps.
    """
    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(directory, runid, doVerbose=False)

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

    # Compute the volume-integrated magnetic pressure.
    if verbose:
        print("Computing volume-integrated magnetic pressure.")
    PbV1, PbV2 = compute_volume_integrated_magnetic_pressure(directory, runid)
    if verbose:
        print("Volume-integrated magnetic pressure (SUM(Pb*dV), code units):")
        print("At start: %s" % PbV1)
        print("At end: %s" % PbV2)
