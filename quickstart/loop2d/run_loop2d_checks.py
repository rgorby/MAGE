#!/usr/bin/env python


"""Perform sample computations on the results of the loop2d quickstart case.

Perform sample computations on the results of the loop2d quickstart case.
"""


# Import standard modules.
import argparse
import os

# Import 3rd-party modules.
import numpy as np

# Import project-specific modules.
import kaipy.gamera.gampp as gampp


# Program constants and defaults

# Default identifier for run.
runid = "loop2d"

# Program description.
description = "Run sample computations on the loop2d quickstart case."


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
        "--verbose", "-v", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def compute_volume_integrated_magnetic_pressure():
    """Compute the volume-integrated magnetic pressure at start and end.

    Compute the volume-integrated magnetic pressure for the first and
    last steps. This is equivalent to computing the total magnetic energy.

    Parameters
    ----------
    None

    Returns
    -------
    Pb_integrated_first, Pb_integrated_last : float
        Volume-integrated magnetic pressure for first and last steps.

    Raises
    ------
    None
    """
    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(".", runid, doVerbose=False)

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


def main():
    """Begin main program."""
    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    debug = args.debug
    verbose = args.verbose

    # Compute the volume-integrated magnetic pressure.
    if verbose:
        print("Computing volume-integrated magnetic pressure.")
    PbV1, PbV2 = compute_volume_integrated_magnetic_pressure()
    if verbose:
        print("Volume-integrated magnetic pressure (SUM(Pb*dV), code units):")
        print(f"At start: {PbV1}")
        print(f"At end: {PbV2}")


if __name__ == "__main__":
    """Begin main program."""
    main()
