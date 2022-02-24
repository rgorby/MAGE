#!/usr/bin/env python


"""Perform quality checks on the results of running the gamera loop2d case..

Examine the results of running the loop2d example through gamera, and
apply consistency checks.
"""


# Import standard modules.
import argparse
# import os

# Import 3rd-party modules.
import numpy as np

# Import project-specific modules.
import kaipy.gamera.gampp as gampp


# Program constants and defaults

# Program description.
description = "Perform consistency checks on the gamera loop2d test case."


def verify_volume_integrated_magnetic_pressure():
    """Verify the volume-integrated magnetic pressure."""

    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(".", "loop2d")

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

    print(Pb_integrated_first, Pb_integrated_last)

    
def create_command_line_parser():
    """Create the command-line argument parser."""
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
        print("Verifying volume-integrated magnetic pressure.")
    verify_volume_integrated_magnetic_pressure()
