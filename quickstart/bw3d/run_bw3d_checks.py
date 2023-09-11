#!/usr/bin/env python


"""Perform sample computations on the results of the bw3d example.

Perform sample computations on the results of the bw3d example.
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
default_runid = "bw3d"

# Program description.
description = (
    "Perform sample computations on the results of the bw3d quickstart case."
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
        "--debug", "-d", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--directory", type=str, metavar="directory", default=os.getcwd(),
        help="Directory to contain files generated for the run "
             "(default: %(default)s)"
    )
    parser.add_argument(
        "--runid", type=str, metavar="runid", default=default_runid,
        help="ID string of the run (default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def compute_asymmetry_metric(directory, runid):
    """Compute the asymmetry metric at first and last steps.

    Compute the asymmetry metric for the first and last steps:

    Pb = (Bx**2 + By**2 + Bz**2)/2
    asymmetry_metric = SUM(ABS(Pb[i,j] - Pb[j,i])*dV)

    Parameters
    ----------
    directory : str
        Path to directory containing model results.
    runid : str
        ID string for model to run.

    Returns
    -------
    asymmetry_metric_first, asymmetry_metric_last : float
        Asymmetry metric (in code units) for first and last steps.
    """
    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(directory, runid, doVerbose=False)

    # Load the grid cell volumes.
    dV = data_pipe.GetVar("dV", None, doVerb=False)[...]

    # Load the magnetic field components from the first step, and use
    # to compute the magnetic pressure and asymmetry metric.
    Bx = data_pipe.GetVar("Bx", data_pipe.s0, doVerb=False)[...]
    By = data_pipe.GetVar("By", data_pipe.s0, doVerb=False)[...]
    Bz = data_pipe.GetVar("Bz", data_pipe.s0, doVerb=False)[...]
    Pb_first = (Bx**2 + By**2 + Bz**2)/2
    asymmetry_metric_first = np.sum(abs(Pb_first - Pb_first.T)*dV)

    # Load the magnetic field components from the last step, and use
    # to compute the magnetic pressure and the asymmetry metric.
    Bx = data_pipe.GetVar("Bx", data_pipe.sFin, doVerb=False)[...]
    By = data_pipe.GetVar("By", data_pipe.sFin, doVerb=False)[...]
    Bz = data_pipe.GetVar("Bz", data_pipe.sFin, doVerb=False)[...]
    Pb_last = (Bx**2 + By**2 + Bz**2)/2
    asymmetry_metric_last = np.sum(abs(Pb_last - Pb_last.T)*dV)

    return asymmetry_metric_first, asymmetry_metric_last


def main():
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    directory = args.directory
    runid = args.runid
    verbose = args.verbose

    # Compute the asymmetry metric.
    if verbose:
        print("Computing asymmetry metric.")
    asymmetry_metric_first, asymmetry_metric_last = (
        compute_asymmetry_metric(directory, runid)
    )
    print("Asymmetry metric (SUM(ABS(Pb - Pb.T)*dV), code units):")
    print(f"At start: {asymmetry_metric_first}")
    print(f"At end: {asymmetry_metric_last}")


if __name__ == "__main__":
    """Begin main program."""
    main()
