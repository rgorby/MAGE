#!/usr/bin/env python


"""Perform quality checks on the results of running the gamera bw3d case.

Examine the results of running the bw3d example through gamera, and
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
description = "Perform consistency checks on the gamera bw3d test case."

# Location of quick-look plots file.
quicklook_file = "bw3d_quicklook.png"

    
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


def compute_symmetry_metric():
    """Compute the symmetry metric at start and end.
    
    Compute the symmetry metric for the first and last steps:

    Pb = (Bx**2 + By**2 + Bz**2)/2
    symmetry_metric = SUM( (Pb[i,j] - Pb[j,i])*dV)

    Parameters
    ----------
    None

    Returns
    -------
    symmetry_metric_first, symmetry_metric_last : float
        Symmetry metric (in code units) for first and last steps.
    """

    # Open a pipe to the data file.
    data_pipe = gampp.GameraPipe(".", "bw3d", doVerbose=False)

    # Load the grid cell volumes.
    dV = data_pipe.GetVar("dV", None, doVerb=False)[...]

    # Load the magnetic field components from the first step, and use
    # to compute the magnetic pressure.
    Bx = data_pipe.GetVar("Bx", data_pipe.s0, doVerb=False)[...]
    By = data_pipe.GetVar("By", data_pipe.s0, doVerb=False)[...]
    Bz = data_pipe.GetVar("Bz", data_pipe.s0, doVerb=False)[...]
    Pb_first = (Bx**2 + By**2 + Bz**2)/2
    symmetry_metric_first = np.sum((Pb_first - Pb_first.T)*dV)

    # Load the magnetic field components from the last step, and use
    # to compute the magnetic pressure.
    Bx = data_pipe.GetVar("Bx", data_pipe.sFin, doVerb=False)[...]
    By = data_pipe.GetVar("By", data_pipe.sFin, doVerb=False)[...]
    Bz = data_pipe.GetVar("Bz", data_pipe.sFin, doVerb=False)[...]
    Pb_last = (Bx**2 + By**2 + Bz**2)/2
    Pb_integrated_last = np.sum(Pb_last*dV)
    symmetry_metric_last = np.sum((Pb_last - Pb_last.T)*dV)

    return symmetry_metric_first, symmetry_metric_last


def create_quicklook_plot():
    """Create the quick-look plots in a file.
    
    Create the quicklook plots summarizing the bw3d run.

    Parameters
    ----------
    None

    Returns
    -------
    quicklook_file : str
        Path to the file containing the quick-look plots.
    """
    cmd = "ql_bw3d.py"
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
        print("Computing symmetry metric.")
    symmetry_metric_first, symmetry_metric_last = compute_symmetry_metric()
    print("Symmetry metric (SUM((P - P.T)*dV):")
    print("At start: %s (code units)" % symmetry_metric_first)
    print("At end: %s (code units)" % symmetry_metric_last)

    # Generate the quick-look plot.
    if verbose:
        print("Creating quick-look plot.")
    quicklook_file = create_quicklook_plot()
    print("The quicklook plot is in %s." % quicklook_file)
