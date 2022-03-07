#!/usr/bin/env python


"""Make a quick-look plot for the helio_serial example run.

The quick-look plot is created by the msphpic.py script.
"""


# Import standard modules.
import argparse
import os
import subprocess

# Import 3rd-party modules.
import matplotlib
import matplotlib.pyplot as plt

# Import project-specific modules.
import kaipy.gamera.gampp as gampp


# Program constants and defaults

# Default identifier for model to run,
default_runid = "helio_serial"

# Program description.
description = "Create a quick-look plot (Pb at start and end) for the %s test case." % default_runid


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
        "--directory", type=str, metavar="directory", default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "--runid", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def create_quicklook_plot(directory, runid):
    """Create the quicklook plot for the run.
    
    Create the quicklook plot for the run.
    
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
    """
    # Save the starting directory.
    initial_directory = os.getcwd()

    # Move to the directory containing the results.
    os.chdir(directory)

    # cmd = "msphpic.py"
    # args = []
    # subprocess.run([cmd] + args)
    # figure_file_name = os.path.join(directory, "qkpic.png")

    # # Move back to the starting directory.
    # os.chdir(initial_directory)

    # return figure_file_name
    return None


if __name__ == "__main__":
    """Make a quick-look plot for the helio_serial example run."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    directory = args.directory
    runid = args.runid
    verbose = args.verbose

    if verbose:
        print("Creating quicklook plot.")
    quicklook_file = create_quicklook_plot(directory, runid)
    if verbose:
        print("The quicklook plot is in %s." % quicklook_file)
