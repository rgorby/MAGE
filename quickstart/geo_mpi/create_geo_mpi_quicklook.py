#!/usr/bin/env python


"""Make a quick-look plot for the geo_mpi quickstart case.

The quick-look plot is created by the msphpic.py script.
"""


# Import standard modules.
import argparse
import os
import subprocess

# Import 3rd-party modules.

# Import project-specific modules.


# Program constants and defaults

# Default identifier for model to run,
runid = "geo_mpi"

# Program description.
description = "Create a quick-look plot for the geo_mpi quicklook case."


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the command-line argument parser.

    Parameters
    ----------
    None

    Returns
    -------
    parse : argparse.ArgumentParser
        Parser for command-line arguments.

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


def create_quicklook_plot():
    """Create the quicklook plot for the run.

    Create the quicklook plot for the run.

    Parameters
    ----------
    None

    Returns
    -------
    figure_file_name : str
        Path to quicklook plot file.

    Raises
    ------
    None
    """
    # Run the quicklook plot generation script.
    cmd = "msphpic.py"
    args = [cmd, "-id", runid]
    subprocess.run(args, check=True)
    figure_file_name = "qkpic.png"
    return figure_file_name


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

    if verbose:
        print("Creating quicklook plot.")
    quicklook_file = create_quicklook_plot()
    if verbose:
        print(f"The quicklook plot is in {quicklook_file}.")


if __name__ == "__main__":
    """Begin main program."""
    main()
