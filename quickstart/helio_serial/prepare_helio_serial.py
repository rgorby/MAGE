#!/usr/bin/env python


"""Prepare to run serial kaiju on the helio_serial quickstart case.

Prepare to run serial kaiju on the helio_serial quickstart case. Perform any
required preprocessing steps, and create the PBS script to run the code.
"""


# Import standard modules.
import argparse
import os
import shutil
import subprocess

# Import 3rd-party modules.

# Import project-specific modules.


# Program constants and defaults

# Default identifier for run.
runid = "helio_serial"

# Program description.
description = "Prepare to run serial kaiju on the helio_serial quickstart case."

# Location of template .ini file.
ini_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "helio_serial",
    "helio_serial_template.ini"
)

# Location of template PBS script.
pbs_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "helio_serial",
    "helio_serial_template.pbs"
)

# Location of .ini file for wsa2gamera.py.
wsa2gamera_ini_path = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "helio_serial", "startup.config"
)

# Location of FITS file for wsa2gamera.py.
wsa2gamera_fits_path = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "helio_serial",
    "vel_201708132000R002_ahmi.fits"
)


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

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


def run_preprocessing_steps():
    """Run any preprocessing steps needed for this run.

    Perform required preprocessing steps.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Copy the .ini and FITS file for the boundary conditions.
    shutil.copy(wsa2gamera_ini_path, ".")
    shutil.copy(wsa2gamera_fits_path, ".")

    # Create the grid and inner boundary conditions files.
    # wsa2gamera.py must be in PATH.
    cmd = "wsa2gamera.py"
    args = [cmd, "startup.config"]
    subprocess.run(args, check=True)


def create_ini_file():
    """Create the .ini file from a template.

    Create the .ini file from a template.

    Parameters
    ----------
    None

    Returns
    -------
    ini_file : str
        Path to .ini file.

    Raises
    ------
    None
    """
    # Read the file template.
    with open(ini_template) as t:
        lines = t.readlines()

    # Process the template here.

    # Write the processed .ini file.
    ini_file = f"{runid}.ini"
    with open(ini_file, "w") as f:
        f.writelines(lines)
    return ini_file


def convert_ini_to_xml(ini_file):
    """Convert the .ini file to XML.

    Convert the .ini file to a .xml file.

    Parameters
    ----------
    ini_file : str
        Path to .ini file to convert.

    Returns
    -------
    xml_file : str
        Path to .xml just-created XML file.

    Raises
    ------
    None
    """
    cmd = "XMLGenerator.py"  # Must be in PATH.
    xml_file = f"{runid}.xml"
    args = [cmd, ini_file, xml_file]
    subprocess.run(args, check=True)
    return xml_file


def create_pbs_job_script():
    """Create the PBS job script for the run.

    Create the PBS job script from a template.

    Parameters
    ----------
    None

    Returns
    -------
    pbs_file : str
        Path to PBS job script.

    Raises
    ------
    None
    """
    # Read the template.
    with open(pbs_template) as t:
        lines = t.readlines()

    # Process the template here.

    # Write out the processed file.
    pbs_file = f"{runid}.pbs"
    with open(pbs_file, "w") as f:
        f.writelines(lines)
    return pbs_file


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

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps()

    # Create the .ini file.
    if verbose:
        print("Creating .ini file.")
    ini_file = create_ini_file()

    # Convert the .ini file to a .xml file.
    if verbose:
        print("Converting .ini file to .xml file.")
    xml_file = convert_ini_to_xml(ini_file)

    # Create the PBS job script.
    if verbose:
        print("Creating PBS job script.")
    pbs_file = create_pbs_job_script()
    if verbose:
        print(f"The PBS job script {pbs_file} is ready.")
        print("Edit this file as needed for your system "
              "(see comments in the file for more information).")
        print("Submit the job to PBS with the command:")
        print(f"    qsub {pbs_file}")


if __name__ == "__main__":
    """Begin main program."""
    main()
