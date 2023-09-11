#!/usr/bin/env python


"""Prepare to run serial kaiju on the geo_serial quickstart case.

Prepare to run serial kaiju on the geo_serial quickstart case. Perform any
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
runid = "geo_serial"

# Program description.
description = "Prepare to run serial kaiju on the geo_serial quickstart case."

# Location of template .ini file.
ini_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "geo_serial",
    "geo_serial_template.ini"
)

# Location of template PBS script.
pbs_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "geo_serial",
    "geo_serial_template.pbs"
)

# Start and stop date of data to fetch from CDAWeb.
start_date = "2016-08-09T09:00:00"
stop_date = "2016-08-09T10:00:00"


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
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
    # Create the grid file using "double" resolution.
    cmd = "genLFM.py"
    args = [cmd, "-gid", "D"]
    subprocess.run(args, check=True)

    # Create the solar wind file.
    cmd = "cda2wind.py"
    args = [cmd, "-t0", start_date, "-t1", stop_date, "-interp"]
    subprocess.run(args, check=True)

    # Create the RCM configuration file.
    cmd = "genRCM.py"
    args = [cmd]
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
        Path to the .ini file to convert.
    
    Returns
    -------
    xml_file : str
        Path to the resulting XML file.

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
    """Create the PBS job script.

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
