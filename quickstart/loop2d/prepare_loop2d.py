#!/usr/bin/env python


"""Prepare for running serial kaiju on the loop2d example.

Perform the preprocessing required to run the serial kaiju code on the
loop2d example. Create any required data files, and create the PBS
script to run the code.
"""


# Import standard modules.
import argparse
import os
import subprocess

# Import 3rd-party modules.

# Import project-specific modules.


# Program constants and defaults

# Default identifier for model to run,
default_runid = "loop2d"

# Program description.
description = "Prepare to run serial kaiju on the loop2d quickstart case."

# Location of template .ini file file for run.
ini_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "loop2d", "loop2d_template.ini"
)

# Location of template PBS script script for run.
pbs_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "loop2d", "loop2d_template.pbs"
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


def run_preprocessing_steps(directory, runid):
    """Run any preprocessing steps needed for the run.

    Perform required preprocessing steps.

    Parameters
    ----------
    directory : str
        Path to directory to receive preprocessing results.
    runid : str
        ID string for the model to run.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # The loop2d example does not require any preprocessing steps.


def create_ini_file(directory, runid):
    """Create the .ini file from a template.

    Create the .ini file from a template.

    Parameters
    ----------
    directory : str
        Path to directory to receive .ini file.
    runid : str
        ID string for the model to run.

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

    # Write the processed .ini file to the run directory.
    ini_file = os.path.join(directory, f"{runid}.ini")
    with open(ini_file, "w") as f:
        f.writelines(lines)
    return ini_file


def convert_ini_to_xml(ini_file, xml_file):
    """Convert the .ini file to XML.

    Convert the .ini file to a .xml file.

    Parameters
    ----------
    ini_file : str
        Path to .ini file to convert.
    xml_file : str
        Path to .xml file to create.

    Returns
    -------
    None

    Raises
    ------
    None
    """
    cmd = "XMLGenerator.py"  # Must be in PATH.
    args = [cmd, ini_file, xml_file]
    subprocess.run(args, check=True)


def create_pbs_job_script(directory, runid):
    """Create the PBS job script for the run.

    Create the PBS job script from a template.

    Parameters
    ----------
    directory : str
        Path to directory to contain PBS job script.
    runid : str
        ID string for model to run.

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
    pbs_file = os.path.join(directory, f"{runid}.pbs")
    with open(pbs_file, "w") as f:
        f.writelines(lines)
    return pbs_file


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    debug = args.debug
    directory = args.directory
    runid = args.runid
    verbose = args.verbose

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps(directory, runid)

    # Create the .ini file.
    if verbose:
        print("Creating .ini file for run.")
    ini_file = create_ini_file(directory, runid)

    # Convert the .ini file to a .xml file.
    if verbose:
        print("Converting .ini file to .xml file for run.")
    xml_file = os.path.join(directory, f"{runid}.xml")
    convert_ini_to_xml(ini_file, xml_file)

    # Create the PBS job script.
    if verbose:
        print("Creating PBS job script for run.")
    pbs_file = create_pbs_job_script(directory, runid)
    if verbose:
        print(f"The PBS job script {pbs_file} is ready.")
        print("Edit this file as needed for your system "
              "(see comments in the file for more information).")
        print("Submit the job to PBS with the command:")
        print(f"    qsub {pbs_file}")
  
