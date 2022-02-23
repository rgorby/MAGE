#!/usr/bin/env python


"""Prepare the PBS script for a loop2d job.

Create the PBS script required to submit a gamera loop2d model to PBS.
"""


# Import standard modules.
import argparse
import os

# Import 3rd-party modules.

# Import project-specific modules.


# Program constants and defaults

# Program description.
description = "Prepare and run gamera on the loop2d test case."

# Location of template .ini file and .ini file for run.
ini_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "loop2d", "loop2d.ini.template"
)
ini_file = "loop2d.ini"

# Location of template XML file and XML file for run.
xml_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "loop2d", "loop2d.xml.template"
)
xml_file = "loop2d.xml"

# Location of template PBS script and PBS script for run.
pbs_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", "loop2d", "loop2d.pbs.template"
)
pbs_file = "loop2d.pbs"


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


def run_preprocessing_steps():
    """Run any preprocessing steps needed for the run.
    
    Run any required preprocessing steps to prepare for the model run.

    There are no preprocessing steps for the loop2d model.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # The loop2d example does not require any preprocessing steps.


def create_ini_file():
    """Create the .ini file from a template.
    
    Create the .ini file describing the loop2d model run.

    For now, we simply make a copy of the .ini template.

    Parameters
    ----------
    None

    Returns
    -------
    ini_file : str
        Path to the .ini file for the loop2d model run.
    """
    # Just use the template for now.
    with open(ini_template) as t:
        lines = t.readlines()
    with open(ini_file, "w") as f:
        f.writelines(lines)
    return ini_file


def convert_ini_to_xml(ini_file):
    """Convert the .ini file to XML.
    
    Convert the .ini file describing the loop2d run to the corresponding
    XML file.

    Parameters
    ----------
    ini_file : str
        Path to the .ini file to convert.
    
    Returns
    -------
    xml_file : str
        Path to the resulting XML file.
    """
    # No conversion is performed yet. Just print the template.
    with open(xml_template) as t:
        lines = t.readlines()
    with open(xml_file, "w") as f:
        f.writelines(lines)
    return xml_file


def create_pbs_job_script():
    """Create the PBS job script for the run.
    
    Create the PBS job script which can be submitted to PBS to perform
    the loop2d model run.

    Parameters
    ----------
    None

    Returns
    -------
    pbs_file : str
        Path to PBS job description file.
    """
    # No changes yet. Just print the template.
    with open(pbs_template) as t:
        lines = t.readlines()
    with open(pbs_file, "w") as f:
        f.writelines(lines)
    return pbs_file


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    verbose = args.verbose

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps()

    # Create the .ini file.
    if verbose:
        print("Creating .ini file for run.")
    ini_file = create_ini_file()

    # Convert the .ini file to a .xml file.
    if verbose:
        print("Converting .ini file to .xml file for run.")
    xml_file = convert_ini_to_xml(ini_file)

    # Create the PBS job script.
    if verbose:
        print("Creating PBS job script for run.")
    pbs_file = create_pbs_job_script()

    print("""
The PBS job script is ready for submission.
You may submit it with the following command:

    qsub %s""" % pbs_file)