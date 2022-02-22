#!/usr/bin/env python


"""Run the gamera code on the loop2d test case and perform quality checks.

Run the loop2d example case using the gamera code. Generate associated data
products, and perform QA by comparing the generated results to the expected
values.
"""


# Import standard modules.
import argparse

# Import 3rd-party modules.

# Import project-specific modules.


# Program constants and defaults

# Program description.
description = "Prepare and run gamera on the loop2d test case."

# Location of template .ini file and .ini file for run.
ini_template = "loop2d.ini.template"
ini_file = "loop2d.ini"

# Location of template XML file and XML file for run.
xml_template = "loop2d.xml.template"
xml_file = "loop2d.xml"

# Location of template PBS script and PBS script for run.
pbs_template = "loop2d.pbs.template"
pbs_file = "loop2d.pbs"


def create_command_line_parser():
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-d", "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--startdate", type=str, default="2016-08-09T09:00:00",
        help="Specify the start date in ISO 8601 format (default: %(default)s)."
    )
    parser.add_argument(
        "--stopdate", type=str, default="2016-08-09T10:00:00",
        help="Specify the stop date in ISO 8601 format (default: %(default)s)."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def run_preprocessing_steps():
    """Run any preprocessing steps needed for the run."""
    # The loop2d example does not require any preprocessing steps.


def create_ini_file():
    """Create the .ini file from a template."""
    # Just use the template for now.
    with open(ini_template) as t:
        lines = t.readlines()
    with open(ini_file, "w") as f:
        f.writelines(lines)


def convert_ini_to_xml():
    """Convert the .ini file to XML."""
    # No conversion is performed yet. Just print the template.
    with open(xml_template) as t:
        lines = t.readlines()
    with open(xml_file, "w") as f:
        f.writelines(lines)


def create_pbs_job_script():
    """Create the PBS job script for the run."""
    # No changes yet. Just print the template.
    with open(pbs_template) as t:
        lines = t.readlines()
    with open(pbs_file, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    start_date = args.startdate
    stop_date = args.stopdate
    verbose = args.verbose

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps()

    # Create the .ini file.
    if verbose:
        print("Creating .ini file for run.")
    create_ini_file()

    # Convert the .ini file to a .xml file.
    if verbose:
        print("Converting .ini file to .xml file for run.")
    convert_ini_to_xml()

    # Create the PBS job script.
    if verbose:
        print("Creating PBS job script for run.")
    create_pbs_job_script()
