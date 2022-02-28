#!/usr/bin/env python


"""Prepare for running MPI gamera on the bw3d example.

Perform the preprocessing required to run the MPI gamera code on the
bw3d example. Create any required data files, and create the PBS
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
default_runid = "bw3d"

# Program description.
description = "Prepare to run MPI gamera on the %s model." % default_runid

# Location of template .ini file.
ini_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", default_runid, "%s.ini.template"
    % default_runid
)

# Location of template XML file.
xml_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", default_runid, "%s.xml.template"
    % default_runid
)

# Location of template PBS script.
pbs_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", default_runid, "%s.pbs.template"
    % default_runid
)


def create_command_line_parser():
    """Create the command-line argument parser.

    Prepare the command-line parser.

    Parameters
    ----------
    None

    Returns
    -------
    parse : argparse.ArgumentParser
        Parser for command-line arguments.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-d", "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--directory", type=str, metavar="directory", default=os.getcwd(),
        help="Directory to contain files generated for the run (default: %(default)s)"
    )
    parser.add_argument(
        "--runid", type=str, metavar="runid", default=default_runid,
        help="ID string of the run (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def run_preprocessing_steps(directory, runid):
    """Run any preprocessing steps needed for the run.

    Run any required preprocessing steps to prepare for the model run.

    There are no preprocessing steps for the bw3d model.

    Parameters
    ----------
    directory : str
        Path to directory to receive preprocessing results.
    runid : str
        ID string for the model to run.

    Returns
    -------
    None
    """
    # The bw3d example does not require any preprocessing steps.


def create_ini_file(directory, runid):
    """Create the .ini file from a template.
    
    Create the .ini file describing the bw3d model run.

    For now, we simply make a copy of the .ini template.

    Parameters
    ----------
    directory : str
        Path to directory to receive .ini file.
    runid : str
        ID string for the model to run.

    Returns
    -------
    ini_file : str
        Path to the .ini file.
    """
    with open(ini_template) as t:
        lines = t.readlines()
    # Process the template here.
    ini_file = os.path.join(directory, "%s.ini" % runid)
    with open(ini_file, "w") as f:
        f.writelines(lines)
    return ini_file


def convert_ini_to_xml(ini_file, xml_file):
    """Convert the .ini file to XML.

    Convert the .ini file to a .xml file.

    NOTE: This conversion does not work yet, so just copy the template.

    Parameters
    ----------
    ini_file : str
        Path to .ini file to convert.
    xml_file : str
        Path to .xml file to create.

    Returns
    -------
    None
    """
    # cmd = "XMLGenerator.py"
    # args = [ini_file, xml_file]
    # subprocess.run([cmd] + args)

    # No conversion is performed yet. Just process the XML template.
    with open(xml_template) as t:
        lines = t.readlines()
    # Process the template here.
    with open(xml_file, "w") as f:
        f.writelines(lines)


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
    """
    with open(pbs_template) as t:
        lines = t.readlines()
    # Process the template here.
    pbs_file = os.path.join(directory, "%s.pbs" % runid)
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
    xml_file = os.path.join(directory, "%s.xml" % runid)
    convert_ini_to_xml(ini_file, xml_file)

    # Create the PBS job script.
    if verbose:
        print("Creating PBS job script for run.")
    pbs_file = create_pbs_job_script(directory, runid)
    if verbose:
        print("The PBS job script %s is ready." % pbs_file)
        print("Submit the job to PBS with the command:")
        print("    qsub %s" % pbs_file)
