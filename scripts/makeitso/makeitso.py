#!/usr/bin/env python


"""makeitso for the kaiju software.

This master script is used to perform all of the steps needed to prepare
and run a kaiju job. This script is interactive - the user is prompted for
each decision that must be made to prepare for the run.
"""


# Import standard modules.
import argparse
import os
import subprocess

# Import 3rd-party modules.

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = "Interactive script to prepare and run a kaiju job."

# Default run ID string.
DEFAULT_RUNID = "helio"

# Default directory to run in.
DEFAULT_RUN_DIRECTORY = "."

# Default HPC system.
DEFAULT_HPC_SYSTEM = "pleiades"

# Default path to configuration file for wsa2gamera.py.
DEFAULT_CONFIG_FILE = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "makeitso", "wsa2gamera.ini"
)

# Default run type.
DEFAULT_RUN_TYPE = "serial"

# Location of template .ini file.
INI_TEMPLATE = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "makeitso",
    f"{DEFAULT_RUNID}_template.ini"
)

# Location of template PBS script.
pbs_template = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "makeitso",
    f"{DEFAULT_RUNID}_template.pbs"
)


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "-d", "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def get_run_options():
    """Prompt the user for run options.

    Prompt the user for run options.

    Parameters
    ----------
    None

    Returns
    -------
    options : dict
        Dictionary of program options, each entry maps str to str.
    """

    # Initialize the dictionary of program options.
    options = {}

    # Get the run ID.
    runid = input(f"Enter the name of the run ({DEFAULT_RUNID}): ")
    if runid == "":
        runid = DEFAULT_RUNID
    options["runid"] = runid

    # Enter the working directory for the run.
    run_directory = input(
        f"Enter the path to the run directory ({DEFAULT_RUN_DIRECTORY}): "
    )
    if run_directory == "":
        run_directory = DEFAULT_RUN_DIRECTORY
    options["run_directory"] = run_directory

    # Specify the HPC system to use.
    hpc_system = input(
        f"Enter the name of the HPC system to use ({DEFAULT_HPC_SYSTEM}): "
    )
    if hpc_system == "":
        hpc_system = DEFAULT_HPC_SYSTEM
    options["hpc_system"] = hpc_system

    # Specify the path to the location of the configuration file for
    # wsa2gamera.py.
    config_file = input(
        f"Enter the path to the wsa2gamera.py configuration file to use "
        f"({DEFAULT_CONFIG_FILE}): "
    )
    if config_file == "":
        config_file = DEFAULT_CONFIG_FILE
    options["config_file"] = config_file

    # Specify the run type (MPI or serial).
    run_type = input(
        f"Specify the run type (serial or mpi) ({DEFAULT_RUN_TYPE}): "
    )
    if run_type == "":
        run_type = DEFAULT_RUN_TYPE
    options["run_type"] = run_type

    # Return the options dictionary.
    return options


def run_preprocessing_steps(options):
    """Execute any preprocessing steps required for the run.

    Execute any preprocessing steps required for the run.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Returns
    -------
    None
    """
    # Save the current directory.
    original_directory = os.getcwd()

    # Move to the output directory.
    os.chdir(options["run_directory"])

    # Create the grid and inner boundary conditions files.
    # NOTE: Assumes wsa2gamera.py is in PATH.
    cmd = "wsa2gamera.py"
    args = [cmd, options["config_file"]]
    subprocess.run(
        args, check=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    # Print captured output if needed.

    # Move back to the originaldirectory.
    os.chdir(original_directory)


def create_ini_file(options):
    """Create the gamhelio .ini file from a template.

    Create the gamhelio .ini file from a template.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Returns
    -------
    ini_file : str
        Path to the .ini file for the gamhelio run.
    """
    # Read the template.
    with open(INI_TEMPLATE) as t:
        lines = t.readlines()

    # Process the template here.

    # Write out the .ini file.
    ini_file = os.path.join(
        options["run_directory"], f"{options['runid']}.ini"
    )
    with open(ini_file, "w") as f:
        f.writelines(lines)

    # Return the path to the .ini file.
    return ini_file


def convert_ini_to_xml(options, ini_file):
    """Convert the .ini file to XML.

    Convert the .ini file describing the run to an XML file.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.
    ini_file : str
        Path to the .ini file to convert.

    Returns
    -------
    xml_file : str
        Path to the resulting XML file.
    """
    # Put the XML file in the same directory as the .ini file.
    xml_file = os.path.join(
        options["run_directory"], f"{options['runid']}.xml"
    )

    # Convert the .ini file to .xml.
    # NOTE: assumes XMLGenerator.py is in PATH.
    cmd = "XMLGenerator.py"
    args = [cmd, ini_file, xml_file]
    subprocess.run(
        args, check=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    # Print captured output if needed.

    # Return the path to the XML file.
    return xml_file


def create_pbs_job_script(options):
    """Create the PBS job script for the run.

    Create the PBS job script from a template.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Returns
    -------
    pbs_script : str
        Path to PBS job script.
    """
    # Put the PBS script in the run directory.
    pbs_script = os.path.join(
        options["run_directory"], f"{options['runid']}.pbs"
    )

    # Read the PBS script template.
    with open(pbs_template) as t:
        lines = t.readlines()

    # Process the template here.

    # Write the PBS job script.
    with open(pbs_script, "w") as f:
        f.writelines(lines)

    # Return the path to the PBS script.
    return pbs_script


def main():
    """Main program code for makeitso.

    This is the main program code for makeitso.

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
    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    verbose = args.verbose
    if debug:
        print(f"args = {args}")

    # Fetch the run options.
    options = get_run_options()
    if debug:
        print(f"options = {options}")

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps(options)

    # Create the .ini file.
    if verbose:
        print("Creating .ini file for run.")
    ini_file = create_ini_file(options)
    if debug:
        print(f"ini_file = {ini_file}")

    # Convert the .ini file to a .xml file.
    if verbose:
        print("Converting .ini file to .xml file.")
    xml_file = convert_ini_to_xml(options, ini_file)
    if debug:
        print(f"xml_file = {xml_file}")

    # Create the PBS job script.
    if verbose:
        print("Creating PBS job script for run.")
    pbs_script = create_pbs_job_script(options)
    if verbose:
        print(f"The PBS job script {pbs_script} is ready.")


if __name__ == "__main__":
    main()
