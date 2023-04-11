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
    os.environ["KAIJUHOME"], "kaipy", "gamhelio", "ConfigScripts",
    "startup.config"
)

# Location of template .ini file.
ini_template = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "makeitso",
    f"{DEFAULT_RUNID}.ini.template"
)

# Location of template PBS script.
# pbs_template = os.path.join(
#     os.environ["KAIJUHOME"], "quickstart", default_runid, "%s.pbs.template"
#     % default_runid
# )


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
    args = [options["config_file"]]
    subprocess.run([cmd] + args)

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
    # Just use the template for now.
    with open(ini_template) as t:
        lines = t.readlines()
    # Process the template here.
    ini_file = os.path.join(directory, "%s.ini" % runid)
    with open(ini_file, "w") as f:
        f.writelines(lines)
    return ini_file


# def convert_ini_to_xml(ini_file, xml_file):
#     """Convert the .ini file to XML.
    
#     Convert the .ini file describing the helio_mpi run to the corresponding
#     XML file.

#     Parameters
#     ----------
#     ini_file : str
#         Path to the .ini file to convert.
#     xml_file : str
#         Path to the resulting XML file.
    
#     Returns
#     -------
#     None
#     """
#     # cmd = "XMLGenerator.py"
#     # args = [ini_file, xml_file]
#     # subprocess.run([cmd] + args)

#     # No conversion is performed yet. Just process the XML template.
#     with open(xml_template) as t:
#         lines = t.readlines()
#     # Process the template here.
#     with open(xml_file, "w") as f:
#         f.writelines(lines)


# def create_pbs_job_script(directory, runid):
#     """Create the PBS job script for the run.

#     Create the PBS job script from a template.

#     Parameters
#     ----------
#     directory : str
#         Path to directory to contain PBS job script.
#     runid : str
#         ID string for model to run.

#     Returns
#     -------
#     pbs_file : str
#         Path to PBS job script.
#     """
#     with open(pbs_template) as t:
#         lines = t.readlines()
#     # Process the template here.
#     pbs_file = os.path.join(directory, "%s.pbs" % runid)
#     with open(pbs_file, "w") as f:
#         f.writelines(lines)
#     return pbs_file


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

    # Initialize the dictionary of program options.
    options = {}

    # Get the run ID.
    runid = input(f"Enter the name of the run ({DEFAULT_RUNID}): ")
    if runid is "":
        runid = DEFAULT_RUNID
    if debug:
        print(f"runid = {runid}")
    options["runid"] = runid

    # Enter the working directory for the run.
    run_directory = input(
        f"Enter the path to the run directory ({DEFAULT_RUN_DIRECTORY}): "
    )
    if run_directory is "":
        run_directory = DEFAULT_RUN_DIRECTORY
    if debug:
        print(f"run_directory = {run_directory}")
    options["run_directory"] = run_directory

    # Specify the HPC system to use.
    hpc_system = input(
        f"Enter the name of the HPC system to use ({DEFAULT_HPC_SYSTEM}): "
    )
    if hpc_system is "":
        hpc_system = DEFAULT_HPC_SYSTEM
    if debug:
        print(f"hpc_system = {hpc_system}")
    options["hpc_system"] = hpc_system

    # Specify the path to the location of the configuration file for
    # wsa2gamera.py.
    config_file = input(
        f"Enter the path to the wsa2gamera.py configuration file to use "
        f"({DEFAULT_CONFIG_FILE}): "
    )
    if config_file is "":
        config_file = DEFAULT_CONFIG_FILE
    if debug:
        print(f"config_file = {config_file}")
    options["config_file"] = config_file

    # Summarize the inputs.
    if verbose:
        print("Summary of input:")
        print(f"  options = {options}")

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
    # if verbose:
    #     print("Converting .ini file to .xml file for run.")
    # xml_file = os.path.join(directory, "%s.xml" % runid)
    # convert_ini_to_xml(ini_file, xml_file)

    # # Create the PBS job script.
    # if verbose:
    #     print("Creating PBS job script for run.")
    # pbs_file = create_pbs_job_script(directory, runid)
    # if verbose:
    #     print("The PBS job script %s is ready." % pbs_file)
    #     print("Submit the job to PBS with the command:")
    #     print("    qsub %s" % pbs_file)


if __name__ == "__main__":
    main()
