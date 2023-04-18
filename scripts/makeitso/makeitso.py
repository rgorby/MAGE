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
from jinja2 import Template

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = "Interactive script to prepare and run a kaiju job."

# Default run ID string.
DEFAULT_RUNID = "helio"

# Default directory to run in.
# DEFAULT_RUN_DIRECTORY = "."
DEFAULT_RUN_DIRECTORY = "/Users/winteel1/cgs/runs/test/makeitso"

# Default HPC system.
DEFAULT_HPC_SYSTEM = "pleiades"

# Default run type.
DEFAULT_RUN_TYPE = "serial"

# Default number of hours after spinup to simulate.
DEFAULT_TFIN = "200.0"

# Default screen output interval (in timesteps).
DEFAULT_TSOUT = "50"

# Default simulated time interval between outputs of results as "Step#nn"
# groups in the output HDF5 files. For gamhelio, units are simulated hours.
DEFAULT_DTOUT = "10.0"

# Default MPI decomposition in i-dimension.
DEFAULT_IPDIR_N = 2

# Default MPI decomposition in j-dimension.
DEFAULT_JPDIR_N = 2

# Default MPI decomposition in k-dimension.
DEFAULT_KPDIR_N = 2

# Default spinup time.
DEFAULT_TSPIN = 200.0

# Path to WSA FITS file for creating boundary conditions.
DEFAULT_WSAFILE = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "makeitso",
    "vel_201708132000R002_ahmi.fits"
)

# Location of template wsa2gamera.py .ini file.
WSA2GAMERA_INI_TEMPLATE = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "makeitso", "wsa2gamera_template.ini"
)

# Location of template gamhelio.x .ini file.
GAMHELIO_INI_TEMPLATE = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "makeitso",
    f"{DEFAULT_RUNID}_template.ini"
)

# Location of template PBS script.
GAMHELIO_PBS_TEMPLATE = os.path.join(
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

    # The following parameters are HPC-system-specific.

    # Specify the HPC system to use.
    hpc_system = input(
        f"Enter the name of the HPC system to use ({DEFAULT_HPC_SYSTEM}): "
    )
    if hpc_system == "":
        hpc_system = DEFAULT_HPC_SYSTEM
    options["hpc_system"] = hpc_system

    # Specify the run type (MPI or serial).
    run_type = input(
        f"Specify the run type (serial or mpi) ({DEFAULT_RUN_TYPE}): "
    )
    if run_type == "":
        run_type = DEFAULT_RUN_TYPE
    options["run_type"] = run_type

    # Enter the working directory for the run.
    run_directory = input(
        f"Enter the path to the run directory ({DEFAULT_RUN_DIRECTORY}): "
    )
    if run_directory == "":
        run_directory = DEFAULT_RUN_DIRECTORY
    options["run_directory"] = run_directory

    #-------------------------------------------------------------------------

    # The following parameters are run-specific.

    #-------------------------------------------------------------------------

    # Parameters for the [sim] section of the .ini file.

    # Get the run ID.
    runid = input(f"Enter the name of the run ({DEFAULT_RUNID}): ")
    if runid == "":
        runid = DEFAULT_RUNID
    options["runid"] = runid

    #-------------------------------------------------------------------------

    # Parameters for the [time] section of the .ini file.

    # Specify the number of hours to simulate after spinup.
    tFin = input(
        f"Specify number of hours after spinup to be simulated ({DEFAULT_TFIN}): "
    )
    if tFin == "":
        tFin = DEFAULT_TFIN
    options["tFin"] = tFin

    #-------------------------------------------------------------------------

    # Parameters for the [spinup] section of the .ini file.

    # Specify the number of simulated hours to spinup the simulation.
    tSpin = input(
        f"Specify number of hours of spinup to be simulated ({DEFAULT_TSPIN}): "
    )
    if tSpin == "":
        tSpin = DEFAULT_TSPIN
    options["tSpin"] = tSpin

    # Start output at spinup start.
    options["tIO"] = f"-{tSpin}"

    #-------------------------------------------------------------------------

    # Parameters for the [output] section of the .ini file.

    # Specify the screen output interval. This is the interval (in timesteps)
    # between screen dumps of information during the simulation.
    tsOut = input(
        f"Specify screen output interval (in timesteps) ({DEFAULT_TSOUT}): "
    )
    if tsOut == "":
        tsOut = DEFAULT_TSOUT
    options["tsOut"] = tsOut

    # Specify the simulation step output interval. This is the time interval
    # between outputs of data slices in "Step#nn" groups in the output HDF5
    # file. For gamhelio, the units of dtOut are *hours*. For gamera, the
    # units of dtOut are *seconds*.
    dtOut = input(
        f"Specify step output interval (in simulated hours) ({DEFAULT_DTOUT}): "
    )
    if dtOut == "":
        dtOut = DEFAULT_DTOUT
    options["dtOut"] = dtOut

    #-------------------------------------------------------------------------

    # Parameters for the [physics] section of the .ini file.

    #-------------------------------------------------------------------------

    # Parameters for the [prob] section of the .ini file.

    #-------------------------------------------------------------------------

    # Parameters for the [restart] section of the .ini file.

    #-------------------------------------------------------------------------

    # Parameters for the [iPdir] section of the .ini file.
    iPdir_N = input(
        f"Specify the nuber of MPI chunks in the i-dimension: ({DEFAULT_IPDIR_N}): "
    )
    if iPdir_N == "":
        iPdir_N = DEFAULT_IPDIR_N
    options["iPdir_N"] = iPdir_N

    #-------------------------------------------------------------------------

    # Parameters for the [jPdir] section of the .ini file.
    jPdir_N = input(
        f"Specify the nuber of MPI chunks in the k-dimension: ({DEFAULT_KPDIR_N}): "
    )
    if jPdir_N == "":
        jPdir_N = DEFAULT_KPDIR_N
    options["jPdir_N"] = jPdir_N

    #-------------------------------------------------------------------------

    # Parameters for the [kPdir] section of the .ini file.
    kPdir_N = input(
        f"Specify the nuber of MPI chunks in the j-dimension: ({DEFAULT_KPDIR_N}): "
    )
    if kPdir_N == "":
        kPdir_N = DEFAULT_KPDIR_N
    options["kPdir_N"] = kPdir_N

    #-------------------------------------------------------------------------

    # Parameters for the .ini file for wsa2gamera.py..

    # Specify the path to the WSA FITS file to use for initial conditions.
    wsafile = input(
        f"Path to WSA FITS file for initial conditions ({DEFAULT_WSAFILE}): "
    )
    if wsafile == "":
        wsafile = DEFAULT_WSAFILE
    options["wsafile"] = wsafile

    #-------------------------------------------------------------------------

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

    # Read and create the template, then render and write it.
    with open(WSA2GAMERA_INI_TEMPLATE) as f:
        template_content = f.read()
    template = Template(template_content)
    ini_content = template.render(options)
    ini_file = os.path.join(options["run_directory"], "wsa2gamera.ini")
    with open(ini_file, "w") as f:
        f.write(ini_content)

    # Create the grid and inner boundary conditions files.
    # NOTE: Assumes wsa2gamera.py is in PATH.
    cmd = "wsa2gamera.py"
    args = [cmd, "wsa2gamera.ini"]
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
    # Read and create the template, then render and write it.
    with open(GAMHELIO_INI_TEMPLATE) as f:
        template_content = f.read()
    template = Template(template_content)
    ini_content = template.render(options)
    ini_file = os.path.join(
        options["run_directory"], f"{options['runid']}.ini"
    )
    with open(ini_file, "w") as f:
        f.write(ini_content)

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
    # Read and create the template, then render and write it.
    with open(GAMHELIO_PBS_TEMPLATE) as f:
        template_content = f.read()
    template = Template(template_content)
    ini_content = template.render(options)
    pbs_script = os.path.join(
        options["run_directory"], f"{options['runid']}.pbs"
    )
    with open(pbs_script, "w") as f:
        f.write(ini_content)

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
