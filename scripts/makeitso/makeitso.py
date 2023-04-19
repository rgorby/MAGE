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


# Program defaults

# Default and valid HPC systems.
DEFAULT_HPC_SYSTEM = "cheyenne"
VALIDS_HPC_SYSTEM = ["cheyenne", "pleiades"]

# Default run type.
DEFAULT_RUN_TYPE = "mpi"
VALIDS_RUN_TYPE = ["mpi", "serial"]

# Global defaults across all HPC systems and run types.
defaults_global = {
    "run_directory": ".",
    "runid": "gamhelio",
    "tFin": "200.0",
    "tSpin": "200.0",
    "dtOut": "10.0",
    "tsOut": "50",
    "iPdir_N": "2",
    "jPdir_N": "1",
    "kPdir_N": "2",
}

# Defaults for serial runs on cheyenne.
defaults_cheyenne_serial = {
    # "pbs_walltime": "12:00:00",
    # "pbs_queue": "regular",
    # "pbs_account": "UJHB0015",
    # "select": "1",
    # "ncpus": "36",
    # "ompthreads": "36",
}

# Defaults for MPI runs on cheyenne.
defaults_cheyenne_mpi = {
    "pbs_walltime": "12:00:00",
    # "pbs_queue": "regular",
    # "pbs_account": "UJHB0015",
    # "select": "4",
    # "ncpus": "28",
    # "ompthreads": "18",
    # "dtOut": "10.0",
}

# Defaults for serial runs on pleiades.
defaults_pleiades_serial = {
    # "pbs_walltime": "12:00:00",
    # "pbs_queue": "normal",

    # "pbs_account": "UNKNOWN",
    # "select": "1",
    # "ncpus": "28",
    # "ompthreads": "28",
}

# Defaults for MPI runs on pleiades.
defaults_pleiades_mpi = {
    # "pbs_walltime": "12:00:00",
    # "pbs_queue": "normal",
    # "pbs_account": "UNKNOWN",
    # "select": "4",
    # "ncpus": "28",
    # "ompthreads": "14",
}

# Gather all defaults in one dictionary.
all_defaults = {
    "cheyenne": {
        "mpi": defaults_cheyenne_mpi,
        "serial": defaults_cheyenne_serial,
    },
    "pleiades": {
        "mpi": defaults_pleiades_mpi,
        "serial": defaults_pleiades_serial,
    },
}

# Location of templates of .ini and .pbs files for gamhelio.x.
# template_root = os.path.join(os.environ["KAIJUHOME"], "scripts", "makeitso")
# GAMHELIO_INI_TEMPLATES = {
#     "cheyenne": {
#         "mpi": os.path.join(template_root, "cheyenne_mpi_gamhelio_template.ini"),
#         "serial": os.path.join(template_root, "cheyenne_serial_gamhelio_template.ini"),
#     },
#     "pleiades": {
#         "mpi": os.path.join(template_root, "pleiades_mpi_gamhelio_template.ini"),
#         "serial": os.path.join(template_root, "pleiades_serial_gamhelio_template.ini"),
#     },
# }
# GAMHELIO_PBS_TEMPLATES = {
#     "cheyenne": {
#         "mpi": os.path.join(template_root, "cheyenne_mpi_gamhelio_template.pbs"),
#         "serial": os.path.join(template_root, "cheyenne_serial_gamhelio_template.pbs"),
#     },
#     "pleiades": {
#         "mpi": os.path.join(template_root, "pleiades_mpi_gamhelio_template.pbs"),
#         "serial": os.path.join(template_root, "pleiades_serial_gamhelio_template.pbs"),
#     },
# }

# # Path to WSA FITS file for creating boundary conditions.
# DEFAULT_WSAFILE = os.path.join(
#     os.environ["KAIJUHOME"], "scripts", "makeitso",
#     "vel_201708132000R002_ahmi.fits"
# )

# # Location of template wsa2gamera.py .ini file.
# WSA2GAMERA_INI_TEMPLATE = os.path.join(
#     os.environ["KAIJUHOME"], "scripts", "makeitso", "wsa2gamera_template.ini"
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


def get_run_option(
        prompt="Enter value for option:",
        name=None, valids=None, default=None
):
    """Prompt the user for a single run option.

    Prompt the user for a single run option.

    Parameters
    ----------
    prompt : str
        Prompt string for option
    name : str, default None
        Name of option
    default : str, default None
        Default value for option, as a string
    valids : list of str, default None
        Valid values for option

    Returns
    -------
    option : str
        Value of option as a string.
    """
    s = prompt
    if valids is not None:
        vs = "|".join(valids)
        s += f" ({vs})"
    if default is not None:
        s += f" (default {default})"
    option = input(f"{s}: ")
    if option == "":
        option = default
    if valids is not None and option not in valids:
        raise TypeError(f"Invalid value for option {name}: {option}!")
    return option


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

    # Fetch the options that determine the sets of defaults and templates to
    # use.

    # Specify the HPC system.
    options["hpc_system"] = get_run_option(
        name="hpc_system",
        prompt="Name of HPC system",
        valids=VALIDS_HPC_SYSTEM,
        default=DEFAULT_HPC_SYSTEM
    )

    # Specify the run type.
    options["run_type"] = get_run_option(
        name="run_type",
        prompt="Run type",
        valids=VALIDS_RUN_TYPE,
        default=DEFAULT_RUN_TYPE
    )

    # Start with the global defaults.
    defaults = defaults_global

    # Select the remaining defaults based on HPC system and run type.
    defaults.update(all_defaults[options["hpc_system"]][options["run_type"]])

    #-------------------------------------------------------------------------

    # The following parameters are run-specific.

    # Working directory for the run
    options["run_directory"] = get_run_option(
        name="run_directory",
        prompt="Run directory",
        default=defaults["run_directory"]
    )

    #-------------------------------------------------------------------------

    # Strings [A]B are the names of sections (A) and parameters (B) in the
    # .ini file that will be created from the template for this HPC system and
    # run type.

    # [sim] runid: Run ID string
    options["runid"] = get_run_option(
        name="runid",
        prompt="Run ID",
        default=defaults["runid"]
    )

    # [time] tFin: Number of hours to simulate after spinup
    options["tFin"] = get_run_option(
        name="tFin",
        prompt="Time duration to simulate (hours)",
        default=defaults["tFin"]
    )

    # [spinup] tSpin: Number of simulated hours for spinup time
    options["tSpin"] = get_run_option(
        name="tSpin",
        prompt="Spinup time duration (simulated hours)",
        default=defaults["tSpin"]
    )

    # [spinup] tIO: Simulated time (hours) to start screen output.
    options["tIO"] = f"-{options['tSpin']}"

    # [output] dtOut: Screen output interval (timesteps)
    options["dtOut"] = get_run_option(
        name="dtOut",
        prompt="Screen output interval (timesteps)",
        default=defaults["dtOut"]
    )

    # [output] tsOut: Time step output interval to HDF5 (simulated hours)
    options["tsOut"] = get_run_option(
        name="tsOut",
        prompt="Timestep slice output interval (simulated hours)",
        default=defaults["tsOut"]
    )

    #-------------------------------------------------------------------------

    # Parameters specific to MPI runs.

    if options["run_type"] == "mpi":
        # [iPdir] N: Number of MPI chunks in i-dimension
        options["iPdir_N"] = get_run_option(
            name="iPdir_N",
            prompt="Number of MPI chunks in i-dimension",
            default=defaults["iPdir_N"]
        )
        # [jPdir] N: Number of MPI chunks in j-dimension
        options["jPdir_N"] = get_run_option(
            name="jPdir_N",
            prompt="Number of MPI chunks in j-dimension",
            default=defaults["jPdir_N"]
        )
        # [kPdir] N: Number of MPI chunks in k-dimension
        options["kPdir_N"] = get_run_option(
            name="kPdir_N",
            prompt="Number of MPI chunks in k-dimension",
            default=defaults["kPdir_N"]
        )
    elif options["run_type"] == "serial":
        pass
    else:
        raise TypeError(f"Invalid run type: {options['run_type']}!")

    #-------------------------------------------------------------------------

    # PBS job parameters

    # Requested wall time as hh:mm:ss.
    options["pbs_walltime"] = get_run_option(
        name="pbs_walltime",
        prompt="PBS walltime request (hh:mm:ss)",
        default=defaults["pbs_walltime"]
    )

    # pbs_walltime = input(
    #     f"Specify the walltime request (hh:mm:ss): ({defaults['pbs_walltime']}): "
    # )
    # if pbs_walltime == "":
    #     pbs_walltime = defaults['pbs_walltime']
    # options["pbs_walltime"] = pbs_walltime

    # # PBS queue name.
    # pbs_queue = input(
    #     f"Specify the PBS queue: ({defaults['pbs_queue']}): "
    # )
    # if pbs_queue == "":
    #     pbs_queue = defaults['pbs_queue']
    # options["pbs_queue"] = pbs_queue

    # # PBS account name.
    # pbs_account = input(
    #     f"Specify PBS account: ({defaults['pbs_account']}): "
    # )
    # if pbs_account == "":
    #     pbs_account = defaults['pbs_account']
    # options["pbs_account"] = pbs_account

    # # Number of nodes to use.
    # pbs_select = input(
    #     f"Specify the number of nodes to use: ({defaults['select']}): "
    # )
    # if pbs_select == "":
    #     pbs_select = defaults['select']
    # options["pbs_select"] = pbs_select

    # # Number of cores per node.
    # pbs_ncpus = input(
    #     f"Specify the number of cores per node: ({defaults['ncpus']}): "
    # )
    # if pbs_ncpus == "":
    #     pbs_ncpus = defaults['ncpus']
    # options["pbs_ncpus"] = pbs_ncpus

    # # Number of OMP threads per MPI rank.
    # pbs_ompthreads = input(
    #     f"Specify the number of OMP threads per MPI rank: ({defaults['ompthreads']}): "
    # )
    # if pbs_ompthreads == "":
    #     pbs_ompthreads = defaults['ompthreads']
    # options["pbs_ompthreads"] = pbs_ompthreads

    # #-------------------------------------------------------------------------

    # # Parameters for the .ini file for wsa2gamera.py.

    # # Specify the path to the WSA FITS file to use for initial conditions.
    # wsafile = input(
    #     f"Path to WSA FITS file for initial conditions ({DEFAULT_WSAFILE}): "
    # )
    # if wsafile == "":
    #     wsafile = DEFAULT_WSAFILE
    # options["wsafile"] = wsafile

    #-------------------------------------------------------------------------

    # Return the options dictionary.
    return options


# def run_preprocessing_steps(options):
#     """Execute any preprocessing steps required for the run.

#     Execute any preprocessing steps required for the run.

#     Parameters
#     ----------
#     options : dict
#         Dictionary of program options, each entry maps str to str.

#     Returns
#     -------
#     None
#     """
#     # Save the current directory.
#     original_directory = os.getcwd()

#     # Move to the output directory.
#     os.chdir(options["run_directory"])

#     # Read and create the template, then render and write it.
#     with open(WSA2GAMERA_INI_TEMPLATE) as f:
#         template_content = f.read()
#     template = Template(template_content)
#     ini_content = template.render(options)
#     ini_file = os.path.join(options["run_directory"], "wsa2gamera.ini")
#     with open(ini_file, "w") as f:
#         f.write(ini_content)

#     # Create the grid and inner boundary conditions files.
#     # NOTE: Assumes wsa2gamera.py is in PATH.
#     cmd = "wsa2gamera.py"
#     args = [cmd, "wsa2gamera.ini"]
#     subprocess.run(
#         args, check=True,
#         stdout=subprocess.PIPE, stderr=subprocess.STDOUT
#     )
#     # Print captured output if needed.

#     # Move back to the originaldirectory.
#     os.chdir(original_directory)


# def create_ini_file(options):
#     """Create the gamhelio .ini file from a template.

#     Create the gamhelio .ini file from a template.

#     Parameters
#     ----------
#     options : dict
#         Dictionary of program options, each entry maps str to str.

#     Returns
#     -------
#     ini_file : str
#         Path to the .ini file for the gamhelio run.
#     """
#     # Read and create the template, then render and write it.
#     with open(GAMHELIO_INI_TEMPLATE) as f:
#         template_content = f.read()
#     template = Template(template_content)
#     ini_content = template.render(options)
#     ini_file = os.path.join(
#         options["run_directory"], f"{options['runid']}.ini"
#     )
#     with open(ini_file, "w") as f:
#         f.write(ini_content)

#     # Return the path to the .ini file.
#     return ini_file


# def convert_ini_to_xml(options, ini_file):
#     """Convert the .ini file to XML.

#     Convert the .ini file describing the run to an XML file.

#     Parameters
#     ----------
#     options : dict
#         Dictionary of program options, each entry maps str to str.
#     ini_file : str
#         Path to the .ini file to convert.

#     Returns
#     -------
#     xml_file : str
#         Path to the resulting XML file.
#     """
#     # Put the XML file in the same directory as the .ini file.
#     xml_file = os.path.join(
#         options["run_directory"], f"{options['runid']}.xml"
#     )

#     # Convert the .ini file to .xml.
#     # NOTE: assumes XMLGenerator.py is in PATH.
#     cmd = "XMLGenerator.py"
#     args = [cmd, ini_file, xml_file]
#     subprocess.run(
#         args, check=True,
#         stdout=subprocess.PIPE, stderr=subprocess.STDOUT
#     )
#     # Print captured output if needed.

#     # Return the path to the XML file.
#     return xml_file


# def create_pbs_job_script(options):
#     """Create the PBS job script for the run.

#     Create the PBS job script from a template.

#     Parameters
#     ----------
#     options : dict
#         Dictionary of program options, each entry maps str to str.

#     Returns
#     -------
#     pbs_script : str
#         Path to PBS job script.
#     """
#     # Read and create the template, then render and write it.
#     with open(GAMHELIO_PBS_TEMPLATE) as f:
#         template_content = f.read()
#     template = Template(template_content)
#     ini_content = template.render(options)
#     pbs_script = os.path.join(
#         options["run_directory"], f"{options['runid']}.pbs"
#     )
#     with open(pbs_script, "w") as f:
#         f.write(ini_content)

#     # Return the path to the PBS script.
#     return pbs_script


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

    # # Run the preprocessing steps.
    # if verbose:
    #     print("Running preprocessing steps.")
    # run_preprocessing_steps(options)

    # # Create the .ini file.
    # if verbose:
    #     print("Creating .ini file for run.")
    # ini_file = create_ini_file(options)
    # if debug:
    #     print(f"ini_file = {ini_file}")

    # # Convert the .ini file to a .xml file.
    # if verbose:
    #     print("Converting .ini file to .xml file.")
    # xml_file = convert_ini_to_xml(options, ini_file)
    # if debug:
    #     print(f"xml_file = {xml_file}")

    # # Create the PBS job script.
    # if verbose:
    #     print("Creating PBS job script for run.")
    # pbs_script = create_pbs_job_script(options)
    # if verbose:
    #     print(f"The PBS job script {pbs_script} is ready.")


if __name__ == "__main__":
    main()
