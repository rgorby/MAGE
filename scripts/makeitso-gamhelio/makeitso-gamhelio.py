#!/usr/bin/env python


"""makeitso for the MAGE heliosphere software.

This script performs all of the steps needed to prepare to run a MAGE
heliosphere simulation run. By default, this script is interactive - the user
is prompted for each decision that must be made to prepare for the run, based
on the current "--mode" setting.

The modes are:

"BASIC" (the default) - the user is prompted to set only a small subset of MAGE
parameters. All "INTERMEDIATE"- and "EXPERT"-mode parameters are automatically
set to default values.

"INTERMEDIATE" - The user is prompted for "BASIC" and "INTERMEDIATE"
parameters, with "EXPERT" parameters set to defaults.

"EXPERT" - The user is prompted for *all* adjustable parameters.
"""


# Import standard modules.
import argparse
import datetime
import json
import os
import subprocess

# Import 3rd-party modules.
from astropy.io import fits
from jinja2 import Template

# Import project modules.
from kaipy.kdefs import JD2MJD


# Program constants

# Program description.
DESCRIPTION = "Interactive script to prepare and run a MAGE heliosphere job."

# Indent level for JSON output.
JSON_INDENT = 4

# Path to current kaiju installation
KAIJUHOME = os.environ["KAIJUHOME"]

# Path to directory containing support files for makeitso.
SUPPORT_FILES_DIRECTORY = os.path.join(KAIJUHOME, "scripts",
                                       "makeitso-gamhelio")

# Path to option descriptions file.
OPTION_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_descriptions.json"
)

# Location of template wsa2gamera.py .ini file.
WSA2GAMERA_INI_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY,
                                       "wsa2gamera_template.ini")

# # Global defaults across all HPC systems and run types.
# defaults_global = {
#     "run_directory": ".",
#     "wsafile": DEFAULT_WSAFILE,
#     "kaiju_home": DEFAULT_KAIJU_HOME,
#     "sim_runid": "gamhelio",
#     "time_tFin": "200.0",
#     "spinup_tSpin": "200.0",
#     "output_dtOut": "10.0",
#     "output_tsOut": "50",
#     "iPdir_N": "2",
#     "jPdir_N": "2",
#     "kPdir_N": "2",
#     "wsa2gamera_Grid_Ni": "128",
#     "wsa2gamera_Grid_Nj": "64",
#     "wsa2gamera_Grid_Nk": "128",
# }

# # Defaults for MPI runs on cheyenne.
# defaults_cheyenne_mpi = {
#     "pbs_account": "",
#     "pbs_queue": "regular",
#     "pbs_walltime": "00:30:00",
#     "pbs_select": "4",
#     "pbs_ncpus": "36",
#     "pbs_mpiprocs": "2",
#     "pbs_ompthreads": "18",
# }

# # Defaults for serial runs on cheyenne.
# defaults_cheyenne_serial = {
#     "pbs_account": "",
#     "pbs_queue": "regular",
#     "pbs_walltime": "02:00:00",
#     "pbs_select": "1",
#     "pbs_ncpus": "36",
#     "pbs_ompthreads": "36",
# }

# # Defaults for MPI runs on pleiades.
# defaults_pleiades_mpi = {
#     "pbs_queue": "normal",
#     "pbs_walltime": "00:30:00",
#     "pbs_select": "4",
#     "pbs_ncpus": "28",
#     "pbs_mpiprocs": "2",
#     "pbs_ompthreads": "14",
# }

# # Defaults for serial runs on pleiades.
# defaults_pleiades_serial = {
#     "pbs_queue": "normal",
#     "pbs_walltime": "02:00:00",
#     "pbs_select": "1",
#     "pbs_ncpus": "28",
#     "pbs_ompthreads": "28",
# }

# # Gather all defaults in one dictionary.
# all_defaults = {
#     "cheyenne": {
#         "mpi": defaults_cheyenne_mpi,
#         "serial": defaults_cheyenne_serial,
#     },
#     "pleiades": {
#         "mpi": defaults_pleiades_mpi,
#         "serial": defaults_pleiades_serial,
#     },
# }

# # Location of templates of .ini and .pbs files for gamhelio.x.
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

    Raises
    ------
    None
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--clobber", action="store_true",
        help="Overwrite existing options file (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--mode", default="BASIC",
        help="User mode (BASIC|INTERMEDIATE|EXPERT) (default: %(default)s)."
    )
    parser.add_argument(
        "--options_path", "-o", default=None,
        help="Path to JSON file of options (default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def get_run_option(name, description, mode="BASIC"):
    """Prompt the user for a single run option.

    Prompt the user for a single run option. If no user input is provided,
    the default value is returned for the option. If valids are provided, the
    new value is compared against the valids, and rejected if not in the
    valids list.

    Parameters
    ----------
    name : str, default None
        Name of option
    description : dict, default None
        Dictionary of metadata for the option.
    mode : str
        User experience mode: "BASIC", "INTERMEDIATE", or "ADVANCED".

    Returns
    -------
    value_str : str
        Value of option as a string.

    Raises
    ------
    None
    """
    # Extract prompt, default, and valids.
    level = description["LEVEL"]
    prompt = description.get("prompt", "")
    default = description.get("default", None)
    valids = description.get("valids", None)

    # Compare the current mode to the parameter level setting. If the variable
    # level is higher than the user mode, just use the default.
    if mode == "BASIC" and level in ["INTERMEDIATE", "EXPERT"]:
        return default
    if mode == "INTERMEDIATE" and level in ["EXPERT"]:
        return default

    # If provided, add the valid values in val1|val2 format to the prompt.
    if valids is not None:
        vs = "|".join(valids)
        prompt += f" ({vs})"

    # If provided, add the default to the prompt.
    if default is not None:
        prompt += f" [{default}]"

    # Prompt the user and fetch the input until a good value is provided.
    ok = False
    while not ok:

        # Fetch input from the user.
        option_value = input(f"{prompt}: ")

        # Use the default if no user input provided.
        if option_value == "":
            option_value = default

        # Validate the result. If bad, start over.
        if valids is not None and option_value not in valids:
            print(f"Invalid value for option {name}: {option_value}!")
            continue

        # Keep this result.
        ok = True

    # Return the option as a string.
    return str(option_value)


def read_mjdc_from_fits(wsa_file):
    """Read the MJDc value from a WSA FITS file.

    Read the MJDc value from a WSA FITS file.

    Parameters
    ----------
    wsa_file : str
        Path to WSA FITS file

    Returns
    -------
    mjdc : str
        MJDc value from FITS file, in 'YYYY-MM-DDTHH:MM:SS' format

    Raises
    ------
    None
    """
    # Open the WSA FITS file.
    hdu = fits.open(wsa_file)

    # Read the JD for the center of the map.
    jdc = hdu[0].header['JULDATE']

    # Convert JD to MJD.
    mjdc = jdc - JD2MJD

    # Return the MJDc value.
    return mjdc


def prompt_user_for_run_options(args):
    """Prompt the user for run options.

    Prompt the user for run options.

    NOTE: In this function, the complete set of parameters is split up
    into logical groups. This is done partly to make the organization of the
    parameters more obvious, and partly to allow the values of options to
    depend upon previously-specified options.

    Parameters
    ----------
    args : dict
        Dictionary of command-line options

    Returns
    -------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Raises
    ------
    None
    """
    # Save the user mode.
    mode = args.mode

    # Read the dictionary of option descriptions.
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)

    # Initialize the dictionary of program options.
    options = {}

    # -------------------------------------------------------------------------

    # General options for the simulation
    o = options["simulation"] = {}
    od = option_descriptions["simulation"]

    # Prompt for the name of the job.
    for on in ["job_name"]:
        o[on] = get_run_option(on, od[on], mode)

    # Path to the WSA FITS file to use for initial conditions
    for on in ["wsa_file"]:
        o[on] = get_run_option(on, od[on], mode)

    # Fetch the MJDc value from the WSA FITS file.
    mjdc = read_mjdc_from_fits(o[on])
    o["mjdc"] = mjdc

    # Prompt for the start and stop date of the run.
    for on in ["start_date", "stop_date"]:
        o[on] = get_run_option(on, od[on], mode)

    # Compute the total simulation time in seconds, use as segment duration
    # default.
    date_format = '%Y-%m-%dT%H:%M:%S'
    start_date = o["start_date"]
    stop_date = o["stop_date"]
    t1 = datetime.datetime.strptime(start_date, date_format)
    t2 = datetime.datetime.strptime(stop_date, date_format)
    simulation_duration = (t2 - t1).total_seconds()
    od["segment_duration"]["default"] = str(simulation_duration)

    # Ask if the user wants to split the run into multiple segments.
    # If so, prompt for the segment duration. If not, use the default
    # for the segment duration (which is the simulation duration).
    for on in ["use_segments"]:
        o[on] = get_run_option(on, od[on], mode)
    if o["use_segments"].upper() == "Y":
        for on in ["segment_duration"]:
            o[on] = get_run_option(on, od[on], mode)
    else:
        o["segment_duration"] = od["segment_duration"]["default"]

    # Compute the number of segments based on the simulation duration and
    # segment duration, with 1 additional segment just for spinup. Add 1 if
    # there is a remainder.
    num_segments = simulation_duration/float(o["segment_duration"])
    if num_segments > int(num_segments):
        num_segments += 1
    num_segments = int(num_segments) + 1

    # Prompt for the remaining parameters.
    for on in ["hpc_system"]:
        o[on] = get_run_option(on, od[on], mode)

    # -------------------------------------------------------------------------

    # PBS options
    o = options["pbs"] = {}

    # Common (HPC platform-independent) options
    od = option_descriptions["pbs"]["_common"]
    od["account_name"]["default"] = os.getlogin()
    od["kaiju_install_directory"]["default"] = KAIJUHOME
    od["kaiju_build_directory"]["default"] = os.path.join(KAIJUHOME, "build_mpi")
    od["num_segments"]["default"] = str(num_segments)
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # HPC platform-specific options
    hpc_platform = options["simulation"]["hpc_system"]
    od = option_descriptions["pbs"][hpc_platform]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # -------------------------------------------------------------------------

    # wsa2gamera.py options
    options["wsa2gamera"] = {}

    # [Gamera] options
    o = options["wsa2gamera"]["Gamera"] = {}
    od = option_descriptions["wsa2gamera"]["Gamera"]
    od["GridDir"]["default"] = options["pbs"]["run_directory"]
    od["IbcDir"]["default"] = options["pbs"]["run_directory"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [Grid] options
    o = options["wsa2gamera"]["Grid"] = {}
    od = option_descriptions["wsa2gamera"]["Grid"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [WSA] options
    o = options["wsa2gamera"]["WSA"] = {}
    od = option_descriptions["wsa2gamera"]["WSA"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [Constants] options
    o = options["wsa2gamera"]["Constants"] = {}
    od = option_descriptions["wsa2gamera"]["Constants"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # [Normalization] options
    o = options["wsa2gamera"]["Normalization"] = {}
    od = option_descriptions["wsa2gamera"]["Normalization"]
    for on in od:
        o[on] = get_run_option(on, od[on], mode)

    # -------------------------------------------------------------------------

    print(options)

    # Return the options dictionary.
    return options


# def get_run_options():
#     """Prompt the user for run options.

#     Prompt the user for run options.

#     Parameters
#     ----------
#     None

#     Returns
#     -------
#     options : dict
#         Dictionary of program options, each entry maps str to str.
#     """

#     # Initialize the dictionary of program options.
#     options = {}

#     # Fetch the options that determine the sets of defaults and templates to
#     # use.

#     # Specify the HPC system.
#     options["hpc_system"] = get_run_option(
#         name="hpc_system",
#         prompt="Name of HPC system",
#         valids=VALIDS_HPC_SYSTEM,
#         default=DEFAULT_HPC_SYSTEM
#     )

#     # Specify the run type.
#     options["run_type"] = get_run_option(
#         name="run_type",
#         prompt="Run type",
#         valids=VALIDS_RUN_TYPE,
#         default=DEFAULT_RUN_TYPE
#     )

#     # Start with the global defaults.
#     defaults = defaults_global

#     # Select the remaining defaults based on HPC system and run type.
#     defaults.update(all_defaults[options["hpc_system"]][options["run_type"]])

#     # -------------------------------------------------------------------------

#     # The following parameters are run-specific.

#     # Working directory for the run
#     options["run_directory"] = get_run_option(
#         name="run_directory",
#         prompt="Run directory",
#         default=defaults["run_directory"]
#     )

#     # Path to the WSA FITS file to use for initial conditions
#     options["wsafile"] = get_run_option(
#         name="wsafile",
#         prompt="Path to WSA FITS file for initial conditions",
#         default=defaults["wsafile"]
#     )

#     # Number of grid points in i-dimension
#     options["wsa2gamera_Grid_Ni"] = get_run_option(
#         name="wsa2gamera_Grid_Ni",
#         prompt="Number of grid points in i-dimension",
#         default=defaults["wsa2gamera_Grid_Ni"]
#     )

#     # Number of grid points in j-dimension
#     options["wsa2gamera_Grid_Nj"] = get_run_option(
#         name="wsa2gamera_Grid_Nj",
#         prompt="Number of grid points in j-dimension",
#         default=defaults["wsa2gamera_Grid_Nj"]
#     )

#     # Number of grid points in k-dimension
#     options["wsa2gamera_Grid_Nk"] = get_run_option(
#         name="wsa2gamera_Grid_Nk",
#         prompt="Number of grid points in k-dimension",
#         default=defaults["wsa2gamera_Grid_Nk"]
#     )

#     # Path to kaiju installation to use
#     options["kaiju_home"] = get_run_option(
#         name="kaiju_home",
#         prompt="Path to kaiju installation",
#         default=defaults["kaiju_home"],
#     )

#     # Path to kaiju binaries
#     if options["run_type"] == "mpi":
#         defaults["kaiju_build_bin"] = os.path.join(options["kaiju_home"], "build_mpi", "bin")
#     elif options["run_type"] == "serial":
#         defaults["kaiju_build_bin"] = os.path.join(options["kaiju_home"], "build_serial", "bin")
#     else:
#         raise TypeError(f"Invalid run type: {options['run_type']}!")
#     options["kaiju_build_bin"] = get_run_option(
#         name="kaiju_build_bin",
#         prompt="Path to kaiju build bin/ directory",
#         default=defaults["kaiju_build_bin"]
#     )

#     # -------------------------------------------------------------------------

#     # Strings [A]B are the names of sections (A) and parameters (B) in the
#     # .ini file that will be created from the template for this HPC system and
#     # run type.

#     # [sim] runid: Run ID string
#     options["sim_runid"] = get_run_option(
#         name="sim_runid",
#         prompt="Run ID",
#         default=defaults["sim_runid"]
#     )

#     # [time] tFin: Number of hours to simulate after spinup
#     options["time_tFin"] = get_run_option(
#         name="time_tFin",
#         prompt="Time duration to simulate (hours)",
#         default=defaults["time_tFin"]
#     )

#     # [spinup] tSpin: Number of simulated hours for spinup time
#     options["spinup_tSpin"] = get_run_option(
#         name="spinup_tSpin",
#         prompt="Spinup time duration (simulated hours)",
#         default=defaults["spinup_tSpin"]
#     )

#     # [spinup] tIO: Simulated time (hours) to start screen output.
#     options["spinup_tIO"] = f"-{options['spinup_tSpin']}"

#     # [output] dtOut: Screen output interval (timesteps)
#     options["output_dtOut"] = get_run_option(
#         name="output_dtOut",
#         prompt="Screen output interval (timesteps)",
#         default=defaults["output_dtOut"]
#     )

#     # [output] tsOut: Time step output interval to HDF5 (simulated hours)
#     options["output_tsOut"] = get_run_option(
#         name="output_tsOut",
#         prompt="Timestep slice output interval (simulated hours)",
#         default=defaults["output_tsOut"]
#     )

#     # -------------------------------------------------------------------------

#     # Parameters specific to MPI runs.

#     if options["run_type"] == "mpi":
#         # [iPdir] N: Number of MPI chunks in i-dimension
#         options["iPdir_N"] = get_run_option(
#             name="iPdir_N",
#             prompt="Number of MPI chunks in i-dimension",
#             default=defaults["iPdir_N"]
#         )
#         # [jPdir] N: Number of MPI chunks in j-dimension
#         options["jPdir_N"] = get_run_option(
#             name="jPdir_N",
#             prompt="Number of MPI chunks in j-dimension",
#             default=defaults["jPdir_N"]
#         )
#         # [kPdir] N: Number of MPI chunks in k-dimension
#         options["kPdir_N"] = get_run_option(
#             name="kPdir_N",
#             prompt="Number of MPI chunks in k-dimension",
#             default=defaults["kPdir_N"]
#         )

#     # -------------------------------------------------------------------------

#     # PBS job parameters

#     # PBS account name
#     if "pbs_account" in defaults:
#         options["pbs_account"] = get_run_option(
#             name="pbs_account",
#             prompt="PBS account name",
#             default=defaults["pbs_account"]
#         )

#     # PBS queue name
#     options["pbs_queue"] = get_run_option(
#         name="pbs_queue",
#         prompt="PBS queue name",
#         default=defaults["pbs_queue"]
#     )

#     # Requested wall time as hh:mm:ss
#     options["pbs_walltime"] = get_run_option(
#         name="pbs_walltime",
#         prompt="PBS walltime request (hh:mm:ss)",
#         default=defaults["pbs_walltime"]
#     )

#     # Number of compute nodes to use
#     options["pbs_select"] = get_run_option(
#         name="pbs_select",
#         prompt="Number of compute nodes to use",
#         default=defaults["pbs_select"]
#     )

#     # Number of cores per compute node
#     options["pbs_ncpus"] = get_run_option(
#         name="pbs_ncpus",
#         prompt="Number of cores per compute node",
#         default=defaults["pbs_ncpus"]
#     )

#     # Number of MPI ranks to run on each compute node
#     # Should be the same as the number of CPU sockets in the node.
#     if options["run_type"] == "mpi":
#         options["pbs_mpiprocs"] = get_run_option(
#             name="pbs_mpiprocs",
#             prompt="Number of MPI ranks per compute node",
#             default=defaults["pbs_mpiprocs"]
#         )

#     # Number of OMP threads per MPI rank
#     if options["run_type"] == "mpi":
#         options["pbs_ompthreads"] = get_run_option(
#             name="pbs_ompthreads",
#             prompt="Number of OMP threads per MPI rank",
#             default=defaults["pbs_ompthreads"]
#         )
#     elif options["run_type"] == "serial":
#         options["pbs_ompthreads"] = get_run_option(
#             name="pbs_ompthreads",
#             prompt="Number of OMP threads per node",
#             default=defaults["pbs_ompthreads"]
#         )
#         pass
#     else:
#         raise TypeError(f"Invalid run type: {options['run_type']}!")

#     # -------------------------------------------------------------------------

#     # Return the options dictionary.
#     return options


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
    # Read and create the template for the ini file for wsa2gamera.py, then
    # render and write it.
    with open(WSA2GAMERA_INI_TEMPLATE) as f:
        template_content = f.read()
    template = Template(template_content)
    ini_content = template.render(options)
    ini_file = os.path.join(options["pbs"]["run_directory"], "wsa2gamera.ini")
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


def create_ini_files(options):
    """Create the MAGE .ini files from a template.

    Create the MAGE .ini files from a template.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Returns
    -------
    ini_files : list of str
        Paths to the .ini files for the gamera run.

    Raises
    ------
    None
    """
#     # Read and create the template, then render and write it.
#     template_file = GAMHELIO_INI_TEMPLATES[options["hpc_system"]][options["run_type"]]
#     with open(template_file) as f:
#         template_content = f.read()
#     template = Template(template_content)
#     ini_content = template.render(options)
#     ini_file = os.path.join(
#         options["run_directory"], f"{options['sim_runid']}.ini"
#     )
#     with open(ini_file, "w") as f:
#         f.write(ini_content)

#     # Save the options dictionary as a JSON file.
#     with open("options.json", "w") as f:
#         json.dump(options, f)

#     # Return the path to the .ini file.
#     return ini_file


def convert_ini_to_xml(ini_files):
    """Convert the .ini files to XML.

    Convert the .ini files describing the run to XML files. The intermediate
    .ini files are then deleted.

    Parameters
    ----------
    ini_files : list of str
        Paths to the .ini files to convert.

    Returns
    -------
    xml_files : str
        Paths to the XML files.

    Raises
    ------
    None
    """
    # Convert each .ini file to an .xml file.
    xml_files = []
    for ini_file in ini_files:

        # Put the XML file in the same directory as the .ini file.
        xml_file = ini_file.replace(".ini", ".xml")

        # Convert the .ini file to .xml.
        # NOTE: assumes XMLGenerator.py is in PATH.
        cmd = "XMLGenerator.py"
        args = [cmd, ini_file, xml_file]
        subprocess.run(args, check=True)

        # Add this file to the list of XML files.
        xml_files.append(xml_file)

        # Remove the .ini file.
        os.remove(ini_file)

    # Return the paths to the XML files.
    return xml_files


def create_pbs_scripts(options):
    """Create the PBS job scripts for the run.

    Create the PBS job scripts from a template.

    Parameters
    ----------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Returns
    -------
    pbs_scripts : list of str
        Paths to PBS job script.
    submit_all_jobs_script : str
        Path to script which submits all PBS jobs.

    Raises
    ------
    TypeError:
        For a non-integral of nodes requested
    """
    # # Compute the number of nodes to request based on the MPI decomposition
    # # and the MPI ranks per node.
    # ni = int(options["gamera"]["iPdir"]["N"])
    # nj = int(options["gamera"]["jPdir"]["N"])
    # nk = int(options["gamera"]["kPdir"]["N"])
    # ranks_per_node = int(options["pbs"]["mpiprocs"])
    # select_nodes = ni*nj*nk/ranks_per_node
    # if int(select_nodes) != select_nodes:
    #     raise TypeError(f"Requested non-integral node count ({select_nodes})!")
    # options["pbs"]["select"] = str(int(select_nodes))

    # # Read the template.
    # with open(PBS_TEMPLATE, "r", encoding="utf-8") as f:
    #     template_content = f.read()
    # template = Template(template_content)

    # # Create a PBS script for each segment.
    # pbs_scripts = []
    # for job in range(int(options["pbs"]["num_segments"])):
    #     opt = copy.deepcopy(options)  # Need a copy of options
    #     runid = opt["simulation"]["job_name"]
    #     segment_id = f"{runid}-{job:02d}"
    #     opt["simulation"]["segment_id"] = segment_id
    #     pbs_content = template.render(opt)
    #     pbs_script = os.path.join(
    #         opt["pbs"]["run_directory"],
    #         f"{opt['simulation']['segment_id']}.pbs"
    #     )
    #     pbs_scripts.append(pbs_script)
    #     with open(pbs_script, "w", encoding="utf-8") as f:
    #         f.write(pbs_content)

    # # Create a single script which will submit all of the PBS jobs in order.
    # submit_all_jobs_script = f"{options['simulation']['job_name']}_pbs.sh"
    # with open(submit_all_jobs_script, "w", encoding="utf-8") as f:
    #     s = pbs_scripts[0]
    #     cmd = f"job_id=`qsub {s}`\n"
    #     f.write(cmd)
    #     cmd = f"echo $job_id\n"
    #     f.write(cmd)
    #     for s in pbs_scripts[1:]:
    #         cmd = "old_job_id=$job_id\n"
    #         f.write(cmd)
    #         cmd = f"job_id=`qsub -W depend=afterok:$old_job_id {s}`\n"
    #         f.write(cmd)
    #         cmd = f"echo $job_id\n"
    #         f.write(cmd)

    # # Return the paths to the PBS scripts.
    # return pbs_scripts, submit_all_jobs_script

#     # Read and create the template, then render and write it.
#     template_file = GAMHELIO_PBS_TEMPLATES[options["hpc_system"]][options["run_type"]]
#     with open(template_file) as f:
#         template_content = f.read()
#     template = Template(template_content)
#     ini_content = template.render(options)
#     pbs_script = os.path.join(
#         options["run_directory"], f"{options['sim_runid']}.pbs"
#     )
#     with open(pbs_script, "w") as f:
#         f.write(ini_content)

#     # Return the path to the PBS script.
#     return pbs_script


def main():
    """Main program code for makeitso-gamhelio.

    This is the main program code for makeitso-gamhelio.

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
    if args.debug:
        print(f"args = {args}")
    clobber = args.clobber
    debug = args.debug
    options_path = args.options_path
    verbose = args.verbose

    # Fetch the run options.
    if options_path:
        # Read the run options from a JSON file.
        with open(options_path, "r", encoding="utf-8") as f:
            options = json.load(f)
    else:
        # Prompt the user for the run options.
        options = prompt_user_for_run_options(args)
    if debug:
        print(f"options = {options}")

    # Move to the run directory.
    os.chdir(options["pbs"]["run_directory"])

    # Save the options dictionary as a JSON file in the current directory.
    path = f"{options['simulation']['job_name']}.json"
    if os.path.exists(path):
        if not clobber:
            raise FileExistsError(f"Options file {path} exists!")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps(options)

    # # Create the .ini file(s).
    # if verbose:
    #     print("Creating .ini file(s) for run.")
    # ini_files = create_ini_files(options)
    # if debug:
    #     print(f"ini_files = {ini_files}")

    # # Convert the .ini file(s) to .xml files(s).
    # if verbose:
    #     print("Converting .ini file(s) to .xml file(s).")
    # xml_files = convert_ini_to_xml(ini_files)
    # if debug:
    #     print(f"xml_files = {xml_files}")

    # # Create the PBS job script(s).
    # if verbose:
    #     print("Creating PBS job script(s) for run.")
    # pbs_scripts, all_jobs_script = create_pbs_scripts(options)
    # if verbose:
    #     print(f"The PBS job scripts {pbs_scripts} are ready.")
    # print(f"The PBS scripts {pbs_scripts} have been created, each with a "
    #       "corresponding XML file. To submit the jobs with the proper "
    #       "dependency (to ensure each segment runs in order), please run the "
    #       f"script {all_jobs_script} like this:\n"
    #       f"bash {all_jobs_script}")


if __name__ == "__main__":
    """Begin main program."""
    main()

