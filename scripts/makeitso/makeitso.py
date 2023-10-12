#!/usr/bin/env python


"""makeitso for the MAGE magnetosphere software.

This master script is used to perform all of the steps needed to prepare
and run a kaiju job. This script is interactive - the user is prompted for
each decision that must be made to prepare for the run.
"""


# Import standard modules.
import argparse
import json
import os
import subprocess

# Import 3rd-party modules.
from jinja2 import Template

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = "Interactive script to prepare a MAGE magnetosphere model run."

# Indent level for JSON output.
JSON_INDENT = 4


# Program defaults

# Module sets by platform and run type.
modules = {
    "cheyenne": [
        "cmake/3.22.0",
        "git/2.33.1",
        "ncarenv/1.3",
        "intel/2022.1",
        "geos/3.10.1",
        "ncarcompilers/0.5.0",
        "mpt/2.25",
        "hdf5-mpi/1.12.2",
    ],
    "derecho": [
        "ncarenv/23.06",
        "cmake/3.26.3",
        "craype/2.7.20",
        "intel/2023.0.0",
        "geos/3.9.1",
        "ncarcompilers/1.0.0",
        "cray-mpich/8.1.25",
        "hdf5-mpi/1.12.2",
    ],
    "pleiades": [
        "nas"
        "pkgsrc/2022Q1-rome",
        "comp-intel/2020.4.304",
        "mpi-hpe/mpt.2.23",
        "hdf5/1.8.18_mpt",
    ],
}

# Path to directory containing support files.
SUPPORT_FILES_DIRECTORY = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "makeitso"
)

# Path to option descriptions file.
OPTION_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_descriptions.json"
)

# Path to template .ini file.
INI_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.ini")

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.pbs")


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
        "--debug", "-d", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--options", "-o", default=None,
        help="Path to JSON file of options (default: %(default)s)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def get_run_option(name, description):
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

    Returns
    -------
    value_str : str
        Value of option as a string.

    Raises
    ------
    None
    """
    # Extract prompt, default, and valids.
    prompt = description.get("prompt", "")
    default = description.get("default", None)
    valids = description.get("valids", None)

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


def get_run_options():
    """Prompt the user for run options.

    Prompt the user for run options.

    NOTE: In this function, the complete set of parameters is split up
    into logical groups. This is done partly to make the organization of the
    parameters more obvious, and partly to allow the values of options to
    depend upon previously-specified options.

    Parameters
    ----------
    None

    Returns
    -------
    options : dict
        Dictionary of program options, each entry maps str to str.

    Raises
    ------
    None
    """
    # Read the dictionary of option descriptions.
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        options_yaml = json.load(f)
        option_descriptions = options_yaml["option_descriptions"]

    # Initialize the dictionary of program options.
    options = {}

    # Fetch the HPC system and run type.
    for name in ["hpc_system"]:
        options[name] = get_run_option(name, option_descriptions[name])

    # Specify the location of the kaiju installation to use.
    # The installation is assumed to be the installation which contains
    # this running script. Note that this code requires that the KAIJUHOME
    # environment variable is set. This environment variable is set
    # automatically when the user runs the setupEnvironment.[c]sh script
    # or its equivalent.
    options["kaijuhome"] = os.environ["KAIJUHOME"]

    # Compute the default build directory based on KAIJUHOME and the run
    # type. The build directory is the directory where cmake and make were
    # run during the build process. Typically, this is done in the kaiju
    # subdirectory "build_mpi" (for MPI builds).
    option_descriptions["kaiju_build_directory"]["default"] = (
        f"{options['kaijuhome']}/build_mpi"
    )

    # Fetch options about the kaiju software to use.
    option_names = [
        "kaiju_build_directory"
    ]
    for name in option_names:
        options[name] = get_run_option(name, option_descriptions[name])

    # Fetch high-level options about the run.
    # NOTE: All files must be in the run directory.
    option_names = [
        "run_directory", "runid", "LFM_grid_type", "solar_wind_file"
    ]
    for name in option_names:
        options[name] = get_run_option(name, option_descriptions[name])

    # PBS job options
    # Assumes all nodes have same architecture.
    hpc_system = options["hpc_system"]
    option_names = [
        "pbs_account_name", "pbs_queue", "pbs_walltime",
        "pbs_select", "pbs_ncpus",
        "pbs_mpiprocs", "pbs_ompthreads",
        "pbs_num_helpers", "pbs_numjobs",
    ]
    for name in option_names:
        options[name] = get_run_option(name, option_descriptions[name])

    # The helper nodes use one OMP thread per core.
    options["pbs_helper_ompthreads"] = options["pbs_ncpus"]

    # Specify the command to launch an MPI program.
    options["pbs_mpiexec_command"] = "mpiexec omplace"

    # Assemble the module load statements into a single string.
    hpc_system = options["hpc_system"]
    module_names = modules[hpc_system]
    pbs_module_load = ""
    for module_name in module_names:
        pbs_module_load += f"module load {module_name}\n"
    pbs_module_load = pbs_module_load.rstrip()
    options["pbs_module_load"] = pbs_module_load

    # Specify other environment variables to set.
    options["pbs_additional_environment_variables"] = ""

    # GAMERA section
    option_names = [
        "gamera_sim_doH5g", "gamera_sim_H5Grid", "gamera_sim_icType",
        "gamera_sim_pdmb", "gamera_sim_rmeth",
        "gamera_floors_dFloor", "gamera_floors_pFloor",
        "gamera_timestep_doCPR", "gamera_timestep_limCPR",
        "gamera_restart_doRes", "gamera_restart_resID", "gamera_restart_nRes",
        "gamera_physics_doMHD", "gamera_physics_doBoris", "gamera_physics_Ca",
        "gamera_ring_gid", "gamera_ring_doRing",
        "gamera_ringknobs_doVClean",
        "gamera_wind_tsfile",
        "gamera_source_doSource", "gamera_source_doWolfLim",
        "gamera_source_doBounceDT", "gamera_source_nBounce",
        "gamera_iPdir_N", "gamera_iPdir_bcPeriodic",
        "gamera_jPdir_N", "gamera_jPdir_bcPeriodic",
        "gamera_kPdir_N", "gamera_kPdir_bcPeriodic",
    ]
    for name in option_names:
        options[name] = get_run_option(name, option_descriptions[name])

    # VOLTRON section
    option_names = [
        "voltron_time_tFin",
        "voltron_spinup_doSpin", "voltron_spinup_tSpin",
        "voltron_output_dtOut", "voltron_output_tsOut",
        "voltron_coupling_dtCouple", "voltron_coupling_rTrc",
        "voltron_coupling_imType", "voltron_coupling_doQkSquish",
        "voltron_coupling_qkSquishStride",
        "voltron_restart_dtRes",
        "voltron_imag_doInit",
        "voltron_ebsquish_epsSquish",
    ]
    for name in option_names:
        options[name] = get_run_option(name, option_descriptions[name])

    # CHIMP section
    option_names = [
        "chimp_units_uid",
        "chimp_fields_grType",
        "chimp_domain_dtype",
        "chimp_tracer_epsds",
    ]
    for name in option_names:
        options[name] = get_run_option(name, option_descriptions[name])

    # REMIX section
    option_names = [
        "remix_conductance_doStarlight", "remix_conductance_doRamp",
        "remix_precipitation_aurora_model_type", "remix_precipitation_alpha",
        "remix_precipitation_beta",
    ]
    for name in option_names:
        options[name] = get_run_option(name, option_descriptions[name])

    # RCM section
    option_names = [
        "rcm_rcmdomain_domType",
        "rcm_ellipse_xSun", "rcm_ellipse_yDD", "rcm_ellipse_xTail",
        "rcm_ellipse_isDynamic", "rcm_grid_LowLat", "rcm_grid_HiLat",
        "rcm_grid_doLatStretch",
        "rcm_plasmasphere_isDynamic", "rcm_plasmasphere_initKp",
        "rcm_plasmasphere_doRefill", "rcm_plasmasphere_DenPP0",
        "rcm_plasmasphere_tAvg",
    ]
    for name in option_names:
        options[name] = get_run_option(name, option_descriptions[name])

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

    Raises
    ------
    None
    """
    # Create the LFM grid file.
    # NOTE: Assumes genLFM.py is in PATH.
    cmd = "genLFM.py"
    args = [cmd, "-gid", options["LFM_grid_type"]]
    subprocess.run(args, check=True)

    # Create the RCM configuration file.
    # NOTE: Assumes genRCM.py is in PATH.
    cmd = "genRCM.py"
    args = [cmd]
    subprocess.run(args, check=True)


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
    # Read and create the template.
    template_file = INI_TEMPLATE
    with open(template_file, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)

    # Initialize the list of file paths.
    ini_files = []

    # If a single segment was requested, create a single file.
    # If multiple segments were requested, create an .ini file for each
    # segment.
    if int(options["pbs_numjobs"]) == 1:
        ini_content = template.render(options)
        ini_file = os.path.join(

            options["run_directory"], f"{options['runid']}.ini"
        )
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(ini_content)
        ini_files.append(ini_file)
    else:
        for job in range(int(options["pbs_numjobs"])):
            opt = options  # Need a copy of options
            if job > 0:
                opt["gamera_restart_doRes"] = "T"
            ini_content = template.render(options)
            ini_file = os.path.join(
                options["run_directory"], f"{options['runid']}-{job:02d}.ini"
            )
            ini_files.append(ini_file)
            with open(ini_file, "w", encoding="utf-8") as f:
                f.write(ini_content)

    # Return the paths to the .ini files.
    return ini_files


def convert_ini_to_xml(ini_files):
    """Convert the .ini files to XML.

    Convert the .ini files describing the run to XML files.

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

    Raises
    ------
    None
    """
    # Read and create the template.
    with open(PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)

    # Create a PBS script for each segment.
    pbs_scripts = []
    if int(options["pbs_numjobs"]) == 1:
        pbs_content = template.render(options)
        pbs_script = os.path.join(
            options["run_directory"], f"{options['runid']}.pbs"
        )
        with open(pbs_script, "w", encoding="utf-8") as f:
            f.write(pbs_content)
            pbs_scripts.append(pbs_script)
    else:
        for segment in range(int(options["pbs_numjobs"])):
            pbs_content = template.render(options)
            pbs_script = os.path.join(
                options["run_directory"],
                f"{options['runid']}-{segment:02d}.pbs"
            )
            pbs_scripts.append(pbs_script)
            with open(pbs_script, "w", encoding="utf-8") as f:
                f.write(pbs_content)

    # Create a single script which will submit all of the PBS jobs in order.
    path = "submit_pbs.sh"
    with open(path, "w", encoding="utf-8") as f:
        s = pbs_scripts[0]
        cmd = f"job_id=`qsub {s}`\n"
        f.write(cmd)
        cmd = f"echo $job_id\n"
        f.write(cmd)
        for s in pbs_scripts[1:]:
            cmd = "old_job_id=$job_id\n"
            f.write(cmd)
            cmd = f"job_id=`qsub -W depend=afterok:$old_job_id {s}`\n"
            f.write(cmd)
            cmd = f"echo $job_id\n"
            f.write(cmd)
    
    # Return the paths to the PBS scripts.
    return pbs_scripts


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
    options_path = args.options
    verbose = args.verbose
    if debug:
        print(f"args = {args}")

    # Fetch the run options.
    if options_path:
        # Read the run options from a JSON file.
        with open(options_path, "r", encoding="utf-8") as f:
            options = json.load(f)
    else:
        # Prompt the user for the run options.
        options = get_run_options()
    if debug:
        print(f"options = {options}")

    # Save the options dictionary as a JSON file.
    path = os.path.join(options["run_directory"], "options.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)

    # Save the current directory.
    original_directory = os.getcwd()

    # Move to the output directory.
    os.chdir(options["run_directory"])

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps(options)

    # Create the .ini file(s).
    if verbose:
        print("Creating .ini file(s) for run.")
    ini_files = create_ini_files(options)
    if debug:
        print(f"ini_files = {ini_files}")

    # Convert the .ini file(s) to .xml files(s).
    if verbose:
        print("Converting .ini file(s) to .xml file(s).")
    xml_files = convert_ini_to_xml(ini_files)
    if debug:
        print(f"xml_files = {xml_files}")

    # Create the PBS job script(s).
    if verbose:
        print("Creating PBS job script(s) for run.")
    pbs_scripts = create_pbs_scripts(options)
    if verbose:
        print(f"The PBS job scripts {pbs_scripts} are ready.")

    # Move back to the original directory.
    os.chdir(original_directory)


if __name__ == "__main__":
    """Begin main program."""
    main()
