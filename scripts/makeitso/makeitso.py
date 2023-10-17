#!/usr/bin/env python


"""makeitso for the MAGE magnetosphere software.

This master script is used to perform all of the steps needed to prepare
and run a kaiju job. This script is interactive - the user is prompted for
each decision that must be made to prepare for the run.
"""


# Import standard modules.
import argparse
import copy
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

# Path to current kaiju installation
KAIJUHOME = os.environ["KAIJUHOME"]

# Path to directory containing support files for makeitso.
SUPPORT_FILES_DIRECTORY = os.path.join(KAIJUHOME, "scripts", "makeitso")

# Path to option descriptions file.
OPTION_DESCRIPTIONS_FILE = os.path.join(
    SUPPORT_FILES_DIRECTORY, "option_descriptions.json"
)

# Path to template .ini file.
INI_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.ini")

# Path to template .pbs file.
PBS_TEMPLATE = os.path.join(SUPPORT_FILES_DIRECTORY, "template.pbs")

# Program defaults

# # Module sets by platform and run type.
# modules = {
#     "cheyenne": [
#         "cmake/3.22.0",
#         "git/2.33.1",
#         "ncarenv/1.3",
#         "intel/2022.1",
#         "geos/3.10.1",
#         "ncarcompilers/0.5.0",
#         "mpt/2.25",
#         "hdf5-mpi/1.12.2",
#     ],
#     "derecho": [
#         "ncarenv/23.06",
#         "cmake/3.26.3",
#         "craype/2.7.20",
#         "intel/2023.0.0",
#         "geos/3.9.1",
#         "ncarcompilers/1.0.0",
#         "cray-mpich/8.1.25",
#         "hdf5-mpi/1.12.2",
#     ],
#     "pleiades": [
#         "nas"
#         "pkgsrc/2022Q1-rome",
#         "comp-intel/2020.4.304",
#         "mpi-hpe/mpt.2.23",
#         "hdf5/1.8.18_mpt",
#     ],
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
        "--advanced", action="store_true",
        help="Configure advanced parameters (default: %(default)s)."
    )
    parser.add_argument(
        "--clobber", action="store_true",
        help="Overwrite existing options file (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true",
        help="Print debugging output (default: %(default)s)."
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


def get_run_option(name, description, advanced=False):
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
    advanced : bool
        True if parameter is from advanced configuration

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

    # If not in advanced mode, and this is an advanced option, take the
    # default.
    if not advanced and "advanced" in description:
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
    # Read the dictionary of option descriptions.
    with open(OPTION_DESCRIPTIONS_FILE, "r", encoding="utf-8") as f:
        option_descriptions = json.load(f)

    # Initialize the dictionary of program options.
    options = {}

    #-------------------------------------------------------------------------

    # General options for the simulation
    options["simulation"] = {}
    o = options["simulation"]
    od = option_descriptions["simulation"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    #-------------------------------------------------------------------------

    # PBS options
    options["pbs"] = {}
    o = options["pbs"]

    # Common (HPC platform-independent) options
    od = option_descriptions["pbs"]["_common"]
    od["account_name"]["default"] = os.getlogin()
    od["kaiju_install_directory"]["default"] = KAIJUHOME
    od["kaiju_build_directory"]["default"] = os.path.join(KAIJUHOME, "build_mpi")
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # Platform-specific options
    hpc_platform = o["hpc_system"]
    od = option_descriptions["pbs"][hpc_platform]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    #-------------------------------------------------------------------------

    # GAMERA options
    options["gamera"] = {}

    # Common options
    options["gamera"]["_common"] = {}
    o = options["gamera"]["_common"]
    od = option_descriptions["gamera"]["_common"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <sim> options
    options["gamera"]["sim"] = {}
    o = options["gamera"]["sim"]
    od = option_descriptions["gamera"]["sim"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <floors> options
    options["gamera"]["floors"] = {}
    o = options["gamera"]["floors"]
    od = option_descriptions["gamera"]["floors"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <timestep> options
    options["gamera"]["timestep"] = {}
    o = options["gamera"]["timestep"]
    od = option_descriptions["gamera"]["timestep"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <restart> options
    options["gamera"]["restart"] = {}
    o = options["gamera"]["restart"]
    od = option_descriptions["gamera"]["restart"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <physics> options
    options["gamera"]["physics"] = {}
    o = options["gamera"]["physics"]
    od = option_descriptions["gamera"]["physics"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <ring> options
    options["gamera"]["ring"] = {}
    o = options["gamera"]["ring"]
    od = option_descriptions["gamera"]["ring"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <ringknobs> options
    options["gamera"]["ringknobs"] = {}
    o = options["gamera"]["ringknobs"]
    od = option_descriptions["gamera"]["ringknobs"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <wind> options
    options["gamera"]["wind"] = {}
    o = options["gamera"]["wind"]
    od = option_descriptions["gamera"]["wind"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <source> options
    options["gamera"]["source"] = {}
    o = options["gamera"]["source"]
    od = option_descriptions["gamera"]["source"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <iPdir> options
    options["gamera"]["iPdir"] = {}
    o = options["gamera"]["iPdir"]
    od = option_descriptions["gamera"]["iPdir"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <jPdir> options
    options["gamera"]["jPdir"] = {}
    o = options["gamera"]["jPdir"]
    od = option_descriptions["gamera"]["jPdir"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <kPdir> options
    options["gamera"]["kPdir"] = {}
    o = options["gamera"]["kPdir"]
    od = option_descriptions["gamera"]["kPdir"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    #-------------------------------------------------------------------------

    # VOLTRON options
    options["voltron"] = {}

    # <time> options
    options["voltron"]["time"] = {}
    o = options["voltron"]["time"]
    od = option_descriptions["voltron"]["time"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <spinup> options
    options["voltron"]["spinup"] = {}
    o = options["voltron"]["spinup"]
    od = option_descriptions["voltron"]["spinup"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <output> options
    options["voltron"]["output"] = {}
    o = options["voltron"]["output"]
    od = option_descriptions["voltron"]["output"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <coupling> options
    options["voltron"]["coupling"] = {}
    o = options["voltron"]["coupling"]
    od = option_descriptions["voltron"]["coupling"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <restart> options
    options["voltron"]["restart"] = {}
    o = options["voltron"]["restart"]
    od = option_descriptions["voltron"]["restart"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <imag> options
    options["voltron"]["imag"] = {}
    o = options["voltron"]["imag"]
    od = option_descriptions["voltron"]["imag"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <ebsquish> options
    options["voltron"]["ebsquish"] = {}
    o = options["voltron"]["ebsquish"]
    od = option_descriptions["voltron"]["ebsquish"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    #-------------------------------------------------------------------------

    # CHIMP options
    options["chimp"] = {}

    # <units> options
    options["chimp"]["units"] = {}
    o = options["chimp"]["units"]
    od = option_descriptions["chimp"]["units"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <fields> options
    options["chimp"]["fields"] = {}
    o = options["chimp"]["fields"]
    od = option_descriptions["chimp"]["fields"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <domain> options
    options["chimp"]["domain"] = {}
    o = options["chimp"]["domain"]
    od = option_descriptions["chimp"]["domain"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <tracer> options
    options["chimp"]["tracer"] = {}
    o = options["chimp"]["tracer"]
    od = option_descriptions["chimp"]["tracer"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    #-------------------------------------------------------------------------

    # REMIX options
    options["remix"] = {}

    # <conductance> options
    options["remix"]["conductance"] = {}
    o = options["remix"]["conductance"]
    od = option_descriptions["remix"]["conductance"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <precipitation> options
    options["remix"]["precipitation"] = {}
    o = options["remix"]["precipitation"]
    od = option_descriptions["remix"]["precipitation"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    #-------------------------------------------------------------------------

    # RCM options
    options["rcm"] = {}

    # <conductance> options
    options["rcm"]["rcmdomain"] = {}
    o = options["rcm"]["rcmdomain"]
    od = option_descriptions["rcm"]["rcmdomain"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <ellipse> options
    options["rcm"]["ellipse"] = {}
    o = options["rcm"]["ellipse"]
    od = option_descriptions["rcm"]["ellipse"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <grid> options
    options["rcm"]["grid"] = {}
    o = options["rcm"]["grid"]
    od = option_descriptions["rcm"]["grid"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

    # <plasmasphere> options
    options["rcm"]["plasmasphere"] = {}
    o = options["rcm"]["plasmasphere"]
    od = option_descriptions["rcm"]["plasmasphere"]
    for on in od:
        o[on] = get_run_option(on, od[on], args.advanced)

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

    Raises
    ------
    None
    """
    # Create the LFM grid file.
    # NOTE: Assumes genLFM.py is in PATH.
    cmd = "genLFM.py"
    args = [cmd, "-gid", options["gamera"]["_common"]["LFM_grid_type"]]
    subprocess.run(args, check=True)

    # Create the solar wind file by fetching data from CDAWeb.
    # NOTE: Assumes cda2wind.py is in PATH.
    # cmd = "cda2wind.py"
    # args = [cmd, "-t0", start_date, "-t1", stop_date, "-interp"]
    # subprocess.run(args, check=True)

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
    if int(options["pbs"]["num_jobs"]) == 1:
        ini_content = template.render(options)
        ini_file = os.path.join(
            options["pbs"]["run_directory"],
            f"{options['simulation']['runid']}.ini"
        )
        with open(ini_file, "w", encoding="utf-8") as f:
            f.write(ini_content)
        ini_files.append(ini_file)
    else:
        for job in range(int(options["pbs"]["num_jobs"])):
            opt = copy.deepcopy(options)  # Need a copy of options
            if job > 0:
                opt["gamera_restart_doRes"] = "T"
            ini_content = template.render(options)
            ini_file = os.path.join(
                options["pbs"]["run_directory"], f"{options['simulation']['runid']}-{job:02d}.ini"
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
    # Read the template.
    with open(PBS_TEMPLATE, "r", encoding="utf-8") as f:
        template_content = f.read()
    template = Template(template_content)

    # Create a PBS script for each segment.
    # options_pbs = options["pbs"]
    # hpc_system = options_pbs["hpc_system"]
    pbs_scripts = []
    if int(options["pbs"]["num_jobs"]) == 1:
        pbs_content = template.render(options)
        pbs_script = os.path.join(
            options["pbs"]["run_directory"],
            f"{options['simulation']['runid']}.pbs"
        )
        with open(pbs_script, "w", encoding="utf-8") as f:
            f.write(pbs_content)
            pbs_scripts.append(pbs_script)
    else:
        for segment in range(int(options["pbs"]["num_jobs"])):
            pbs_content = template.render(options)
            pbs_script = os.path.join(
                options["pbs"]["run_directory"],
                f"{options['simulation']['runid']}-{segment:02d}.pbs"
            )
            pbs_scripts.append(pbs_script)
            with open(pbs_script, "w", encoding="utf-8") as f:
                f.write(pbs_content)

    # Create a single script which will submit all of the PBS jobs in order.
    path = f"{options['simulation']['runid']}_pbs.sh"
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
    path = f"{options['simulation']['runid']}.json"
    if os.path.exists(path):
        if not clobber:
            raise FileExistsError(f"Options file {path} exists!")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(options, f, indent=JSON_INDENT)

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


if __name__ == "__main__":
    """Begin main program."""
    main()
