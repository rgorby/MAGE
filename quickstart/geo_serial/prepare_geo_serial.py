#!/usr/bin/env python


"""Prepare for running serial kaiju on the geo_serial example.

Perform the preprocessing required to run the serial kaiju code on the
geo_serial example. Create any required data files, and create the PBS
script to run the code.
"""


# Import standard modules.
import argparse
import os
import shutil
import subprocess

# Import 3rd-party modules.

# Import project-specific modules.


# Program constants and defaults

# Default identifier for model to run.
default_runid = "geo_serial"

# Program description.
description = "Prepare to run serial kaiju on the %s quickstart case." % default_runid

# Location of template .ini file for run.
ini_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", default_runid, "%s_template.ini"
    % default_runid
)

# Location of template PBS script for run.
pbs_template = os.path.join(
    os.environ["KAIJUHOME"], "quickstart", default_runid, "%s_template.pbs"
    % default_runid
)

# Name of HDF5 file containing solar wind data for initial conditions.
sw_file_name = "bcwind.h5"


def create_command_line_parser():
    """Create the command-line argument parser.
    
    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
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
        "--startdate", type=str, default="2016-08-09T09:00:00",
        help="Specify the start date in ISO 8601 format (default: %(default)s)."
    )
    parser.add_argument(
        "--stopdate", type=str, default="2016-08-09T10:00:00",
        help="Specify the stop date in ISO 8601 format (default: %(default)s)."
    )
    parser.add_argument(
        "--swfile", type=str,
        help="Specify an existing file of solar wind data for initialization. "
            "If used, --startdate and --stopdate must be set to the time "
            "limits of this file."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def run_preprocessing_steps(directory, runid, startdate, stopdate, swfile=None):
    """Run any preprocessing steps needed for this run.

    Perform required preprocessing steps.

    Parameters
    ----------
    directory : str
        Path to directory to receive preprocessing results.
    runid : str
        ID string for the model to run.
    startdate, stopdate : str
        Start and stop date & time for run in ISO 8601 format.
    swfile : str, default None
        Path to optional file of solar wind data for initial conditions. If
        provided, data will not be fetched from CDAWeb for the specified
        time range.

    Returns
    -------
    None
    """
    # Save the current directory.
    original_directory = os.getcwd()

    # Move to the output directory.
    os.chdir(directory)

    # Create the grid file using "double" resolution.
    cmd = "genLFM.py"
    args = ["-gid", "D"]
    subprocess.run([cmd] + args)

    # Create or copy the solar wind file.
    if swfile is not None:
        # Use an existing solar wind data file.
        shutil.copyfile(swfile, sw_file_name)
    else:
        # Fetch data from CDAWeb.
        cmd = "cda2wind.py"
        args = ["-t0", startdate, "-t1", stopdate, "-interp"]
        subprocess.run([cmd] + args)

    # Create the RCM configuration file.
    cmd = "genRCM.py"
    args = []
    subprocess.run([cmd] + args)

    # Move back to the originaldirectory.
    os.chdir(original_directory)


def create_ini_file(directory, runid):
    """Create the .ini file from a template.

    Create the .ini file from a template.

    Parameters
    ----------
    directory : str
        Path to directory to receive .ini file.
    runid : str
        ID string for the model to run.

    Returns
    -------
    ini_file : str
        Path to .ini file.
    """
    # Read the file template.
    with open(ini_template) as t:
        lines = t.readlines()

    # Process the template here.

    # Write the processed .ini file to the run directory.
    ini_file = os.path.join(directory, "%s.ini" % runid)
    with open(ini_file, "w") as f:
        f.writelines(lines)
    return ini_file


def convert_ini_to_xml(ini_file, xml_file):
    """Convert the .ini file to XML.
    
    Convert the .ini file to a .xml file.

    Parameters
    ----------
    ini_file : str
        Path to the .ini file to convert.
    xml_file : str
        Path to the resulting XML file.
    
    Returns
    -------
    None
    """
    cmd = "XMLGenerator.py"
    args = [ini_file, xml_file]
    subprocess.run([cmd] + args)


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
    # Read the template.
    with open(pbs_template) as t:
        lines = t.readlines()

    # Process the template here.

    # Write out the processed file.
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
    startdate = args.startdate
    stopdate = args.stopdate
    swfile = args.swfile
    verbose = args.verbose

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps(directory, runid, startdate, stopdate, swfile)

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
        print("Edit this file as needed for your system (see comments in %s"
              " for more information)." % pbs_file)
        print("Submit the job to PBS with the command:")
        print("    qsub %s" % pbs_file)
