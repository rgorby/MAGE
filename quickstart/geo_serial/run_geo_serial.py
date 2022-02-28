#!/usr/bin/env python


"""Prepare the PBS script for a geo_serial job.

Create the PBS script required to submit a gamera geo_serial model to PBS.
"""


# Import standard modules.
import argparse
import subprocess

# Import 3rd-party modules.

# Import project-specific modules.


# Program constants and defaults

# Identifier for run
runid = "geo_serial"

# Program description.
description = "Prepare and run gamera on the %s test case." % runid

# # Location of template .ini file and .ini file for run.
# ini_template = os.path.join(
#     os.environ["KAIJUHOME"], "quickstart", runid, "%s.ini.template" % runid
# )
# ini_file = "%s.ini" % runid

# # Location of template XML file and XML file for run.
# xml_template = os.path.join(
#     os.environ["KAIJUHOME"], "quickstart", runid, "%s.xml.template" % runid
# )
# xml_file = "%s.xml" % runid

# # Location of template PBS script and PBS script for run.
# pbs_template = os.path.join(
#     os.environ["KAIJUHOME"], "quickstart", runid, "%s.pbs.template" % runid
# )
# pbs_file = "%s.pbs" % runid


def create_command_line_parser():
    """Create the command-line argument parser.
    
    Ceate the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.
    """
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


def run_preprocessing_steps(startdate, stopdate):
    """Run any preprocessing steps needed for the geo_serial run.

    Run any required preprocessing steps to prepare for the geo_serial run.

    Parameters
    ----------
    startdate, stopdate : str
        Start and stop date & time for run in ISO 8601 format.

    Returns
    -------
    None
    """
    # Create the grid file.
    cmd = "genLFM.py"
    subprocess.run([cmd])

    # Create the solar wind file.
    cmd = "omni2wind.py"
    args = ["-t0", startdate, "-t1", "stopdate", "-interp"]
    subprocess.run([cmd] + args)

    # Create the RCM configuration file.
    cmd = "genRCM.py"
    subprocess.run([cmd])


# def create_ini_file():
#     """Create the .ini file from a template.
    
#     Create the .ini file describing the geo_serial model run.

#     For now, we simply make a copy of the .ini template.

#     Parameters
#     ----------
#     None

#     Returns
#     -------
#     ini_file : str
#         Path to the .ini file for the geo_serial model run.
#     """
#     # Just use the template for now.
#     with open(ini_template) as t:
#         lines = t.readlines()
#     with open(ini_file, "w") as f:
#         f.writelines(lines)
#     return ini_file


# def convert_ini_to_xml(ini_file):
#     """Convert the .ini file to XML.
    
#     Convert the .ini file describing the geo_serial run to the corresponding
#     XML file.

#     Parameters
#     ----------
#     ini_file : str
#         Path to the .ini file to convert.
    
#     Returns
#     -------
#     xml_file : str
#         Path to the resulting XML file.
#     """
#     # No conversion is performed yet. Just print the template.
#     with open(xml_template) as t:
#         lines = t.readlines()
#     with open(xml_file, "w") as f:
#         f.writelines(lines)
#     return xml_file


# def create_pbs_job_script():
#     """Create the PBS job script for the run.
    
#     Create the PBS job script which can be submitted to PBS to perform
#     the geo_serial model run.

#     Parameters
#     ----------
#     None

#     Returns
#     -------
#     pbs_file : str
#         Path to PBS job description file.
#     """
#     # No changes yet. Just print the template.
#     with open(pbs_template) as t:
#         lines = t.readlines()
#     with open(pbs_file, "w") as f:
#         f.writelines(lines)
#     return pbs_file


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    startdate = args.startdate
    stopdate = args.stopdate
    verbose = args.verbose

    # Run the preprocessing steps.
    if verbose:
        print("Running preprocessing steps.")
    run_preprocessing_steps(startdate, stopdate)

#     # Create the .ini file.
#     if verbose:
#         print("Creating .ini file for run.")
#     ini_file = create_ini_file()

#     # Convert the .ini file to a .xml file.
#     if verbose:
#         print("Converting .ini file to .xml file for run.")
#     xml_file = convert_ini_to_xml(ini_file)

#     # Create the PBS job script.
#     if verbose:
#         print("Creating PBS job script for run.")
#     pbs_file = create_pbs_job_script()

#     print("""
# The PBS job script is ready for submission.
# You may submit it with the following command:

#     qsub %s""" % pbs_file)