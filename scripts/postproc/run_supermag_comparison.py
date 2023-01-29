#!/usr/bin/env python

"""Run a SuperMag comparison for a MAGE magnetosphere run.

Perform a comparison of ground magnetic field perturbations computed for a
MAGE magnetosphere simulation with measured data from SuperMag.

This code caches SuperMAG data in your home directory, in a subdirectory
called "supermag". Make sure this directory exists (even if it is empty)
before running this script. Note that this caching can be very time-consuming
the first time it is run.

Author
------
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import math
import os
import subprocess

# Import 3rd-party modules.
import matplotlib as mpl
import matplotlib.pyplot as plt

# Import project-specific modules.
import kaipy.supermage as sm


# Program constants and defaults

# Program description.
description = "Compare MAGE ground delta-B to SuperMag measurements."

# Default identifier for results to read.
default_runid = "msphere"

# Location of template XML file.
xml_template_path = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "postproc", "calcdb.xml.template"
)

# Name of XML file read by calcdb.x.
xml_filename_template = "calcdb_RUNID.xml"

# Command to run for calcdb.x.
calcdb_cmd = "calcdb.x"

# calcdb.x file to capture stdout, stderr.
calcdb_output_file = "calcdb.out"

# Tempalte for filename for delta-B values.
deltab_filename_template = "%s.deltab.h5"

# User ID to use when fetching SupreMAG data.
supermag_user_id = "ewinter"

# Number of seconds in a day.
SECONDS_PER_DAY = 86400

# Local cache folder for SuperMAG data.
SUPERMAG_CACHE_FOLDER = os.path.join(os.environ["HOME"], "supermag")

# Limit of fraction of bad data to keep SuperMAG results.
SUPERMAG_BAD_FRAC = 0.1

# Template for file containing index plot.
INDEX_PLOT_FILENAME_TEMPLATE = "%s_indices.png"

# Template for file containing contour plots.
CONTOUR_PLOT_FILENAME_TEMPLATE = "%s_contours.png"


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
        "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-d", type=str, metavar="directory", default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "-id", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def create_calcdb_xml_file(runid):
    """Create the XML input file for calcdb.x from a template.

    Create the XML input file for calcdb.x from a template. The file is
    created in the current directory.

    Parameters
    ----------
    runid : str
        runid for MAGE results file.
    
    Returns
    -------
    xml_file : str
        Name of XML file.
    """
    # Read the template file.
    with open(xml_template_path) as t:
        lines = t.readlines()

    # Process the template here.
    # <HACK>
    # This should be done with a proper templating package.
    lines[3] = lines[3].replace('RUNID', runid)
    lines[5] = lines[5].replace('EBFILE', runid)
    # </HACK>

    # Write out the processed XML.
    xml_file = xml_filename_template.replace('RUNID', runid)
    with open(xml_file, "w") as f:
        f.writelines(lines)

    # Return the name of the XML file.
    return xml_file


def compute_ground_delta_B(runid):
    """Compute ground delta B values for a MAGE run.

    Compute ground delta B values for a MAGE run. The computation is done with
    the program calcdb.x.

    Parameters
    ----------
    runid : str
        runid for MAGE results file.

    Returns
    -------
    delta_B_file : str
        Name of file in cirrent directory containing calcdb.x results.
    """
    # Create the XML file for calcdb.x from the template.
    xml_file = create_calcdb_xml_file(runid)

    # Run the command to compute ground delta B values.
    # Capture stdout and stderr to a file, raise exception on non-zero exit.
    # result = subprocess.run([calcdb_cmd, xml_file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # with open(calcdb_output_file, "w") as f:
    #     f.write(result.stdout.decode('utf-8'))

    # Compute the name of the file containing the delta B values.
    delta_B_file = deltab_filename_template % runid

    # Verify that the file was created.
    assert os.path.isfile(delta_B_file)

    # Return the name of the file containing the delta B results.
    return delta_B_file


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    results_directory = args.d
    runid = args.id
    verbose = args.verbose
    if debug:
        print("args = %s" % args)

    # Move to the results directory.
    os.chdir(results_directory)

    # Compute the ground delta B values for this run.
    delta_B_file = compute_ground_delta_B(runid)
    if debug:
        print("delt_B_file = %s" % delta_B_file)

    # Read the delta B values.
    SIM = sm.ReadSimData(delta_B_file)

    # Fetch the SuperMag indices for the desired time range.
    supermag_user = supermag_user_id    # username used with SuperMag
    datetime_start = SIM['td'][0]       # start time of simulation data
    datetime_end = SIM['td'][-1]        # end time of simulation data
    # Fetch data in 1-day chunks to ensure proper caching?
    num_days = int(math.ceil((datetime_end - datetime_start).seconds/SECONDS_PER_DAY))
    if debug:
        print("supermag_user = %s" % supermag_user)
        print("datetime_start = %s" % datetime_start)
        print("datetime_end = %s" % datetime_end)
        print("num_days = %s" % num_days)
    SMI = sm.FetchSMIndices(
        user=supermag_user,
        start=datetime_start,
        numofdays=num_days
    )
    if debug:
        print("SMI = %s" % SMI)

    # Fetch the SuperMag station data for the desired time range.
    local_supermag_cache_path = SUPERMAG_CACHE_FOLDER
    if debug:
        print("local_supermag_cache_path = %s" % local_supermag_cache_path)
    SM = sm.FetchSMData(
        user=supermag_user,
        start=datetime_start,
        numofdays=num_days,
        savefolder=local_supermag_cache_path,
        badfrac=SUPERMAG_BAD_FRAC
    )
    if debug:
        print("SM = %s" % SM)

    # Interpolate the simulated delta B to the measurement times from SuperMag.
    SMinterp = sm.InterpolateSimData(SIM=SIM, SM=SM)
    if debug:
        print("SMinterp = %s" % SMinterp)

    # Create the plots in memory.
    mpl.use('Agg')

    # Make the indices plot.
    sm.MakeIndicesPlot(SMI=SMI, SMinterp=SMinterp, fignumber=1)
    plt.savefig(INDEX_PLOT_FILENAME_TEMPLATE % runid)

    # Make the contour plots.
    sm.MakeContourPlots(SM=SM, SMinterp=SMinterp, fignumber=2)
    plt.savefig(CONTOUR_PLOT_FILENAME_TEMPLATE % runid)
