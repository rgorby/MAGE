#!/usr/bin/env python

"""Run a SuperMag comparison for a MAGE magnetosphere run.

Perform a comparison of ground magnetic field perturbations computed for a
MAGE magnetosphere simulation with measured data from SuperMag.

Author
------
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import os
# import shutil
import subprocess

# Import 3rd-party modules.
import matplotlib as mpl
import matplotlib.pyplot as plt

# Import project-specific modules.
import kaipy.supermage as sm


# Program constants and defaults

# Program description.
description = "Compare MAGE ground delta-B to SuperMag measurements."

# Location of template XML file.
xml_template = os.path.join(
    os.environ["KAIJUHOME"], "scripts", "postproc", "calcdb.xml.template"
)

# Name of XML file read by calcdb.x.
xml_filename_template = "calcdb_RUNID.xml"


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
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        "mage_results_path",
        help='Path to a result file for a MAGE magnetosphere run.'
    )
    return parser


def filename_to_runid(filename):
    """Parse the runid from a MAGE results file name.

    Parse the runid from a MAGE results file name.

    The runid is all text before the first period in the name.

    Parameters
    ----------
    filename : str
        Name of MAGE results file.

    Returns
    -------
    runid : str
        The MAGE runid for the file.
    """
    parts = filename.split('.')
    runid = parts[0]
    return runid


def create_xml_file(runid):
    """Create the XML input file for calcdb.x from a template.

    Create the XML input file for calcdb.x from a template.

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
    with open(xml_template) as t:
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
        Name of file containing calcdb.x results.
    """
    # Create the XML file for calcdb.x from the template.
    xml_file = create_xml_file(runid)

    # Run the command to compute ground delta B values.
    cmd = 'calcdb.x'
    args = [xml_file]
    subprocess.run([cmd] + args)

    # Compute the name of the file containing the delta B values.
    delta_B_file = runid + 'deltab.h5'
    return delta_B_file


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    verbose = args.verbose
    mage_results_path = args.mage_results_path

    # Split the MAGE results path into a directory and a file.
    (mage_results_dir, mage_results_file) = os.path.split(mage_results_path)

    # Compute the runid from the file name.
    runid = filename_to_runid(mage_results_file)

    # Move to the results directory.
    os.chdir(mage_results_dir)

    # Compute the ground delta B values for this run.
    # delta_B_file = compute_ground_delta_B(runid)
    delta_B_file = 'msphere.deltab.h5'

    # Read the delta B values.
    SIM = sm.ReadSimData(delta_B_file)

    # Fetch the SuperMag indices for the desired time range.
    user = 'ewinter'    # username used with SuperMag
    start = SIM['td'][0] # start time of simulation data
    numofdays = 3        # retrieve 4 days of data 
    SMI  = sm.FetchSMIndices(user, start, numofdays)

    # Fetch the SuperMag data for the desired time range.
    savefolder = '/Users/winteel1/supermag'
    SM = sm.FetchSMData(user, start, numofdays, savefolder, badfrac = 0.1)

    # Interpolate the simulated delta B to the measurement times from SuperMag.
    SMinterp = sm.InterpolateSimData(SIM, SM)

    # Create the plots in memory.
    mpl.use('Agg')

    # Make the indices plot.
    sm.MakeIndicesPlot(SMI, SMinterp, fignumber = 1)
    plt.savefig(runid + '_indices.png')

    # Make the contour plots.
    sm.MakeContourPlots(SM, SMinterp, maxx = 1000, fignumber = 2)
    plt.savefig(runid + '_contours.png')
