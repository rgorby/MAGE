#!/usr/bin/env python

"""Compare gamhelio results with spacecraft data.

Compare heliospher model results from gamhelio with data measured by spacecraft.

Authors
-------
M. Wiltberger
E. Winter (eric.winter@jhuapl.edu)
"""


# Include standard modules.
import argparse
from argparse import RawTextHelpFormatter
import os
import sys

# Include 3rd-party modules.
import numpy as np
# from   astropy.time import Time
# import h5py
# import spacepy.datamodel as dm

# Include project modules.
import kaipy.kaiH5 as kaiH5
# import kaipy.kaiViz as kv
import kaipy.kaiTools as kaiTools
# import kaipy.kaijson as kj
import kaipy.satcomp.scutils as scutils


# Program constants.

# Program description string.
description = """Extracts information from satellite trajectory for various
spacecraft. Spacecraft data is pulled from CDAWeb. Output CDF files
contain data pulled from CDAWeb along with data extracted from gamhelio.
Image files of satellite comparisons are also produced."""

# Defaut number of segments to process.
default_numseg = 1

# Default run ID string.
default_runid = "hsphere"

# Path to file of spacecraft data.
package_directory = os.path.dirname(os.path.abspath(__file__))
spacecraft_data_file = os.path.join(
    package_directory, "..", "..", "kaipy", "satcomp", "sc_helio.json"
)


def create_command_line_parser():
    """Create the command-line argument parser.

    Prepare the command-line parser.

    Parameters
    ----------
    None

    Returns
    -------
    parse : argparse.ArgumentParser
        Parser for command-line arguments.
    """
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-c", "--cmd", type=str, metavar="command", default=None,
        help="Full path to sctrack.x command (default: %(default)s)."
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-i", "--id", type=str, metavar="runid", default=default_runid,
        help="ID string of the run (default: %(default)s)"
    )
    parser.add_argument(
        "-k", "--keep", action="store_true", default=False,
        help="Keep intermediate files (default: %(default)s).")
    parser.add_argument(
        "-n", "--numSeg", type=int, metavar="number_segments", default=default_numseg,
        help="Number of segments to simultaneously process (default: %(default)s).")
    parser.add_argument(
        "-p", "--path", type=str, metavar="path", default=os.getcwd(),
        help="Path to directory containing gamhelio results (default: %(default)s)"
    )
    parser.add_argument(
        "-s", "--satId", type=str, metavar="satellite_id", default=None,
        help="Name of Satellite to compare"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


if __name__ == "__main__":
    """Begin main program."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    cmd = args.cmd
    debug = args.debug
    ftag = args.id
    keep = args.keep
    numSegments = args.numSeg
    fdir = args.path
    scRequested = args.satId
    verbose = args.verbose

    # If the current directory is used, get the complete path.
    if fdir == ".":
        fdir = os.getcwd()

    # If the path to sctrack.x was not provided, construct it based on the
    # installation directory.
    if cmd == None:
        cmd = os.path.join(os.getenv("KAIJUHOME"), "build", "bin", "sctrack.x")
    if not (os.path.isfile(cmd) and os.access(cmd, os.X_OK)):
        print(cmd, "either not found or not executable")
        sys.exit()

    # Read the list of available spacecraft from the YAML configuration file.
    scIds = scutils.getScIds(spacecraft_data_file, doPrint=verbose)

    # Compute the path to the gamhelio output file to examine.
    (fname, isMPI, Ri, Rj, Rk) = kaiTools.getRunInfo(fdir, ftag)

    # Determine the number of steps in the gamhelio output file, and get a
    # list of the step indices.
    nsteps, sIds = kaiH5.cntSteps(fname)

	# Pull the timestep information (the "MJD" attribute) from the gamhelio
    # output files. This fetches the "MJD" attribute from each of the top-
    # level groups called "Step#[\d]+". These MJD values are floats.
    gamMJD = kaiH5.getTs(fname, sIds, aID="MJD")

    # Now get the "time" attribute from each step. This represents the time
    # in seconds since the simulation start time.
    gamT = kaiH5.getTs(fname, sIds, aID="time")

    # Convert the MJD values to Universal Time datetime objects.
    gamUT = kaiTools.MJD2UT(gamMJD)

    # Use the first non zero time as the inital MJD.
    # N.B. THIS SKIPS THE FIRST SIMULATION STEP SINCE IT HAS gamT[0] = 0.
    loc = np.argwhere(gamT > 0.0)[0][0]
    t0 = gamUT[loc]  # First non-0 time
    t1 = gamUT[-1]  # Last time

    # Compute the time interval (seconds) between step 1 and step 2.
    deltaT = np.round(gamT[loc + 1] - gamT[loc])

    # Save the (float) MJD of the first step.
    mjdFileStart = gamMJD[loc]

    # Save the elapsed simulation time in seconds of the first step.
    secFileStart = gamT[loc]

    # Determine the list of spacecraft to fetch data from.
    if scRequested is None:
        # This is a list of dictionaries read from the JSON spacecraft file.
        # The keys are the spacecraft ID strings.
        scToDo = scIds
    else:
        # This is a list of spacecraft ID strings from the command line. 
        scToDo = [scRequested]

    # Fetch the data for each spacecraft in the list.
    for scId in scToDo:
        print("Getting spacecraft data for %s." % scId)

        # Try to fetch the data for the current spacecraft. The data will be
        # between the first and last simulation times, and are assumed to have
        # a time interval equal to the simulation interval deltaT (?).
        status, data = scutils.getSatData(
            scIds[scId],
            t0.strftime("%Y-%m-%dT%H:%M:%SZ"),
            t1.strftime("%Y-%m-%dT%H:%M:%SZ"),
            deltaT
        )

        # If no data was found for the spacecraft, go to the next.
        if status["http"]["status_code"] != 200 or data is None:
            print("No data available for %s." % scId)
        else:
            # Extract the spacecraft data coprresponding to the gamera simulation times.
            print("Extracting GAMERA data.")
            scutils.extractGAMERA(
                data, scIds[scId], scId, mjdFileStart, secFileStart, fdir, ftag,
                cmd, numSegments, keep
            )
            # scutils.matchUnits(data)
        # 		cdfname = os.path.join(fdir, scId + '.comp.cdf')
        # 		if os.path.exists(cdfname):
        # 			print('Deleting %s' % cdfname)
        # 			os.system('rm %s' % cdfname)
        # 		print('Creating CDF file',cdfname,'with',scId,'and GAMERA data')
        # 		dm.toCDF(cdfname,data)
        # 		plotname = os.path.join(fdir,scId+'.png')
        # 		print('Plotting results to',plotname)
        # 		kv.compPlot(plotname,scId,data)
        # 		print('Computing Errors')
        # 		errname = os.path.join(fdir,scId+'-error.txt')
        # 		scutils.errorReport(errname,scId,data)
        # 		plotname = os.path.join(fdir,scId+'-traj.png')
        # 		print('Plotting trajectory to',plotname)
        # 		kv.trajPlot(plotname,scId,data)
    pass
