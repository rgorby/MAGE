#!/usr/bin/env python

"""Compare gamhelio results with spacecraft data.

Compare heliospheric model results from gamhelio with data measured by
spacecraft.

Authors
-------
Eric Winter (eric.winter@jhuapl.edu)
Mike Wiltberger
"""


# Include standard modules.
import argparse
from argparse import RawTextHelpFormatter
import os

# Include 3rd-party modules.
import numpy as np
import spacepy.datamodel as dm

# Include project modules.
import kaipy.kaiH5 as kaiH5
import kaipy.kaiViz as kv
import kaipy.kaiTools as kaiTools
import kaipy.satcomp.scutils as scutils


# Program constants.

# Program description string.
description = """Extract satellite trajectory and observations for various
heliospheric spacecraft from CDAWeb. Produce comparisons between the
observations and corresponding gamhelio model results."""

# Default path to sctrack.x
default_cmd = os.path.join(
    os.environ["KAIJUHOME"], "build", "bin", "sctrack.x"
)

# Default time interval for ephemeris data returned from CDAWeb (seconds).s
default_deltaT = 3600.00  # 1 hour

# Default run ID string.
default_runid = "hsphere"

# Defaut number of segments to process.
default_numSeg = 1

# Default path to model results directory.
default_path=os.getcwd()

# Path to file of heliospheric spacecraft data.
spacecraft_data_file = os.path.join(
    os.environ["KAIJUHOME"], "kaipy", "satcomp", "sc_helio.json"
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
        "-c", "--cmd", type=str, metavar="command", default=default_cmd,
        help="Full path to sctrack.x command (default: %(default)s)."
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--deltaT", type=float, metavar="deltaT", default=default_deltaT,
        help="Time interval (seconds) for ephemeris points from CDAWeb (default: %(default)s)."
    )
    parser.add_argument(
        "-i", "--id", type=str, metavar="runid", default=default_runid,
        help="ID string of the run (default: %(default)s)"
    )
    parser.add_argument(
        "-k", "--keep", action="store_true", default=False,
        help="Keep intermediate files (default: %(default)s).")
    parser.add_argument(
        "-n", "--numSeg", type=int, metavar="number_segments",default=default_numSeg,
        help="Number of segments to simultaneously process (default: %(default)s).")
    parser.add_argument(
        "-p", "--path", type=str, metavar="path", default=default_path,
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
    deltaT = args.deltaT
    ftag = args.id
    keep = args.keep
    numSegments = args.numSeg
    fdir = args.path
    scRequested = args.satId
    verbose = args.verbose
    if debug:
        print("args = %s" % args)

    # Read the list of available spacecraft from the YAML configuration file.
    scIds = scutils.getScIds(spacecraft_data_file, doPrint=verbose)
    if debug:
        print("scIds = %s" % scIds)

    # Compute the path to the gamhelio output file to examine.
    (fname, isMPI, Ri, Rj, Rk) = kaiTools.getRunInfo(fdir, ftag)
    if debug:
        print("fname = %s" % fname)
        print("isMPI = %s" % isMPI)
        print("(Ri, Rj, Rk) = (%s, %s, %s)" % (Ri, Rj, Rk))

    # Determine the number of steps in the gamhelio output file, and get a
    # list of the step indices.
    (nsteps, sIds) = kaiH5.cntSteps(fname)
    if debug:
        print("nsteps = %s" % nsteps)
        print("sIds = %s" % sIds)

	# Pull the timestep information (the "MJD" attribute) from the gamhelio
    # output files. This fetches the "MJD" attribute from each of the top-
    # level groups called "Step#[\d]+". These MJD values are floats.
    gamMJD = kaiH5.getTs(fname, sIds, aID="MJD")
    if debug:
        print("gamMJD = %s" % gamMJD)

    # Now get the "time" attribute from each step. For gamhelio, these elapsed
    # times are in seconds since the start of the simulation.
    gamT = kaiH5.getTs(fname, sIds, aID="time")
    if debug:
        print("gamT = %s" % gamT)

    # Convert the MJD values to Universal Time datetime objects.
    gamUT = kaiTools.MJD2UT(gamMJD)
    if debug:
        print("gamUT = %s" % gamUT)

    # Use the first (non zero?) time as the inital MJD.
    # N.B. THIS SKIPS THE FIRST SIMULATION STEP SINCE IT TYPICALLY HAS
    # gamT[0] = 0.
    # Use the last time as the last MJD.
    loc = np.argwhere(gamT > 0.0)[0][0]
    t0 = gamUT[loc]  # First non-0 time
    t1 = gamUT[-1]
    if debug:
        print("t0 = %s" % t0)
        print("t1 = %s" % t1)

    # Save the (float) MJD of the first used step.
    mjdFileStart = gamMJD[loc]
    if debug:
        print("mjdFileStart = %s" % mjdFileStart)

    # Save the elapsed simulation time (seconds) of the first step.
    secFileStart = gamT[loc]
    if debug:
        print("secFileStart = %s" % secFileStart)

    # Determine the list of IDs of spacecraft to fetch data from.
    if scRequested:
        scToDo = [scRequested]
    else:
        scToDo = list(scIds.keys())
    if debug:
        print("scToDo = %s" % scToDo)

    # Fetch the ephemeris and observed data for each spacecraft in the list.
    for scId in scToDo:

        # Fetch the ephemeris and observed data for the current spacecraft.
        if verbose:
            print("Getting ephemeris and instrument data from CDAWeb for %s." % scId)
        status, data = scutils.getSatData(
            scIds[scId],  # ID string of spacecraft
            t0.strftime("%Y-%m-%dT%H:%M:%SZ"),  # Start time for gamhelio results
            t1.strftime("%Y-%m-%dT%H:%M:%SZ"),  # Stop time for gamhelio results
            deltaT  # Time interval (seconds) for data returned from CDAWeb
        )
        if debug:
            print("status = %s" % status)
            print("data = %s" % data)

        # If no data was found for the spacecraft, go to the next.
        if status["http"]["status_code"] != 200 or data is None:
            print("No data available for %s." % scId)
        else:
            # Use the spacecraft trajectory to interpolate simulated
            # observations from the gamhelio output.
            if verbose:
                print("Interpolating simulated observations along trajectory.")
            scutils.extractGAMHELIO(
                data, scIds[scId], scId, mjdFileStart, secFileStart, fdir, ftag,
                cmd, numSegments, keep
            )
            cdfname = os.path.join(fdir, scId + ".comp.cdf")
            if os.path.exists(cdfname):
                if verbose:
                    print("Deleting existing CDF comparison file %s" % cdfname)
                os.system("rm %s" % cdfname)
            if verbose:
                print("Creating CDF file %s with %s and GAMERA data" % (cdfname, scId))
            dm.toCDF(cdfname, data)
            plotname = os.path.join(fdir, scId + ".png")
            if verbose:
                print("Plotting results to %s." % plotname)
            kv.compPlot(plotname, scId, data)
            if verbose:
                print("Computing errors.")
            errname = os.path.join(fdir, scId + "-error.txt")
            if verbose:
                print("Writing errors to %s." % errname)
            scutils.errorReport(errname, scId, data)
            plotname = os.path.join(fdir, scId + "-traj.png")
            if verbose:
                print("Plotting trajectory to %s." % plotname)
            kv.trajPlot(plotname, scId, data)
