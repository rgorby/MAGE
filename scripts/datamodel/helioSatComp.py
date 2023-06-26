#!/usr/bin/env python

"""Compare gamhelio results with spacecraft data.

Compare heliospheric model results from gamhelio with data measured by
spacecraft.

Note that the terms "ephemeris" and "trajectory" are used interchangeably.

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

# Default time interval for ephemeris data returned from CDAWeb (seconds).
default_deltaT = 3600.00  # 1 hour

# Default run ID string.
default_runid = "wsa"

# Defaut number of segments to process.
default_numSeg = 1

# Default path to model results directory.
default_path = os.getcwd()

# Path to file of heliospheric spacecraft metadata.
helio_sc_metadata_path = os.path.join(
    os.environ["KAIJUHOME"], "kaipy", "satcomp", "sc_helio.json"
)


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the command-line parser.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Parser for command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=RawTextHelpFormatter
    )
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
        help="Time interval (seconds) for ephemeris points from CDAWeb " +
             "(default: %(default)s)."
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
    cmd_sctrack = args.cmd
    debug = args.debug
    cdaweb_data_interval = args.deltaT
    gh_run_id = args.id
    keep = args.keep
    num_segments = args.numSeg
    gh_result_directory = args.path
    sc_to_compare = args.satId
    verbose = args.verbose
    if debug:
        print("args = %s" % args)

    # Read the list of available spacecraft from the YAML configuration file.
    if verbose:
        print("Reading heliosphere spacecraft metadata from %s." %
              helio_sc_metadata_path)
    sc_metadata = scutils.getScIds(helio_sc_metadata_path, doPrint=verbose)
    if debug:
        print("sc_metadata = %s" % sc_metadata)

    # Compute the path to the gamhelio output files to examine.
    if verbose:
        print("Looking for gamhelio results for run %s in %s." %
              (gh_run_id, gh_result_directory))
    (gh_result_path, is_MPI, Ri, Rj, Rk) = kaiTools.getRunInfo(
        gh_result_directory, gh_run_id
    )
    if debug:
        print("gh_result_path = %s" % gh_result_path)
        print("is_MPI = %s" % is_MPI)
        print("(Ri, Rj, Rk) = (%s, %s, %s)" % (Ri, Rj, Rk))

    # Determine the number of timesteps in the gamhelio output files, and get
    # a list of the timestep indices.
    if verbose:
        print("Counting timesteps in %s." % gh_result_path)
    (gh_n_timesteps, gh_timestep_indices) = kaiH5.cntSteps(gh_result_path)
    if verbose:
        print("  Found %s timesteps." % gh_n_timesteps)
    if debug:
        print("gh_n_timesteps = %s" % gh_n_timesteps)
        print("gh_timestep_indices = %s" % gh_timestep_indices)

	# Pull the timestep information from the gamhelio output files. This
    # fetches the "MJD" attribute from each of the top-level groups called
    # "Step#[\d]+". These MJD values are floats, and are assumed to be in
    # increasing order.
    if verbose:
        print("Reading timestep times from %s." % gh_result_path)
    gh_timestep_MJDs = kaiH5.getTs(
        gh_result_path, gh_timestep_indices, aID="MJD"
    )
    if debug:
        print("gh_timestep_MJDs = %s" % gh_timestep_MJDs)

    # Get the MJDc value for use in computing the gamhelio frame. This value
    # was specified in the WSA file of initial conditions.
    if verbose:
        print("Reading MJDc value to use for constructing gamhelio coordinate frame.")
    gh_MJDc = scutils.read_MJDc(gh_result_path)
    if verbose:
        print("  Found MJDc = %s" % gh_MJDc)
    if debug:
        print("gh_MJDc = %s" % gh_MJDc)

    # Now get the "time" attribute from each step. For gamhelio, these elapsed
    # times are in seconds since the start of the simulation, and are assumed
    # to be in increasing order
    if verbose:
        print("Reading timestep elapsed seconds from %s." % gh_result_path)
    gh_timestep_elapsed_seconds = kaiH5.getTs(
        gh_result_path, gh_timestep_indices, aID="time"
    )
    if debug:
        print("gh_timestep_elapsed_seconds = %s" %
              gh_timestep_elapsed_seconds)

    # Convert the MJD values to Universal Time datetime objects.
    if verbose:
        print("Converting timestep MJDs to UTC datetimes.")
    gh_timestep_UT = kaiTools.MJD2UT(gh_timestep_MJDs)
    if debug:
        print("gh_timestep_UT = %s" % gh_timestep_UT)

    # Use the first (positive) elapsed time as the initial time.
    # N.B. THIS SKIPS THE FIRST SIMULATION STEP SINCE IT TYPICALLY HAS
    # gh_timestep_elapsed_seconds[0] = 0.
    # Use the last time as the last MJD.
    if verbose:
        print("Identifying first timestep with elapsed seconds > 0.")
    gh_first_step_used = np.argwhere(gh_timestep_elapsed_seconds > 0.0)[0][0]
    if verbose:
        print("  Found first non-zero elapsed seconds at timestep %s." %
              gh_first_step_used)
    gh_t0 = gh_timestep_UT[gh_first_step_used]
    gh_t1 = gh_timestep_UT[-1]
    if verbose:
        print("Using %s as the starting datetime for data comparison." % gh_t0)
        print("Using %s as the ending datetime for data comparison." % gh_t1)
    if debug:
        print("gh_t0 = %s" % gh_t0)
        print("gh_t1 = %s" % gh_t1)

    # Construct the string versions of the first and last times.
    datestr_start = gh_t0.strftime("%Y-%m-%dT%H:%M:%SZ")
    datestr_end = gh_t1.strftime("%Y-%m-%dT%H:%M:%SZ")
    if debug:
        print("datestr_start = %s" % datestr_start)
        print("datestr_end = %s" % datestr_end)

    # Save the (float) MJD of the first used step.
    gh_first_MJD = gh_timestep_MJDs[gh_first_step_used]
    if verbose:
        print("Using %s as the starting MJD for data comparison." % gh_first_MJD)
    if debug:
        print("gh_first_MJD = %s" % gh_first_MJD)

    # Save the elapsed simulation time (seconds) of the first used step.
    gh_first_elapsed_seconds = gh_timestep_elapsed_seconds[gh_first_step_used]
    if verbose:
        print("Using %s as the starting elapsed seconds for data comparison." %
              gh_first_elapsed_seconds)
    if debug:
        print("gh_first_elapsed_seconds = %s" % gh_first_elapsed_seconds)

    # Determine the list of IDs of spacecraft to fetch data from. If no
    # spacecraft were specified on the command line, use all spacecraft
    # listed in the heliosphere spacecraft metadata file.
    if sc_to_compare:
        sc_to_compare = sc_to_compare.split(",")
    else:
        sc_to_compare = list(sc_metadata.keys())
    if verbose:
        print("Comparing gamhelio results to data measured by:")
        for sc_id in sc_to_compare:
            print("  %s" % sc_id)
    if debug:
        print("sc_to_compare = %s" % sc_to_compare)

    # Fetch the ephemeris and observed data for each spacecraft in the list.
    for sc_id in sc_to_compare:

        # Fetch the ephemeris and observed data for the current spacecraft.
        if verbose:
            print("Fetching ephemeris and instrument data for %s from CDAWeb." % sc_id)
        sc_data = scutils.get_helio_cdaweb_data(
            sc_id, sc_metadata[sc_id],
            datestr_start, datestr_end, cdaweb_data_interval,
            verbose=verbose, debug=debug
        )
        if debug:
            print("sc_data = %s" % sc_data)

        # If no data was found for the spacecraft, go to the next.
        if sc_data is None:
            print("No data found for %s." % sc_id)
            continue

        # At this point, the data object contains only the raw spacecraft
        # ephemeris, and the raw spacecraft-measured data, as returned from
        # CDAWeb.

        # Ingest the CDAWeb data. This means convert the data as originally
        # retrieved from CDAWeb to the units and coordinate systems used by
        # gamhelio so that comparisons may be made. This will create new
        # variables in the data object, with names analogous to the
        # corresponding gamhelio variables.
        if verbose:
            print("Converting CDAWeb data for %s into gamhelio format." % sc_id)
        scutils.ingest_helio_cdaweb_data(
            sc_id, sc_data, sc_metadata[sc_id], gh_MJDc,
            verbose=verbose, debug=debug
        )

        # Use the spacecraft trajectory to interpolate simulated observations
        # from the gamhelio output.
        if verbose:
            print("Interpolating gamhelio results along %s trajectory." % sc_id)
        scutils.interpolate_gamhelio_results_to_trajectory(
            sc_data, sc_metadata[sc_id], sc_id,
            gh_first_MJD, gh_first_elapsed_seconds,
            gh_result_directory, gh_run_id,
            cmd_sctrack, num_segments, keep, gh_MJDc
        )

        # Save the important measured and simulated data as a CDF file for
        # comparison.
        cdf_path = os.path.join(gh_result_directory, sc_id + ".comp.cdf")
        if verbose:
            print("Saving comparison data to %s." % cdf_path)
        if os.path.exists(cdf_path):
            if verbose:
                print("Deleting existing CDF comparison file %s" % cdf_path)
            os.system("rm %s" % cdf_path)
        if verbose:
            print("Creating CDF file %s with %s and GAMERA data" % (cdf_path, sc_id))
        dm.toCDF(cdf_path, sc_data)

        # Compute the errors in the simulated data relative to the measured
        # data and save in a file.
        error_file_path = os.path.join(gh_result_directory, sc_id + "-error.txt")
        if verbose:
            print("Computing gamhelio-%s errors and saving to %s." %
                  (sc_id, error_file_path))
        scutils.write_helio_error_report(error_file_path, sc_id, sc_data)

        # Make a comparison plot of the measured and simulated results.
        plot_file_path = os.path.join(gh_result_directory, sc_id + ".png")
        if verbose:
            print("Saving gamhelio-%s comparison plots in %s." %
                  (sc_id, plot_file_path))
        kv.helioCompPlot_new(plot_file_path, sc_id, sc_data)

        # Plot the spacecraft trajectory.
        # plot_file_path = os.path.join(gh_result_directory, sc_id + "-traj.png")
        # if verbose:
        #     print("Plotting %s trajectory in spacecraft frame to %s." % (sc_id, plot_file_path))
        # kv.helioTrajPlot(plot_file_path, sc_id, sc_data)
