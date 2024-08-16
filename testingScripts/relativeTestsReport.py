#!/usr/bin/env python

"""Create the MAGE comparative test report.

This script creates the MAGE comparative test report.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import glob
import os
import platform
import re
import subprocess
import sys

# Import 3rd-party modules.
from astropy.time import Time
import matplotlib as mpl
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

# Import project modules.
import common
import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv


# Program constants

# Program description.
DESCRIPTION = 'Create the MAGE comparative test report.'

# Root of directory tree for all tests.
MAGE_TEST_ROOT = os.environ['MAGE_TEST_ROOT']

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for unit tests
TEST_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'compTest')

# Glob pattern for individual weekly dash directories
COMP_TESTS_DIRECTORY_GLOB_PATTERN = 'compTest_*'

# Regular expression for git hash read from weekly dash output log.
GIT_HASH_PATTERN = 'Git hash   = (.{8})'

# Name of subdirectory containing binaries and test results.
BIN_DIR = 'bin'

# Name of file containg PBS job IDs.
JOB_LIST_FILE = 'jobs.txt'

# String naming branch or commit used in this test.
BRANCH_OR_COMMIT = os.environ['BRANCH_OR_COMMIT']

# List of basic test videos to generate
BASIC_TESTS = [
    {
        'case1':'relSerialRelease',
        'case2':'relMpi44Release',
        'ts':-60,
        'te':120,
        'videoName':'serialVsMpi',
        'postMessage':'Comparing Non-MPI case to Mpi4x4 Case\n\n'
    },
    {
        'case1':'relMpi44Release',
        'case2':'relMpi44ResRelease',
        'ts':-60,
        'te':120,
        'videoName':'MpiRestart',
        'postMessage':'Comparing Mpi4x4 Case to a Similar Case Restarted at 60 Minutes\n\n'
    }
]

# List of additional videos to generate if performing full test suite
EXPANDED_TESTS = [
]

def processCases(testDict):
    if verbose:
        print(f"Creating video comparing case {testDict['case1']} to {testDict['case2'}")
    
    # Create the video
    cmd = f"gamerrVid.py -d1 {testDict['case1']} -d2 {testDict['case2']} -o {testDict['videoName']} "
        f" -ts {testDict['ts']} -te {testDict['te']} -dt 60 -Nth 36"
    if debug:
        print(f"cmd = {cmd}")
    try:
        _ = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'ERROR: Unable to create comparative video.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f'See log for output.\n',
            file=sys.stderr
        )


def main():
    """Begin main program.

    This is the main program code.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Raises
    ------
    subprocess.CalledProcessError
        If an exception occurs in subprocess.run()
    """
    # Set up the command-line parser.
    parser = common.create_command_line_parser(DESCRIPTION)

    # Add additional arguments
    parser.add_argument(
        '--allTests', '-a', action='store_true',
        help='Perform all comparative tests (default: %(default)s).'
    )

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    debug = args.debug
    be_loud = args.loud
    # slack_on_fail = args.slack_on_fail
    is_test = args.test
    verbose = args.verbose
    allTests = args.allTests

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Move to the top-level weekly dash directory.
    os.chdir(TEST_DIRECTORY)

    # ------------------------------------------------------------------------

    # Get list of weekly dash directories.
    test_directories = glob.glob(COMP_TESTS_DIRECTORY_GLOB_PATTERN)
    if debug:
        print(f"comp_test_directories = {test_directories}")

    # <HACK>
    # Use only the first mdirectory for now.
    test_directory = test_directories[0]
    # </HACK>

    # ------------------------------------------------------------------------

    # Read results from the latest run.
    if verbose:
        print(f"Reading results for latest run in {test_directory}.")

    # Go to test folder
    os.chdir(test_directory)

    # Move down to the build directory
    os.chdir(BIN_DIR)

    # ------------------------------------------------------------------------

    # Generate the comparative videos

    for case in BASIC_TESTS:
        processCases(case)

    if allTests:
        for case in EXPANDED_TESTS:
            processCases(case)

    # -----------------------------------------------------------------------

    # If loud mode is on, post results to Slack.
    if be_loud:
        slack_client = common.slack_create_client()
        if debug:
            print(f"slack_client = {slack_client}")
        message = (
            'Comparative test result plots complete on branch '
            f"{BRANCH_OR_COMMIT}.\n"
            ' Latest  results attached as replies to this '
            'message.\n'
        )
        message += (
            f"Test results are in {os.getcwd()}.\n"
        )
        slack_response = common.slack_send_message(
            slack_client, message, is_test=is_test)
        if slack_response['ok']:
            parent_ts = slack_response['ts']
            message = (
                'All results are from a Quad Resolution Run using'
                ' various Gamera and Voltron settings.'
            )
            slack_response = common.slack_send_message(
                slack_client, message, thread_ts=parent_ts, is_test=is_test)
            for case in BASIC_TESTS:
                slack_response = common.slack_send_image(
                    slack_client, case['videoName']+'.mp4', initial_comment=case['postMessage'],
                    thread_ts=parent_ts, is_test=is_test
                )
            if allTests:
                for case in EXPANDED_TESTS:
                    slack_response = common.slack_send_image(
                        slack_client, case['videoName']+'.mp4', initial_comment=case['postMessage'],
                        thread_ts=parent_ts, is_test=is_test
                    )

        else:
            print('Failed to post parent message and images to Slack.')

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}")


if __name__ == '__main__':
    main()
