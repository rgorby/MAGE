#!/usr/bin/env python

"""Report on the MAGE Fortran unit test results.

This script generates a report on the results of the most recent run of
the MAGE Fortran unit tests.

This script *assumes* that enough time has passed that the unit test runs
(which were submitted as PBS jobs) have completed. Those *should* take
around 20 minutes. Add that to the build time for the unitTest.py script,
and a conservative guess is to run this script 2 hours after running the
unitTest.py script.

Authors
-------
Jeff Garretson
Eric Winter
"""


# Import standard modules.
import datetime
import glob
import os
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Report on the MAGE Fortran unit test results.'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# glob pattern for naming unit test directories
UNIT_TEST_DIRECTORY_GLOB_PATTERN = 'unitTest_*'

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = 'bin'

# Name of file containing job IDs for each unit test directory.
JOB_ID_LIST_FILE = 'jobs.txt'

# Name of file to receive unit test report.
UNIT_TEST_REPORT_FILE = 'Report.txt'


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
    None
    """
    # Set up the command-line parser.
    parser = common.create_command_line_parser(DESCRIPTION)

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    debug = args.debug
    be_loud = args.loud
    is_test = args.test
    verbose = args.verbose

    # -------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # -------------------------------------------------------------------------

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # -------------------------------------------------------------------------

    # Move to the MAGE installation directory.
    os.chdir(KAIJUHOME)

    # -------------------------------------------------------------------------

    # Find the current branch.
    git_branch_name = common.git_get_branch_name()
    if debug:
        print(f"git_branch_name = {git_branch_name}")

    # -------------------------------------------------------------------------

    # Get list of unit test directories.
    unit_test_directories = glob.glob(UNIT_TEST_DIRECTORY_GLOB_PATTERN)
    if debug:
        print(f"unit_test_directories = {unit_test_directories}")

    # -------------------------------------------------------------------------

    # Initialize result flags.
    myError = False
    jobKilled = False
    okFailure = False
    okCount = 0

    # Check the results in each unit test directory.
    for unit_test_directory in unit_test_directories:
        if verbose:
            print(f"Checking unit test results in {unit_test_directory}.")

        # Move back to the kaiju home.
        os.chdir(KAIJUHOME)

        # Move to the directory containing the unit test results.
        path = os.path.join(unit_test_directory, BUILD_BIN_DIR)
        if debug:
            print(f"path = {path}")
        os.chdir(path)

        # Check for the job ID list file. If not found, skip the rest of this
        # loop.
        if not os.path.isfile(JOB_ID_LIST_FILE):
            print(
                f"Job list file {JOB_ID_LIST_FILE} in {unit_test_directory}"
                ' not found, skipping report for this directory.'
            )
            continue

        # Read the job IDs from the job ID list file.
        with open(JOB_ID_LIST_FILE, 'r', encoding='utf-8') as f:
            lines = f. readlines()
        job_ids = [_.rstrip() for _ in lines]
        if debug:
            print(f"job_ids = {job_ids}")

        # <HACK>
        # This needs to be reorganized.

        # Compute the names of the job log files.
        job_file_0 = f"testResGen.o{job_ids[0]}"
        job_file_1 = f"caseTests.o{job_ids[1]}"
        job_file_2 = f"nonCaseTests1.o{job_ids[2]}"
        # job_file_3 = f"nonCaseTests2.o{job_ids[3]}"  # SKIP FOR NOW
        if debug:
            print(f"job_file_o = {job_file_0}")
            print(f"job_file_1 = {job_file_1}")
            print(f"job_file_2 = {job_file_2}")
            # print(f"job_file_3 = {job_file_3}")

        # Combine the results of each test log file.
        bigFile = []
        # job_files = [job_file_0, job_file_1, job_file_2, job_file_3]
        job_files = [job_file_0, job_file_1, job_file_2]
        for job_file in job_files:
            with open(job_file, 'r', encoding='utf-8') as f:
                bigFile += f.readlines()
            bigFile.append('\n\n\n')

        # Scan through for some key things like "error" and "job killed"
        for line in bigFile:
            line = line.rstrip()
            if line == ' OK':
                okCount += 1
            if 'error' in line:
                myError = True
            elif 'job killed' in line:
                jobKilled = True

        # There should be exactly 6 OKs (8 if job_file_3 is used).
        if okCount != 6:
            okFailure = True
        else:
            okFailure = False

        if debug:
            print(f"myError = {myError}")
            print(f"jobKilled = {jobKilled}")
            print(f"okFailure = {okFailure}")
            print(f"okCount = {okCount}")

        # End of unit test directory loop

    # Move back to the kaiju home.
    os.chdir(KAIJUHOME)

    # If no tests were complete, say so.
    if len(unit_test_directories) == 0:
        print(
            'Results of unit test report `unitTestReport.py`:\n'
            f"kaiju code branch: `{git_branch_name}`\n"
            'No Fortran unit test results were available.'
        )

    # ------------------------------------------------------------------------

    # NOTE: Assumes only 1 module set was used.

    # Detail the test results
    test_details_message = ''
    if myError:
        test_details_message += 'Errors occurred during testing.\n'
    if jobKilled:
        test_details_message += 'The testing job was killed.\n'
    if okFailure:
        test_details_message += 'There was not the correct OK count.\n'
    else:
        test_details_message += 'Individual tests: *ALL PASSED*\n'

    # Summarize the test results.
    test_summary_message = (
        'Summary of Fortran unit test results from `unitTestReport.py`: '
    )
    if myError or jobKilled or okFailure:
        test_summary_message += '*FAILED*\n'
    else:
        test_summary_message += '*ALL PASSED*\n'

    # Print the test results summary and details.
    print(test_summary_message)
    print(test_details_message)

    # If loud mode is on, post report to Slack.
    if be_loud:
        test_summary_message += 'Details in thread for this messsage.\n'
        slack_response_summary = common.slack_send_message(
            slack_client, test_summary_message, is_test=is_test
        )
        if slack_response_summary['ok']:
            thread_ts = slack_response_summary['ts']
            slack_response_details = common.slack_send_message(
                slack_client, test_details_message, thread_ts=thread_ts,
                is_test=is_test
            )
        else:
            print('*ERROR* Unable to post test result summary to Slack.')

    # -------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    main()
