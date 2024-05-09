#!/usr/bin/env python

"""Create the report for MAGE tests run using Intel Inspector.

This script creates the report for MAGE tests run using Intel Inspector.

This script assumes it is run in the directory containing the test results.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import re
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Create report for Intel Inspector tests.'


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

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    debug = args.debug
    be_loud = args.loud
    is_test = args.test

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Check for for the job list file.
    if not os.path.exists('jobs.txt'):
        print('Nothing to report on.')
        sys.exit(0)

    # Read in the jobs.txt file to get the job numbers
    with open('jobs.txt', 'r', encoding='utf-8') as f:
        job1 = f.readline().rstrip()
        job2 = f.readline().rstrip()
    if debug:
        print(f"job1 = {job1}")
        print(f"job2 = {job2}")

    # Combine the job output files.
    jobFile1 = f"memCheck.o{job1}"
    jobFile2 = f"threadCheck.o{job2}"
    if debug:
        print(f"jobFile1 = {jobFile1}")
        print(f"jobFile2 = {jobFile2}")

    # Check the the logs exist.
    if not os.path.exists(jobFile1):
        print(f"{jobFile1} not found, tests are incomplete.")
        sys.exit(0)
    if not os.path.exists(jobFile2):
        print(f"{jobFile2} not found, tests are incomplete.")
        sys.exit(0)

    # Set up regex for counting problems found.
    problemPattern = '(\d+) new problem\(s\) found'

    # ------------------------------------------------------------------------

    # Search the output logs for errors.

    # Memory
    memErrs = False
    memErrsFile = 'combinedMemErrs.txt'
    for root, dirs, files in os.walk('.'):
        for d in dirs:
            if 'memResults' in d:
                try:
                    memOut = subprocess.check_output(
                        ['inspxe-cl', '-report summary', '-result-dir ' + d,
                         '-s-f memSuppress.sup'],
                        stderr=subprocess.STDOUT, universal_newlines=True)
                except subprocess.CalledProcessError as memProcErr:
                    # we need to handle non-zero error code
                    memOut = memProcErr.output
                problemMatch = re.search(problemPattern, memOut)
                if not problemMatch or int(problemMatch.group(1)) > 0:
                    memErrs = True
                    try:
                        memOut = subprocess.check_output(
                            ['inspxe-cl', '-report problems',
                             '-result-dir ' + d, '-s-f memSuppress.sup',
                             '-report-all'],
                            stderr=subprocess.STDOUT, universal_newlines=True)
                    except subprocess.CalledProcessError as memProcErr:
                        # we need to handle non-zero error code
                        memOut = memProcErr.output
                    with open(memErrsFile, 'a', encoding='utf-8') as memFile:
                        memFile.write(memOut)
                        memFile.write('\n')

    # Thread
    threadErrs = False
    threadErrsFile = 'combinedThreadErrs.txt'
    for root, dirs, files in os.walk('.'):
        for d in dirs:
            if 'threadResults' in d:
                try:
                    threadOut = subprocess.check_output(
                        ['inspxe-cl', '-report summary', '-result-dir ' + d,
                         '-s-f threadSuppress.sup'],
                        stderr=subprocess.STDOUT, universal_newlines=True
                    )
                except subprocess.CalledProcessError as threadProcErr:
                    # we need to handle non-zero error code
                    threadOut = threadProcErr.output
                problemMatch = re.search(problemPattern, threadOut)
                if not problemMatch or int(problemMatch.group(1)) > 0:
                    threadErrs = True
                    try:
                        threadOut = subprocess.check_output([
                            'inspxe-cl', '-report problems',
                            '-result-dir ' + d, '-s-f threadSuppress.sup',
                            '-report-all'],
                            stderr=subprocess.STDOUT, universal_newlines=True)
                    except subprocess.CalledProcessError as threadProcErr:
                        # we need to handle non-zero error code
                        threadOut = threadProcErr.output
                    with open(
                        threadErrsFile, 'a', encoding='utf-8'
                    ) as threadFile:
                        threadFile.write(threadOut)
                        threadFile.write('\n')

    # -----------------------------------------------------------------------

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    test_report_details_string += 'Results of memory tests:\n'
    with open(memErrsFile, 'r', encoding='utf-8') as f:
        test_report_details_string += f.readlines()
    test_report_details_string += '\n'
    test_report_details_string += 'Results of thread tests:\n'
    with open(threadErrsFile, 'r', encoding='utf-8') as f:
        test_report_details_string += f.readlines()

    # Summarize the test results
    test_report_summary_string = (
        'Intel Inspector tests complete (`intelChecksReport.py`).'
    )

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If loud mode is on, post report to Slack.
    if be_loud:
        test_report_summary_string += 'Details in thread for this messsage.\n'
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_summary_string, is_test=is_test
        )
        if slack_response_summary['ok']:
            thread_ts = slack_response_summary['ts']
            slack_response_details = common.slack_send_message(
                slack_client, test_report_details_string, thread_ts=thread_ts,
                is_test=is_test
            )
            if 'ok' not in slack_response_details:
                print('*ERROR* Unable to post test details to Slack.')
        else:
            print('*ERROR* Unable to post test summary to Slack.')

    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    main()
