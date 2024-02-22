#!/usr/bin/env python

"""Create a report for the results of the MAGE pythin unit tests.

Create a report for the results of the MAGE pythin unit tests.

Authors
-------
Jeff Garretson (jeffrey.garretson@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import datetime
import logging
import os
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Generate a report for MAGE python unit test results.'

# Name of directory for running python unit tests
PYTHON_UNIT_TEST_DIRECTORY = 'pytests'

# Name of python unit test log file
PYTHON_UNIT_TEST_LOG_FILE = 'kaiju-pyunit.txt'


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
    # Turn on logging.
    logging.basicConfig(level=logging.DEBUG, filename='pyunitReport.log')

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

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    #--------------------------------------------------------------------------

    # Set up for communication with Slack.

    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    #--------------------------------------------------------------------------

    # Move to the MAGE installation directory.
    kaiju_home = os.environ['KAIJUHOME']
    os.chdir(kaiju_home)

    #--------------------------------------------------------------------------

    # Move to the pytests folder.
    work  = os.path.join(kaiju_home, PYTHON_UNIT_TEST_DIRECTORY)
    os.chdir(work)

    # Check that the python unit tests completed.
    if os.path.exists(PYTHON_UNIT_TEST_LOG_FILE):
        with open(PYTHON_UNIT_TEST_LOG_FILE, encoding='utf-8') as f:
            lines = f.readlines()
        last_line = lines[-1]
    else:
        message = 'Python Unit Tests job did not complete\n\n'
        if not is_test and be_loud:
            common.slack_send_message(slack_client, message)
        print("Python Unit Tests job did not complete")
        sys.exit(0)

    # Examine the tests results for failures.
    has_fail = False
    has_error = False
    has_pass = False
    if 'fail' in last_line:
        has_fail = True
    if 'error' in last_line:
        has_error = True
    if 'passed' in last_line:
        has_pass = True

    # Construct message based on results.
    if is_test:
        # This is a test run so don't post to Slack.
        if has_error:
            print("Python Unit Tests did not run sucessfully")
        if has_fail:
            print("Python Unit Tests Failed")
        if has_pass and not has_error and not has_fail:
            print("Python Unit Tests Passed")
        if not has_pass and not has_error and not has_fail:
            print('Unexpected error occured during python unit tests\n\n')
    elif be_loud: # Send all messages to Slack
        # Send all messages to Slack.
        if has_error:
            common.slack_send_message(
                slack_client, 'Python Unit Tests did not run succesfully\n\n')
        if has_fail:
            common.slack_send_message(
                slack_client, 'Python Unit Tests Failed\n\n')
        if has_pass and not has_error and not has_fail:
            common.slack_send_message(slack_client, 'Python Unit Tests Pass')
        if not has_pass and not has_error and not has_fail:
            common.slack_send_message(
                slack_client, 'Unexpected error occured during python unit tests\n\n')
    else:
        # Only report errors to Slack
        if has_error:
            common.slack_post_message(
                slack_client, 'Python Unit Tests did not run succesfully\n\n')
        if has_fail:
            common.slack_send_message(
                slack_client, 'Python Unit Tests Failed\n\n')
        if not has_pass and not has_error and not has_fail:
            common.slack_send_message(
                slack_client, 'Unexpected error occured during python unit tests\n\n')

    #--------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
