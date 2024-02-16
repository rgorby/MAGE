#!/usr/bin/env python

"""Run MAGE initial condition tests.

This script runs a series of initial condition tests of the MAGE software.

Authors
-------
Jeff Garretson (jeffrey.garretson@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import datetime
import glob
import os
import subprocess
import sys

# Import 3rd-party modules.
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE initial condition testing'


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.

    Raises
    ------
    None
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--account', default=None,
        help='PBS account to use for testing (default: %(default)s)'
    )
    parser.add_argument(
        '--all', '-a', action='store_true',
        help='Run all tests (default: %(default)s).'
    )
    parser.add_argument(
        '--debug', '-d', action='store_true',
        help='Print debugging output (default: %(default)s).'
    )
    parser.add_argument(
        '--force', '-f', action='store_true',
        help='Force all tests to run (default: %(default)s).'
    )
    parser.add_argument(
        '--loud', '-l', action='store_true',
        help='Enable loud mode (post results to Slack) (default: %(default)s).'
    )
    parser.add_argument(
        '--test', '-t', action='store_true',
        help='Enable testing mode (default: %(default)s).'
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Print verbose output (default: %(default)s).'
    )
    return parser


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
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    account = args.account
    doAll = args.all
    debug = args.debug
    forceRun = args.force
    beLoud = args.loud
    isTest = args.test
    verbose = args.verbose

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")

    #--------------------------------------------------------------------------

    # Set up for communication with Slack.

    # Get the Slack API token
    slack_token = os.environ['SLACK_BOT_TOKEN']
    if debug:
        print(f"slack_token = {slack_token}")

    # Create the Slack client.
    slack_client = WebClient(token=slack_token)
    if debug:
        print(f"slack_client = {slack_client}")

    #--------------------------------------------------------------------------

    # Determine the path to the MAGE installation to use for testing.

    # Fetch the path to this running script.
    called_from = os.path.dirname(os.path.abspath(__file__))
    if debug:
        print(f"called_from = {called_from}")

    # Assume this script is in a subdirectory of the kaiju directory.
    os.chdir(called_from)
    os.chdir('..')

    # Use this directory as the home directory for testing.
    home = os.getcwd()
    if debug:
        print(f"home = {home}")
    if verbose:
        print('I am the initial condition test script. This is my current '
              'home directory:')
        print(home)

    #--------------------------------------------------------------------------

    # Clean up from previous tests.
    if verbose:
        print('Cleaning up from previous initial condition test builds.')
    os.chdir(home)
    os.system('rm -rf ICBuilds')

    #--------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Go to the module sets folder.
    path = os.path.join(home, 'testingScripts', 'mage_build_test_modules')
    os.chdir(path)
    if debug:
        print(f"cwd = {os.getcwd()}")

    # Read the list of  module sets to use for Intel checks.
    with open('initial_condition_checks.lst', encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [s.rstrip() for s in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    #--------------------------------------------------------------------------
    
    # Get a list of initial conditions to try, ignoring files in the
    # "deprecated" folder.
    # GAMERA ONLY FOR NOW
    os.chdir(home)
    if debug:
        print(f"Now in directory {os.getcwd()}.")
    os.chdir('src/gamera/ICs')
    if debug:
        print(f"Now in directory {os.getcwd()}.")
    initial_condition_directory = os.getcwd() + '/'
    if debug:
        print(f"initial_condition_directory = {initial_condition_directory}")
    initial_condition_paths = []
    for root, directories, filenames in os.walk(initial_condition_directory):
        if 'deprecated' not in root and 'underdev' not in root:
            for filename in filenames:
                initial_condition_paths.append(os.path.join(root, filename))
    if debug:
        print(f"initial_condition_paths = {initial_condition_paths}")

    #--------------------------------------------------------------------------

    # Run the initial condition tests with each set of modules.

    # Create the root directory for the initial condition tests.
    initial_condition_test_root = os.path.join(home, 'ICBuilds')
    if debug:
        print(f"initial_condition_test_root = {initial_condition_test_root}")
    os.mkdir(initial_condition_test_root)

    # Run initial condition tests with each set of modules.
    for module_list_file in module_list_files:
        if verbose:
            print('Performing initial condition checks with module set '
                  f"{module_list_file}.")

        # Extract the name of the list.
        module_list_name = module_list_file.replace('.lst', '')
        if debug:
            print(f"module_list_name = {module_list_name}")

        # Read this module list file, extracting cmake environment and
        # options, if any.
        path = os.path.join(home, 'testingScripts', 'mage_build_test_modules',
                            module_list_file)
        with open(path, encoding='utf-8') as f:
            lines = f.readlines()
        cmake_env = ''
        label = 'CMAKE_ENV='
        if lines[0].startswith(label):
            cmake_env = lines[0][len(label):].rstrip()
            lines.pop(0)  # Remove cmake environment line.
        cmake_options = ''
        label = 'CMAKE_OPTIONS='
        if lines[0].startswith(label):
            cmake_options = lines[0][len(label):].rstrip()
            lines.pop(0)  # Remove cmake options line.
        module_names = [line.rstrip() for line in lines]
        if debug:
            print(f"module_names = {module_names}")

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Test each initial condition.
        for initial_condition_path in initial_condition_paths:

            # Extract the initial condition name.
            initial_condition_name = os.path.basename(initial_condition_path)
            if debug:
                print(f"initial_condition_name = {initial_condition_name}")

            # Make a directory for this test, and go there.
            build_directory = os.path.join(
                home, 'ICBuilds',
                f"gamera_{initial_condition_name}_{module_list_name}"
            )
            if debug:
                print(f"build_directory = {build_directory}")
            os.mkdir(build_directory)
            os.chdir(build_directory)
            if debug:
                print(f"Now in directory {os.getcwd()}.")

            # Assemble the cmake options for this initial condition.
            cmake_options = (
                f"{cmake_options} -DGAMIC:FILEPATH={initial_condition_path}"
            )
            if debug:
                print(f"cmake_options = {cmake_options}")

            # Run cmake to build the Makefile.
            cmd = (
                f"{module_cmd}; {cmake_env} cmake {cmake_options} ../.."
            )
            if debug:
                print(f"cmd = {cmd}")
            cproc = subprocess.run(cmd, shell=True, check=True, text=True)

            # Run the build.
            cmd = f"{module_cmd}; make gamera.x"
            if debug:
                print(f"cmd = {cmd}")
            cproc = subprocess.run(cmd, shell=True, check=True, text=True)

            # <HACK>
            message = (
                f"Build using module set {module_list_name} for "
                f"initial conditions {initial_condition_name} "
                f"returned {make_return_code}."
            )
            # </HACK>

            # If this is a test run, don't post to Slack.
            if isTest:
                pass
            else:
                # If loud, or an error occurred, send Slack message.
                if beLoud:
                    try:
                        response = slack_client.chat_postMessage(
                            channel="#kaijudev",
                            text=message,
                        )
                    except SlackApiError as e:
                        # You will get a SlackApiError if "ok" is False
                        # str like 'invalid_auth', 'channel_not_found'
                        assert e.response["error"]

            # Send message to stdout.
            print(message)

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
