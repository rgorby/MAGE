#!/usr/bin/env python

"""Run MAGE Intel tests.

This script runs a series of tests of the MAGE software using Intel tools.

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
import sys
import subprocess
import time

# Import 3rd-party modules.
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = "Script for MAGE Ontel checks"


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
        "--account", default=None,
        help="PBS account to use for testing (default: %(default)s)"
    )
    parser.add_argument(
        "--all", "-a", action="store_true",
        help="Run all tests (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--force", "-f", action="store_true",
        help="Force all tests to run (default: %(default)s)."
    )
    parser.add_argument(
        "--loud", "-l", action="store_true",
        help="Enable loud mode (default: %(default)s)."
    )
    parser.add_argument(
        "--test", "-t", action="store_true",
        help="Enable testing mode (default: %(default)s)."
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
        help="Print verbose output (default: %(default)s)."
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
    slack_token = os.environ["SLACK_BOT_TOKEN"]
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
        print('I am the unit test script. This is my current home directory:')
        print(home)

    #--------------------------------------------------------------------------

    # Clean up from previous tests.
    os.chdir(home)
    os.system('rm -rf intelChecks*')

    #--------------------------------------------------------------------------

    # Run the Intel checks with each set of modules.
    
    # Go back to scripts folder
    path = os.path.join(home, 'testingScripts', 'mage_build_test_modules')
    if debug:
        print(f"path = {path}")
    os.chdir(path)

    # Get a list of build module sets.
    module_list_files = glob.glob('*.lst')
    if debug:
        print(f"module_list_files = {module_list_files}")

    # <HACK>
    # Just use first module set for now.
    module_list_files = [module_list_files[0]]
    # </HACK>

    # Run unit tests with each set of modules.
    for module_list_file in module_list_files:
        if debug:
            print(f"module_list_file = {module_list_file}")

        # Extract the name of the list.
        module_list_name = module_list_file.replace('.lst', '')
        if debug:
            print(f"module_list_name = {module_list_name}")

        # Read this module list file, extracting cmake environment and options,
        # if any.
        path = os.path.join(home, 'testingScripts', 'mage_build_test_modules',
                            module_list_file)
        with open(path, encoding="utf-8") as f:
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

        # Add the additional flags needed for intel checks.
        cmake_options += ' -DDISABLE_DEBUG_BOUNDS_CHECKS=ON'
        cmake_options += ' -DCMAKE_BUILD_TYPE=DEBUG'

        # Make a directory for this test, and go there.
        dir_name = f"intelChecks_{module_list_name}"
        build_directory = os.path.join(home, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Assemble the commands to load the listed modules.
        module_cmd = f"module --force purge; module load {' '.join(module_names)}"
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Run cmake to build the Makefile.
        cmake_cmd = f"{module_cmd}; {cmake_env} cmake {cmake_options} .."
        if debug:
            print(f"cmake_cmd = {cmake_cmd}")
        cmake_process = subprocess.Popen(cmake_cmd, shell=True)
        if debug:
            print(f"cmake_process = {cmake_process}")
        cmake_process.wait()

        # Run the build.
        make_cmd = module_cmd + '; make gamera_mpi voltron_mpi'
        if debug:
            print(f"make_cmd = {make_cmd}")
        make_process = subprocess.Popen(make_cmd, shell=True)
        if debug:
            print(f"make_process = {make_process}")
        make_process.wait()

        # Copy in the test PBS scripts and files.
        testing_files_home = '../testingScripts'
        subprocess.call(
            f"cp {testing_files_home}/tinyCase.xml bin/",
            shell=True
        )
        subprocess.call(
            f"cp {testing_files_home}/intelCheckSubmitMem.pbs bin/",
            shell=True
        )
        subprocess.call(
            f"cp {testing_files_home}/intelCheckSubmitThread.pbs bin/",
            shell=True
        )
        subprocess.call(
            f"cp {testing_files_home}/bcwind.h5 ./bin",
            shell=True
        )
        subprocess.call(
            f"cp {testing_files_home}/lfmD.h5 ./bin",
            shell=True
        )
        subprocess.call(
            f"cp {testing_files_home}/rcmconfig.h5 ./bin",
            shell=True
        )
        subprocess.call(
            f"cp {testing_files_home}/memSuppress.sup ./bin",
            shell=True
        )
        subprocess.call(
            f"cp {testing_files_home}/threadSuppress.sup ./bin",
            shell=True
        )

        # Go to the bin directory for testing.
        path = os.path.join(home, dir_name, 'bin')
        os.chdir(path)

        # Submit the job to run the Intel memory checks.
        module_list = ' '.join(module_names)
        qsub_cmd = f"qsub -A {account} -v MODULE_LIST='{module_list}',KAIJUROOTDIR={home} intelCheckSubmitMem.pbs"
        if debug:
            print(f"qsub_cmd = {qsub_cmd}")
        submission = subprocess.Popen(
            qsub_cmd, shell=True, stdout=subprocess.PIPE
        )
        submission.wait()
        readString = submission.stdout.read().decode('ascii')
        if debug:
            print(f"readString = {readString}")
        memCheckJob = readString.split('.')[0]
        if debug:
            print(f"memCheckJob = {memCheckJob}")

        # Submit the job to run the Intel thread checks.
        module_list = ' '.join(module_names)
        qsub_cmd = f"qsub -A {account} -v MODULE_LIST='{module_list}',KAIJUROOTDIR={home} intelCheckSubmitThread.pbs"
        if debug:
            print(f"qsub_cmd = {qsub_cmd}")
        submission = subprocess.Popen(
            qsub_cmd, shell=True, stdout=subprocess.PIPE
        )
        submission.wait()
        readString = submission.stdout.read().decode('ascii')
        if debug:
            print(f"readString = {readString}")
        threadCheckJob = readString.split('.')[0]
        if debug:
            print(f"threadCheckJob = {threadCheckJob}")

        # Record the job IDs.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            f.write(f"{memCheckJob}\n")
            f.write(f"{threadCheckJob}\n")

        # <HACK>
        message = f"Intel tests submitted in jobs {memCheckJob}, {threadCheckJob}"
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
                    assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

        # Send message to stdout.
        print(message)

# # SUBMIT FOLLOW-UP JOB FOR SLACK POSTING
# #os.chdir(home)
# #os.chdir('kaiju/testingScripts')
# #arguments = 'qsub intelCheckReportSubmit.pbs -W depend=after:'
# #arguments = arguments + numberString
# #print(arguments)

# # WAIT ABOUT 1 MINUTE
# #time.sleep(60)

# #report = subprocess.call(arguments, shell=True, stdout=subprocess.PIPE)

# # FINISHED

# # If not a test, send message to Slack
# #if (not isTest):
#     # Try to send Slack message
# #    try:
# #        response = client.chat_postMessage(
# #            channel="#kaijudev",
# #            text=myText,
# #        )
# #    except SlackApiError as e:
#         # You will get a SlackApiError if "ok" is False
# #        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
