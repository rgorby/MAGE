#!/usr/bin/env python


"""Run MAGE build regression tests.

This script runs a series of builds of the MAGE software using sets of modules
listed in files under:

$KAIJUHOME/testingScripts/mage_build_modules

Any file in this directory which ends in .lst will be treated as a module
list file, and used for a build test. The results of each build test are
posted to Slack.

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

# Import 3rd-party modules.
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = "Script for MAGE build testing"


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
        print('I am the build script. This is my current home directory:')
        print(home)

    #--------------------------------------------------------------------------

    # Clean up the results from previous builds.
    if verbose:
        print("Cleaning up from previous build tests.")
    os.system("rm -rf build*/ testFolder")
    if verbose:
        print('The current test directory contents are:')
        os.system('ls')

    #--------------------------------------------------------------------------

    # Do a preliminary cmake run to generate the list of executables.

    # Find the current branch.
    os.chdir(home)
    p = subprocess.Popen('git symbolic-ref --short HEAD', shell=True, stdout=subprocess.PIPE)
    git_branch = p.stdout.read().decode('ascii').rstrip()
    if debug:
        print(f"git_branch = {git_branch}")

    # Make and move to the preliminary build folder.
    os.system('mkdir testFolder')
    os.chdir('testFolder')

    # Read this module list file, extracting cmake environment and options, if any.
    path = os.path.join(home, 'testingScripts', 'mage_build_test_modules',
                        '01.lst')
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    if debug:
        print(f"lines = {lines}")
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
    module_cmd = 'module --force purge; module load'
    for module_name in module_names:
        module_cmd += f" {module_name}"
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

    # Build the list of executable targets.
    make_cmd = f"{module_cmd}; make help | grep '\.x'"
    if debug:
        print(f"make_cmd = {make_cmd}")
    listProcess = subprocess.Popen(make_cmd, shell=True, stdout=subprocess.PIPE)
    if debug:
        print(f"listProcess = {listProcess}")
    listProcess.wait()
    executableString = listProcess.stdout.read().decode('ascii')
    if debug:
        print(f"executableString = {executableString}")
    executableList = executableString.splitlines()
    if debug:
        print(f"executableList = {executableList}")

    # Remove the first four characters (dots and spaces).
    executableList = [e[4:] for e in executableList]
    if debug:
        print(f"executableList = {executableList}")

    #--------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Go to the module sets folder.
    path = os.path.join(home, 'testingScripts', 'mage_build_test_modules')
    os.chdir(path)
    if debug:
        print(f"cwd = {os.getcwd()}")

    # Read the list of  module sets to use for build tests.
    with open('build_test.lst', encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [s.rstrip() for s in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    #--------------------------------------------------------------------------

    # Do a build with each set of modules.

    # Perform a test build with each set of modules.
    for module_list_file in module_list_files:
        if debug:
            print(f"module_list_file = {module_list_file}")

        # Extract the name of the list.
        module_list_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_list_name = {module_list_name}")

        # Read this module list file, extracting cmake environment and options, if any.
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

        # Make a directory for this build, and go there.
        dir_name = f"build_{module_list_name}"
        build_directory = os.path.join(home, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Assemble the commands to load the listed modules.
        module_cmd = 'module --force purge; module load'
        for module_name in module_names:
            module_cmd += f" {module_name}"
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
        make_cmd = f"{module_cmd}; make"
        if debug:
            print(f"make_cmd = {make_cmd}")
        make_process = subprocess.Popen(make_cmd, shell=True)
        if debug:
            print(f"make_process = {make_process}")
        make_process.wait()

        # Create a test result message.
        message = '*IGNORE - TESTING*\n'
        message = f"Building MAGE branch {git_branch} with module set {module_list_file}: {module_names}\n"

        # Check for all executables
        os.chdir('bin')
        missing = []
        for executable in executableList:
            if not os.path.isfile(executable):
                missing.append(executable)
        if debug:
            print(f"missing = {missing}")
        os.chdir('..')
        isPerfect = True
        if len(missing) > 0:
            isPerfect = False
            for executable in missing:
                message += f"I couldn't build {executable}.\n"
        else:
            message += f"Everything built properly on branch {git_branch} with module set {module_list_file}!"
        if debug:
            print(f"message = {message}")

        # If this is a test run, don't post to Slack.
        if isTest:
            pass
        else:
            # If loud, or an error occurred, send Slack message.
            if beLoud or not isPerfect:
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

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
