#!/usr/bin/env python


"""Run MAGE build regression tests.

This script runs a series of builds of the MAGE software using sets of modules
listed in files under:

$KAIJUHOME/testingScripts/mage_build_modules

This script reads the file build_test.lst from this directory, and
uses the contents as a list of module list files to use for MAGE build
tests.

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
import shutil
import subprocess
import sys

# Import 3rd-party modules.
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE build testing'


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
        '--account', default='UJHB0019',
        help="PBS account to use for testing (default: %(default)s)"
    )
    parser.add_argument(
        '--all', '-a', action='store_true',
        help="Run all tests (default: %(default)s)."
    )
    parser.add_argument(
        '--debug', '-d', action='store_true',
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        '--force', '-f', action='store_true',
        help="Force all tests to run (default: %(default)s)."
    )
    parser.add_argument(
        '--loud', '-l', action='store_true',
        help="Enable loud mode (post results to Slack) (default: %(default)s)."
    )
    parser.add_argument(
        '--test', '-t', action='store_true',
        help="Enable testing mode (default: %(default)s)."
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true',
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

    # Assume this script is in a subdirectory of the kaiju directory,
    # and go there.
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
        print('Cleaning up from previous build tests.')
    directories = glob.glob('build_*')
    directories.append('testFolder')
    for directory in directories:
        try:
            shutil.rmtree(directory)
        except:
            pass
    # <HACK>
    # Remove the pFUnit compiled code to prevent using it during the
    # build test. If PFUNIT-4.2 is in kaiju/external during a build,
    # make will try to build the unit test code even if it is not
    # requested, which causes fatal errors when building with a module
    # set that uses a non-Intel compioler, since pFUnit was built with
    # the Intel compiler.
    pfunit_binary_directories = [
        'FARGPARSE-1.1',
        'GFTL-1.3',
        'GFTL_SHARED-1.2',
        'PFUNIT-4.2',
    ]
    for directory in pfunit_binary_directories:
        path = os.path.join(home, 'external', directory)
        try:
            shutil.rmtree(path)
        except:
            pass
    # </HACK>
    if verbose:
        print('The current test directory contents are:')
        os.system('ls')

    #--------------------------------------------------------------------------

    # Do a preliminary cmake run to generate the list of executables.

    # Find the current branch.
    os.chdir(home)
    cmd = 'git symbolic-ref --short HEAD'
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    git_branch = cproc.stdout.rstrip()
    if debug:
        print(f"git_branch = {git_branch}")

    # Make and move to the preliminary build folder.
    os.mkdir('testFolder')
    os.chdir('testFolder')

    # Read the first module list file, extracting cmake environment
    # and options, if any.
    path = os.path.join(home, 'testingScripts', 'mage_build_test_modules',
                        '01.lst')
    with open(path, encoding='utf-8') as f:
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
    module_cmd = 'module --force purge; module load '
    module_cmd += ' '.join(module_names)
    if debug:
        print(f"module_cmd = {module_cmd}")

    # Run cmake to build the Makefile.
    cmd = f"{module_cmd}; {cmake_env} cmake {cmake_options} .."
    if debug:
        print(f"cmd = {cmd}")
    # <HACK> To ignore cmake error on bcwind.h5 for now.
    try:
        cproc = subprocess.run(cmd, shell=True, check=True, text=True)
    except:
        pass
    # </HACK>

    # Build the list of executable targets.
    cmd = f"{module_cmd}; make help | grep '\.x'"
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    executableList = cproc.stdout.splitlines()
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
    for module_list_file in module_list_files:
        if debug:
            print(f"module_list_file = {module_list_file}")

        # Extract the name of the list.
        module_list_name = module_list_file.rstrip('.lst')
        if verbose:
            print(f"Starting build with module list {module_list_name}.")

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

        # Make a directory for this build, and go there.
        dir_name = f"build_{module_list_name}"
        build_directory = os.path.join(home, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Assemble the commands to load the listed modules.
        module_cmd = 'module --force purge; module load '
        module_cmd += ' '.join(module_names)
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Run cmake to build the Makefile.
        cmd = f"{module_cmd}; {cmake_env} cmake {cmake_options} .."
        if debug:
            print(f"cmd = {cmd}")
        # <HACK> To ignore cmake error on bcwind.h5 for now.
        try:
            cproc = subprocess.run(cmd, shell=True, check=True, text=True)
        except:
            pass
        # </HACK>

        # Run the build.
        cmd = f"{module_cmd}; make {' '.join(executableList)}"
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, check=True, text=True)

        # Create a test result message.
        message = (
            f"Built MAGE branch {git_branch} with module set "
            f"{module_list_file}: {module_names}\n"
        )

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
            message += (
                f"Everything built properly on branch {git_branch} "
                f"with module set {module_list_file}!"
            )
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
                        channel='#kaijudev',
                        text=message,
                    )
                except SlackApiError as e:
                    # You will get a SlackApiError if "ok" is False
                    assert e.response['error']  # str like 'invalid_auth', 'channel_not_found'

        # Send message to stdout.
        print(message)

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
