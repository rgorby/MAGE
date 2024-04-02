#!/usr/bin/env python

"""Run MAGE initial condition build tests.

This script runs a series of initial condition build tests of the MAGE
software. Each set of initial conditions provided in the kaiju source
tree is tested with each set of modules used for these tests. The module
sets are listed in files under:

$KAIJUHOME/testingScripts/mage_build_test_modules

This script reads the file initial_condition_build_test.lst from this
directory, and uses the contents as a list of module list files to use for
MAGE iinitial condition build tests.

NOTE: These tests are performed on a load-balance-assigned login node on
derecho. No PBS job is submitted.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import shutil
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE initial condition build testing'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Path directory for initial condition build tests
IC_BUILD_TEST_DIRECTORY = os.path.join(KAIJUHOME, 'ICBuilds')

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Path to file containing names of module lists to use for initial
# condition build tests
# Path to module list file to use when generating the list of executables
IC_BUILD_TEST_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY,
                                       'initial_condition_build_test.lst')

# Path to directory containing initial condition source code
IC_SRC_DIRECTORY = os.path.join(KAIJUHOME, 'src', 'gamera', 'ICs')

# Prefix for subdirectory names for individual builds.
IC_BUILD_DIR_PREFIX = 'gamera_'

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = 'bin'


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

    # Clean up from previous tests.
    if verbose:
        print('Cleaning up from previous tests.')
    directories = [IC_BUILD_TEST_DIRECTORY]
    for directory in directories:
        try:
            shutil.rmtree(directory)
        except FileNotFoundError:
            pass  # These directories may not exist.

    # <HACK>
    # Remove the pFUnit compiled code to prevent using it during the
    # build test. If PFUNIT-4.2 is in kaiju/external during a build,
    # make will try to build the unit test code even if it is not
    # requested, which causes fatal errors when building with a module
    # set that uses a non-Intel compioler, since pFUnit was built with
    # the Intel compiler.
    directories = [
        'FARGPARSE-1.1',
        'GFTL-1.3',
        'GFTL_SHARED-1.2',
        'PFUNIT-4.2',
    ]
    for directory in directories:
        path = os.path.join(KAIJUHOME, 'external', directory)
        try:
            shutil.rmtree(path)
        except FileNotFoundError:
            pass  # These directories may not exist.
    # </HACK>

    # -------------------------------------------------------------------------

    # Find the current branch.
    git_branch_name = common.git_get_branch_name()
    if debug:
        print(f"git_branch_name = {git_branch_name}")

    # -------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    with open(IC_BUILD_TEST_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # -------------------------------------------------------------------------

    # Create the top-level directory for the initial condition tests.
    os.mkdir(IC_BUILD_TEST_DIRECTORY)

    # -------------------------------------------------------------------------

    # Get a list of initial conditions to try, ignoring files in the
    # "deprecated" folder. GAMERA ONLY FOR NOW.
    initial_condition_paths = []
    for root, directories, filenames in os.walk(IC_SRC_DIRECTORY):
        if 'deprecated' not in root and 'underdev' not in root:
            for filename in filenames:
                initial_condition_paths.append(os.path.join(root, filename))
    if debug:
        print(f"initial_condition_paths = {initial_condition_paths}")

    # -------------------------------------------------------------------------

    # Build using each set of modules and each initial condition.

    # Initalize test results for all module sets and initial conditions to
    # False (failed).
    test_passed = []
    for _ in module_list_files:
        test_passed.append([False]*len(initial_condition_paths))

    # Define the make command for each build.
    make_cmd = 'make gamera.x'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Build with each initial condition with each set of modules.
    for (i_test, module_list_file) in enumerate(module_list_files):
        if verbose:
            print('Performing initial condition build tests with module set '
                  f"{module_list_file}.")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_set_name = {module_set_name}.")

        # Read this module list file, extracting cmake environment and
        # options, if any.
        path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
        if debug:
            print(f"path = {path}")
        module_names, cmake_environment, cmake_options = (
            common.read_build_module_list_file(path)
        )
        if debug:
            print(f"module_names = {module_names}")
            print(f"cmake_environment = {cmake_environment}")
            print(f"cmake_options = {cmake_options}")

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge"
            f"; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Build with each initial condition.
        for (j_ic, initial_condition_path) in enumerate(initial_condition_paths):
            if verbose:
                print(f"Building with module set {module_set_name} and IC "
                      f"file {initial_condition_path}")

            # Extract the initial condition name.
            initial_condition_name = os.path.basename(initial_condition_path)
            initial_condition_name = initial_condition_name.rstrip('.F90')
            if debug:
                print(f"initial_condition_name = {initial_condition_name}")

            # Make a directory for this test, and go there.
            build_directory = os.path.join(
                IC_BUILD_TEST_DIRECTORY,
                f"{IC_BUILD_DIR_PREFIX}{initial_condition_name}"
                f"_{module_set_name}"
            )
            if debug:
                print(f"build_directory = {build_directory}")
            os.mkdir(build_directory)
            os.chdir(build_directory)

            # Add cmake options for initial condition test builds.
            IC_cmake_options = (
                f"{cmake_options} -DGAMIC:FILEPATH={initial_condition_path}"
            )
            if debug:
                print(f"IC_cmake_options = {IC_cmake_options}")

            # Run cmake to build the Makefile.
            if verbose:
                print(
                    'Running cmake to create Makefile for module set'
                    f" {module_set_name},"
                    f" initial condition {initial_condition_name}."
                )
            cmd = (
                f"{module_cmd}"
                f"; {cmake_environment} cmake {IC_cmake_options}"
                f" {KAIJUHOME}"
                '>& cmake.out'
            )
            if debug:
                print(f"cmd = {cmd}")
            try:
                _ = subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(
                    f"ERROR: cmake for module set {module_set_name}, "
                    f"initial condition {initial_condition_name} failed.\n"
                    f"e.cmd = {e.cmd}\n"
                    f"e.returncode = {e.returncode}\n"
                    f"See {os.path.join(build_directory, 'cmake.out')}"
                    ' for output from cmake.\n'
                    "Skipping remaining steps for module set "
                    f"{module_set_name}, initial condition "
                    "{initial_condition_name}.",
                    file=sys.stderr
                )
                continue

            # Run the build.
            if verbose:
                print(
                    'Running make to build kaiju for module set'
                    f" {module_set_name}, initial condition"
                    f" {initial_condition_name}."
                )
            cmd = f"{module_cmd}; {make_cmd} >& make.out"
            if debug:
                print(f"cmd = {cmd}")
            try:
                _ = subprocess.run(cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(
                    f"ERROR: make for module set {module_set_name}, "
                    f"initial condition {initial_condition_name} failed.\n"
                    f"e.cmd = {e.cmd}\n"
                    f"e.returncode = {e.returncode}\n"
                    f"See {os.path.join(build_directory, 'make.out')}"
                    ' for output from make.\n'
                    "Skipping remaining steps for module set "
                    f"{module_set_name}, initial condition "
                    "{initial_condition_name}.",
                    file=sys.stderr
                )
                continue

            # Check for gamera.x.
            executable_list = ['gamera.x']
            missing = []
            for executable in executable_list:
                path = os.path.join(build_directory, BUILD_BIN_DIR, executable)
                if not os.path.isfile(path):
                    missing.append(executable)
            if len(missing) > 0:
                for executable in missing:
                    print(f"ERROR: Did not build {executable}.")
            else:
                test_passed[i_test][j_ic] = True

        # End loop over initial conditions
    # End loop over module sets

    # -------------------------------------------------------------------------

    # Summarize the test results
    test_summary_message = (
        'Results of initial condition build tests (`ICtest.py`):\n'
    )
    for (i_test, module_list_file) in enumerate(module_list_files):
        for (j_ic, initial_condition_path) in enumerate(initial_condition_paths):
            module_set_name = module_list_file.rstrip('.lst')
            initial_condition_name = os.path.basename(initial_condition_path)
            initial_condition_name = initial_condition_name.rstrip('.F90')
            test_summary_message += (
                f"Module set `{module_set_name}`, initial condition "
                f"`{initial_condition_name}`: "
            )
            if test_passed[i_test][j_ic]:
                test_summary_message += 'PASSED\n'
            else:
                test_summary_message += '*FAILED*\n'
    print(test_summary_message)

    # If loud mode is on, post report to Slack.
    if be_loud:
        message = 'Results of initial condition build tests (`ICtest.py`): '
        if 'FAILED' in test_summary_message:
            message += '*FAILED*\n'
        else:
            message += '*ALL PASSED*\n'
        message += 'Details in thread for this messsage.\n'
        slack_response = common.slack_send_message(
            slack_client, message, is_test=is_test
        )
        if slack_response['ok']:
            thread_ts = slack_response['ts']
            slack_response = common.slack_send_message(
                slack_client, test_summary_message, thread_ts=thread_ts,
                is_test=is_test
            )
        else:
            print('*ERROR* Unable to post test summary to Slack.')

    # -------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
