#!/usr/bin/env python

"""Run MAGE build regression tests.

This script runs a series of builds of the MAGE software using sets of
modules listed in files under:

$KAIJUHOME/testingScripts/mage_build_modules

This script reads the file build_test.lst from this directory, and
uses the contents as a list of module list files to use for MAGE build
tests.

NOTE: These tests are performed on a load-balance-assigned login node on
derecho. No PBS job is submitted.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import glob
import os
import shutil
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE build testing'

# Prefix for naming build test directories
BUILD_TEST_DIRECTORY_PREFIX = 'build_'

# glob pattern for naming build test directories
BUILD_TEST_DIRECTORY_GLOB_PATTERN = 'build_*'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Path to directory to use for building executable list
EXECUTABLE_LIST_BUILD_DIRECTORY = os.path.join(KAIJUHOME, 'testFolder')

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Path to module list file to use when generating the list of executables
EXECUTABLE_LIST_MODULE_LIST = os.path.join(MODULE_LIST_DIRECTORY, '01.lst')

# Path to file containing list of module sets to use for build tests
BUILD_TEST_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'build_test.lst')

# Name of subdirectory of current build subdirectory containing binaries
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
    directories = glob.glob(BUILD_TEST_DIRECTORY_GLOB_PATTERN)
    directories.append(EXECUTABLE_LIST_BUILD_DIRECTORY)
    for directory in directories:
        try:
            shutil.rmtree(directory)
        except FileNotFoundError:
            pass  # EXECUTABLE_LIST_BUILD_DIRECTORY may not exist.

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

    # Do a preliminary cmake run to generate the list of executables.

    if verbose:
        print(f"Generating list of executables for branch {git_branch_name}.")

    # Create and move to the preliminary build folder.
    os.mkdir(EXECUTABLE_LIST_BUILD_DIRECTORY)
    os.chdir(EXECUTABLE_LIST_BUILD_DIRECTORY)

    # Read the module list file for building the executable list,
    if verbose:
        print(f"Reading {EXECUTABLE_LIST_MODULE_LIST} to determine"
              ' modules and build options for generating list of executables'
              ' to build.')
    module_names, cmake_environment, cmake_options = (
        common.read_build_module_list_file(EXECUTABLE_LIST_MODULE_LIST)
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

    # Run cmake to build the Makefile.
    if verbose:
        print('Running cmake for executable list generation.')
    cmd = (
        f"{module_cmd}"
        f"; {cmake_environment} cmake {cmake_options}"
        f" {KAIJUHOME}")
    if debug:
        print(f"cmd = {cmd}")
    try:
        _ = subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(
            'cmake for building executable list failed.\n'
            f"e.cmd = {e.cmd}"
            f"e.returncode = {e.returncode}\n"
            'See test log for output.\n'
            'Unable to generate executable list.',
            file=sys.stderr
        )
        raise

    # Run make to build the list of executable targets.
    if verbose:
        print('Running make for executable list generation.')
    cmd = f"{module_cmd}; make help | grep '\.x'"
    if debug:
        print(f"cmd = {cmd}")
    try:
        cproc = subprocess.run(cmd, shell=True, check=True,
                               text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(
            'make for building executable list failed.\n'
            f"e.cmd = {e.cmd}\n"
            f"e.returncode = {e.returncode}\n"
            f"e.stdout = {e.stdout}\n"
            f"e.stderr = {e.stderr}\n"
            'Unable to generate executable list.',
            file=sys.stderr
        )
        raise
    executable_list = cproc.stdout.splitlines()

    # Remove the first four characters (dots and spaces).
    executable_list = [_[4:] for _ in executable_list]
    if debug:
        print(f"executable_list = {executable_list}")

    # Create the make command to build the executable list.
    make_cmd = f"make {' '.join(executable_list)}"
    if debug:
        print(f"make_cmd = {make_cmd}")

    # -------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    with open(BUILD_TEST_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # -------------------------------------------------------------------------

    # Initalize test results for all module sets to False (failed).
    test_passed = [False]*len(module_list_files)

    # Do a build with each set of modules.
    for (i_test, module_list_file) in enumerate(module_list_files):
        if verbose:
            print('Performing build test with module_list_file '
                  f"{module_list_file}")

        # Extract the name of the list.
        module_list_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_list_name = {module_list_name}.")

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

        # Make a directory for this build, and go there.
        dir_name = f"{BUILD_TEST_DIRECTORY_PREFIX}{module_list_name}"
        build_directory = os.path.join(KAIJUHOME, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Run cmake to build the Makefile.
        if verbose:
            print(
                'Running cmake to create Makefile for module set'
                f" {module_list_name}."
            )
        cmd = (
            f"{module_cmd}"
            f"; {cmake_environment} cmake {cmake_options}"
            f" {KAIJUHOME}")
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"cmake for module set {module_list_name} failed.\n"
                f"e.cmd = {e.cmd}"
                f"e.returncode = {e.returncode}\n"
                'See test log for output.\n'
                f"Skipping remaining steps for module set {module_list_name}",
                file=sys.stderr
            )
            continue

        # Run the build.
        if verbose:
            print(
                'Running make to build kaiju for module set'
                f" {module_list_name}."
            )
        cmd = f"{module_cmd}; {make_cmd}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"make for module set {module_list_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See test log for output.\n'
                f"Skipping remaining steps for module set {module_list_name}",
                file=sys.stderr
            )
            continue

        # Check for all executables
        if verbose:
            print(f"Checking build results for module set {module_list_name}.")
        missing = []
        for executable in executable_list:
            path = os.path.join(build_directory, BUILD_BIN_DIR, executable)
            if not os.path.isfile(path):
                missing.append(executable)
        if len(missing) > 0:
            for executable in missing:
                print(f"Did not build {executable}.")
        else:
            test_passed[i_test] = True

    # End of loop over module sets.

    # -------------------------------------------------------------------------

    # Summarize the test results
    test_summary_message = 'Results of build tests:\n'
    for (i_test, module_list_file) in enumerate(module_list_files):
        module_list_name = module_list_file.rstrip('.lst')
        test_summary_message += f"Module set {module_list_name}: "
        if test_passed[i_test]:
            test_summary_message += 'PASSED\n'
        else:
            test_summary_message += '*FAILED*\n'

    # If loud mode is on, post report to Slack.
    if be_loud:
        common.slack_send_message(slack_client, test_summary_message,
                                  is_test=is_test)

    # Send message to stdout.
    print(test_summary_message)

    # -------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
