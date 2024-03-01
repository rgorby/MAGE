#!/usr/bin/env python

"""Run MAGE initial condition build tests.

This script runs a series of initial condition build tests of the MAGE
software.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import argparse
import datetime
import glob
import os
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE initial condition build testing'

# Name of top-level directory for initial condition build tests
IC_BUILD_TEST_DIRECTORY = 'ICBuilds'

# Subdirectory of KAIJUHOME containing the test scripts
KAIJU_TEST_SCRIPTS_DIRECTORY = 'testingScripts'

# Subdirectory of KAIJU_TEST_SCRIPTS_DIRECTORY containing module lists
MODULE_LIST_DIRECTORY = 'mage_build_test_modules'

# Name of file containing names of modules lists to use for initial
# condition build tests
IC_BUILD_TEST_LIST_FILE = 'initial_condition_build_test.lst'

# Path to directory containing initial condition source code
IC_SRC_DIRECTORY = os.path.join('src', 'gamera', 'ICs')


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

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")

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

    # Clean up from previous tests.
    if verbose:
        print('Cleaning up from previous tests.')
    # shutil.rmtree() does not remove non-empty directories.
    cmd = f"rm -rf {IC_BUILD_TEST_DIRECTORY}"
    cproc = subprocess.run(cmd, shell=True, check=True, text=True)

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
        path = os.path.join(kaiju_home, 'external', directory)
        try:
            shutil.rmtree(path)
        except:
            pass
    # </HACK>

    #--------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
                        MODULE_LIST_DIRECTORY, IC_BUILD_TEST_LIST_FILE)
    with open(path, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [s.rstrip() for s in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    #--------------------------------------------------------------------------
    
    # Get a list of initial conditions to try, ignoring files in the
    # "deprecated" folder. GAMERA ONLY FOR NOW
    initial_condition_directory = os.path.join(kaiju_home, IC_SRC_DIRECTORY)
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

    # Build using each set of modules and each initial condition.

    # Create the root directory for the initial condition tests.
    initial_condition_test_root = os.path.join(kaiju_home, IC_BUILD_TEST_DIRECTORY)
    if debug:
        print(f"initial_condition_test_root = {initial_condition_test_root}")
    os.mkdir(initial_condition_test_root)

    # Run initial condition tests with each set of modules.
    for module_list_file in module_list_files:
        if verbose:
            print('Performing initial condition build tests with module set '
                  f"{module_list_file}.")

        # Extract the name of the list.
        module_list_name = module_list_file.replace('.lst', '')
        if debug:
            print(f"module_list_name = {module_list_name}")

        # Read this module list file, extracting cmake environment and
        # options, if any.
        path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
                            MODULE_LIST_DIRECTORY, module_list_file)
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
            f"module --force purge; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Build with each initial condition.
        for initial_condition_path in initial_condition_paths:

            # Extract the initial condition name.
            initial_condition_name = os.path.basename(initial_condition_path)
            if debug:
                print(f"initial_condition_name = {initial_condition_name}")
            if verbose:
                print(f"Building with module set {module_list_name} and initial condition {initial_condition_name}.")

            # Make a directory for this test, and go there.
            build_directory = os.path.join(
                kaiju_home, IC_BUILD_TEST_DIRECTORY,
                f"gamera_{initial_condition_name}_{module_list_name}"
            )
            if debug:
                print(f"build_directory = {build_directory}")
            os.mkdir(build_directory)
            os.chdir(build_directory)

            # <HACK>
            # Add cmake options for initial condition test builds.
            cmake_options += f" -DGAMIC:FILEPATH={initial_condition_path}"
            if debug:
                print(f"cmake_options = {cmake_options}")
            # </HACK>

            # Run cmake to build the Makefile.
            cmd = f"{module_cmd}; {cmake_environment} cmake {cmake_options} {kaiju_home}"
            if debug:
                print(f"cmd = {cmd}")
            # <HACK>
            # Ignore cmake error on bcwind,h5.
            try:
                cproc = subprocess.run(cmd, shell=True, check=True, text=True)
            except:
                pass
            # </HACK>

            # Run the build.
            cmd = f"{module_cmd}; make gamera.x"
            if debug:
                print(f"cmd = {cmd}")
            cproc = subprocess.run(cmd, shell=True, check=True, text=True)
            make_return_code = cproc.returncode

            # <HACK>
            message = (
                f"Build using module set {module_list_name} for "
                f"initial conditions {initial_condition_name} "
                f"returned {make_return_code}."
            )
            # </HACK>

            # If this is a test run, don't post to Slack. Otherwise,
            # if loud, send Slack message.
            if not is_test and be_loud:
                common.slack_send_message(slack_client, message)

            # Send message to stdout.
            print(message)

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
