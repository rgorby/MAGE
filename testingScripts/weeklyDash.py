#!/usr/bin/env python


"""Run the MAGE weekly dash tests.

This script runs the MAGE weekly dash tests.

Authors
-------
Jeff Garretson
Eric Winter
"""


# Import standard modules.
import datetime
import glob
import os
import platform
import shutil
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Run the MAGE weekly dash tests.'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Top-level directory for testing
KAIJU_TESTING_HOME = '/glade/work/ewinter/mage_testing/derecho'

# Glob pattern for weekly dash test direectories
WEEKLY_DASH_DIRECTORY_GLOB_PATTERN = 'weeklyDash_*'

# Name of working directory containing dash restart files.
WORKING_DASH_RESTART_DIRECTORY = 'dashRestarts'

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Name of file containing names of modules lists to use for weekly dash
WEEKLY_DASH_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'weekly_dash.lst')

# Prefix for weekly dash directory name
WEEKLY_DASH_DIRECTORY_PREFIX = 'weeklyDash_'

# Subdirectory of build directory containing compiled products to use in tests
BIN_DIR = 'bin'

# Source directory for weekly dash restart files
WEEKLY_DASH_RESTART_SRC_DIRECTORY = os.path.join(
    os.environ['MAGE_TEST_ROOT'], 'dashRestarts'
)

# List of weekly dash test files to copy
WEEKLY_DASH_TEST_FILES = [
    'weeklyDashGo.xml',
    'weeklyDashGo.pbs',
]

# Path to directory containing master-branch reference results.
REFERENCE_RESULTS_DIRECTORY_MASTER = os.path.join(
    KAIJU_TESTING_HOME, 'weekly_dash_files', 'reference_results',
    'master'
)

# Path to lfmQ.h5 for master-branch reference results.
REFERENCE_RESULTS_LFM_GRID_FILE_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'lfmQ.h5'
)

# Path to bcwind.h5 for master-branch reference results.
REFERENCE_RESULTS_BCWIND_FILE_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'bcwind.h5'
)

# Path to rcmconfig.h5 for master-branch reference results.
REFERENCE_RESULTS_RCMCONFIG_FILE_MASTER = os.path.join(
    REFERENCE_RESULTS_DIRECTORY_MASTER, 'rcmconfig.h5'
)


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
    account = args.account
    debug = args.debug
    be_loud = args.loud
    is_test = args.test
    verbose = args.verbose

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}"
              f" on {platform.node()}")
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

    # Clean up from previous builds.
    if verbose:
        print('Cleaning up from previous tests.')
    directories = glob.glob(WEEKLY_DASH_DIRECTORY_GLOB_PATTERN)
    directories.append(WORKING_DASH_RESTART_DIRECTORY)
    for directory in directories:
        if debug:
            print(f"Trying to remove {directory}.")
        try:
            shutil.rmtree(directory)
        except FileNotFoundError:
            pass  # directory may not exist.

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

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    with open(WEEKLY_DASH_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # -------------------------------------------------------------------------

    # Create the make command to build the code.
    make_cmd = 'make voltron_mpi.x'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Run the tests with each set of modules.
    for module_list_file in module_list_files:
        if verbose:
            print('Performing weekly dash with module set '
                  f"{module_list_file}")

        # Extract the name of the list.
        module_set_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_set_name = {module_set_name}.")

        # Read this module list file, extracting cmake environment and
        # options, if any.
        path = os.path.join(MODULE_LIST_DIRECTORY, module_list_file)
        if verbose:
            print(f"Reading {path}.")
        module_names, cmake_environment, cmake_options = (
            common.read_build_module_list_file(path)
        )
        if debug:
            print(f"module_names = {module_names}")
            print(f"cmake_environment = {cmake_environment}")
            print(f"cmake_options = {cmake_options}")

        # Add the cmake option for the weekly dash build.
        cmake_options += ' -DCMAKE_BUILD_TYPE=Release'
        if debug:
            print(f"cmake_options = {cmake_options}")

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge"
            f"; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Make a directory for this build, and go there.
        dir_name = f"{WEEKLY_DASH_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(KAIJUHOME, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Run cmake to build the Makefile.
        if verbose:
            print(
                'Running cmake to create Makefile for module set'
                f" {module_set_name}."
            )
        cmd = (
            f"{module_cmd}"
            f"; {cmake_environment} cmake {cmake_options}"
            f" {KAIJUHOME}"
            '>& cmake.out'
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                f"ERROR: cmake for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'cmake.out')}"
                ' for output from cmake.\n'
                f"Skipping remaining steps for module set {module_set_name}\n"
            )
            continue

        # Run the build.
        if verbose:
            print(
                'Running make to build kaiju for module set'
                f" {module_set_name}."
            )
        cmd = f"{module_cmd}; {make_cmd} >& make.out"
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                f"ERROR: make for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'make.out')}"
                ' for output from make.\n'
                f"Skipping remaining steps for module set {module_set_name}\n"
            )
            continue

        # Move into the bin directory to run the tests.
        os.chdir(BIN_DIR)

        # Generate the LFM grid file.
        if verbose:
            print('Creating LFM grid file.')
        cmd = 'genLFM.py -gid Q'
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                'ERROR: Unable to create LFM grid file for module set '
                f"{module_set_name}.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See testing log for output from genLFM.py.\n'
                f"Skipping remaining steps for module set {module_set_name}\n"
            )
            continue

        # Compare the new LFM grid file to the reference version.
        if verbose:
            print('Comparing new LFM grid file to reference version.')
        cmd = f"h5diff {REFERENCE_RESULTS_LFM_GRID_FILE_MASTER} lfmQ.h5"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                'ERROR: Error comparing LFM grid file for module set '
                f"{module_set_name} to reference version.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See testing log for output from h5diff.\n'
                f"Skipping remaining steps for module set {module_set_name}\n"
            )
            continue
        grid_diff = cproc.stdout.rstrip()
        if debug:
            print(f"grid_diff = {grid_diff}")
        if grid_diff != '':
            test_summary_message += ('LFM grid for weekly dash has changed on branch'
                  f" {git_branch_name}. Case cannot be run. Please"
                  ' re-generate restart data, and ensure the grid change'
                  ' was intentional.\n')

        # Generate the solar wind boundary condition file.
        if verbose:
            print('Creating solar wind initial conditions file.')
        cmd = (
            'cda2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00')
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                'ERROR: Unable to create solar wind boundary conditions file'
                f" for module set {module_set_name}.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See testing log for output from cda2wind.py.\n'
                f"Skipping remaining steps for module set {module_set_name}\n"
            )
            continue

        # Compare the new solar wind boundary condition file to the original.
        if verbose:
            print('Comparing new solar wind file to reference version.')
        cmd = f"h5diff {REFERENCE_RESULTS_BCWIND_FILE_MASTER} bcwind.h5"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                'ERROR: Error comparing solar wind file for module set '
                f"{module_set_name} to reference version.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See testing log for output from h5diff.\n'
                f"Skipping remaining steps for module set {module_set_name}\n"
            )
            continue
        wind_diff = cproc.stdout.rstrip()
        if debug:
            print(f"wind_diff = {wind_diff}")
        if wind_diff != '':
            test_summary_message += ('Solar wind data for weekly dash has changed on branch'
                  f" {git_branch_name}. Case cannot be run. Please"
                  ' re-generate restart data, and ensure the change'
                  ' was intentional.\n')
            continue

        # Generate the RCM configuration file.
        if verbose:
            print('Creating RCM configuration file.')
        cmd = 'genRCM.py'
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                'ERROR: Unable to create RCM configuration file'
                f" for module set {module_set_name}.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See testing log for output from genRCM.py.\n'
                f"Skipping remaining steps for module set {module_set_name}\n"
            )
            continue

        # Compare the new RCM configuration file to the original.
        if verbose:
            print('Comparing new RCM configuration file to reference version.')
        cmd = f"h5diff {REFERENCE_RESULTS_RCMCONFIG_FILE_MASTER} rcmconfig.h5"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                'ERROR (expected): Error comparing RCM configuration file for'
                f" module set {module_set_name} to reference version.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See testing log for output from h5diff.\n'
            )
            # continue  ERROR EXPECTED FOR NOW
        rcm_diff = cproc.stdout.rstrip()
        if debug:
            print(f"rcm_diff = {rcm_diff}")
        if rcm_diff != '':
            test_summary_message += ('RCM configuration for weekly dash has changed on branch'
                  f" {git_branch_name}. Case cannot be run. Please"
                  ' re-generate restart data, and ensure the change'
                  ' was intentional.\n')
            # continue  ERROR EXPECTED FOR NOW

        # Copy files needed for the weekly dash job.
        if verbose:
            print('Copying files needed for weekly dash.')
        for filename in WEEKLY_DASH_TEST_FILES:
            from_file = os.path.join(TEST_SCRIPTS_DIRECTORY, filename)
            to_file = os.path.join('.', filename)
            shutil.copyfile(from_file, to_file)

        # Create and submit the qsub command for a model run with these files.
        if verbose:
            print('Preparing to submit weekly dash model run.')
        cmd = (
            f"qsub -A {account} "
            f"-v MODULE_LIST='{' '.join(module_names)}'"
            f",KAIJUROOTDIR={KAIJUHOME} weeklyDashGo.pbs"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True,
                                   text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            test_summary_message += (
                'ERROR: Unable to submit job request for module set '
                f"{module_set_name}.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                'See testing log for output.\n'
                f"Skipping remaining steps for module set {module_set_name}\n"
            )
            continue
        firstJobNumber = cproc.stdout.split('.')[0]
        if debug:
            print(f"firstJobNumber = {firstJobNumber}")

        # Save the job number in a file.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            f.write(f"{firstJobNumber}\n")

        # Summarize the test results
        test_summary_message = (
            f"Weekly dash run with module set {module_set_name} "
            f"(`weeklyDash.py`) submitted as job {firstJobNumber}.\n"
        )
        print(test_summary_message)

        # If loud mode is on, post report to Slack.
        if be_loud:
            message = 'Results of weekly dash submission (`weeklyDash.py`): '
            if 'ERROR' in test_summary_message:
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
    main()
