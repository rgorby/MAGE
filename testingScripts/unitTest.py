#!/usr/bin/env python

"""Run MAGE Fortran unit tests.

This script runs a series of unit tests of the MAGE Fortran software. These
tests are run as PBS jobs on derecho. There will be one job which generates
the data for testing, then 1 or more dependent jobs that use the newly-
generated data for unit testing.

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
DESCRIPTION = 'Script for MAGE Fortran unit testing'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Prefix for naming unit test directories
UNIT_TEST_DIRECTORY_PREFIX = 'unitTest_'

# glob pattern for naming unit test directories
UNIT_TEST_DIRECTORY_GLOB_PATTERN = 'unitTest_*'

# Top-level direectory for testing on derecho.
DERECHO_TESTING_HOME = '/glade/u/home/ewinter/work/mage_testing/derecho'

# Home directory for pFUnit compiled code
PFUNIT_HOME = os.path.join(
    DERECHO_TESTING_HOME, 'pfunit', 'pFUnit-4.2.0', 'ifort-23-mpich-derecho'
)

# List of pFUnit directories to copy
PFUNIT_BINARY_DIRECTORIES = [
    os.path.join(PFUNIT_HOME, 'FARGPARSE-1.1'),
    os.path.join(PFUNIT_HOME, 'GFTL-1.3'),
    os.path.join(PFUNIT_HOME, 'GFTL_SHARED-1.2'),
    os.path.join(PFUNIT_HOME, 'PFUNIT-4.2'),
]

# Path to kaiju subdirectory for external code
KAIJU_EXTERNAL_DIRECTORY = os.path.join(KAIJUHOME, 'external')

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Name of file containing names of modules lists to use for unit tests
UNIT_TEST_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'unit_test.lst')

# Path to directory containing unit test files.
TESTS_DIRECTORY = os.path.join(KAIJUHOME, 'tests')

# PBS scripts for unit test jobs.
UNIT_TEST_PBS_SCRIPTS = [
    'genTestData.pbs',
    'runCaseTests.pbs',
    'runNonCaseTests1.pbs',
    # 'runNonCaseTests2.pbs',  # Hangs for 12 hours
]

# Input files for unit tests
UNIT_TEST_DATA_INPUT_DIRECTORY = os.path.join(os.environ['MAGE_TEST_ROOT'],
                                              'unit_test_inputs')
UNIT_TEST_DATA_INPUT_FILES = [
    'bcwind.h5',
    'geo_mpi.xml',
    'lfmD.h5',
    'rcmconfig.h5',
]

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
    subprocess.CalledProcessError
        If an exception occurs in subprocess.run()
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

    # Clean up the results from previous tests.
    if verbose:
        print('Cleaning up from previous tests.')
    directories = glob.glob(UNIT_TEST_DIRECTORY_GLOB_PATTERN)
    for directory in directories:
        shutil.rmtree(directory)

    # -------------------------------------------------------------------------

    # Make a copy of the pFUnit code under kaiju/external. Delete
    # existing copies.
    if verbose:
        print('Copying compiled pFUnit binaries.')
    for directory in PFUNIT_BINARY_DIRECTORIES:
        from_path = directory
        dir_name = os.path.split(from_path)[-1]
        to_path = os.path.join(KAIJU_EXTERNAL_DIRECTORY, dir_name)
        try:
            shutil.rmtree(to_path)
        except FileNotFoundError:
            pass  # Might not exist.
        shutil.copytree(from_path, to_path)

    # -------------------------------------------------------------------------

    # Find the current branch.
    git_branch_name = common.git_get_branch_name()
    if debug:
        print(f"git_branch_name = {git_branch_name}")

    # -------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    with open(UNIT_TEST_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # Compute the list of module set names.
    module_set_names = [_.rstrip('.lst') for _ in module_list_files]
    if debug:
        print(f"module_set_names = {module_set_names}")

    # -------------------------------------------------------------------------

    # Create the common make command for all module sets.
    make_cmd = 'make gamera_mpi voltron_mpi allTests'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Initalize job ID to None for all module set/PBS script combinations.
    job_ids = []
    for _ in module_list_files:
        job_ids.append([None]*len(UNIT_TEST_PBS_SCRIPTS))

    # Run the tests with each set of modules.
    for (i_set, module_list_file) in enumerate(module_list_files):
        module_set_name = module_set_names[i_set]
        if verbose:
            print('Performing Fortran unit tests with module set '
                  f"{module_set_name}.")

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

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge"
            f"; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # <HACK>
        # Extra argument needed for unit test build.
        cmake_options += ' -DCMAKE_BUILD_TYPE=RELWITHDEBINFO'
        # </HACK>

        # Make a directory for this test, and go there.
        dir_name = f"{UNIT_TEST_DIRECTORY_PREFIX}{module_set_name}"
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
            f"{module_cmd}; {cmake_environment} cmake {cmake_options}"
            f" {KAIJUHOME} >& cmake.out"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            path = os.path.join(build_directory, 'cmake.out')
            print(
                f"ERROR: cmake for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {path} for output from cmake.\n"
                f"Skipping remaining steps for module set {module_set_name}",
                file=sys.stderr
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
            path = os.path.join(build_directory, 'make.out')
            print(
                f"ERROR: make for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {path} for output from make.\n"
                f"Skipping remaining steps for module set {module_set_name}",
                file=sys.stderr
            )
            continue

        # Copy in the PBS scripts for unit testing.
        for filename in UNIT_TEST_PBS_SCRIPTS:
            from_path = os.path.join(TESTS_DIRECTORY, filename)
            to_path = os.path.join(BUILD_BIN_DIR, filename)
            shutil.copyfile(from_path, to_path)

        # Copy in inputs for unit test data generation.
        for filename in UNIT_TEST_DATA_INPUT_FILES:
            from_path = os.path.join(UNIT_TEST_DATA_INPUT_DIRECTORY, filename)
            to_path = os.path.join(BUILD_BIN_DIR, filename)
            shutil.copyfile(from_path, to_path)

        # Go to the bin directory for testing.
        os.chdir(BUILD_BIN_DIR)

        # Submit the jobs to create the test data and run the unit
        # tests. Note that the unit test jobs will only run if the
        # data generation job completes successfully.
        for (j_pbs, pbs_file) in enumerate(UNIT_TEST_PBS_SCRIPTS):
            job_id = None
            cmd = (
                f"qsub -A {account} -v MODULE_LIST='{' '.join(module_names)}',"
                f"KAIJUROOTDIR={KAIJUHOME}"
            )
            # <HACK>
            # Assumes data generation job is first.
            if j_pbs > 0:
                cmd += f" -W depend=afterok:{job_ids[i_set][0]}"
            # </HACK>
            cmd += f" {pbs_file}"
            if debug:
                print(f"cmd = {cmd}")
            try:
                cproc = subprocess.run(cmd, shell=True, check=True,
                                       text=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                path = os.path.join(build_directory, BUILD_BIN_DIR,
                                    f"qsub_{j_pbs}.out")
                with open(path, 'w', encoding='utf-8') as f:
                    f.write(e.stdout)
                print(
                    'ERROR: Job submission failed.\n'
                    f"e.cmd = {e.cmd}\n"
                    f"e.returncode = {e.returncode}\n"
                    f"See {path} for output from qsub.\n"
                    'Skipping remaining steps for module set '
                    f"{module_set_name}.",
                    file=sys.stderr
                )
                continue

            # Save the job ID.
            job_id = cproc.stdout.split('.')[0]
            if debug:
                print(f"job_id = {job_id}")
            job_ids[i_set][j_pbs] = job_id
        # End of loop over PBS scripts

        # Record the job IDs in a text file.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            for job_id in job_ids[i_set]:
                f.write(f"{job_id}\n")

    # End of loop over module sets

    # -------------------------------------------------------------------------

    # Summarize the test results
    test_summary_message = (
        'Results of Fortran unit test submission `unitTest.py`):\n'
    )
    for (i_set, module_set_name) in enumerate(module_set_names):
        for (j_pbs, pbs_file) in enumerate(UNIT_TEST_PBS_SCRIPTS):
            test_summary_message += (
                f"Module set `{module_set_name}`, submit `{pbs_file}`: "
            )
            if job_ids[i_set][j_pbs] is not None:
                test_summary_message += f"{job_ids[i_set][j_pbs]}\n"
            else:
                test_summary_message += '*FAILED*\n'
    print(test_summary_message)

    # If loud mode is on, post report to Slack.
    if be_loud:
        common.slack_send_message(slack_client, test_summary_message,
                                  is_test=is_test)

    # -------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
