#!/usr/bin/env python

"""Run MAGE Fortran unit tests.

This script runs a series of unit tests of the MAGE Fortran software.

Authors
-------
Jeff Garretson (jeffrey.garretson@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
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

# Prefix for naming unit test directories
UNIT_TEST_DIRECTORY_PREFIX = 'unitTest_'

# glob pattern for naming unit test directories
UNIT_TEST_DIRECTORY_GLOB_PATTERN = 'unitTest_*'

# Home directory for pFUnit compiled code
PFUNIT_HOME = '/glade/u/home/ewinter/work/mage_testing/derecho/pfunit/pFUnit-4.2.0/ifort-23-mpich-derecho'

# Name of kaiju subdirectory for external code
KAIJU_EXTERNAL_DIRECTORY = 'external'

# Subdirectory of KAIJUHOME containing the test scripts
KAIJU_TEST_SCRIPTS_DIRECTORY = 'testingScripts'

# Subdirectory of KAIJU_TEST_SCRIPTS_DIRECTORY containing module lists
MODULE_LIST_DIRECTORY = 'mage_build_test_modules'

# Name of file containing names of modules lists to use for unit tests
UNIT_TEST_LIST_FILE = 'unit_test.lst'

# List of pFUnit directories to copy
PFUNIT_BINARY_DIRECTORIES = [
    'FARGPARSE-1.1',
    'GFTL-1.3',
    'GFTL_SHARED-1.2',
    'PFUNIT-4.2',
]

# PBS scripts for unit test jobs.
UNIT_TEST_PBS_SCRIPTS = [
    'genTestData.pbs',
    'runCaseTests.pbs',
    'runNonCaseTests1.pbs',
    # 'runNonCaseTests2.pbs',  # Hangs for 12 hours
]

# Subdirectory of kaiju containing unit test data files.
TESTS_DIRECTORY = 'tests'

# Subdirectory of kaiju/tests containing unit test data files.
UNIT_TEST_DATA_DIRECTORY = 'voltron_mpi'

# Data files for unit tests
UNIT_TEST_DATA_FILES = [
    'bcwind.h5',
    'lfmD.h5',
    'rcmconfig.h5',
    'geo_mpi.xml',
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

    # Clean up the results from previous tests.
    if verbose:
        print(f'Cleaning up from previous tests.')
    directories = glob.glob(UNIT_TEST_DIRECTORY_GLOB_PATTERN)
    for directory in directories:
        shutil.rmtree(directory)

    #--------------------------------------------------------------------------

    # Make a copy of the pFUnit code under kaiju/external. Delete
    # existing copies.
    for directory in PFUNIT_BINARY_DIRECTORIES:
        from_path = os.path.join(PFUNIT_HOME, directory)
        to_path = os.path.join(KAIJU_EXTERNAL_DIRECTORY, directory)
        try:
            shutil.rmtree(to_path)
        except:
            pass
        shutil.copytree(from_path, to_path)

    #--------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for unit tests.
    path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
                        MODULE_LIST_DIRECTORY, UNIT_TEST_LIST_FILE)
    with open(path, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [s.rstrip() for s in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # <HACK>
    # Just use first module set for now.
    module_list_files = [module_list_files[0]]
    # </HACK>

    #--------------------------------------------------------------------------

    # Run the unit tests with each set of modules.

    # Run unit tests with each set of modules.
    for module_list_file in module_list_files:
        if verbose:
            print('Performing Fortran unit tests with module set '
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

        # Make a directory for this test, and go there.
        dir_name = f"{UNIT_TEST_DIRECTORY_PREFIX}{module_list_name}"
        build_directory = os.path.join(kaiju_home, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Assemble the commands to load the listed modules.
        module_cmd = f"module --force purge; module load {' '.join(module_names)}"
        if debug:
            print(f"module_cmd = {module_cmd}")

        # <HACK>
        # Extra argument needed for unit test build.
        cmake_options += ' -DCMAKE_BUILD_TYPE=RELWITHDEBINFO'
        # </HACK>

        # Run cmake to build the Makefile.
        cmd = f"{module_cmd}; {cmake_environment} cmake {cmake_options} {kaiju_home}"
        if debug:
            print(f"cmd = {cmd}")
        # <HACK> To ignore cmake error on bcwind.h5 for now.
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except:
            pass
        # </HACK>

        # Run the build.
        cmd = f"{module_cmd}; make gamera_mpi voltron_mpi allTests"
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, check=True, text=True)

        # Copy in the test PBS scripts.
        for filename in UNIT_TEST_PBS_SCRIPTS:
            from_path = os.path.join(kaiju_home, TESTS_DIRECTORY, filename)
            to_path = os.path.join(BUILD_BIN_DIR, filename)
            shutil.copyfile(from_path, to_path)
        for filename in UNIT_TEST_DATA_FILES:
            from_path = os.path.join(kaiju_home, TESTS_DIRECTORY, UNIT_TEST_DATA_DIRECTORY, filename)
            to_path = os.path.join(BUILD_BIN_DIR, filename)
            shutil.copyfile(from_path, to_path)
            
        # Go to the bin directory for testing.
        os.chdir(BUILD_BIN_DIR)
    
        # Submit the jobs to create the test data and run the unit
        # tests.  Note that the unit test jobs will only run if the
        # data generation job completes successfully.
        job_ids = []
        for pbs_file in UNIT_TEST_PBS_SCRIPTS:
            cmd = (f"qsub -A {account} -v MODULE_LIST='{' '.join(module_names)}',"
                   f"KAIJUROOTDIR={kaiju_home}")
            if len(job_ids) > 0:
                cmd += f" -W depend=afterok:{job_ids[0]}"
            cmd += f" {pbs_file}"
            if debug:
                print(f"cmd = {cmd}")
            cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                                   capture_output=True)
            job_id = cproc.stdout.split('.')[0]
            if debug:
                print(f"job_id = {job_id}")
            job_ids.append(job_id)

        # Record the job IDs.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            for job_id in job_ids:
                f.write(f"{job_id}\n")

        # <HACK>
        message = f"Unit tests submitted in jobs {', '.join(job_ids)}"
        # </HACK>

        # If this is a test run, don't post to Slack. Otherwise, if
        # loud, send Slack message.
        if not is_test and be_loud:
            common.slack_send_message(slack_client, message)

        # Send message to stdout.
        print(message)

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
