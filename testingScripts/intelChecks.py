#!/usr/bin/env python

"""Run MAGE Intel Inspector tests.

This script runs a series of tests of the MAGE software using Intel tools.

Authors
-------
Jeff Garretson
Eric Winter

"""


# Import standard modules.
import datetime
import os
import sys
import shutil
import subprocess

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE checks with Intel Inspector tools'

# Home directory of kaiju installation
KAIJUHOME = os.environ['KAIJUHOME']

# Root of directory tree for this set of tests.
MAGE_TEST_SET_ROOT = os.environ['MAGE_TEST_SET_ROOT']

# Directory for Intel Inspector checks
INTEL_CHECKS_DIRECTORY = os.path.join(MAGE_TEST_SET_ROOT, 'intelChecks')

# Path to directory containing the test scripts
TEST_SCRIPTS_DIRECTORY = os.path.join(KAIJUHOME, 'testingScripts')

# Path to directory containing module lists
MODULE_LIST_DIRECTORY = os.path.join(TEST_SCRIPTS_DIRECTORY,
                                     'mage_build_test_modules')

# Name of file containing names of modules lists to use for Intel checks
INTEL_CHECKS_LIST_FILE = os.path.join(MODULE_LIST_DIRECTORY, 'intelChecks.lst')

# Prefix for naming Intel Inspector checks directories
INTEL_CHECKS_DIRECTORY_PREFIX = 'intelChecks_'

# Name of build subdirectory containing binaries
BUILD_BIN_DIR = 'bin'

# Name of PBS account to use for testing jobs.
DERECHO_TESTING_ACCOUNT = os.environ['DERECHO_TESTING_ACCOUNT']

# Token string for access to Slack.
SLACK_BOT_TOKEN = os.environ['SLACK_BOT_TOKEN']

# Branch or commit (or tag) used for testing.
BRANCH_OR_COMMIT = os.environ['BRANCH_OR_COMMIT']


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
    # account = args.account
    debug = args.debug
    be_loud = args.loud
    is_test = args.test
    verbose = args.verbose

    # ------------------------------------------------------------------------

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")
        print(f"Current directory is {os.getcwd()}")

    # ------------------------------------------------------------------------

    # Make a directory to hold all of the Intel Inspector tests.
    if verbose:
        print(f"Creating {INTEL_CHECKS_DIRECTORY}.")
    os.mkdir(INTEL_CHECKS_DIRECTORY)

    # ------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of module sets to use for Intel checks.
    with open(INTEL_CHECKS_LIST_FILE, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [_.rstrip() for _ in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # ------------------------------------------------------------------------

    # Run the Intel Inspector checks with each set of modules.

    # Create the common make command for all module sets.
    make_cmd = 'make gamera_mpi voltron_mpi'
    if debug:
        print(f"make_cmd = {make_cmd}")

    # Run Intel checks with each set of modules.
    for (i_set, module_list_file) in enumerate(module_list_files):
        if verbose:
            print('Performing Intel Inspector checks with module set '
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
            f"module --force purge; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Add the additional flags needed for Intel Inspector checks.
        cmake_options += ' -DDISABLE_DEBUG_BOUNDS_CHECKS=ON'
        cmake_options += ' -DCMAKE_BUILD_TYPE=DEBUG'
        if debug:
            print(f"cmake_options = {cmake_options}")

        # Make a directory for this test, and go there.
        dir_name = f"{INTEL_CHECKS_DIRECTORY_PREFIX}{module_set_name}"
        build_directory = os.path.join(INTEL_CHECKS_DIRECTORY, dir_name)
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
            # NOTE: stdout and stderr goes cmake.out.
            _ = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"ERROR: cmake for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'cmake.out')}"
                ' for output from cmake.\n'
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
            # NOTE: stdout and stderr go into make.out.
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(
                f"ERROR: make for module set {module_set_name} failed.\n"
                f"e.cmd = {e.cmd}\n"
                f"e.returncode = {e.returncode}\n"
                f"See {os.path.join(build_directory, 'make.out')}"
                ' for output from make.\n'
                f"Skipping remaining steps for module set {module_set_name}",
                file=sys.stderr
            )
            continue

        # Go to the bin directory for testing.
        os.chdir(BUILD_BIN_DIR)

        # Copy in the test PBS scripts and files.
        test_files = [
            'tinyCase.xml',
            'intelCheckSubmitMem.pbs',
            'intelCheckSubmitThread.pbs',
            'intelCheckReportSubmit.pbs',
            'bcwind.h5',
            'lfmD.h5',
            'rcmconfig.h5',
            'memSuppress.sup',
            'threadSuppress.sup',
        ]
        cwd = os.getcwd()
        for filename in test_files:
            from_path = os.path.join(TEST_SCRIPTS_DIRECTORY, filename)
            to_path = os.path.join(cwd, filename)
            shutil.copyfile(from_path, to_path)

        # Run each check in its own PBS job.
        pbs_files = [
            'intelCheckSubmitMem.pbs',
            'intelCheckSubmitThread.pbs',
            'intelCheckReportSubmit.pbs',
        ]
        job_ids = []
        for pbs_file in pbs_files:
            cmd = (
                f"qsub -A {DERECHO_TESTING_ACCOUNT} "
                f"-v MODULE_LIST='{' '.join(module_names)}',"
                f"KAIJUROOTDIR={KAIJUHOME},"
                f"MAGE_TEST_SET_ROOT={MAGE_TEST_SET_ROOT},"
                f"DERECHO_TESTING_ACCOUNT={DERECHO_TESTING_ACCOUNT},"
                f"SLACK_BOT_TOKEN={SLACK_BOT_TOKEN}"
            )
            if pbs_file == 'intelCheckReportSubmit.pbs':
                cmd += f" -W depend=after:{':'.join(job_ids)}"
            cmd += f" {pbs_file}"
            if debug:
                print(f"cmd = {cmd}")
            try:
                cproc = subprocess.run(cmd, shell=True, check=True,
                                       text=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                print('qsub failed.')
                print(f"e.cmd = {e.cmd}")
                print(f"e.returncode = {e.returncode}")
                print('See test log for output.')
                print('Skipping remaining steps for module set'
                      f" {module_list_file}.\n")
                continue
            job_id = cproc.stdout.split('.')[0]
            if debug:
                print(f"job_id = {job_id}")
            job_ids.append(job_id)

        # Record the job IDs.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            for job_id in job_ids:
                f.write(f"{job_id}\n")

    # -----------------------------------------------------------------------

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # ------------------------------------------------------------------------

    # Detail the test results
    test_report_details_string = ''
    for (pbs_file, job_id) in zip(pbs_files, job_ids):
        test_report_details_string += f"{pbs_file} submitted as job {job_id}."

    # Summarize the test results
    test_report_summary_string = (
        'Intel Inspector tests submitted (`intelChecks.py`).'
    )

    # Print the test results summary and details.
    print(test_report_summary_string)
    print(test_report_details_string)

    # If loud mode is on, post report to Slack.
    if be_loud:
        test_report_summary_string += 'Details in thread for this messsage.\n'
        slack_response_summary = common.slack_send_message(
            slack_client, test_report_summary_string, is_test=is_test
        )
        if slack_response_summary['ok']:
            thread_ts = slack_response_summary['ts']
            slack_response_details = common.slack_send_message(
                slack_client, test_report_details_string, thread_ts=thread_ts,
                is_test=is_test
            )
            if 'ok' not in slack_response_details:
                print('*ERROR* Unable to post test details to Slack.')
        else:
            print('*ERROR* Unable to post test summary to Slack.')


    # ------------------------------------------------------------------------

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    main()
