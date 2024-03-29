#!/usr/bin/env python

"""Run MAGE Intel tests.

This script runs a series of tests of the MAGE software using Intel tools.

Authors
-------
Jeff Garretson
Eric Winter
"""


# Import standard modules.
import datetime
import glob
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

# Prefix for naming Intel Inspector checks directories
INTEL_CHECKS_DIRECTORY_PREFIX = 'intelChecks_'

# glob pattern for naming Intel Inspector checks directories
INTEL_CHECKS_DIRECTORY_GLOB_PATTERN = 'intelChecks_*'

# Subdirectory of KAIJUHOME containing the test scripts
KAIJU_TEST_SCRIPTS_DIRECTORY = 'testingScripts'

# Subdirectory of KAIJU_TEST_SCRIPTS_DIRECTORY containing module lists
MODULE_LIST_DIRECTORY = 'mage_build_test_modules'

# Name of file containing names of modules lists to use for Intel Inspector
# tests
INTEL_CHECKS_LIST_FILE = 'intel_checks.lst'


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
    # doAll = args.all
    debug = args.debug
    # forceRun = args.force
    be_loud = args.loud
    is_test = args.test
    verbose = args.verbose

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")

    # -------------------------------------------------------------------------

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    # -------------------------------------------------------------------------

    # Move to the MAGE installation directory.
    kaiju_home = os.environ['KAIJUHOME']
    os.chdir(kaiju_home)

    # -------------------------------------------------------------------------

    # Clean up from previous tests.
    if verbose:
        print('Cleaning up from previous Intel Inspector checks.')
    os.chdir(kaiju_home)
    directories = glob.glob(INTEL_CHECKS_DIRECTORY_GLOB_PATTERN)
    for directory in directories:
        shutil.rmtree(directory)

    # -------------------------------------------------------------------------

    # Find the current branch.
    git_branch_name = common.git_get_branch_name()
    if debug:
        print(f"git_branch_name = {git_branch_name}")

    # -------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
                        MODULE_LIST_DIRECTORY, INTEL_CHECKS_LIST_FILE)
    with open(path, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [s.rstrip() for s in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    # <HACK>
    # Just use first module set for now.
    module_list_files = [module_list_files[0]]
    # </HACK>

    # -------------------------------------------------------------------------

    # Run the Intel Inspector checks with each set of modules.

    # Initialize the report string.
    message = f"Running `{sys.argv[0]}`.\n"

    # <HACK>
    message += (
        'WARNING: Inspector checks all disabled due to failure to run, or'
        ' timeout issues.'
    )
    # </HACK>

    # Run Intel checks with each set of modules.
    for module_list_file in module_list_files:
        if verbose:
            print('Performing Intel Inspector checks with module set '
                  f"{module_list_file}.")

        # Extract the name of the list.
        module_list_name = module_list_file.rstrip('.lst')
        if debug:
            print(f"module_list_name = {module_list_name}.")

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

        # Add the additional flags needed for Intel Inspector checks.
        cmake_options += ' -DDISABLE_DEBUG_BOUNDS_CHECKS=ON'
        cmake_options += ' -DCMAKE_BUILD_TYPE=DEBUG'
        if debug:
            print(f"cmake_options = {cmake_options}")

        # Make a directory for this test, and go there.
        dir_name = f"intelChecks_{module_list_name}"
        dir_name = f"{INTEL_CHECKS_DIRECTORY_PREFIX}{module_list_name}"
        build_directory = os.path.join(kaiju_home, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge; module load {' '.join(module_names)}"
        )
        module_cmd += '; module list'
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Run cmake to build the Makefile.
        cmd = (
            f"{module_cmd}; {cmake_environment} cmake {cmake_options}"
            f" {kaiju_home}"
        )
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            message += 'cmake failed.\n'
            message += f"e.cmd = {e.cmd}\n"
            message += f"e.returncode = {e.returncode}\n"
            message += 'See test log for output.\n'
            message += (
                'Skipping remaining steps for module set'
                f" {module_list_file}.\n"
            )
            continue

        # Run the build.
        cmd = f"{module_cmd}; make gamera_mpi voltron_mpi"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            message += 'make failed.\n'
            message += f"e.cmd = {e.cmd}\n"
            message += f"e.returncode = {e.returncode}\n"
            message += 'See test log for output.\n'
            message += (
                'Skipping remaining steps for module set'
                f" {module_list_file}.\n"
            )
            continue

        # Go to the bin directory for testing.
        os.chdir('bin')

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
            from_path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
                                     filename)
            to_path = os.path.join(cwd, filename)
            shutil.copyfile(from_path, to_path)

        # Run each check in its own PBS job.
        pbs_files = [
            # 'intelCheckSubmitMem.pbs',
            # 'intelCheckSubmitThread.pbs',
            # 'intelCheckReportSubmit.pbs',
        ]
        job_ids = []
        for pbs_file in pbs_files:
            cmd = (
                f"qsub -A {account}"
                f" -v MODULE_LIST='{' '.join(module_names)}',"
                f"KAIJUROOTDIR={kaiju_home}")
            if pbs_file == 'intelCheckReportSubmit.pbs':
                cmd += f" -W depend=after:{':'.join(job_ids)}"
            cmd += f" {pbs_file}"
            if debug:
                print(f"cmd = {cmd}")
            try:
                cproc = subprocess.run(cmd, shell=True, check=True,
                                       text=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                message += 'qsub failed.\n'
                message += f"e.cmd = {e.cmd}\n"
                message += f"e.returncode = {e.returncode}\n"
                message += 'See test log for output.\n'
                message += (
                    'Skipping remaining steps for module set'
                    f" {module_list_file}.\n"
                )
                continue
            job_id = cproc.stdout.split('.')[0]
            if debug:
                print(f"job_id = {job_id}")
            job_ids.append(job_id)

            # Update the Slack message.
            message += f"Submitted {pbs_file} as job {job_id}."

        # Record the job IDs.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            for job_id in job_ids:
                f.write(f"{job_id}\n")

    # If this is a test run, don't post to Slack. Otherwise, if loud,
    # send Slack message.
    if is_test:
        pass
    elif be_loud:
        if debug:
            print('Sending summary message to Slack.')
        common.slack_send_message(slack_client, message)

    # Send message to stdout.
    print(message)

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
