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
import shutil
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Run the MAGE weekly dash tests.'

# Glob pattern for weekly dash test direectories
WEEKLY_DASH_DIRECTORY_GLOB_PATTERN = 'weeklyDash_*'

# Source directory for weekly dash restart files
WEEKLY_DASH_RESTART_SRC_DIRECTORY = '/glade/work/ewinter/mage_testing/derecho/dashRestarts'

# Prefix for weekly dash directory name
WEEKLY_DASH_DIRECTORY_PREFIX = 'weeklyDash_'

# Name of working directory containing dash restart files.
WORKING_DASH_RESTART_DIRECTORY = 'dashRestarts'

# Subdirectory of KAIJUHOME containing the test scripts
KAIJU_TEST_SCRIPTS_DIRECTORY = 'testingScripts'

# Subdirectory of KAIJU_TEST_SCRIPTS_DIRECTORY containing module lists
MODULE_LIST_DIRECTORY = 'mage_build_test_modules'

# Name of file containing names of modules lists to use for weekly dash
WEEKLY_DASH_LIST_FILE = 'weekly_dash.lst'

# Subdirectory of build directory containing compiled products to use in tests
BIN_DIR = 'bin'

# List of weekly dash test files to copy
WEEKLY_DASH_TEST_FILES = [
    'weeklyDashGo.xml',
    'weeklyDashGo.pbs',
]


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
    if debug:
        print(f"kaiju_home = {kaiju_home}")
    os.chdir(kaiju_home)
    if verbose:
        print(f"Now in directory {os.getcwd()}.")

    #--------------------------------------------------------------------------

    # Find the current branch.
    git_branch_name = common.git_get_branch_name()
    if debug:
        print(f"git_branch_name = {git_branch_name}")

    #--------------------------------------------------------------------------

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
        except:
            pass

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
        if verbose:
            print(f"<HACK> Trying to remove pFUnit code {path}.")
        try:
            shutil.rmtree(path)
        except:
            pass
        if verbose:
            print(f"</HACK> Done trying to remove pFUnit code {path}.")
    # </HACK>

    #--------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
                        MODULE_LIST_DIRECTORY, WEEKLY_DASH_LIST_FILE)
    if verbose:
        print(f"Reading module list from {path}.")
    module_list_files, _, _ = common.read_build_module_list_file(path)
    if debug:
        print(f"module_list_files = {module_list_files}")

    # <HACK>
    # Just use first module set for now.
    module_list_files = [module_list_files[0]]
    if verbose:
        print(f"<HACK> Only using first module set {module_list_files[0]}! </HACK>.")
    # </HACK>

    # #--------------------------------------------------------------------------

    # Run the tests with each set of modules.
    for module_list_file in module_list_files:
        if verbose:
            print('Performing weekly dash with module set '
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

        # Add the cmake option for the weekly dash build.
        cmake_options += ' -DCMAKE_BUILD_TYPE=Release'
        if debug:
            print(f"cmake_options = {cmake_options}")

        # Make a directory for this build, and go there.
        dir_name = f"{WEEKLY_DASH_DIRECTORY_PREFIX}{module_list_name}"
        build_directory = os.path.join(kaiju_home, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)
        if verbose:
            print(f"Now in directory {os.getcwd()}.")

        # Assemble the commands to load the listed modules.
        module_cmd = f"module --force purge; module load {' '.join(module_names)}"
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Run cmake to build the Makefile.
        if verbose:
            print('Running cmake to generate Makefile.')
        cmd = (f"{module_cmd}; {cmake_environment} cmake {cmake_options} "
               f"{kaiju_home}")
        if debug:
            print(f"cmd = {cmd}")
        # <HACK> Ignore cmake error on bcwind.h5 for now.
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except:
            if verbose:
                print('<HACK> Ignoring cmake error for bcwind.h5 rule. </HACK>')
        if debug:
            print(f"cproc = {cproc}")
        # </HACK>

        # Run the build.
        cmd = f"{module_cmd}; make voltron_mpi.x"
        if debug:
            print(f"cmd = {cmd}")
        if verbose:
            print('Running make to build code.')
        cproc = subprocess.run(cmd, shell=True, check=True, text=True)
        if debug:
            print(f"cproc = {cproc}")

        # Move into the bin directory to run the tests.
        os.chdir(BIN_DIR)
        if verbose:
            print(f"Now in directory {os.getcwd()}.")

        # Copy the restart data and reference results.
        # if verbose:
            print('Copying restart data and reference results for weekly dash.')
        from_glob = os.path.join(WEEKLY_DASH_RESTART_SRC_DIRECTORY, '*')
        for from_path in glob.glob(from_glob):
            if debug:
                print(f"from_path = {from_path}")
            filename = os.path.split(from_path)[-1]
            to_path = os.path.join('.', filename)
            shutil.copyfile(from_path, to_path)

        # Generate the LFM grid file.
        if verbose:
            print('Creating LFM grid file.')
        cmd = 'genLFM.py -gid Q'
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, check=True)
        if debug:
            print(f"cproc = {cproc}")
        shutil.move('lfmQ.h5', 'NEWlfmX.h5')

        # Compare the new LFM grid file to the original.
        cmd = 'h5diff lfmX.h5 NEWlfmX.h5'
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, text=True, capture_output=True)
        if debug:
            print(f"cproc = {cproc}")
        grid_diff = cproc.stdout.rstrip()
        if debug:
            print(f"grid_diff = {grid_diff}")
        if grid_diff != '' or cproc.returncode != 0:
            message = (
                f"Quad grid for weekly dash has changed on branch {git_branch_name}."
                'Case cannot be run. Please re-generate restart data, and ensure the grid change was intentional.'
            )
            if not is_test:
                common.slack_send_message(slack_client, message)
            else:
                print(message)

        # Generate the solar wind boundary condition file.
        if verbose:
            print('Creating solar wind initial conditions file.')
        cmd = 'cda2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00 -o NEWbcwind.h5'
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, check=True)
        if debug:
            print(f"cproc = {cproc}")

        # Compare the new solar wind boundary condition file to the original.
        cmd = 'h5diff bcwind.h5 NEWbcwind.h5'
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, text=True, capture_output=True)
        if debug:
            print(f"cproc = {cproc}")
        wind_diff = cproc.stdout.rstrip()
        if debug:
            print(f"wind_diff = {wind_diff}")
        if wind_diff != '' or cproc.returncode != 0:
            message = (
                f"Solar wind file for weekly dash has changed on branch {git_branch_name}."
                'Case cannot be run. Please re-generate restart data, and ensure the wind data change was intentional.'
            )
            if not is_test:
                common.slack_send_message(slack_client, message)
            else:
                print(message)

        # Generate the RCM configuration file.
        cmd = 'genRCM.py -o NEWrcmconfig.h5'
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, check=True)
        if debug:
            print(f"cproc = {cproc}")

        # Compare the new RCM configuration file to the original.
        cmd = 'h5diff rcmconfig.h5 NEWrcmconfig.h5'
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, text=True, capture_output=True)
        if debug:
            print(f"cproc = {cproc}")
        rcm_diff = cproc.stdout.rstrip()
        if debug:
            print(f"rcm_diff = {rcm_diff}")
        if rcm_diff != '' or cproc.returncode != 0:
            message = (
                f"rcmconfig for weekly dash has changed on branch {git_branch_name}."
                'Case cannot be run. Please re-generate restart data, and ensure the rcmconfig change was intentional.'
            )
            if not is_test:
                common.slack_send_message(slack_client, message)
            else:
                print(message)

        # Copy files needed for the weekly dash job.
        if verbose:
            print('Copying files needed for weekly dash.')
        for f in WEEKLY_DASH_TEST_FILES:
            from_file = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY, f)
            to_file = os.path.join('.', f)
            shutil.copyfile(from_file, to_file)

        # Create the qsub command to submit a model run with these files.
        cmd = f"qsub -A {account} -v MODULE_LIST='{' '.join(module_names)}',KAIJUROOTDIR={kaiju_home} weeklyDashGo.pbs"
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                               capture_output=True)
        if debug:
            print(f"cproc = {cproc}")
        firstJobNumber = cproc.stdout.split('.')[0]
        if debug:
            print(f"firstJobNumber = {firstJobNumber}")

        # Save the job number in a file.
        with open('jobs.txt', 'w', encoding='utf-8') as f:
            f.write(f"{firstJobNumber}\n")

        message = f"Run started on branch {git_branch_name} as jobid {firstJobNumber}."
        print(message)

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
