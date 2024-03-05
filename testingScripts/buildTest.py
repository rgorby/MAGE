#!/usr/bin/env python

"""Run MAGE build regression tests.

This script runs a series of builds of the MAGE software using sets of
modules listed in files under:

$KAIJUHOME/testingScripts/mage_build_modules

This script reads the file build_test.lst from this directory, and
uses the contents as a list of module list files to use for MAGE build
tests.

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

# Name of kaiju subdirectory to use for building executable list
EXECUTABLE_LIST_BUILD_DIRECTORY = 'testFolder'

# Subdirectory of KAIJUHOME containing the test scripts
KAIJU_TEST_SCRIPTS_DIRECTORY = 'testingScripts'

# Subdirectory of KAIJU_TEST_SCRIPTS_DIRECTORY containing module lists
MODULE_LIST_DIRECTORY = 'mage_build_test_modules'

# Name of module list file to use when generating the list of executables
EXECUTABLE_LIST_BUILD_MODULE_LIST = '01.lst'

# Name of file containing names of modules lists to use for build tests
BUILD_TEST_LIST_FILE = 'build_test.lst'

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

    # Clean up from previous tests.
    if verbose:
        print('Cleaning up from previous tests.')
    directories = glob.glob(BUILD_TEST_DIRECTORY_GLOB_PATTERN)
    directories.append(EXECUTABLE_LIST_BUILD_DIRECTORY)
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

    # Find the current branch.
    git_branch_name = common.git_get_branch_name()
    if debug:
        print(f"git_branch_name = {git_branch_name}")

    #--------------------------------------------------------------------------

    # Do a preliminary cmake run to generate the list of executables.

    if verbose:
        print(f"Generating list of executables for branch {git_branch_name}.")

    # Make and move to the preliminary build folder.
    os.mkdir(EXECUTABLE_LIST_BUILD_DIRECTORY)
    os.chdir(EXECUTABLE_LIST_BUILD_DIRECTORY)

    # Read the module list file for building the executable list,
    # extracting cmake environment and options, if any.
    path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
                        MODULE_LIST_DIRECTORY,
                        EXECUTABLE_LIST_BUILD_MODULE_LIST)
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
    module_cmd = f"module --force purge; module load {' '.join(module_names)}"
    if debug:
        print(f"module_cmd = {module_cmd}")

    # Run cmake to build the Makefile.
    cmd = (f"{module_cmd}; {cmake_environment} cmake {cmake_options} "
           f"{kaiju_home}")
    if debug:
        print(f"cmd = {cmd}")
    # <HACK> To ignore cmake error on bcwind.h5 for now.
    try:
        cproc = subprocess.run(cmd, shell=True, check=True)
    except:
        pass
    # </HACK>

    # Build the list of executable targets.
    cmd = f"{module_cmd}; make help | grep '\.x'"
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    executable_list = cproc.stdout.splitlines()

    # Remove the first four characters (dots and spaces).
    executable_list = [e[4:] for e in executable_list]
    if debug:
        print(f"executable_list = {executable_list}")

    #--------------------------------------------------------------------------

    # Make a list of module sets to build with.

    # Read the list of  module sets to use for build tests.
    path = os.path.join(kaiju_home, KAIJU_TEST_SCRIPTS_DIRECTORY,
                        MODULE_LIST_DIRECTORY, BUILD_TEST_LIST_FILE)
    with open(path, encoding='utf-8') as f:
        lines = f.readlines()
    module_list_files = [s.rstrip() for s in lines]
    if debug:
        print(f"module_list_files = {module_list_files}")

    #--------------------------------------------------------------------------

    # Initialize the build test report string.
    message = f"Running {sys.argv[0]}.\n"

    # Do a build with each set of modules.
    for module_list_file in module_list_files:
        if verbose:
            print('Performing build test with module_list_file '
                  f"{module_list_file}")

        # Create a test result message.
        message += (
            f"Building MAGE branch {git_branch_name} with module set "
            f"{module_list_file}.\n"
        )

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

        # Make a directory for this build, and go there.
        dir_name = f"{BUILD_TEST_DIRECTORY_PREFIX}{module_list_name}"
        build_directory = os.path.join(kaiju_home, dir_name)
        if debug:
            print(f"build_directory = {build_directory}")
        os.mkdir(build_directory)
        os.chdir(build_directory)

        # Assemble the commands to load the listed modules.
        module_cmd = (
            f"module --force purge; module load {' '.join(module_names)}"
        )
        if debug:
            print(f"module_cmd = {module_cmd}")

        # Run cmake to build the Makefile.
        cmd = f"{module_cmd}; {cmake_environment} cmake {cmake_options} {kaiju_home}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            message += 'cmake failed.\n'
            message += f"e.cmd = {e.cmd}\n"
            message += f"e.returncode = {e.returncode}\n"
            message += 'See test log for output.\n'
            message += f"Skipping remaining steps for module set {module_list_file}.\n"
            continue

        # Run the build.
        cmd = f"{module_cmd}; make {' '.join(executable_list)}"
        if debug:
            print(f"cmd = {cmd}")
        try:
            cproc = subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            message += 'make failed.\n'
            message += f"e.cmd = {e.cmd}\n"
            message += f"e.returncode = {e.returncode}\n"
            message += 'See test log for output.\n'
            message += f"Skipping remaining steps for module set {module_list_file}.\n"
            continue

        # Check for all executables
        missing = []
        for executable in executable_list:
            path = os.path.join(build_directory, BUILD_BIN_DIR, executable)
            if not os.path.isfile(path):
                missing.append(executable)
        if debug:
            print(f"missing = {missing}")
        if len(missing) > 0:
            for executable in missing:
                message += f"I couldn't build {executable}.\n"
        else:
            message += (
                f"Everything built properly on branch {git_branch_name} "
                f"with module set {module_list_file}.\n"
            )

    # If this is a test run, don't post to Slack. Otherwise, if loud,
    # send Slack message.
    if debug:
        print('Sending build test report to Slack.')
    if is_test:
        pass
    elif be_loud:
        common.slack_send_message(slack_client, message)

    # Send message to stdout.
    print(message)

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
