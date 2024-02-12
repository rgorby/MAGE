#!/usr/bin/env python


"""Master script for automated MAGE regression testing.

This script controls automated MAGE regression testing on
derecho. This script should be run in a cron job using the testing
setup script run_mage_test_00.sh.

Authors
-------
Jeff Garretson (jeffrey.garretson@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)

"""


# Import standard modules.
import argparse
import datetime
import os
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = 'Master script for MAGE regression testing'


def create_command_line_parser():
    """Create the command-line argument parser.

    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.

    Raises
    ------
    None
    """
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        '--account', default='UJHB0019',
        help="PBS account to use for testing (default: %(default)s)"
    )
    parser.add_argument(
        '--all', '-a', action='store_true',
        help="Run all tests (default: %(default)s)."
    )
    parser.add_argument(
        '--debug', '-d', action='store_true',
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        '--force', '-f', action='store_true',
        help="Force all tests to run (default: %(default)s)."
    )
    parser.add_argument(
        '--loud', '-l', action='store_true',
        help="Enable loud mode (post results to Slack) (default: %(default)s)."
    )
    parser.add_argument(
        '--test', '-t', action='store_true',
        help="Enable testing mode (default: %(default)s)."
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help="Print verbose output (default: %(default)s)."
    )
    return parser


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
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    account = args.account
    doAll = args.all
    debug = args.debug
    forceRun = args.force
    beLoud = args.loud
    isTest = args.test
    verbose = args.verbose

    if verbose:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")

    #--------------------------------------------------------------------------

    # Move to the MAGE installation directory.

    # Fetch the directory containing this running script.
    called_from = os.path.dirname(os.path.abspath(__file__))
    if debug:
        print(f"called_from = {called_from}")

    # Assume this script is in a subdirectory of the kaiju directory,
    # and go there.
    path = os.path.join(called_from, '..')
    if debug:
        print(f"path = {path}")
    os.chdir(path)

    #--------------------------------------------------------------------------

    # Check to see if any changes have been made in the MAGE
    # repository.
    if verbose:
        print('Attempting git pull ...')
    cmd = 'git pull'
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    git_pull_output = cproc.stdout.rstrip()
    if debug:
        print(git_pull_output)

    # Find the current branch.
    cmd = 'git symbolic-ref --short HEAD'
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    git_branch = cproc.stdout.rstrip()
    if debug:
        print(f"git_branch = {git_branch}")

    # If there are no changes to the code, and the forced-test flag is
    # not set, abort this test.
    if git_pull_output == 'Already up to date.' and not forceRun:
        print(f"No test today. Branch {git_branch} is already up to date!")
        exit(0)

    #--------------------------------------------------------------------------

    # Run the tests.

    # Move to the directory containing the test scripts.
    os.chdir('testingScripts')
    if debug:
        print(f"cwd = {os.getcwd()}")

    # Assemble command-line flags to pass to the individual test
    # scripts.
    subArgString = ''
    if doAll:
        subArgString += ' -a'
    if debug:
        subArgString += ' -d'
    if forceRun:
        subArgString += ' -f'
    if beLoud:
        subArgString += ' -l'
    if isTest:
        subArgString += ' -t'
    if verbose:
        subArgString += ' -v'
    subArgString += f" --account {account}"
    if debug:
        print(f"subArgString = {subArgString}")

    # Run the tests.
    if doAll:

        # Run all tests.
        if verbose:
            print('Running all tests.')

        # if verbose:
        #     print('Running build tests.')
        # cmd = f"python buildTest.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        if verbose:
            print('Running Fortran unit tests.')
        cmd = f"python unitTest.py {subArgString}"
        if debug:
            print(f"cmd = {cmd}")
        cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running Intel Inspector checks.')
        # cmd = f"python intelChecks.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running initial condition tests.')
        # cmd = f"python ICtest.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running initial conditions test report.')
        # cmd = f"python ICtestReport.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running python unit tests.')
        # cmd = f"python pyunitTest.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running weekly dash.')
        # cmd = f"python weeklyDash.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

    else:

        # Run only typical tests.
        if verbose:
            print('Running typical tests.')

        # if verbose:
        #     print('Running build tests.')
        # cmd = f"python buildTest.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running Fortran unit tests.')
        # cmd = f"python unitTest.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running initial condition tests.')
        # cmd = f"python ICtest.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running initial conditions test report.')
        # cmd = f"python ICtestReport.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

        # if verbose:
        #     print('Running python unit tests.')
        # cmd = f"python pyunitTest.py {subArgString}"
        # if debug:
        #     print(f"cmd = {cmd}")
        # cproc = subprocess.run(cmd, shell=True, check=True)

    if verbose:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
