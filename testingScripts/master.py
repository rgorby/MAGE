#!/usr/bin/env python


"""Top-level script for automated MAGE regression testing.

This script controls automated MAGE regression testing on derecho. This script
should be scheduled to run every night at 2300 CT on derecho. The cron job
entry should look something like this:

30 23 * * * bash -c 'ssh derecho /glade/work/ewinter/mage_testing/derecho/scripts/run_mage_test_00.sh' >> /glade/work/ewinter/mage_testing/derecho/logs/test.out 2>&1

Note that thos crontab must be used on the cron host connected to derecho.

The bash script run_mage_test_00.sh sets up the environment for this script,
then runs it.

Authors
-------
Jeff Garretson (jeffrey.garretson@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import datetime
import os
import sys
import subprocess

# Import 3rd-party modules.

# Import project modules.


# Program constants

# Program description.
DESCRIPTION = "Master script for MAGE regression testing"


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
        "--account", default=None,
        help="PBS account to use for testing (default: %(default)s)"
    )
    parser.add_argument(
        "--all", "-a", action="store_true",
        help="Run all tests (default: %(default)s)."
    )
    parser.add_argument(
        "--debug", "-d", action="store_true",
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "--force", "-f", action="store_true",
        help="Force all tests to run (default: %(default)s)."
    )
    parser.add_argument(
        "--loud", "-l", action="store_true",
        help="Enable loud mode (default: %(default)s)."
    )
    parser.add_argument(
        "--test", "-t", action="store_true",
        help="Enable testing mode (default: %(default)s)."
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true",
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

    if debug:
        print(f"Starting {sys.argv[0]} at {datetime.datetime.now()}")

    #--------------------------------------------------------------------------

    # Determine the path to the MAGE installation to use for testing.

    # Fetch the path to this running script.
    called_from = os.path.dirname(os.path.abspath(__file__))
    if debug:
        print(f"called_from = {called_from}")

    # Assume this script is in a subdirectory of the kaiju directory.
    os.chdir(called_from)
    os.chdir('..')

    # Use this directory as the home directory for testing.
    home = os.getcwd()
    if debug:
        print(f"home = {home}")
    if verbose:
        print('I am the master script. This is my current working directory:')
        print(home)

    #--------------------------------------------------------------------------

    # Check to see if any changes have been made in the MAGE repository.
    os.system('git status')
    if verbose:
        print('Attempting git pull via subprocess...')
        p = subprocess.Popen('git pull', shell=True, stdout=subprocess.PIPE)
        text = p.stdout.read().decode('ascii').rstrip()
    if verbose:
        print(text)

    # Find the current branch.
    p = subprocess.Popen('git symbolic-ref --short HEAD', shell=True, stdout=subprocess.PIPE)
    git_branch = p.stdout.read().decode('ascii').rstrip()
    if debug:
        print(f"git_branch = {git_branch}")

    # If there are no changes to the code, and the forced-test flag is not set,
    # abort this test.
    if text == 'Already up to date.' and not forceRun:
        print("No test today. Branch " + git_branch + " is already up to date!")
        exit(0)

    #--------------------------------------------------------------------------

    # Move to the directory containing the test scripts, and run them.
    os.chdir('testingScripts')
    if verbose:
        print('I made it this far!')
        print(os.path.dirname(os.path.abspath(__file__)))

    # Assemble command-line flags to pass to the individual test scripts.
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

    # Run the tests. Send output to stdout. Wait for each test to finish before
    # running the next test.
    if doAll:
        # Run all tests.
        # buildTest = subprocess.Popen(f"python3 buildTest.py {subArgString}",
        #                              shell=True)
        # buildTest.wait()
        if verbose:
            print("Running unit tests.")
        unitTest = subprocess.Popen(f"python3 unitTest.py {subArgString}",
                                    shell=True)
        unitTest.wait()
        # intelTest = subprocess.Popen(f"python3 intelChecks.py {subArgString}",
        #                              shell=True)
        # intelTest.wait()
        # ICTest = subprocess.Popen(f"python3 ICtest.py {subArgString}",
        #                           shell=True)
        # ICTest.wait()
        # ICReport = subprocess.Popen(f"python3 ICtestReport.py {subArgString}",
        #                             shell=True)
        # ICReport.wait()
        # pyunitTest = subprocess.Popen(f"python3 pyunitTest.py {subArgString}",
        #                               shell=True)
        # pyunitTest.wait()
        # weeklyDash = subprocess.Popen(f"python3 weeklyDash.py {subArgString}",
        #                               shell=True)
        # weeklyDash.wait()
    else:
        pass
    # Run only typical tests.
    # buildTest = subprocess.Popen(f"python3 buildTest.py {subArgString}",
    #                              shell=True, check=True)
    # buildTest.wait()
    # unitTest = subprocess.Popen(f"python3 unitTest.py {subArgString}",
    #                             shell = True)
    # unitTest.wait()
    # intelTest = subprocess.Popen(f"python3 intelChecks.py {subArgString}",
    #                              shell=True)
    # intelTest.wait()
    # ICTest = subprocess.Popen(f"python3 ICtest.py {subArgString}",
    #                           shell=True)
    # ICTest.wait()
    # ICReport = subprocess.Popen(f"python3 ICtestReport.py {subArgString}",
    #                             shell=True)
    # ICReport.wait()
    # pyunitTest = subprocess.Popen(f"python3 pyunitTest.py {subArgString}",
    #                               shell=True)
    # pyunitTest.wait()
    # weeklyDash = subprocess.Popen(f"python3 weeklyDash.py {subArgString}",
    #                               shell=True)
    # weeklyDash.wait()

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
