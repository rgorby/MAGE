#!/usr/bin/env python

import argparse
import os
import subprocess
import pathlib
import datetime

def create_command_line_parser():
    """Create the command-line argument parser.
    
    Create the parser for command-line arguments.
    
    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.
    
    """
    parser = argparse.ArgumentParser(description="Script to help setup automated tests within a kaiju repo")
    
    parser.add_argument(
        "-A", required=True,
        help="Charge code to use when running tests."
    )
    parser.add_argument(
        "-ce", required=True,
        help="Conda environment name to load with conda module"
    )
    parser.add_argument(
        "-p", default="economy",
        help="Job priority to use when running tests (default: %(default)s)."
    )
    
    parser.add_argument(
        "--unitTests",  action='store_true',default=False,
        help="Run unit tests (default: %(default)s)."
    )
    parser.add_argument(
        "--weeklyDash",  action='store_true',default=False,
        help="Run weekly dash (default: %(default)s)."
    )
    parser.add_argument(
        "--compTests",  action='store_true',default=False,
        help="Run default subset of comparative tests (default: %(default)s)."
    )
    parser.add_argument(
        "--compTestsFull",  action='store_true',default=False,
        help="Run full suite of comparative tests (over-rides --compTests) (default: %(default)s)."
    )
    parser.add_argument(
        "--buildTest",  action='store_true',default=False,
        help="Run build tests (default: %(default)s)."
    )
    parser.add_argument(
        "--icTests",  action='store_true',default=False,
        help="Run tests to build Initial Condition files (default: %(default)s)."
    )
    parser.add_argument(
        "--memCheck",  action='store_true',default=False,
        help="Run memory check test (default: %(default)s)."
    )
    parser.add_argument(
        "--threadCheck",  action='store_true',default=False,
        help="Run thread check test (default: %(default)s)."
    )
    
    parser.add_argument(
        "--all",  action='store_true',default=False,
        help="Run all tests (default: %(default)s)."
    )
    
    return parser

def main():
    """Helper script to run automated tests locally inside a kaiju repository
    """
    # Set up the command-line parser.
    parser = create_command_line_parser()
    
    # Parse the command-line arguments.
    args = parser.parse_args()
    
    # find repo home directory
    called_from = os.path.dirname(os.path.abspath(__file__))
    os.chdir(called_from)
    os.chdir('..')
    homeDir = os.getcwd()
    
    # Check for necessary environment variables
    if 'KAIJUHOME' not in os.environ:
        print("The setupEnvironment.sh script must be sourced for the repo this script resides in before calling it.")
        exit()
    if pathlib.Path(homeDir).resolve() != pathlib.Path(os.environ['KAIJUHOME']).resolve():
        print("The setupEnvironment.sh script must be sourced for the repo this script resides in before calling it.")
        exit()
    if 'KAIPYHOME' not in os.environ:
        print("The setupEnvironment.sh script for ANY kaipy repo must be sourced before calling this.")
        exit()
    
    # Set environment variables
    os.environ['MAGE_TEST_ROOT'] = homeDir
    os.environ['MAGE_TEST_RUNS_ROOT']=os.path.join(os.environ['MAGE_TEST_ROOT'],"test_runs")
    os.environ['DERECHO_TESTING_ACCOUNT'] = args.A
    os.environ['DERECHO_TESTING_QUEUE'] = 'main'
    os.environ['DERECHO_TESTING_PRIORITY'] = args.p
    os.environ['SLACK_BOT_TOKEN'] = '' # help ensure no accidental reporting to slack
    os.environ['PYTHONUNBUFFERED']='TRUE'
    os.environ['CONDA_ENVIRONMENT']=args.ce
    gitBranch  = subprocess.run(['git','branch','--show-current'], stdout=subprocess.PIPE).stdout.decode('utf-8')
    if 'not a git repository' in gitBranch:
        print("This script must be executed inside a kaiju git repository")
        exit()
    gitBranch = gitBranch.strip()
    os.environ['BRANCH_OR_COMMIT'] = gitBranch
    currenttime = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    test_set_dir = f"{currenttime}-{gitBranch}"
    os.environ['MAGE_TEST_SET_ROOT'] = os.path.join(os.environ['MAGE_TEST_RUNS_ROOT'],test_set_dir)
    os.makedirs(os.environ['MAGE_TEST_SET_ROOT'], exist_ok=True)
    os.chdir(os.environ['MAGE_TEST_SET_ROOT'])
    
    print(f"Running tests on branch {gitBranch}")
    print(f"Using charge code {args.A} with priority {args.p}")
    print(f"Running in folder {test_set_dir}")
    
    # Adjust test options
    if args.all:
        args.unitTests = True
        args.weeklyDash = True
        args.compTests = True
        args.compTestsFull = True
        args.buildTest = True
        args.icTests = True
        args.memCheck = True
        args.threadCheck = True
    
    if args.compTestsFull:
        args.compTests = False
    
    # Run Tests
    if args.unitTests:
        print("Running unit tests")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','unitTest.py'),'-tv'])
    if args.weeklyDash:
        print("Weekly dash not yet supported")
    if args.compTests:
        print("Running default comparative tests subset")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','relativeTests.py'),'-tv'])
    if args.compTestsFull:
        print("Running full comparative tests")
        subprocess.call(['python', os.path.join(os.environ['MAGE_TEST_ROOT'],'testingScripts','relativeTests.py'),'-tva'])
    if args.buildTest:
        print("Build test not yet supported")
    if args.icTests:
        print("Initial Condition tests not yet supported")
    if args.memCheck:
        print("Memory check test not yet supported")
    if args.threadCheck:
        print("Thread check test not yet supported")

if __name__ == "__main__":
    main()

