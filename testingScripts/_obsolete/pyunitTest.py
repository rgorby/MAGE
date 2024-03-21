#!/usr/bin/env python

"""Run MAGE python unit tests.

This script runs a series of unit tests of the MAGE python software.

Authors
-------
Jeff Garretson
Eric Winter
"""


# Import standard modules.
import argparse
import datetime
import os
import subprocess
import sys

# Import 3rd-party modules.

# Import project modules.
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Script for MAGE python unit testing'

# Name of directory for running python unit tests
PYTHON_UNIT_TEST_DIRECTORY = 'pytests'

# Name of python unit test log file
PYTHON_UNIT_TEST_LOG_FILE = 'kaiju-pyunit.txt'


def genPyunitPbsScript(fdir, logname, account):
    headerString = """#!/bin/bash
#PBS -A %s
#PBS -N %s
#PBS -j oe
#PBS -q main
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=128
"""
    moduleString = """
module purge
module load ncarenv/23.06
module load cmake/3.26.3
module load craype/2.7.20
module load intel/2023.0.0
module load ncarcompilers/1.0.0
module load cray-mpich/8.1.25
module load hdf5-mpi/1.12.2
module load mkl/2023.0.0
module list

export MAGE_TEST_ROOT='/glade/work/ewinter/mage_testing/derecho'
export CONDARC="${MAGE_TEST_ROOT}/condarc"
export CONDA_ENVS_PATH="${MAGE_TEST_ROOT}/conda"
export DERECHO_TESTING_ACCOUNT=UJHB0019
export SLACK_BOT_TOKEN='xoxb-1065817665921-1413594823303-gUePq3obrqlPmlCHC5E7rKVP'

mage_miniconda3="${MAGE_TEST_ROOT}/miniconda3"
mage_conda="${mage_miniconda3}/bin/conda"
echo "mage_conda=$mage_conda"
__conda_setup="$($mage_conda 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$mage_miniconda3/etc/profile.d/conda.sh" ]; then
        . "$mage_miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="$mage_miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate kaiju-3.8-testing
"""
    commandString = """
source %s
source $HOME/local/cdf/3.9.0/bin/definitions.B
cd %s
date
echo 'Running pytest'
pytest > %s
date
"""
    pbsFileName = os.path.join(fdir, 'pyunit.pbs')
    with open(pbsFileName, 'w', encoding='utf-8') as f:
        f.write(headerString % (account, 'pyunit'))
        f.write(moduleString)
        f.write(commandString % (os.path.abspath(os.path.join(fdir, '../scripts/setupEnvironment.sh')),
                                 fdir, logname))
    return pbsFileName


def genPyunitReportPbsScript(fdir, args, account):
    headerString = """#!/bin/bash
#PBS -A %s
#PBS -N %s
#PBS -j oe
#PBS -q main
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=128
"""
    moduleString = """
module --force purge
module load ncarenv/23.06
module load cmake/3.26.3
module load craype/2.7.20
module load intel/2023.0.0
module load ncarcompilers/1.0.0
module load cray-mpich/8.1.25
module load hdf5-mpi/1.12.2
module load mkl/2023.0.0
module list

export MAGE_TEST_ROOT='/glade/work/ewinter/mage_testing/derecho'
export CONDARC="${MAGE_TEST_ROOT}/condarc"
export CONDA_ENVS_PATH="${MAGE_TEST_ROOT}/conda"
export DERECHO_TESTING_ACCOUNT=UJHB0019
export SLACK_BOT_TOKEN='xoxb-1065817665921-1413594823303-gUePq3obrqlPmlCHC5E7rKVP'

mage_miniconda3="${MAGE_TEST_ROOT}/miniconda3"
mage_conda="${mage_miniconda3}/bin/conda"
echo "mage_conda=$mage_conda"
__conda_setup="$($mage_conda 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$mage_miniconda3/etc/profile.d/conda.sh" ]; then
        . "$mage_miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="$mage_miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate kaiju-3.8-testing

"""
    commandString = """
source %s
source $HOME/local/cdf/3.9.0/bin/definitions.B
cd %s
date
echo 'Running Report'
../testingScripts/pyunitReport.py %s
date
"""
    pbsFileName = os.path.join(fdir, 'pyunitReport.pbs')
    with open(pbsFileName, 'w', encoding='utf-8') as f:
        f.write(headerString%(account, 'pyunitReport'))
        f.write(moduleString)
        f.write(commandString%(os.path.abspath(os.path.join(fdir, '../scripts/setupEnvironment.sh')),
                               fdir, ' '.join(args)))
    return pbsFileName


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

    #--------------------------------------------------------------------------

    # Move to the MAGE installation directory.
    kaiju_home = os.environ['KAIJUHOME']
    os.chdir(kaiju_home)

    #--------------------------------------------------------------------------

    # Clean up from previous tests.

    # Move to the pytests folder.
    work = os.path.join(kaiju_home, PYTHON_UNIT_TEST_DIRECTORY)
    os.chdir(work)

    # Remove the old log file.
    try:
        os.remove(PYTHON_UNIT_TEST_LOG_FILE)
    except:
        pass

    #--------------------------------------------------------------------------

    # Run the python unit tests as a PBS job.

    # Create the PBS script for the python unit tests.
    pyPbsName = genPyunitPbsScript(work, PYTHON_UNIT_TEST_LOG_FILE,
                                   account=account)
    if debug:
        print(f"pyPbsName = {pyPbsName}")

    # Submit the unit test script for python.
    cmd = f"qsub {pyPbsName}"
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    readString = cproc.stdout.rstrip()
    if debug:
        print(f"readString = {readString}")
    job_name_1 = readString.split('.')[0]
    if debug:
        print(f"job_name_1 = {job_name_1}")

    #--------------------------------------------------------------------------

    # Run the python unit test report as a PBS job.

    # Assemble additional arguments for the job.
    passargs = []
    if is_test:
        passargs.append('-t')
    if be_loud:
        passargs.append('-l')

    # Create the PBS script for the python unit test reports.
    pyReportName = genPyunitReportPbsScript(work, passargs, account=account)
    if debug:
        print(f"pyReportName = {pyReportName}")

    # Submit the unit test report script for python.
    cmd = f"qsub -W depend=afterany:{job_name_1} {pyReportName}"
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True,
                           capture_output=True)
    readString = cproc.stdout.rstrip()
    if debug:
        print(f"readString = {readString}")
    job_name_2 = readString.split('.')[0]
    if debug:
        print(f"job_name_2 = {job_name_2}")

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
