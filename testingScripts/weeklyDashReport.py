#!/usr/bin/env python


"""Create the MAGE weekly dash test report.

This script creates the MAGE weekly dash test report.

Authors
-------
Jeff Garretson (jeffrey.garretson@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import datetime
import os
import re
import shutil
import subprocess
import sys

# Import 3rd-party modules.
from astropy.time import Time
import h5py
import matplotlib as mpl
mpl.use('Agg')
from slack_sdk.errors import SlackApiError

# Import project modules.
import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv
from kaipy.testing import common


# Program constants

# Program description.
DESCRIPTION = 'Create the MAGE weekly dash test report.'

# Default path to wiki directory to receive results.
DEFAULT_WIKI_PATH = os.path.join(
    os.environ['MAGE_TEST_ROOT'], 'workspace', 'kaijuWiki', 'wiki'
)

# Name of subdirectory containing binaries and test results.
BIN_DIR = 'bin'

# Name of file containg PBS job IDs.
JOB_LIST_FILE = 'jobs.txt'

# Name of directory containing weekly dash results.
WEEKLY_DASH_DIRECTORY = 'weeklyDash_01'

# Name of weekly dash log file.
WEEKLY_DASH_LOG_FILE = 'weeklyDashGo.out'

# Subdirectory of wiki repository for weekly dash results.
WEEKLY_DASH_WIKI_DIRECTORY = 'weeklyDash'

# Regular expression for git hash read from weekly dash output log.
GIT_HASH_PATTERN = 'Git hash   = (.{8})'

# Name of voltron output file.
VOLTRON_OUTPUT_FILE = 'msphere.volt.h5'

# Name of file for previous data.
PREVIOUS_DATA_FILE = 'previousData.h5'


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

    # Add an option for the wiki path.
    parser.add_argument('--wiki_path', '-w', type=str, default=DEFAULT_WIKI_PATH,
                        help='Wiki path to receive results')

    # Parse the command-line arguments.
    args = parser.parse_args()
    if args.debug:
        print(f"args = {args}")
    account = args.account
    debug = args.debug
    be_loud = args.loud
    is_test = args.test
    verbose = args.verbose
    wiki_path = args.wiki_path

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
    if debug:
        print(f"Now in directory {os.getcwd()}.")

    #--------------------------------------------------------------------------

    # Find the current branch.
    git_branch_name = common.git_get_branch_name()
    if debug:
        print(f"git_branch_name = {git_branch_name}")

    #--------------------------------------------------------------------------

    # Make sure the path to the wiki results exists.
    if debug:
        print(f"wiki_path = {wiki_path}")
    if not os.path.exists(wiki_path):
        print(f"Wiki folder does not exist at {wiki_path}. Please create it.")
        sys.exit(1)

    #--------------------------------------------------------------------------

    # Get the short hostname.
    cmd = 'hostname -s'
    cproc = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)
    if debug:
        print(f"cproc = {cproc}")
    short_hostname = cproc.stdout.rstrip()
    if verbose:
        print(f"Running weekly dash test on host {short_hostname}.")

    # Go to weekly dash folder
    path = os.path.join(kaiju_home, WEEKLY_DASH_DIRECTORY)
    try:
        os.chdir(path)
    except:
        print(f"Weekly dash results folder {path} does not exist.")
        sys.exit(1)
    if debug:
        print(f"Now in {os.getcwd()}")

    # Move down to the directory containing the dash results.
    os.chdir(BIN_DIR)
    if debug:
        print(f"Now in {os.getcwd()}")

    # Check for jobs.txt
    if not os.path.exists(JOB_LIST_FILE):
        print(f"Job list file {JOB_LIST_FILE} not found.")
        sys.exit(1)

    # Read in the jobs.txt file to get the job numbers.
    with open(JOB_LIST_FILE, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    job1 = lines[0].rstrip()
    if debug:
        print(f"job1 = {job1}")

    # Find the job output log.
    job1_output_log = f"wDashGo.o{job1}"
    if debug:
        print(f"job1_output_log = {job1_output_log}")
    if not os.path.exists(job1_output_log):
        print('The dash job has not finished yet.')
        sys.exit(1)

    # Move to wiki repository folder and ensure it is up to date.
    os.chdir(wiki_path)
    if debug:
        print(f"Now in {os.getcwd()}")
    cmd = 'git pull'
    cproc = subprocess.run(cmd, shell=True, text=True, check=True, capture_output=True)
    if debug:
        print(f"cproc = {cproc}")
    text = cproc.stdout.rstrip()
    if debug:
        print(f"text = {text}")
    if 'not a git repository' in text:
        print(f"Wiki path {wiki_path} is not a git repository.")
        sys.exit(1)
    os.chdir(WEEKLY_DASH_WIKI_DIRECTORY)
    if debug:
        print(f"Now in {os.getcwd()}")

    # Read in the git hash from the output file.
    git_hash = 'XXXXXXXX'
    path = os.path.join(kaiju_home, WEEKLY_DASH_DIRECTORY, BIN_DIR, WEEKLY_DASH_LOG_FILE)
    if debug:
        print(f"path = {path}")
    with open(path, 'r', encoding='utf-8') as weekly_dash_output:
        for line in weekly_dash_output:
            git_hash_match = re.search(GIT_HASH_PATTERN, line)
            if git_hash_match:
                git_hash = git_hash_match.group(1)
                break
    if debug:
        print(f"git_hash = {git_hash}")

    # Read performance data from the test results.
    cmd = 'sed --quiet "s/^ \\+UT \\+= \\+\\([0-9-]\\+ [0-9:]\\+\\).*$/\\1/p" '
    cmd += os.path.join(kaiju_home, WEEKLY_DASH_DIRECTORY, BIN_DIR, WEEKLY_DASH_LOG_FILE)
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)
    # if debug:
    #     print(f"cproc = {cproc}")
    utData = cproc.stdout
    # if debug:
    #     print(f"utData = {utData}")
    cmd = 'sed --quiet "s/^ \\+Running @ *\\([0-9]\\+\\.\\?[0-9]*\\)% of real-time.*$/\\1/p" '
    cmd += os.path.join(kaiju_home, WEEKLY_DASH_DIRECTORY, BIN_DIR, WEEKLY_DASH_LOG_FILE)
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)
    # if debug:
    #     print(f"cproc = {cproc}")
    rtData = cproc.stdout
    # if debug:
    #     print(f"rtData = {rtData}")

    # <HACK>
    # There is always one extra line of UT data in the front, strip it.
    utData = utData[1:]
    # </HACK>

    # Split UT and RT data into lists
    utData = utData.splitlines()
    rtData = rtData.splitlines()
    # if debug:
    #     print(f"utData = {utData}")
    #     print(f"rtData = {rtData}")

    # <HACK>
    # Make sure the lists are of equal length now (console output not
    # always reliable).
    if len(utData) > len(rtData):
        utData = utData[:len(rtData)]
        # if debug:
        #     print(f"utData = {utData}")
    elif len(rtData) > len(utData):
        rtData = rtData[:len(utData)]
        # if debug:
        #     print(f"rtData = {rtData}")

    # Convert the rt data to floats.
    rtData_f = [float(x) for x in rtData]
    # if debug:
    #     print(f"rtData_f = {rtData_f}")

    #--------------------------------------------------------------------------

    # Get simulation data from voltron output file.
    fVolt = os.path.join(kaiju_home, WEEKLY_DASH_DIRECTORY, BIN_DIR, VOLTRON_OUTPUT_FILE)
    # if debug:
    #     print(f"fVolt = {fVolt}")

    # Get Dst and CPCP plot info.
    if verbose:
        print(f"Reading Dst and CPCP data from {fVolt}.")
    nSteps, sIds = kh5.cntSteps(fVolt)
    symh = kh5.getTs(fVolt, sIds, 'SymH')
    MJD = kh5.getTs(fVolt, sIds, 'MJD')
    BSDst = kh5.getTs(fVolt, sIds, 'BSDst')
    nCPCP = kh5.getTs(fVolt, sIds, 'cpcpN')
    sCPCP = kh5.getTs(fVolt, sIds, 'cpcpS')
    UT = Time(MJD, format='mjd').isot
    ut_datetime = [datetime.datetime.strptime(ut, '%Y-%m-%dT%H:%M:%S.%f') for ut in UT] # needed to plot symh
    # if debug:
    #     print(f"nSteps = {nSteps}")
    #     print(f"sIds = {sIds}")
    #     print(f"MJD = {MJD}")
    #     print(f"BSDst = {BSDst}")
    #     print(f"nCPCP = {nCPCP}")
    #     print(f"sCPCP = {sCPCP}")
    #     print(f"UT = {UT}")
    #     print(f"ut_datetime = {ut_datetime}")

    # Load old data from an h5 file.
    masterUT = None
    masterRT = None
    masterUTsim = None
    masterDST = None
    masterCPCPn = None
    materCPCPs = None
    masterGitHash = None
    devpriorUT = None
    devpriorRT = None
    devpriorUTsim = None
    devpriorDST = None
    devpriorCPCPn = None
    devpriorCPCPs = None
    devpriorGitHash = None
    devcurrentUT = None
    devcurrentRT = None
    devcurrentUTsim = None
    devcurrentDST = None
    devcurrentCPCPn = None
    devcurrentCPCPs = None
    devcurrentGitHash = None
    if os.path.exists(PREVIOUS_DATA_FILE):
        with h5py.File(PREVIOUS_DATA_FILE, 'r') as data_object:
            if 'masterUT' in data_object:
                masterUT = [x.decode('utf-8') for x in data_object['masterUT']]
                masterRT = data_object['masterRT'][...]
                masterUTsim = [x.decode('utf-8') for x in data_object['masterUTsim']]
                masterDST = data_object['masterDST'][...]
                masterCPCPn = data_object['masterCPCPn'][...]
                masterCPCPs = data_object['masterCPCPs'][...]
                masterGitHash = data_object['masterGitHash'][()].decode('utf-8')
            if 'devpriorUT' in data_object:
                devpriorUT = [x.decode('utf-8') for x in data_object['devpriorUT']]
                devpriorRT = data_object['devpriorRT'][...]
                devpriorUTsim = [x.decode('utf-8') for x in data_object['devpriorUTsim']]
                devpriorDST = data_object['devpriorDST'][...]
                devpriorCPCPn = data_object['devpriorCPCPn'][...]
                devpriorCPCPs = data_object['devpriorCPCPs'][...]
                devpriorGitHash = data_object['devpriorGitHash'][()].decode('utf-8')
            if 'devcurrentUT' in data_object:
                devcurrentUT = [x.decode('utf-8') for x in data_object['devcurrentUT']]
                devcurrentRT = data_object['devcurrentRT'][...]
                devcurrentUTsim = [x.decode('utf-8') for x in data_object['devcurrentUTsim']]
                devcurrentDST = data_object['devcurrentDST'][...]
                devcurrentCPCPn = data_object['devcurrentCPCPn'][...]
                devcurrentCPCPs = data_object['devcurrentCPCPs'][...]
                devcurrentGitHash = data_object['devcurrentGitHash'][()].decode('utf-8')

    # Update appropriate data with new data.
    if git_branch_name == 'master':
        masterUT = utData
        masterRT = rtData_f
        masterUTsim = UT
        masterDST = BSDst
        masterCPCPn = nCPCP
        masterCPCPs = sCPCP
        masterGitHash = git_hash
    elif git_branch_name == 'development':
        devpriorUT = devcurrentUT
        devpriorRT = devcurrentRT
        devpriorUTsim = devcurrentUTsim
        devpriorDST = devcurrentDST
        devpriorCPCPn = devcurrentCPCPn
        devpriorCPCPs = devcurrentCPCPs
        devpriorGitHash = devcurrentGitHash
        devcurrentUT = utData
        devcurrentRT = rtData_f
        devcurrentUTsim = UT
        devcurrentDST = BSDst
        devcurrentCPCPn = nCPCP
        devcurrentCPCPs = sCPCP
        devcurrentGitHash = git_hash

    # Convert date strings into date-time objects
    if masterUT is not None:
        masterUTdt = [datetime.datetime.strptime(ut, '%Y-%m-%d %H:%M:%S') for ut in masterUT]
        masterUTsimdt = [datetime.datetime.strptime(ut, '%Y-%m-%dT%H:%M:%S.%f') for ut in masterUTsim]
    if devpriorUT is not None:
        devpriorUTdt = [datetime.datetime.strptime(ut, '%Y-%m-%d %H:%M:%S') for ut in devpriorUT]
        devpriorUTsimdt = [datetime.datetime.strptime(ut, '%Y-%m-%dT%H:%M:%S.%f') for ut in devpriorUTsim]
    if devcurrentUT is not None:
        devcurrentUTdt = [datetime.datetime.strptime(ut, '%Y-%m-%d %H:%M:%S') for ut in devcurrentUT]
        devcurrentUTsimdt = [datetime.datetime.strptime(ut, '%Y-%m-%dT%H:%M:%S.%f') for ut in devcurrentUTsim]

    #--------------------------------------------------------------------------

    # Make the real-time performance plot.

    # Plot parameters
    line_width = 0.75 # LineWidth
    alpha = 0.25 # Transparency
    grid_color = 'slategrey' # GridColor
    figsize = (14, 7)

    # Create the figure to hold the plot.
    fig = mpl.pyplot.figure(figsize=figsize)
    gs = mpl.gridspec.GridSpec(1, 1, hspace=0.05, wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])

    # Create the plot.
    if masterRT is not None:
        ax.plot(masterUTdt, masterRT, label=f"master ({masterGitHash})", linewidth=line_width)
    if devpriorRT is not None:
        ax.plot(devpriorUTdt, devpriorRT, label=f"dev prior ({devpriorGitHash})", linewidth=line_width)
    if devcurrentRT is not None:
        ax.plot(devcurrentUTdt, devcurrentRT, label=f"dev current ({devcurrentGitHash})", linewidth=line_width)
    ax.legend(loc='lower right', fontsize='small')
    ax.minorticks_on()
    ax.xaxis_date()
    xfmt = mpl.dates.DateFormatter('%H:%M \n%Y-%m-%d')
    ax.set_ylabel('Percent of Real-Time [%]')
    ax.xaxis.set_major_formatter(xfmt)
    ax.xaxis.set_minor_locator(mpl.dates.HourLocator())
    # mpl.pyplot.grid(True)
    # ax.xaxis.grid(True,which='major',linewidth=line_width  ,alpha=alpha,color=grid_color)
    # ax.xaxis.grid(True,which='minor',linewidth=line_width/4,alpha=alpha,color=grid_color)
    # ax.yaxis.grid(True,which='major',linewidth=line_width  ,alpha=alpha,color=grid_color)
    ax.set_title('Real-Time Performance')
    fOut = 'perfPlots.png'
    kv.savePic(fOut)
    mpl.pyplot.close('all')

    #--------------------------------------------------------------------------

    # Make the DST plot.

    # Create the figure to hold the plot.
    fig = mpl.pyplot.figure(figsize=figsize)
    gs = mpl.gridspec.GridSpec(1, 1, hspace=0.05, wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])

    # Create the plot.
    ax.plot(ut_datetime, symh, label='SYM-H', linewidth=2*line_width)
    if masterDST is not None:
        ax.plot(masterUTsimdt, masterDST, label=f"master ({masterGitHash})", linewidth=line_width)
    if devpriorDST is not None:
        ax.plot(devpriorUTsimdt, devpriorDST, label=f"dev prior ({devpriorGitHash})", linewidth=line_width)
    if devcurrentDST is not None:
        ax.plot(devcurrentUTsimdt, devcurrentDST, label=f"dev current ({devcurrentGitHash})", linewidth=line_width)
    ax.legend(loc='upper right', fontsize='small')
    ax.minorticks_on()
    ax.xaxis_date()
    xfmt = mpl.dates.DateFormatter('%H:%M \n%Y-%m-%d')
    ax.set_ylabel('Dst [nT]')
    ax.xaxis.set_major_formatter(xfmt)
    ax.xaxis.set_minor_locator(mpl.dates.HourLocator())
    # mpl.pyplot.grid(True)
    # ax.xaxis.grid(True,which='major',linewidth=line_width  ,alpha=alpha,color=grid_color)
    # ax.xaxis.grid(True,which='minor',linewidth=line_width/4,alpha=alpha,color=grid_color)
    # ax.yaxis.grid(True,which='major',linewidth=line_width  ,alpha=alpha,color=grid_color)
    ax.set_title('BSDst')
    fOut = 'dstPlots.png'
    kv.savePic(fOut)
    mpl.pyplot.close('all')

    #--------------------------------------------------------------------------

    # Make the CPCP plot.

    # Create the figure to hold the plot.
    fig = mpl.pyplot.figure(figsize=figsize)
    gs = mpl.gridspec.GridSpec(1, 1, hspace=0.05, wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])

    # Create the plot.
    if masterCPCPn is not None:
        ax.plot(masterUTsimdt, masterCPCPn, color='orange', linestyle='dotted',
                label=f"master-N ({masterGitHash})", linewidth=line_width)
        ax.plot(masterUTsimdt, masterCPCPs, color='blue', linestyle='dotted',
                label=f"master-S ({masterGitHash})", linewidth=line_width)
    if devpriorCPCPn is not None:
        ax.plot(devpriorUTsimdt, devpriorCPCPn, color='orange', linestyle='dashed',
                label=f"dev prior-N ({devpriorGitHash})", linewidth=line_width)
        ax.plot(devpriorUTsimdt, devpriorCPCPs, color='blue', linestyle='dashed',
                label=f"dev prior-S ({devpriorGitHash})", linewidth=line_width)
    if devcurrentCPCPn is not None:
        ax.plot(devcurrentUTsimdt, devcurrentCPCPn, color='orange', linestyle='solid',
                label=f"dev current-N ({devcurrentGitHash})", linewidth=line_width)
        ax.plot(devcurrentUTsimdt, devcurrentCPCPs, color='blue', linestyle='solid',
                label=f"dev current-S ({devcurrentGitHash})", linewidth=line_width)
    ax.legend(loc='upper right', fontsize='small')
    ax.minorticks_on()
    ax.xaxis_date()
    xfmt = mpl.dates.DateFormatter('%H:%M \n%Y-%m-%d')
    ax.set_ylabel('CPCP [kV]')
    ax.xaxis.set_major_formatter(xfmt)
    ax.xaxis.set_minor_locator(mpl.dates.HourLocator())
    # mpl.pyplot.grid(True)
    # ax.xaxis.grid(True,which='major',linewidth=line_width  ,alpha=alpha,color=grid_color)
    # ax.xaxis.grid(True,which='minor',linewidth=line_width/4,alpha=alpha,color=grid_color)
    # ax.yaxis.grid(True,which='major',linewidth=line_width  ,alpha=alpha,color=grid_color)
    ax.set_title('CPCP')
    fOut = 'cpcpPlots.png'
    kv.savePic(fOut)
    mpl.pyplot.close('all')

    #--------------------------------------------------------------------------

    # Save the new data as h5.
    with h5py.File('previousData.h5', 'w') as data_object:
        if masterUT is not None:
            data_object.create_dataset('masterUT',  data=[x.encode('utf-8') for x in masterUT])
            data_object.create_dataset('masterRT', data=masterRT)
            data_object.create_dataset('masterUTsim', data=[x.encode('utf-8') for x in masterUTsim])
            data_object.create_dataset('masterDST', data=masterDST)
            data_object.create_dataset('masterCPCPn', data=masterCPCPn)
            data_object.create_dataset('masterCPCPs', data=masterCPCPs)
            data_object.create_dataset('masterGitHash', data=masterGitHash.encode('utf-8'))
        if devpriorUT is not None:
            data_object.create_dataset('devpriorUT', data=[x.encode('utf-8') for x in devpriorUT])
            data_object.create_dataset('devpriorRT', data=devpriorRT)
            data_object.create_dataset('devpriorUTsim', data=[x.encode('utf-8') for x in devpriorUTsim])
            data_object.create_dataset('devpriorDST', data=devpriorDST)
            data_object.create_dataset('devpriorCPCPn', data=devpriorCPCPn)
            data_object.create_dataset('devpriorCPCPs', data=devpriorCPCPs)
            data_object.create_dataset('devpriorGitHash', data=devpriorGitHash.encode('utf-8'))
        if devcurrentUT is not None:
            data_object.create_dataset('devcurrentUT', data=[x.encode('utf-8') for x in devcurrentUT])
            data_object.create_dataset('devcurrentRT', data=devcurrentRT)
            data_object.create_dataset('devcurrentUTsim', data=[x.encode('utf-8') for x in devcurrentUTsim])
            data_object.create_dataset('devcurrentDST', data=devcurrentDST)
            data_object.create_dataset('devcurrentCPCPn', data=devcurrentCPCPn)
            data_object.create_dataset('devcurrentCPCPs', data=devcurrentCPCPs)
            data_object.create_dataset('devcurrentGitHash', data=devcurrentGitHash.encode('utf-8'))

    # If I'm on development, copy latest quick look plots over old ones
    if git_branch_name == 'development':
        shutil.copyfile('development_qk_msph.png', 'development_qk_msph_old.png')
        shutil.copyfile('development_qk_rcm.png', 'development_qk_rcm_old.png')
        shutil.copyfile('development_qk_mix.png', 'development_qk_mix_old.png')

    #--------------------------------------------------------------------------

    # Make quick-look plots

    # Move back to the directory containing the simulation results.
    path = os.path.join(kaiju_home, WEEKLY_DASH_DIRECTORY, BIN_DIR)
    if debug:
        print(f"path = {path}")
    os.chdir(path)
    if debug:
        print(f"Now in {os.getcwd()}")

    # Magnetosphere quicklook
    cmd = 'msphpic.py'
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")
    from_path = 'qkmsphpic.png'
    to_path = os.path.join(wiki_path, WEEKLY_DASH_WIKI_DIRECTORY, f"{git_branch_name}_qk_msph.png")
    shutil.copyfile(from_path, to_path)

    # REMIX quicklook
    cmd = 'mixpic.py'
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")
    from_path = 'remix_n.png'
    to_path = os.path.join(wiki_path, WEEKLY_DASH_WIKI_DIRECTORY, f"{git_branch_name}_qk_mix.png")
    shutil.copyfile(from_path, to_path)

    # RCM quicklook
    cmd = 'rcmpic.py'
    if debug:
        print(f"cmd = {cmd}")
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")
    from_path = 'qkrcmpic.png'
    to_path = os.path.join(wiki_path, WEEKLY_DASH_WIKI_DIRECTORY, f"{git_branch_name}_qk_rcm.png")
    shutil.copyfile(from_path, to_path)

    #--------------------------------------------------------------------------

    # Move back to the wiki repository.
    os.chdir(wiki_path)
    os.chdir(WEEKLY_DASH_WIKI_DIRECTORY)

    # Combine quick looks into larger images
    cmd = f"convert master_qk_msph.png -gravity NorthWest -pointsize 60 -annotate +0+0 'master ({masterGitHash})' mm.png"
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    cmd = f"convert master_qk_mix.png -gravity NorthWest -pointsize 80 -annotate +0+0 'master ({masterGitHash})' mx.png"
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    cmd = f"convert master_qk_rcm.png -gravity NorthWest -pointsize 80 -annotate +0+0 'master ({masterGitHash})' mr.png"
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    cmd = f"convert development_qk_msph_old.png -gravity NorthWest -pointsize 60 -annotate +0+0 'dev prior ({devpriorGitHash})' dmo.png"
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    cmd = f"convert development_qk_mix_old.png -gravity NorthWest -pointsize 80 -annotate +0+0 'dev prior ({devpriorGitHash})' dxo.png"
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    cmd = f"convert development_qk_rcm_old.png -gravity NorthWest -pointsize 80 -annotate +0+0 'dev prior ({devpriorGitHash})' dro.png"
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    # cmd = f"convert development_qk_msph.png -gravity NorthWest -pointsize 60 -annotate +0+0 'dev current ({devcurrentGitHash})' dm.png"
    # cproc = subprocess.run(cmd, shell=True, check=True)
    # if debug:
    #     print(f"cproc = {cproc}")

    # cmd = f"convert development_qk_mix.png -gravity NorthWest -pointsize 80 -annotate +0+0 'dev current ({devcurrentGitHash})' dx.png"
    # cproc = subprocess.run(cmd, shell=True, check=True)
    # if debug:
    #     print(f"cproc = {cproc}")

    # cmd = f"convert development_qk_rcm.png -gravity NorthWest -pointsize 80 -annotate +0+0 'dev current ({devcurrentGitHash})' dr.png"
    # cproc = subprocess.run(cmd, shell=True, check=True)
    # if debug:
    #     print(f"cproc = {cproc}")

    # <HACK>
    # For testing
    shutil.copyfile('dmo.png', 'dm.png')
    shutil.copyfile('dxo.png', 'dx.png')
    shutil.copyfile('dro.png', 'dr.png')
    # </HACK>

    cmd = 'convert mm.png dmo.png dm.png -append combined_qk_msph.png'
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    cmd = 'convert mx.png dxo.png dx.png -append combined_qk_mix.png'
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    cmd = 'convert mr.png dro.png dr.png -append combined_qk_rcm.png'
    cproc = subprocess.run(cmd, shell=True, check=True)
    if debug:
        print(f"cproc = {cproc}")

    # Delete the original images.
    for filename in ['mm.png', 'mx.png', 'mr.png',
                     'dmo.png', 'dxo.png', 'dro.png',
                     'dm.png', 'dx.png', 'dr.png']:
        if os.path.exists(filename):
            os.remove(filename)

    # Push the data to the wiki.
    # <HACK>
    debug = False
    # </HACK>
    if not debug:
        cmd = f"git commit -a -m 'New weekly dash data for branch {git_branch_name}'"
        cproc = subprocess.run(cmd, shell=True, check=True)
        if debug:
            print(f"cproc = {cproc}")
        if debug:
            print(f"text = {text}")
        cmd = 'git push'
        cproc = subprocess.run(cmd, shell=True, check=True)
        if debug:
            print(f"cproc = {cproc}")
        if debug:
            print(f"text = {text}")

    #--------------------------------------------------------------------------

    # Post results to Slack.

    # Set up for communication with Slack.
    slack_client = common.slack_create_client()
    if debug:
        print(f"slack_client = {slack_client}")

    if not is_test and be_loud:
        if verbose:
            print('Posting test results to Slack.')
        message = f"Weekly results complete on branch {git_branch_name}, run on host {short_hostname}. Latest comparative results attached as replies to this message.\nOr up-to-date results can be viewed on the wiki at https://bitbucket.org/aplkaiju/kaiju/wiki/weeklyDash/dashStatus"
        slack_response = common.slack_send_message(slack_client, message)
        if slack_response['ok']:
            parent_ts = slack_response['ts']
            message = 'This was a 4x4x1 (IxJxK) decomposed Quad Resolution Run using 8 nodes for Gamera, 1 for Voltron, and 3 Squish Helper nodes (12 nodes total)'
            slack_response = common.slack_send_message(slack_client, message,
                                                       thread_ts=parent_ts)
            image_files = [
                'perfPlots.png',
                'dstPlots.png',
                'cpcpPlots.png',
                'combined_qk_msph.png',
                'combined_qk_mix.png',
                'combined_qk_rcm.png',
            ]
            image_comments = [
                'Real-Time Performance\n\n',
                'DST Plots\n\n',
                'CPCP Plots\n\n',
                'Quick-Looks Magnetosphere\n\n',
                'Quick-Looks Remix\n\n',
                'Quick-Looks RCM\n\n',
            ]
            for (f, c) in zip(image_files, image_comments):
                try:
                    slack_response = common.slack_send_image(
                        slack_client, f, initial_comment=c,
                        thread_ts=parent_ts,
                    )
                    assert slack_response['ok']
                except SlackApiError as e:
                    # You will get a SlackApiError if "ok" is False
                    assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
        else:
            print("Failed to post parent message to Slack, could not attach replies either.")

    else:
        print(f"weekly run completed successfully on branch {git_branch_name}")

    if debug:
        print(f"Ending {sys.argv[0]} at {datetime.datetime.now()}")


if __name__ == '__main__':
    """Call main program function."""
    main()
