import glob
import os
import sys
import subprocess
# from os.path import expanduser
# sys.path.insert(1, "./python-slackclient")
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError
# import logging
# logging.basicConfig(level=logging.DEBUG)
import argparse

print(f"Starting {sys.argv[0]}", flush=True)

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d',action='store_true',default=False, help='Enables debugging output')
parser.add_argument('-t',action='store_true',default=False, help='Enables testing mode')
parser.add_argument('-l',action='store_true',default=False, help='Enables loud mode')
parser.add_argument('-a',action='store_true',default=False, help='Run all tests')
parser.add_argument('-f',action='store_true',default=False, help='Force the tests to run')
parser.add_argument('--account',type=str, default='', help='qsub account number')

args = parser.parse_args()
debug = args.d
isTest = args.t
beLoud = args.l
doAll = args.a
forceRun = args.f
account = args.account
if debug:
    print(f"args = {args}", flush=True)

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
if debug:
    print(f"slack_token = {slack_token}", flush=True)
client = WebClient(token=slack_token)
if debug:
    print(f"client = {client}", flush=True)

# Get CWD and move to main kaiju folder
calledFrom = os.path.dirname(os.path.abspath(__file__))
if debug:
    print(f"calledFrom = {calledFrom}", flush=True)
os.chdir(calledFrom)
origCWD = os.getcwd()
if debug:
    print(f"origCWD = {origCWD}", flush=True)
os.chdir('..')
home = os.getcwd()
if debug:
    print(f"home = {home}", flush=True)
# print("I am the build script. This is my home directory: ", flush=True)
# print(home, flush=True)

# Delete all build folders
os.system("rm -rf build*/")
os.system('ls')

# Create a test build folder, get the list of executables to be generated and store them
os.chdir(home)
os.system("mkdir testFolder")
os.chdir("testFolder")

# get my current branch
p = subprocess.Popen("git symbolic-ref --short HEAD", shell=True, stdout=subprocess.PIPE)
gBranch = p.stdout.read()
gBranch = gBranch.decode('ascii')
gBranch = gBranch.rstrip()
if debug:
    print(f"gBranch = {gBranch}", flush=True)
# print(gBranch, flush=True)

# Set up some MPI modules in order to ask for the correct set of executables
testModules = 'module --force purge'
testModules += '; module load ncarenv/23.06'
testModules += '; module load cmake/3.26.3'
testModules += '; module load craype/2.7.20'
testModules += '; module load intel/2023.0.0'
testModules += '; module load geos/3.9.1'
testModules += '; module load ncarcompilers/1.0.0'
testModules += '; module load cray-mpich/8.1.25'
testModules += '; module load hdf5-mpi/1.12.2'
testModules += '; module load mkl/2023.0.0'
if debug:
    print(f"testModules = {testModules}", flush=True)

os.system("module load cmake/3.26.3; cmake -DENABLE_MPI=ON ..")
cmd = testModules + "; make help | grep '\.x'"
if debug:
    print(f"cmd = {cmd}", flush=True)
listProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
if debug:
    print(f"listProcess = {listProcess}", flush=True)
listProcess.wait()
executableString = listProcess.stdout.read()
if debug:
    print(f"executableString = {executableString}", flush=True)
executableString = executableString.decode('ascii')
if debug:
    print(f"executableString = {executableString}", flush=True)
executableList = executableString.splitlines()
if debug:
    print(f"executableList = {executableList}", flush=True)

# Loop through each entry of the list and remove the first four characters
for index, element in enumerate(executableList):
    executableList[index] = executableList[index][4:]
if debug:
    print(f"executableList = {executableList}", flush=True)

# Go back to scripts folder
os.chdir(home)
os.system("rm -rf testFolder")
os.chdir("testingScripts")

#***********************

# Perform test builds.

# Get a list of build module sets.
module_list_files = glob.glob("build_test_modules/*.lst")
if debug:
    print(f"module_list_files = {module_list_files}", flush=True)

# Perform a test build with each set of modules.
for module_list_file in module_list_files:
    if debug:
        print(f"module_list_file = {module_list_file}", flush=True)

    # Move to the home directory for testing.
    os.chdir(home)

    # Read this module list file.
    path = os.path.join('testingScripts', module_list_file)
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    module_names = [line.rstrip() for line in lines]
    if debug:
        print(f"module_names = {module_names}", flush=True)

    # Extract the name of the list.
    filename = os.path.split(module_list_file)[-1]
    if debug:
        print(f"filename = {filename}", flush=True)
    module_list_name = filename.replace('.lst', '')
    if debug:
        print(f"module_list_name = {module_list_name}", flush=True)

    # Make a directory for this build, and go there.
    dir_name = f"build_{module_list_name}"
    build_directory = os.path.join(home, dir_name)
    if debug:
        print(f"build_directory = {build_directory}", flush=True)
    os.mkdir(build_directory)
    os.chdir(build_directory)

    # Assemble the commands to load the listed modules.
    module_cmd = 'module --force purge; module load'
    for module_name in module_names:
        module_cmd += f" {module_name}"
    if debug:
        print(f"module_cmd = {module_cmd}", flush=True)

    # Run cmake to build the Makefile.
    cmake_cmd = module_cmd + '; FC=`which ifort1` FFLAGS="-qmkl" cmake -DENABLE_MPI=ON -DENABLE_MKL=ON ..'
    if debug:
        print(f"cmake_cmd = {cmake_cmd}", flush=True)
    cmake_process = subprocess.Popen(cmake_cmd, shell=True)
    if debug:
        print(f"cmake_process = {cmake_process}", flush=True)
    cmake_process.wait()

    # Run the build.
    make_cmd = module_cmd + '; make'
    if debug:
        print(f"make_cmd = {make_cmd}", flush=True)
    make_process = subprocess.Popen(make_cmd, shell=True)
    if debug:
        print(f"make_process = {make_process}", flush=True)
    make_process.wait()

    # Create a test result message.
    message = f"IGNORE - TESTING\n*Trying the following module set: {module_names}*\n"

    # Check for all executables
    os.chdir('bin')
    missing = []
    for element in executableList:
        if os.path.isfile(element) == False:
            missing.append(element)
    if debug:
        print(f"missing = {missing}", flush=True)
    os.chdir('..')
    isPerfect = True
    if len(missing) > 0:
        isPerfect = False
        for element in missing:
            message += f"I couldn't build {element}.\n"
    else:
        message += f"Everything built properly on branch {gBranch}!"
    if debug:
        print(f"message = {message}", flush=True)

    # Don't print if it's a test, otherwise print if force, or there's an error.
    # if (not isTest and (beLoud or not isPerfect) ):
    # Try to send Slack message
    try:
        response = client.chat_postMessage(
            channel="#kaijudev",
            text=message,
        )
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    # else:
    #     # If not slack, just print to command line
    #     print(message)

print(f"Ending {sys.argv[0]}", flush=True)
