import glob
import os
import subprocess
import sys
# import subprocess
# from os.path import expanduser
# sys.path.insert(1, "./python-slackclient")
from slack_sdk import WebClient
# from slack.errors import SlackApiError
# import logging
# logging.basicConfig(level=logging.DEBUG)
import argparse

print(f"Starting {sys.argv[0]}")

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
    print(f"args = {args}")

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
if debug:
    print(f"slack_token = {slack_token}")
client = WebClient(token=slack_token)
if debug:
    print(f"client = {client}")

# Get CWD and move to the main Kaiju folder
calledFrom = os.path.dirname(os.path.abspath(__file__))
if debug:
    print(f"calledFrom = {calledFrom}")
origCWD = os.getcwd()
if debug:
    print(f"origCWD = {origCWD}")
os.chdir(calledFrom)
os.chdir('..')
home = os.getcwd()
if debug:
    print(f"home = {home}")

# Delete everything in the unitTest folder
# os.chdir(home)
# os.system('rm -rf unitTest1')
# os.system('rm -rf unitTest2')
# os.system('mkdir unitTest1')
# os.system('mkdir unitTest2')

# Copy pFUnit stuff into Kaiju External
os.chdir(home)
os.system('cp -r /glade/p/hao/msphere/gamshare/pFUnit-4.2.0/ifort-23-mpich-derecho/FARGPARSE-1.1 external')
os.system('cp -r /glade/p/hao/msphere/gamshare/pFUnit-4.2.0/ifort-23-mpich-derecho/GFTL-1.3 external')
os.system('cp -r /glade/p/hao/msphere/gamshare/pFUnit-4.2.0/ifort-23-mpich-derecho/GFTL_SHARED-1.2 external')
os.system('cp -r /glade/p/hao/msphere/gamshare/pFUnit-4.2.0/ifort-23-mpich-derecho/PFUNIT-4.2 external')

# Go back to scripts folder
os.chdir(home)
os.chdir("testingScripts")

#***********************

# Perform unit tests with each module set..

# Get a list of build module sets.
module_list_files = glob.glob("mage_unit_test_modules/*.lst")
if debug:
    print(f"module_list_files = {module_list_files}")

# Perform a test build with each set of modules.
for module_list_file in module_list_files:
    if debug:
        print(f"module_list_file = {module_list_file}")

    # Move to the home directory for testing.
    os.chdir(home)

    # Read this module list file.
    path = os.path.join('testingScripts', module_list_file)
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    module_names = [line.rstrip() for line in lines]
    if debug:
        print(f"module_names = {module_names}")

    # Extract the name of the list.
    filename = os.path.split(module_list_file)[-1]
    if debug:
        print(f"filename = {filename}")
    module_list_name = filename.replace('.lst', '')
    if debug:
        print(f"module_list_name = {module_list_name}")

    # Make a directory for this build, and go there.
    dir_name = f"unit_{module_list_name}"
    build_directory = os.path.join(home, dir_name)
    if debug:
        print(f"build_directory = {build_directory}")
    os.system(f"rm -rf {build_directory}")
    os.mkdir(build_directory)
    os.chdir(build_directory)

    # Assemble the commands to load the listed modules.
    module_cmd = 'module --force purge; module load'
    for module_name in module_names:
        module_cmd += f" {module_name}"
    if debug:
        print(f"module_cmd = {module_cmd}")

    # Run cmake to build the Makefile.
    cmake_cmd = module_cmd + '; FC=`which ifort` FFLAGS="-qmkl" cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO -DENABLE_MPI=ON -DENABLE_MKL=ON ..'
    if debug:
        print(f"cmake_cmd = {cmake_cmd}")
    cmake_process = subprocess.Popen(cmake_cmd, shell=True)
    if debug:
        print(f"cmake_process = {cmake_process}")
    cmake_process.wait()

    # Run the build.
    make_cmd = module_cmd + '; make gamera_mpi; make voltron_mpi; make allTests'
    if debug:
        print(f"make_cmd = {make_cmd}")
    make_process = subprocess.Popen(make_cmd, shell=True)
    if debug:
        print(f"make_process = {make_process}")
    make_process.wait()

    # Copy in the test PBS scripts.
    subprocess.call("cp ../tests/genTestData.pbs ./bin", shell=True)
    subprocess.call("cp ../tests/runNonCaseTests1.pbs ./bin", shell=True)
    subprocess.call("cp ../tests/runNonCaseTests2.pbs ./bin", shell=True)
    subprocess.call("cp ../tests/runCaseTests.pbs ./bin", shell=True)
    subprocess.call("cp ../tests/voltron_mpi/bcwind.h5 ./bin", shell=True)
    subprocess.call("cp ../tests/voltron_mpi/lfmD.h5 ./bin", shell=True)
    subprocess.call("cp ../tests/voltron_mpi/rcmconfig.h5 ./bin", shell=True)

    # Go to the bin directory for testing.
    path = os.path.join(home, dir_name, 'bin')
    os.chdir(path)

    # Make a list of modules to use for this unit test.
    modset = ''
    for line in module_names:
        modset = modset + line + " "

    # Submit job to generate data needed for automated tests.
    arguments = 'qsub -A ' + account + ' -v MODULE_LIST="' + modset + '",KAIJUROOTDIR=' + home + ' genTestData.pbs'
    print(arguments)
    submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
    readString = submission.stdout.read()
    readString = readString.decode('ascii')
    print(readString)
    dataGenJob = readString.split('.')[0]
    print(dataGenJob)

    # now submit the three automated testing jobs, all contingent on
    # the data gen job.
    arguments = 'qsub -A ' + account + ' -v MODULE_LIST="' + modset + '",KAIJUROOTDIR=' + home + ' -W depend=afterok:' + dataGenJob + ' runCaseTests.pbs'
    print(arguments)
    submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
    readString = submission.stdout.read()
    readString = readString.decode('ascii')
    print(readString)
    finalString = readString
    firstJob = readString.split('.')[0]
    print(firstJob)

    arguments = 'qsub -A ' + account + ' -v MODULE_LIST="' + modset + '",KAIJUROOTDIR=' + home + ' -W depend=afterok:' + dataGenJob + ' runNonCaseTests1.pbs'
    print(arguments)
    submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
    readString = submission.stdout.read()
    readString = readString.decode('ascii')
    print(readString)
    finalString = finalString + readString
    secondJob = readString.split('.')[0]
    print (secondJob)

    arguments = 'qsub -A ' + account + ' -v MODULE_LIST="' + modset + '",KAIJUROOTDIR=' + home + ' -W depend=afterok:' + dataGenJob + ' runNonCaseTests2.pbs'
    print(arguments)
    submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
    readString = submission.stdout.read()
    readString = readString.decode('ascii')
    print(readString)
    finalString = finalString + readString
    thirdJob = readString.split('.')[0]
    print (thirddJob)

    file = open("jobs.txt", 'w+')
    file.write(firstJob + "\n")
    file.write(secondJob + "\n")
    file.write(thirdJob)

# # SUBMIT JOB THAT WILL FOLLOW UP ONCE PREVIOUS JOBS HAVE FINISHED

# # HERE IS THE STUFF FOR MOVING TO THE CORRECT FOLDER
# # Move to the correct unitTest folder
# #        arguments = arguments + "cd unitTest" + str(iteration) + "; "
#         # Invoke cmake
# #        arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON;"


# # Change directory to Kaiju repo
# #os.chdir(home)
# #os.chdir("kaiju")

# # Check build directories for good executables
# #myText = ""
# #i = 1
# #isPerfect = True
# #while i <= iteration:
# #    isGamera = False
# #    isVoltron = False
    
#     # Move to next build folder
# #    os.chdir("build" + str(i) + "/bin")

# #    print(os.getcwd())

# #    anyWrong = False
# #    missing = []

#     # Check for all executables
# #    for element in executableList:
# #        isThere = os.path.isfile(element)
# #        if (isThere == False):
# #            anyWrong = True
# #            isPerfect = False
# #            missing.append(element)
    
#     # If any are missing, report. Otherwise, skip
# #    if (not anyWrong):
# #        i = i + 1
#         # Move back out into kaiju folder
# #        os.chdir(home)
# #        os.chdir("kaiju")
# #        continue
    
# #    else:
# #        myText = myText + "*Trying the following module set:*\n"
# #        myText = myText + ModuleList[i - 1]

#         # Which executables failed?
# #        for element in missing:
# #            myText = myText + "I couldn't build " + element + "\n"
            
#         # Move back out into kaiju folder
# #        os.chdir(home)
# #        os.chdir("kaiju")
# #        i = i + 1

# # If nothing was wrong, change myText
# #if (isPerfect == True):
# #    myText = ""
# #    myText = "Everything built properly!"

# # If not a test, send message to Slack
# #if (not isTest):
#     # Try to send Slack message
# #    try:
# #        response = client.chat_postMessage(
# #            channel="#kaijudev",
# #            text=myText,
# #        )
# #    except SlackApiError as e:
#         # You will get a SlackApiError if "ok" is False
# #        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

print(f"Ending {sys.argv[0]}")
