import os
import sys
import subprocess
from os.path import expanduser
sys.path.insert(1, "./python-slackclient")
from slack import WebClient
from slack.errors import SlackApiError
import logging
logging.basicConfig(level=logging.DEBUG)
import time

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
print(slack_token)
client = WebClient(token=slack_token)

# Get CWD and set kaiju to "home"
calledFrom = os.path.dirname(os.path.abspath(__file__))
os.chdir(calledFrom)
orig = os.getcwd()
os.chdir('..')
home = os.getcwd()

isTest = False
beLoud = False

# Check argument flags
if (len(sys.argv) >= 2):
    for i in range(1,len(sys.argv)):
        if(str(sys.argv[i]) == '-t'):
            print("Test Mode: On")
            isTest = True
        elif(str(sys.argv[i]) == '-l'):
            print("Being Loud")
            beLoud = True
        else:
            print("Unrecognized argument: ", sys.argv[i])

os.chdir(home)
# If the weekly dash base folder doesn't exist, need to generate the restart
os.system('rm -r intelChecks')
os.system('mkdir intelChecks')

# Build voltron_mpi.x

# Copy the restart data

# Generate supporting files and compare to originals
grid
bcwind
rcmconfig
omni2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00

# Submit the run

# Go back to scripts folder
os.chdir(home)
os.chdir("testingScripts")

# Read in modules.txt and load only the requested modules
file = open('dashModules.txt', 'r')
modules = file.readlines()
#print(modules)

myModules = []
tempString = ""

# Create List from separate modules
for line in modules:
    myModules.append(line.strip())

for line in myModules:
	print(line)

# Create the list of arguments for the first set
arguments = "module purge; module list;"

for line in myModules:
	arguments = arguments + "module load " + line + ";"

# BUILD EXECUTABLES AND TESTS
# Move to the correct test folder
os.chdir(home)
os.chdir('intelChecks')
#arguments = arguments + "cd" + home + ";"
#arguments = arguments + "cd kaiju/unitTest1;"
# Invoke cmake
arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON -DENABLE_MPI=ON -DENABLE_MKL=ON -DCMAKE_BUILD_TYPE=DEBUG;"
# Make gamera, voltron and allTests
arguments = arguments + "make gamera_mpi; make voltron_mpi;"
print(arguments)
subprocess.call(arguments, shell=True)

os.chdir(home)
os.chdir('testingScripts')
subprocess.call("cp tinyCase.xml ../intelChecks/bin", shell=True)
subprocess.call("cp lfmD.h5 ../intelChecks/bin", shell=True)
subprocess.call("cp bcwind.h5 ../intelChecks/bin", shell=True)
subprocess.call("cp rcmconfig.h5 ../intelChecks/bin", shell=True)
subprocess.call("cp intelCheckSubmitMem.pbs ../intelChecks/bin", shell=True)
subprocess.call("cp intelCheckSubmitThread.pbs ../intelChecks/bin", shell=True)
subprocess.call("cp memSuppress.sup ../intelChecks/bin", shell=True)
subprocess.call("cp threadSuppress.sup ../intelChecks/bin", shell=True)

# SUBMIT INTEL CHECK JOBS
os.chdir(home)
os.chdir('intelChecks/bin')

# list all modules with spaces between them, to be loaded in the qsub scripts
modset = ""
for line in ModuleList[0]:
    modset = modset + line + " "

# submit memory checker
arguments = 'qsub -v MODULE_LIST="' + modset + '" intelCheckSubmitMem.pbs'
print(arguments)
submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
readString = submission.stdout.read()
readString = readString.decode('ascii')
print(readString)

firstJobNumber = readString.split('.')[0]
print(firstJobNumber)

# submit thread checker
arguments = 'qsub -v MODULE_LIST="' + modset + '" intelCheckSubmitThread.pbs'
print(arguments)
submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
readString = submission.stdout.read()
readString = readString.decode('ascii')
print(readString)

secondJobNumber = readString.split('.')[0]
print(secondJobNumber)

file = open("jobs.txt", 'w+')
file.write(firstJobNumber + "\n")
file.write(secondJobNumber)

# SUBMIT FOLLOW-UP JOB FOR SLACK POSTING
#os.chdir(home)
#os.chdir('kaiju/testingScripts')
#arguments = 'qsub intelCheckReportSubmit.pbs -W depend=after:'
#arguments = arguments + numberString
#print(arguments)

# WAIT ABOUT 1 MINUTE
#time.sleep(60)

#report = subprocess.call(arguments, shell=True, stdout=subprocess.PIPE)

# FINISHED

# If not a test, send message to Slack
#if (not isTest):
    # Try to send Slack message
#    try:
#        response = client.chat_postMessage(
#            channel="#kaijudev",
#            text=myText,
#        )
#    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
#        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
