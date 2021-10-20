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

# Delete everything in the unitTest folder
os.chdir(home)
os.system('rm -r intelChecks')
os.system('mkdir intelChecks')


# Go back to scripts folder
os.chdir(home)
os.chdir("testingScripts")

iteration = 1

# Read in modules.txt and load only the requested modules
file = open('intelModules.txt', 'r')
modules = file.readlines()
#print(modules)

ModuleList = []
myModules = []
tempString = ""

# Create List from separate modules
for line in modules:
    if (line.strip() == "##NEW ENVIRONMENT##"):
        # Set aside what we have already
        ModuleList.append(myModules)
        # Reset
        myModules = []
        iteration += 1
    else:
        myModules.append(line.strip())

# Add the last module set
ModuleList.append(myModules)

for setOfModules in ModuleList:
	for line in setOfModules:
		print(line)

# Create the list of arguments for the first set
arguments = "module purge; module list;"

for line in ModuleList[0]:
	arguments = arguments + "module load " + line + ";"

# BUILD EXECUTABLES AND TESTS
# Move to the correct test folder
os.chdir(home)
os.chdir('intelChecks')
#arguments = arguments + "cd" + home + ";"
#arguments = arguments + "cd kaiju/unitTest1;"
# Invoke cmake
arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON -DENABLE_MPI=ON;"
# Make gamera, voltron and allTests
arguments = arguments + "make gamera_mpi; make voltron_mpi;"
print(arguments)
subprocess.call(arguments, shell=True)

os.chdir(home)
os.chdir('testingScripts')
subprocess.call("cp tinyCase.xml ../intelChecks/bin", shell=True)
subprocess.call("cp intelCheckSubmit.pbs ../intelChecks/bin", shell=True)

# SUBMIT INTEL CHECK JOBS
os.chdir(home)
os.chdir('intelChecks/bin')
arguments = 'qsub -V intelCheckSubmit.pbs'
print(arguments)
submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
readString = submission.stdout.read()
readString = readString.decode('ascii')
print(readString)

jobNumber = readString.split('.')[0]
print(jobNumber)

numberString = str(jobNumber)

file = open("jobs.txt", 'w+')
file.write(jobNumber)

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
