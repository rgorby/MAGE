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
from os import path
import re

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
print(slack_token)
client = WebClient(token=slack_token, timeout=120)

# Get CWD and set main kaiju folder to "home"
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

# Go back to scripts folder
os.chdir(home)
os.chdir("testingScripts")

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

print(arguments)
subprocess.call(arguments, shell=True)

# Go to IntelChecks folder
os.chdir(home)
os.chdir('intelChecks/bin')

# Check for jobs.txt
jobsExists = path.exists("jobs.txt")

# If not, end. Otherwise, continue
if (not jobsExists):
    print("Nothing to Test.\n")
    exit()

# Read in the jobs.txt file to get the job numbers
file = open("jobs.txt", 'r')
job1 = file.readline()
job1 = job1.strip()
job2 = file.readline()
job2 = job2.strip()
file.close()

# Take the output files and slap them together
jobFile1 = "memCheck.o" + job1
jobFile2 = "threadCheck.o" + job2

if (not path.exists(jobFile1) or not path.exists(jobFile2)):
    print("One of the jobs isn't complete yet.\n")
    exit()

# Go to IntelChecks folder
os.chdir(home)
os.chdir('intelChecks/bin')

# Setup regex for counting problems found
problemPattern = "(\d+) new problem\(s\) found"

# Memory
memErrs = False
memErrsFile = "combinedMemErrs.txt"
for root, dirs, files in os.walk("."):
    for d in dirs:
        if "memResults" in d:
            try:
                memOut = subprocess.check_output(["inspxe-cl","-report summary","-result-dir " + d,"-s-f memSuppress.sup"], \
                    stderr=subprocess.STDOUT,universal_newlines=True)
            except subprocess.CalledProcessError as memProcErr:
                # we need to handle non-zero error code
                memOut = memProcErr.output
            problemMatch = re.search(problemPattern, memOut)
            if(not problemMatch or int(problemMatch.group(1)) > 0):
                memErrs = True
                try:
                    memOut = subprocess.check_output(["inspxe-cl","-report problems","-result-dir " + d,"-s-f memSuppress.sup","-report-all"], \
                        stderr=subprocess.STDOUT,universal_newlines=True)
                except subprocess.CalledProcessError as memProcErr:
                    # we need to handle non-zero error code
                    memOut = memProcErr.output
                with open(memErrsFile, "a") as memFile:
                    memFile.write(memOut)
                    memFile.write("\n")

if(not isTest and (beLoud or memErrs)):
    if(memErrs):
        try:
            response = client.files_upload(
                file=memErrsFile,
                initial_comment='Memory Access Problems:\n\n',
                channels="#kaijudev",
                )
            assert response['ok']
            slack_file = response['file']
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
        try:
            response = client.chat_postMessage(
                channel="#kaijudev",
                text="Memory Access Problems Detected\nPlease check the errors sent above",
            )
        except SlackApiError as e:
           # You will get a SlackApiError if "ok" is False
           assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    else:
        try:
            response = client.chat_postMessage(
                channel="#kaijudev",
                text="No Memory Access Problems Detected",
            )
        except SlackApiError as e:
           # You will get a SlackApiError if "ok" is False
           assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
else:
    if(memErrs):
        print("*** Memory Access Problems Detected ***\nPlease check the results in " + memErrsFile + "\n")
    else:
        print("No Memory Access Problems Detected")

# Thread
threadErrs = False
threadErrsFile = "combinedThreadErrs.txt"
for root, dirs, files in os.walk("."):
    for d in dirs:
        if "threadResults" in d:
            try:
                threadOut = subprocess.check_output(["inspxe-cl","-report summary","-result-dir " + d,"-s-f threadSuppress.sup"], \
                    stderr=subprocess.STDOUT,universal_newlines=True)
            except subprocess.CalledProcessError as threadProcErr:
                # we need to handle non-zero error code
                threadOut = threadProcErr.output
            problemMatch = re.search(problemPattern, threadOut)
            if(not problemMatch or int(problemMatch.group(1)) > 0):
                threadErrs = True
                try:
                    threadOut = subprocess.check_output(["inspxe-cl","-report problems","-result-dir " + d,"-s-f threadSuppress.sup","-report-all"], \
                        stderr=subprocess.STDOUT,universal_newlines=True)
                except subprocess.CalledProcessError as threadProcErr:
                    # we need to handle non-zero error code
                    threadOut = threadProcErr.output
                with open(threadErrsFile, "a") as threadFile:
                    threadFile.write(threadOut)
                    threadFile.write("\n")

if(not isTest and (beLoud or threadErrs)):
    if(threadErrs):
        try:
            response = client.files_upload(
                file=threadErrsFile,
                initial_comment='Threading Data Race Problems:\n\n',
                channels="#kaijudev",
                )
            assert response['ok']
            slack_file = response['file']
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
        try:
            response = client.chat_postMessage(
                channel="#kaijudev",
                text="Threading Data Race Problems Detected\nPlease check the errors sent above",
            )
        except SlackApiError as e:
           # You will get a SlackApiError if "ok" is False
           assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    else:
        try:
            response = client.chat_postMessage(
                channel="#kaijudev",
                text="No Threading Data Race Problems Detected",
            )
        except SlackApiError as e:
           # You will get a SlackApiError if "ok" is False
           assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
else:
    if(threadErrs):
        print("*** Threading Data Race Problems Detected ***\nPlease check the results in " + threadErrsFile + "\n")
    else:
        print("No Threading Data Race Problems Detected")

# Go to IntelChecks and delete jobs.txt
os.chdir(home)
os.chdir('intelChecks/bin')
os.remove("jobs.txt")

