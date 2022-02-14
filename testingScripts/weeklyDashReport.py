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

# get my current branch
p = subprocess.Popen("git symbolic-ref --short HEAD", shell=True, stdout=subprocess.PIPE)
gBranch = p.stdout.read()
gBranch = gBranch.decode('ascii')
gBranch = gBranch.rstrip()
print(gBranch)

# Go to weekly dash folder
os.chdir(home)
os.chdir('weeklyDash')

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
file.close()

# Take the output files and slap them together
jobFile1 = "wDashGo.o" + job1

if (not path.exists(jobFile1)):
    print("The dash job isn't complete yet.\n")
    exit()

# Announce results
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
                text="Threading Data Race Problems Detected on Branch " + gBranch + "\nPlease check the errors sent above",
            )
        except SlackApiError as e:
           # You will get a SlackApiError if "ok" is False
           assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    else:
        try:
            response = client.chat_postMessage(
                channel="#kaijudev",
                text="No Threading Data Race Problems Detected on Branch " + gBranch,
            )
        except SlackApiError as e:
           # You will get a SlackApiError if "ok" is False
           assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
else:
    if(threadErrs):
        print("*** Threading Data Race Problems Detected on Branch " + gBranch + "***\nPlease check the results in " + threadErrsFile + "\n")
    else:
        print("No Threading Data Race Problems Detected on Branch " + gBranch)

# Delete jobs.txt
os.remove("jobs.txt")

