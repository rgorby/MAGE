import os
import sys
import subprocess
from os.path import expanduser
sys.path.insert(1, "./python-slackclient")
from slack import WebClient
from slack.errors import SlackApiError
import logging
logging.basicConfig(level=logging.DEBUG)
from os import path
import argparse

# read arguments
parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
parser.add_argument('-t',action='store_true',default=False, help='Enables testing mode')
parser.add_argument('-l',action='store_true',default=False, help='Enables loud mode')
parser.add_argument('-a',action='store_true',default=False, help='Run all tests')
parser.add_argument('-f',action='store_true',default=False, help='Force the tests to run')
parser.add_argument('--account',type=str, default='', help='qsub account number')

args = parser.parse_args()
isTest = args.t
beLoud = args.l
doAll = args.a
forceRun = args.f
account = args.account

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
print(slack_token)
client = WebClient(token=slack_token)

# Get CWD and set main kaiju folder to "home"
calledFrom = os.path.dirname(os.path.abspath(__file__))
os.chdir(calledFrom)
orig = os.getcwd()
os.chdir('..')
home = os.getcwd()

# Go to the unit Test directory
os.chdir(home)
os.chdir("unitTest1/bin")

# get my current branch
p = subprocess.Popen("git symbolic-ref --short HEAD", shell=True, stdout=subprocess.PIPE)
gBranch = p.stdout.read()
gBranch = gBranch.decode('ascii')
gBranch = gBranch.rstrip()
print(gBranch)

# Check for jobs.txt
jobsExists = path.exists("jobs.txt")

# If not, end. Otherwise, continue
if (not jobsExists):
    print("Nothing to Test.\n")
    exit()

# Read in the jobs.txt file to get the job numbers
file = open('jobs.txt', 'r')
job1 = file.readline()
job1 = job1.strip()
job2 = file.readline()
job2 = job2.strip()
job3 = file.readline()
job3 = job3.strip()
file.close()

# Take the output files and slap them together
jobFile1 = "caseTests.o" + job1
jobFile2 = "nonCaseTests1.o" + job2
jobFile3 = "nonCaseTests2.o" + job3

if (not path.exists(jobFile1) or not path.exists(jobFile2) or not path.exists(jobFile3)):
    print("One of the jobs isn't complete yet.\n")
    exit()

# Case Tests
file = open(jobFile1, 'r')
bigFile = file.readlines()
file.close()
bigFile.append("\n\n\n")

# Non Case Tests 1
file = open(jobFile2, 'r')
nextFile = file.readlines()
file.close()
bigFile = bigFile + nextFile
bigFile.append("\n\n\n")

#Non Case Tests 2
file = open(jobFile3, 'r')
finalFile = file.readlines()
file.close()
bigFile = bigFile + finalFile

# Scan through for some key things like "error" and "job killed"
myError = False
jobKilled = False
okFailure = False
okCount = 0

for line in bigFile:
    if 'OK' in line:
        okCount += 1

    if 'error' in line:
        myError = True
    
    elif 'job killed' in line:
        jobKilled = True

if okCount is not 8:
    okFailure = True

# delete jobs.txt
os.remove("jobs.txt")

if not okFailure and not myError and not jobKilled:
    if not isTest and beLoud:
        try:
            response = client.chat_postMessage(
           channel="#kaijudev",
           text="Fortran Unit Tests Passed on Branch " + gBranch,
           )
        except SlackApiError as e:
           # You will get a SlackApiError if "ok" is False
           assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    else:
        print("Fortran Unit Tests Passed on Branch " + gBranch)
    exit()

# Write to a file
file = open('Results.txt', 'w+')
file.writelines(bigFile)
file.close()

# Post the file to slack
if not isTest:
    try:
        response = client.files_upload(
            file='Results.txt',
            initial_comment='Unit Test Results:\n\n',
            channels="#kaijudev",
            )
    
        assert response['ok']
        slack_file = response['file']
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

# If there were any issues, also send a message saying "Uh oh" or something
myText = ""

if (myError):
    myText = "There were errors!\n"

if (jobKilled):
    myText = myText + "The job was killed early!\n"

if (okFailure):
    myText = myText + "There were not the correct amount of OKs!\n"

myText = myText + "On Branch " + gBranch + "\n"

# If not a test, send message to Slack
# Try to send Slack message
if not isTest:
    try:
        response = client.chat_postMessage(
       channel="#kaijudev",
       text=myText + "https://tenor.com/Yx4u.gif",
       )
    except SlackApiError as e:
       # You will get a SlackApiError if "ok" is False
       assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

