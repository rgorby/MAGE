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

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
print(slack_token)
client = WebClient(token=slack_token)

# Get CWD and set main kaiju folder to "home"
orig = os.getcwd()
os.chdir('..')
home = os.getcwd()

# Go to the unit Test directory
os.chdir(home)
os.chdir("unitTest1/bin")

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

# Take the two output files and slap them together
extension1 = "o" + job1
extension2 = "o" + job2
extension3 = "o" + job3

# Case Tests
file = open('caseTests.' + extension1, 'r')
bigFile = file.readlines()
file.close()
bigFile.append("\n\n\n")

# Non Case Tests 1
file = open('nonCaseTests1.' + extension2, 'r')
nextFile = file.readlines()
file.close()
bigFile = bigFile + nextFile
bigFile.append("\n\n\n")

#Non Case Tests 2
file = open('nonCaseTests2.' + extension3, 'r')
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

# Write to a file
file = open('Results.txt', 'w+')
file.writelines(bigFile)
file.close()

# Post the file to slack
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

if (myText == ""):
    exit()

# If not a test, send message to Slack
# Try to send Slack message
try:
    response = client.chat_postMessage(
   channel="#kaijudev",
   text=myText + "https://tenor.com/Yx4u.gif",
   )
except SlackApiError as e:
   # You will get a SlackApiError if "ok" is False
   assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

# Go to unitTests1 and delete jobs.txt
os.chdir(home)
os.chdir('unitTest1/bin')
os.remove("jobs.txt")
