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

# Get the home directory
home = expanduser("~")

# Go to the unit Test directory
os.chdir(home)
os.chdir("kaiju/unitTest1/bin")

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
file.close()

# Take the two output files and slap them together
extension1 = "o" + job1
extension2 = "o" + job2

# Case Tests
file = open('caseTests.' + extension1, 'r')
bigFile = file.readlines()
file.close()
bigFile.append("\n\n\n")

# Non Case Tests
file = open('nonCaseTests.' + extension2, 'r')
nextFile = file.readlines()
file.close()
bigFile = bigFile + nextFile

# Scan through for some key things like "error" and "job killed"
myError = False
jobKilled = False

for line in bigFile:
    if 'error' in line:
        myError = True
    
    elif 'job killed' in line:
        jobKilled = True

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

if (myText == ""):
    exit()

# If not a test, send message to Slack
# Try to send Slack message
try:
    response = client.chat_postMessage(
   channel="#kaijudev",
   text=myText,
   )
except SlackApiError as e:
   # You will get a SlackApiError if "ok" is False
   assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

# Go to unitTests1 and delete jobs.txt
os.chdir(home)
os.chdir('kaiju/unitTest1/bin')
os.remove("jobs.txt")
