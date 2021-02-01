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

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
print(slack_token)
client = WebClient(token=slack_token, timeout=120)

# Get the home directory
home = expanduser("~")

# Go to IntelChecks folder
os.chdir(home)
os.chdir('kaiju/intelChecks/bin')

# Check for jobs.txt
jobsExists = path.exists("jobs.txt")

# If not, end. Otherwise, continue
if (not jobsExists):
    exit()

# Go to IntelChecks folder
os.chdir(home)
os.chdir('kaiju/intelChecks/bin/r000mi3')
process = subprocess.call("ls", shell=True)

# Find output files

# Memory
time.sleep(60)

# Open them and record results
file = open("inspxe-cl.txt", 'r')
results = file.readlines()

# Turn list into big string
str1 = ""

for element in results:
	str1 = str1 + element

# Set it as a message for the slack bot
myText = "Intel Memory Test:\nLocated at: "
myText = myText + os.getcwd() + "\n"
myText = myText + str1

print(myText)

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

# Threading
os.chdir(home)
os.chdir('kaiju/intelChecks/bin/r001ti3')
process = subprocess.call("ls", shell=True)
time.sleep(60)

# Open file and record results
file = open("inspxe-cl.txt", 'r')
results = file.readlines()

# Turn list into big string
str1 = ""

for element in results:
	str1 = str1 + element

# Set it as a message for the slack bot
myText = "Intel Thread Checker:\nLocated at: "
myText = myText + os.getcwd() + "\n"
myText = myText + str1

print(myText)

# Try to send Slack message
try:
    response = client.chat_postMessage(
    	channel="#kaijudev",
    	text=myText,
    )
except SlackApiError as e:
   # You will get a SlackApiError if "ok" is False
   assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'


# Go to IntelChecks and delete jobs.txt
os.chdir(home)
os.chdir('kaiju/intelChecks/bin')
os.remove("jobs.txt")
