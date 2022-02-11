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

# Get CWD and set Kaiju to "home"
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

# Go to ICBuilds folder
os.chdir(home)
os.chdir('ICBuilds')

# get my current branch
p = subprocess.Popen("git symbolic-ref --short HEAD", shell=True, stdout=subprocess.PIPE)
gBranch = p.stdout.read()
gBranch = gBranch.decode('ascii')
gBranch = gBranch.rstrip()
print(gBranch)

# Get List of subdirectories
directoryList = os.listdir(os.getcwd())

print(directoryList)

# Check for Gamera

incorrectList = []

# Loop through the list
for element in directoryList:
    # Check for gamera.x
    exists = path.exists(element + "/bin/gamera.x")
    # If not, record an error
    if (not exists):
        print("Found a bad one!")
        temp = element[0:6]
        temp = temp + " couldn't be made with module set "
        moduleSet = element[-1]
        temp = temp + moduleSet
        temp = temp + " and "
        ICfile = element[6:]
        ICfile = ICfile[:-1]
        temp = temp + ICfile + "."
        incorrectList.append(temp)
        temp = ""

print(incorrectList)

# Get the module list
os.chdir(home)
os.chdir("testingScripts")

file = open('modules1.txt', 'r')
modules = file.readlines()
moduleList = []
tempString = ""

for line in modules:
    if (line.strip() == "##NEW ENVIRONMENT##"):
        moduleList.append(tempString)
        tempString = ""
    else:
        tempString += line

moduleList.append(tempString)

for moduleSet in moduleList:
    print(moduleSet)

myText = ""

# Create a long string to send to slack
for element in incorrectList:
    myText = myText + element + "\n"

myText = myText + "\nModule Set 1:\n" + moduleList[0]
myText = myText + "\nModule Set 2:\n" + moduleList[1]
myText = myText + "\nModule Set 3:\n" + moduleList[2] + "\n"

hadErrors = True
if not incorrectList:
    hadErrors = False
    myText = "All ICs built OK on branch " + gBranch

# Try to send Slack message
if(not isTest and (beLoud or hadErrors)):
    try:
        response = client.chat_postMessage(
        	channel="#kaijudev",
        	text=myText,
        )
    except SlackApiError as e:
       # You will get a SlackApiError if "ok" is False
       assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
else:
    print(myText)

