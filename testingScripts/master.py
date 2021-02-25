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

# Get the home directory
home = expanduser("~")

# Get CWD and move to the main Kaiju folder
origCWD = os.getcwd()
os.chdir("..")

# Delete all build folders
os.system("rm -rf build*/")
os.system('ls')

# Git Status and then attempt to pull
os.system('git status')
print('Attempting git pull via subprocess...')
p = subprocess.Popen("git pull", shell=True, stdout=subprocess.PIPE)
text = p.stdout.read()
text = text.decode('ascii')
text = text.rstrip()
print(text)
isTest = False

# Check argument flags
if (len(sys.argv) < 2):
    # If no arguments, check for update
    if (text == 'Already up to date.'):
        # Try to send Slack message
        try:
            response = client.chat_postMessage(
                channel="#kaijudev",
                text='No test today. It is already up to date!',
            )
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
        
        exit()
# Else check for force flag
elif(str(sys.argv[1]) == '-f'):
    print("Buuuuut you forced me to do it anyway...")
# Else check for testing flag
elif(str(sys.argv[1]) == '-t'):
    print("Test Mode: On")
    isTest = True

os.chdir("testingScripts")

if (isTest == True):
    buildTest = subprocess.Popen("python3 buildTest.py -t", shell = True)

elif (not len(sys.argv) < 2):
    buildTest = subprocess.Popen("python3 buildTest.py -f", shell = True)
    unitTest = subprocess.Popen("python3 unitTest.py", shell = True)
    intelTest = subprocess.Popen("python3 intelChecks.py", shell=True)
    ICTest = subprocess.Popen("python3 ICtest.py", shell=True)
    ICReport = subprocess.Popen("python3 ICTestReport.py", shell=True)

else:
    buildTest = subprocess.Popen("python3 buildTest.py", shell = True)
    unitTest = subprocess.Popen("python3 unitTest.py", shell = True)
    intelTest = subprocess.Popen("python3 intelChecks.py", shell=True)
    ICTest = subprocess.Popen("python3 ICtest.py", shell=True)
    ICReport = subprocess.Popen("python3 ICTestReport.py", shell=True)

