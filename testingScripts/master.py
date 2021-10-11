import os
import sys
import subprocess
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

# Get CWD and move to the main Kaiju folder
calledFrom = os.path.dirname(os.path.abspath(__file__))
origCWD = os.getcwd()
os.chdir(calledFrom)
os.chdir('..')
home = os.getcwd()
print("I am the master script. This is my current working directory: ")
print(home)

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
doAll = False
forceRun = False

# Check argument flags
if (len(sys.argv) >= 2):
    if(str(sys.argv[1]) == '-f'):
        print("Buuuuut you forced me to do it anyway...")
        forceRun = True
    elif(str(sys.argv[1]) == '-t'):
        print("Test Mode: On")
        isTest = True
    elif(str(sys.argv[1]) == '-a'):
        print("Running All Tests")
        doAll = True

if(len(sys.argv) < 2 or doAll == True):
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

os.chdir("testingScripts")

print("I made it this far!")
print(os.path.dirname(os.path.abspath(__file__)))

if (isTest == True):
    buildTest = subprocess.Popen("python3 buildTest.py -t", shell = True)

elif (doAll == True):
    buildTest = subprocess.Popen("python3 buildTest.py", shell = True)
    buildTest.wait()
    unitTest = subprocess.Popen("python3 unitTest.py", shell = True)
    unitTest.wait()
    intelTest = subprocess.Popen("python3 intelChecks.py", shell=True)
    intelTest.wait()
    ICTest = subprocess.Popen("python3 ICtest.py", shell=True)
    ICTest.wait()
    ICReport = subprocess.Popen("python3 ICtestReport.py", shell=True)
    ICReport.wait()
    pyunitTest = subprocess.Popen("python3 pyunitTest.py", shell=True)
    pyunitTest.wait()

elif (forceRun == True):
    buildTest = subprocess.Popen("python3 buildTest.py -f", shell = True)
    buildTest.wait()
    unitTest = subprocess.Popen("python3 unitTest.py", shell = True)
    buildTest.wait()
    intelTest = subprocess.Popen("python3 intelChecks.py", shell=True)
    intelTest.wait()
    ICTest = subprocess.Popen("python3 ICtest.py", shell=True)
    ICTest.wait()
    ICReport = subprocess.Popen("python3 ICtestReport.py", shell=True)
    ICReport.wait()
    pyunitTest = subprocess.Popen("python3 pyunitTest.py", shell=True)
    pyunitTest.wait()

else:
    buildTest = subprocess.Popen("python3 buildTest.py", shell = True)
    buildTest.wait()
    #unitTest = subprocess.Popen("python3 unitTest.py", shell = True)
    #unitTest.wait()
    #intelTest = subprocess.Popen("python3 intelChecks.py", shell=True)
    #intelTest.wait()
    ICTest = subprocess.Popen("python3 ICtest.py", shell=True)
    ICTest.wait()
    ICReport = subprocess.Popen("python3 ICtestReport.py", shell=True)
    ICReport.wait()
    pyunitTest = subprocess.Popen("python3 pyunitTest.py", shell=True)
    pyunitTest.wait()

