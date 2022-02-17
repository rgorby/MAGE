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

# get my current branch
p = subprocess.Popen("git symbolic-ref --short HEAD", shell=True, stdout=subprocess.PIPE)
gBranch = p.stdout.read()
gBranch = gBranch.decode('ascii')
gBranch = gBranch.rstrip()
print(gBranch)

isTest = False
doAll = False
forceRun = False
beLoud = False

# Check argument flags
if (len(sys.argv) >= 2):
    for i in range(1,len(sys.argv)):
        if(str(sys.argv[i]) == '-f'):
            print("Buuuuut you forced me to do it anyway...")
            forceRun = True
        elif(str(sys.argv[i]) == '-t'):
            print("Test Mode: On")
            isTest = True
        elif(str(sys.argv[i]) == '-a'):
            print("Running All Tests")
            doAll = True
        elif(str(sys.argv[i]) == '-l'):
            print("Being Loud")
            beLoud = True
        else:
            print("Unrecognized argument: ", sys.argv[i])

if(forceRun == False):
    # If not forced, check for update
    if (text == 'Already up to date.'):
        print("No test today. Branch " + gBranch + " is already up to date!")
        exit()

os.chdir("testingScripts")

print("I made it this far!")
print(os.path.dirname(os.path.abspath(__file__)))

subArgString = ""
if isTest:
    subArgString = subArgString + " -t"
if beLoud:
    subArgString = subArgString + " -l"

if (doAll == True):
    buildTest = subprocess.Popen("python3 buildTest.py"+subArgString, shell = True)
    buildTest.wait()
    unitTest = subprocess.Popen("python3 unitTest.py"+subArgString, shell = True)
    unitTest.wait()
    intelTest = subprocess.Popen("python3 intelChecks.py"+subArgString, shell=True)
    intelTest.wait()
    ICTest = subprocess.Popen("python3 ICtest.py"+subArgString, shell=True)
    ICTest.wait()
    ICReport = subprocess.Popen("python3 ICtestReport.py"+subArgString, shell=True)
    ICReport.wait()
    pyunitTest = subprocess.Popen("python3 pyunitTest.py"+subArgString, shell=True)
    pyunitTest.wait()

else:
    buildTest = subprocess.Popen("python3 buildTest.py"+subArgString, shell = True)
    buildTest.wait()
    #unitTest = subprocess.Popen("python3 unitTest.py"+subArgString, shell = True)
    #unitTest.wait()
    #intelTest = subprocess.Popen("python3 intelChecks.py"+subArgString, shell=True)
    #intelTest.wait()
    ICTest = subprocess.Popen("python3 ICtest.py"+subArgString, shell=True)
    ICTest.wait()
    ICReport = subprocess.Popen("python3 ICtestReport.py"+subArgString, shell=True)
    ICReport.wait()
    pyunitTest = subprocess.Popen("python3 pyunitTest.py"+subArgString, shell=True)
    pyunitTest.wait()

