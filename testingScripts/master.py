import os
import sys
import subprocess
# sys.path.insert(1, "./python-slackclient")
from slack_sdk import WebClient
# from slack.errors import SlackApiError
# import logging
# logging.basicConfig(level=logging.DEBUG)
# import time
import argparse

print(f"Starting {sys.argv[0]}")

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d',action='store_true',default=False, help='Enables debugging output')
parser.add_argument('-t',action='store_true',default=False, help='Enables testing mode')
parser.add_argument('-l',action='store_true',default=False, help='Enables loud mode')
parser.add_argument('-a',action='store_true',default=False, help='Run all tests')
parser.add_argument('-f',action='store_true',default=False, help='Force the tests to run')
parser.add_argument('--account',type=str, default='', help='qsub account number')

args = parser.parse_args()
debug = args.d
isTest = args.t
beLoud = args.l
doAll = args.a
forceRun = args.f
account = args.account
if debug:
    print(f"args = {args}")

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
if debug:
    print(f"slack_token = {slack_token}")
client = WebClient(token=slack_token)
if debug:
    print(f"client = {client}")

# Get CWD and move to the main Kaiju folder
calledFrom = os.path.dirname(os.path.abspath(__file__))
if debug:
    print(f"calledFrom = {calledFrom}")
origCWD = os.getcwd()
if debug:
    print(f"origCWD = {origCWD}")
os.chdir(calledFrom)
os.chdir('..')
home = os.getcwd()
if debug:
    print(f"home = {home}")
# print("I am the master script. This is my current working directory: ")
# print(home)

# Delete all build folders
os.system("rm -rf build*/")
os.system('ls')

# Git Status and then attempt to pull
os.system('git status')
#print('Attempting git pull via subprocess...')
p = subprocess.Popen("git pull", shell=True, stdout=subprocess.PIPE)
text = p.stdout.read()
text = text.decode('ascii')
text = text.rstrip()
#print(text)

# get my current branch
p = subprocess.Popen("git symbolic-ref --short HEAD", shell=True, stdout=subprocess.PIPE)
gBranch = p.stdout.read()
gBranch = gBranch.decode('ascii')
gBranch = gBranch.rstrip()
#print(gBranch)

if(forceRun == False):
    # If not forced, check for update
    if (text == 'Already up to date.'):
        print("No test today. Branch " + gBranch + " is already up to date!")
        exit()

os.chdir("testingScripts")

print("I made it this far!")
print(os.path.dirname(os.path.abspath(__file__)))

subArgString = ""
if debug:
    subArgString += ' -d'
if isTest:
    subArgString = subArgString + " -t"
if beLoud:
    subArgString = subArgString + " -l"
subArgString = subArgString + " --account " + account
if debug:
    print(f"subArgString = {subArgString}")
if (doAll == True):
    buildTest = subprocess.Popen("python3 buildTest.py -d -l"+subArgString, shell = True)
    buildTest.wait()
    unitTest = subprocess.Popen("python3 unitTest.py"+subArgString, shell = True)
    unitTest.wait()
    # intelTest = subprocess.Popen("python3 intelChecks.py"+subArgString, shell=True)
    # intelTest.wait()
    # ICTest = subprocess.Popen("python3 ICtest.py"+subArgString, shell=True)
    # ICTest.wait()
    # ICReport = subprocess.Popen("python3 ICtestReport.py"+subArgString, shell=True)
    # ICReport.wait()
    # pyunitTest = subprocess.Popen("python3 pyunitTest.py"+subArgString, shell=True)
    # pyunitTest.wait()
    # weeklyDash = subprocess.Popen("python3 weeklyDash.py"+subArgString, shell=True)
    # weeklyDash.wait()

# else:
#     buildTest = subprocess.Popen("python3 buildTest.py"+subArgString, shell = True)
#     buildTest.wait()
#     unitTest = subprocess.Popen("python3 unitTest.py"+subArgString, shell = True)
#     unitTest.wait()
#     #intelTest = subprocess.Popen("python3 intelChecks.py"+subArgString, shell=True)
#     #intelTest.wait()
#     ICTest = subprocess.Popen("python3 ICtest.py"+subArgString, shell=True)
#     ICTest.wait()
#     ICReport = subprocess.Popen("python3 ICtestReport.py"+subArgString, shell=True)
#     ICReport.wait()
#     pyunitTest = subprocess.Popen("python3 pyunitTest.py"+subArgString, shell=True)
#     pyunitTest.wait()
#     #weeklyDash = subprocess.Popen("python3 weeklyDash.py"+subArgString, shell=True)
#     #weeklyDash.wait()

print(f"Ending {sys.argv[0]}")
