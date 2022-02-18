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
wikiPath = ""

# Check argument flags
if (len(sys.argv) >= 2):
    i=1
    while(i < len(sys.argv)):
        if(str(sys.argv[i]) == '-t'):
            print("Test Mode: On")
            isTest = True
        elif(str(sys.argv[i]) == '-l'):
            print("Being Loud")
            beLoud = True
        elif(str(sys.argv[i]) == "-w"):
            wikiPath = sys.argv[i+1]
            i=i+1
        else:
            print("Unrecognized argument: ", sys.argv[i])
        i=i+1

print(wikiPath)
if(not os.path.exists(wikiPath)):
    print("Wiki folder does not exist at " + wikiPath)
    exit()

# ensure wikiPath ends in a slash, for later string manipulation
if(wikiPath[-1] != '/'):
    wikiPath = wikiPath + '/'
print(wikiPath)

# Go back to scripts folder
os.chdir(home)
os.chdir("testingScripts")

# get my current branch
p = subprocess.Popen("git symbolic-ref --short HEAD", shell=True, stdout=subprocess.PIPE)
gBranch = p.stdout.read()
gBranch = gBranch.decode('ascii')
gBranch = gBranch.rstrip()
print(gBranch)

if(gBranch != "master" and gBranch != "development"):
    print("storm dash only reported for master and development branches, but this branch is " + gBranch)
    exit()

# Go to weekly dash folder
os.chdir(home)

if(not os.path.exists("weeklyDash")):
    print("dash results folder does not exist")
    exit()

os.chdir('weeklyDash')
os.chdir('bin')

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

# Find out if the run is done
jobFile1 = "wDashGo.o" + job1

if (not path.exists(jobFile1)):
    print("The dash job isn't complete yet.\n")
    exit()

# move to wiki folder and ensure it is up to date
os.chdir(wikiPath)
p = subprocess.Popen("git pull", shell=True, stdout=subprocess.PIPE)
text = p.stdout.read().decode('ascii').rstrip()
print(text)
if('not a git repository' in text):
    print("wiki path is not a git repository")
    exit()

# If this is development branch, move the old data
os.chdir("weeklyDash")
if(gBranch == "development"):
    os.system("mv development_ut.txt      development_ut_old.txt")
    os.system("mv development_rt.txt      development_rt_old.txt")
    os.system("mv development_qk_msph.png development_qk_msph_old.png")
    os.system("mv development_qk_rcm.png  development_qk_rcm_old.png")
    os.system("mv development_qk_mix.png  development_qk_mix_old.png")

# Move back to simulation data
os.chdir(home)
os.chdir('weeklyDash')
os.chdir('bin')

# Get performance data
p = subprocess.Popen('sed --quiet "s/^ \\+UT \\+= \\+2016-08-09 \\+\\([0-9:]\\+\\).*$/\\1/p" weeklyDashGo.out', shell=True, stdout=subprocess.PIPE)
utData = p.stdout.read().decode('ascii')
p = subprocess.Popen('sed --quiet "s/^ \\+Running @ *\\([0-9]\\+\\.\\?[0-9]*\\)% of real-time.*$/\\1/p" weeklyDashGo.out', shell=True, stdout=subprocess.PIPE)
rtData = p.stdout.read().decode('ascii')

#There is always one extra line of UT data in the front, strip it, then write the data to files
utData = ''.join(utData.splitlines(keepends=True)[1:])

utFile = open(wikiPath + "weeklyDash/" + gBranch + '_ut.txt', "w")
utFile.write(utData)
utFile.close()

rtFile = open(wikiPath + "weeklyDash/" + gBranch + '_rt.txt', "w")
rtFile.write(rtData)
rtFile.close()

# Make quick-look plots
subprocess.call('msphpic.py', shell=True)
os.system('cp qkpic.png ' + wikiPath + "weeklyDash/" + gBranch + '_qk_msph.png')
subprocess.call('mixpic.py', shell=True)
os.system('cp remix_n.png ' + wikiPath + "weeklyDash/" + gBranch + '_qk_mix.png')
subprocess.call('rcmpic.py', shell=True)
os.system('cp qkrcmpic.png ' + wikiPath + "weeklyDash/" + gBranch + '_qk_rcm.png')

# Move back to wiki
os.chdir(wikiPath)
os.chdir("weeklyDash")

# Make new performance plot. Could make this prettier, gnuplot for now
subprocess.call('gnuplot ' + home + '/testingScripts/perfPlot.plg', shell=True)

# Combine quick looks into larger images
subprocess.call("convert master_qk_msph.png -gravity NorthWest -pointsize 60 -annotate +0+0 'master' mm.png", shell=True)
subprocess.call("convert master_qk_mix.png  -gravity NorthWest -pointsize 80 -annotate +0+0 'master' mx.png", shell=True)
subprocess.call("convert master_qk_rcm.png  -gravity NorthWest -pointsize 80 -annotate +0+0 'master' mr.png", shell=True)
subprocess.call("convert development_qk_msph_old.png -gravity NorthWest -pointsize 60 -annotate +0+0 'development prior' dmo.png", shell=True)
subprocess.call("convert development_qk_mix_old.png  -gravity NorthWest -pointsize 80 -annotate +0+0 'development prior' dxo.png", shell=True)
subprocess.call("convert development_qk_rcm_old.png  -gravity NorthWest -pointsize 80 -annotate +0+0 'development prior' dro.png", shell=True)
subprocess.call("convert development_qk_msph.png -gravity NorthWest -pointsize 60 -annotate +0+0 'development latest' dm.png", shell=True)
subprocess.call("convert development_qk_mix.png  -gravity NorthWest -pointsize 80 -annotate +0+0 'development latest' dx.png", shell=True)
subprocess.call("convert development_qk_rcm.png  -gravity NorthWest -pointsize 80 -annotate +0+0 'development latest' dr.png", shell=True)
subprocess.call('convert mm.png dmo.png dm.png +append combined_qk_msph.png', shell=True)
subprocess.call('convert mx.png dxo.png dx.png +append combined_qk_mix.png',  shell=True)
subprocess.call('convert mr.png dro.png dr.png +append combined_qk_rcm.png',  shell=True)
if(os.path.exists("mm.png")) : os.remove("mm.png")
if(os.path.exists("mx.png")) : os.remove("mx.png")
if(os.path.exists("mr.png")) : os.remove("mr.png")
if(os.path.exists("dmo.png")) : os.remove("dmo.png")
if(os.path.exists("dxo.png")) : os.remove("dxo.png")
if(os.path.exists("dro.png")) : os.remove("dro.png")
if(os.path.exists("dm.png")) : os.remove("dm.png")
if(os.path.exists("dx.png")) : os.remove("dx.png")
if(os.path.exists("dr.png")) : os.remove("dr.png")

# Push the data to the wiki
p = subprocess.Popen('git commit -a -m "New weekly dash data for branch ' + gBranch + '"', shell=True, stdout=subprocess.PIPE)
text = p.stdout.read().decode('ascii').rstrip()
print(text)
p = subprocess.Popen("git push", shell=True, stdout=subprocess.PIPE)
text = p.stdout.read().decode('ascii').rstrip()
print(text)

# Announce results
if(not isTest and beLoud):
    try:
        response = client.chat_postMessage(
            channel="#kaijudev",
            text="Weekly simulation complete on branch " + gBranch + ". Updated results follow."
        )
    except SlackApiError as e:
       # You will get a SlackApiError if "ok" is False
       assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

    try:
        response = client.files_upload(
            file='perfPlots.png',
            initial_comment='Real-Time Performance\n\n',
            channels="#kaijudev",
            )
        assert response['ok']
        slack_file = response['file']
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

    try:
        response = client.files_upload(
            file='combined_qk_msph.png',
            initial_comment='Quick-Looks Magnetosphere\n\n',
            channels="#kaijudev",
            )
        assert response['ok']
        slack_file = response['file']
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

    try:
        response = client.files_upload(
            file='combined_qk_mix.png',
            initial_comment='Quick-Looks Remix\n\n',
            channels="#kaijudev",
            )
        assert response['ok']
        slack_file = response['file']
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

    try:
        response = client.files_upload(
            file='combined_qk_rcm.png',
            initial_comment='Quick-Looks RCM\n\n',
            channels="#kaijudev",
            )
        assert response['ok']
        slack_file = response['file']
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

else:
    print("weekly run completed successfully on branch " + gBranch)

# Delete jobs.txt
os.chdir(home)
os.chdir('weeklyDash')
os.chdir('bin')
os.remove("jobs.txt")

