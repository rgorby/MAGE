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

# Get CWD and set kaiju to "home"
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

os.chdir(home)
# If the weekly dash base folder doesn't exist, need to generate the restart
if( not path.exists("dashRestarts")):
    os.system('mkdir dashRestarts')

# Build voltron_mpi.x
os.system('rm -r weeklyDash')
os.chdir("weeklyDash")

# Read in modules.txt and load only the requested modules
file = open('../testingScripts/dashModules.txt', 'r')
modules = file.readlines()
#print(modules)

myModules = []
tempString = ""

# Create List from separate modules
for line in modules:
    myModules.append(line.strip())

for line in myModules:
	print(line)

# Create the list of arguments for the first set
arguments = "module purge; module list;"

for line in myModules:
	arguments = arguments + "module load " + line + ";"

# BUILD EXECUTABLES
# Invoke cmake
arguments = arguments + "cmake ../ -DENABLE_MPI=ON -DENABLE_MKL=ON -DCMAKE_BUILD_TYPE=Release;"
# Make voltron_mpi
arguments = arguments + "make voltron_mpi;"
print(arguments)
subprocess.call(arguments, shell=True)

os.chdir("bin")

# copy supporting files from testing folder
subprocess.call("cp ../../testingScripts/weeklyDashGo.xml .", shell=True)
subprocess.call("cp ../../testingScripts/weeklyDashGo.pbs .", shell=True)

# Copy the restart data
subprocess.call("cp ../../dashRestarts/msphere* .", shell=True)

# Generate supporting files and compare to originals
subprocess.call("python genGridLFM -o newGrid.h5", shell=True)
grid
bcwind
rcmconfig
omni2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00

# Submit the run

# list all modules with spaces between them, to be loaded in the qsub scripts
modset = ""
for line in myModules:
    modset = modset + line + " "

# submit weekly dash
arguments = 'qsub -v MODULE_LIST="' + modset + '" weeklyDashGo.pbs'
print(arguments)
submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
readString = submission.stdout.read()
readString = readString.decode('ascii')
print(readString)

firstJobNumber = readString.split('.')[0]
print(firstJobNumber)

file = open("jobs.txt", 'w+')
file.write(firstJobNumber)

