import os
import sys
import subprocess
from os.path import expanduser
sys.path.insert(1, "./python-slackclient")
from slack import WebClient
from slack.errors import SlackApiError
import logging
logging.basicConfig(level=logging.DEBUG)

# Get Slack API token
slack_token = os.environ["SLACK_BOT_TOKEN"]
print(slack_token)
client = WebClient(token=slack_token)

# Get the home directory
home = expanduser("~")

# Change directory to Kaiju repo
os.chdir(home)
os.chdir("kaiju")
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

# If there is no update, then skip the test for the day
if (text == 'Already up to date.'):
    print("Oh crap, it is already up to date!")
    if (len(sys.argv) < 2):
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
    elif(str(sys.argv[1] == 'f')):
        print("Buuuuut you forced me to do it anyway....")

# Go back to scripts folder
os.chdir(home)
os.chdir("kaiju/testingScripts")

iteration = 1

# Read in modules.txt and load only the requested modules
file = open('modules.txt', 'r')
modules = file.readlines()
arguments = "module purge; module list;"

ModuleList = []
tempString = ""

# Add modules to the list to be loaded
for line in modules:
    if (line.strip() == "##NEW ENVIRONMENT##"):
        arguments = arguments + " module list;" # List modules for this build
        # Create the build folder
        arguments = arguments + "cd " + home + "; cd kaiju; mkdir build" + str(iteration) + ";"
        # Move to the build folder
        arguments = arguments + "cd build" + str(iteration) + "; "
        # Make gamera and voltron
        arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON; make gamera_mpi.x; make voltron_mpi.x;"
        print(arguments)
        subprocess.call(arguments, shell=True)
        ModuleList.append('>' + tempString)
        tempString = ""
        iteration += 1
        arguments = "module purge; module list;"
        #print(arguments)
    else:
        tempString += line
        arguments = arguments + " module load " + line.strip() + ";" # Strip off newline characters

arguments = arguments + " module list;" # List modules for this build
# Create the build folder
arguments = arguments + "cd " + home + "; cd kaiju; mkdir build" + str(iteration) + ";"
# Move to the build folder
arguments = arguments + "cd build" + str(iteration) + "; "
# Make gamera and voltron
arguments = arguments + "cmake ..; make gamera_mpi.x; make voltron_mpi.x;"
print(arguments)
subprocess.call(arguments, shell=True)
ModuleList.append(tempString)
arguments = "module purge; module list;"
#print(arguments)
subprocess.call(arguments, shell=True)

# Change directory to Kaiju repo
os.chdir(home)
os.chdir("kaiju")

# Check build directories for good executables
myText = ""
i = 1
isPerfect = True
while i <= iteration:
    isGamera = False
    isVoltron = False
    
    # Move to next build folder
    os.chdir("build" + str(i) + "/bin")

    print(os.getcwd())
    
    # Check for Gamera and Voltron
    isGamera = os.path.isfile('gamera_mpi.x')
    isVoltron = os.path.isfile('voltron_mpi.x')
    
    # Check if both worked. If so, skip.
    if (isGamera and isVoltron):
        i = i + 1
        continue

    else:
        isPerfect = False

    myText = myText + "*Trying the following module set:*\n"
    myText = myText + ModuleList[i - 1]

    if (not isGamera):
        myText = myText + "Whoops! I couldn't build Gamera using that module set.\n"

    if (not isVoltron):
        myText = myText + "Whoops! I couldn't build Voltron using that module set.\n"

    # Move back out into kaiju folder
    os.chdir(home)
    os.chdir("kaiju")

    i = i + 1

# If nothing was wrong, change myText
myText = ""
myText = "Everything built properly!"

# Try to send Slack message
try:
    response = client.chat_postMessage(
        channel="#kaijudev",
        text=myText,
    )
except SlackApiError as e:
    # You will get a SlackApiError if "ok" is False
    assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
