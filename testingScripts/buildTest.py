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

# NOT NEEDED. HANDLED IN MASTER SCRIPT
# Git Status and then attempt to pull
#os.system('git status')
#print('Attempting git pull via subprocess...')
#p = subprocess.Popen("git pull", shell=True, stdout=subprocess.PIPE)
#text = p.stdout.read()
#text = text.decode('ascii')
#text = text.rstrip()
#print(text)
isTest = False
# print(str(sys.argv[1]))

# Check argument flags
if (len(sys.argv) < 2):
    print("No Arguments")
# Else check for testing flag
elif(str(sys.argv[1]) == '-t'):
    print("Test Mode: On")
    isTest = True

# Create a test build folder, get the list of executables to be generated and store them
os.chdir(home)
os.chdir("kaiju")
os.system("mkdir testFolder")
os.chdir("testFolder")
os.system("cmake ..")
listProcess = subprocess.Popen("make help | grep '\.x'", shell=True, stdout=subprocess.PIPE)
executableString = listProcess.stdout.read()
executableString = executableString.decode('ascii')
executableList = executableString.splitlines()

# Loop through each entry of the list and remove the first four characters
for index, element in enumerate(executableList):
    executableList[index] = executableList[index][4:]
print(executableList)

# Go back to scripts folder
os.chdir(home)
os.chdir("kaiju")
os.system("rm -rf testFolder")
os.chdir("testingScripts")

iteration = 1

# Read in modules.txt and load only the requested modules
file = open('modules1.txt', 'r')
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
        # Invoke cmake
        arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON;"
        # Create list of executables
        for element in executableList:
            arguments = arguments + " make " + element + ";"
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
# Invoke cmake
arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON;"
# Create list of executables
for element in executableList:
    arguments = arguments + " make " + element + ";"
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

    anyWrong = False
    missing = []

    # Check for all executables
    for element in executableList:
        isThere = os.path.isfile(element)
        if (isThere == False):
            anyWrong = True
            isPerfect = False
            missing.append(element)
    
    # If any are missing, report. Otherwise, skip
    if (not anyWrong):
        i = i + 1
        # Move back out into kaiju folder
        os.chdir(home)
        os.chdir("kaiju")
        continue
    
    else:
        myText = myText + "*Trying the following module set:*\n"
        myText = myText + ModuleList[i - 1]

        # Which executables failed?
        for element in missing:
            myText = myText + "I couldn't build " + element + "\n"
            
        # Move back out into kaiju folder
        os.chdir(home)
        os.chdir("kaiju")
        i = i + 1

# If nothing was wrong, change myText
if (isPerfect == True):
    myText = ""
    myText = "Everything built properly!"

# If not a test, send message to Slack
if (not isTest):
    # Try to send Slack message
    try:
        response = client.chat_postMessage(
            channel="#kaijudev",
            text=myText,
        )
    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

# If it is a test, just print to the command line
else:
    print(myText)
