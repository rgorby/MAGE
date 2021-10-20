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

# Get CWD and set Kaiju as "home"
calledFrom = os.path.dirname(os.path.abspath(__file__))
os.chdir(calledFrom)
orig = os.getcwd()
os.chdir('..')
home = os.getcwd()

# Check argument flags
if (len(sys.argv) < 2):
    print("No Arguments")
# Else check for testing flag
elif(str(sys.argv[1]) == '-t'):
    print("Test Mode: On")
    isTest = True

# Go back to scripts folder
os.chdir(home)
os.chdir("testingScripts")

iteration = 1

# Read in modules.txt and load only the requested modules
file = open('modules1.txt', 'r')
modules = file.readlines()
arguments = "module purge; module list;"

ModuleList = []
tempString = ""

# Get list of IC files to try (GAMERA ONLY FOR NOW)
os.chdir(home)
os.chdir('src/gamera/ICs')
directory = os.getcwd()
directory = directory + "/"
listFiles = subprocess.Popen("ls", shell=True, stdout=subprocess.PIPE)
fileList = listFiles.stdout.read()
fileList = fileList.decode('ascii')
print(fileList)
print(directory)

# Add modules to the list to be loaded
os.chdir(home)
os.system('mkdir ICBuilds')
os.chdir('ICBuilds')

# Split file list into an actual list
fileList = fileList.splitlines()

# Ignore old IC files in the "deprecated" directory
fileList.remove('deprecated')
print(fileList)


for line in modules:
    if (line.strip() == "##NEW ENVIRONMENT##"):
        for element in fileList:
            arguments = arguments + " module list;" # List modules for this build
            # Create the build folder
            arguments = arguments + "cd " + home + "; cd ICBuilds/;"
            arguments = arguments + "mkdir gamera" + element.strip() + str(iteration) + ";"
            # Move to the build folder
            arguments = arguments + "cd gamera" + element.strip() + str(iteration) + "; "
            # Invoke cmake
            arguments = arguments + "cmake ../../ -DALLOW_INVALID_COMPILERS=ON -DGAMIC:FILEPATH="
            # Add the correct path for GAMIC
            arguments = arguments + directory + element.strip() + "; "
            # Create Gamera
            arguments = arguments + " make gamera.x; "
            print(arguments)
            subprocess.call(arguments, shell=True)
        iteration += 1
        arguments = "module purge; module list;"
        #print(arguments)
    else:
        tempString += line
        arguments = arguments + " module load " + line.strip() + ";" # Strip off newline characters

for element in fileList:
    arguments = arguments + " module list;" # List modules for this build
    # Create the build folder
    arguments = arguments + "cd " + home + "; cd ICBuilds;"
    arguments = arguments + "mkdir gamera" + element.strip() + str(iteration) + ";"
    # Move to the build folder
    arguments = arguments + "cd gamera" + element.strip() + str(iteration) + "; "
    # Invoke cmake
    arguments = arguments + "cmake ../../ -DALLOW_INVALID_COMPILERS=ON -DGAMIC:FILEPATH="
    # Add the correct path for GAMIC
    arguments = arguments + directory + element.strip() + "; "
    # Create Gamera
    arguments = arguments + " make gamera.x; "
    print(arguments)
    subprocess.call(arguments, shell=True)
iteration += 1
arguments = "module purge; module list;"
#print(arguments)

#subprocess.call(arguments, shell=True)
#
## Change directory to Kaiju repo
#os.chdir(home)
#os.chdir("kaiju")
#
## Check build directories for good executables
#myText = ""
#i = 1
#isPerfect = True
#while i <= iteration:
#    isGamera = False
#    isVoltron = False
#    
#    # Move to next build folder
#    os.chdir("build" + str(i) + "/bin")
#
#    print(os.getcwd())
#
#    anyWrong = False
#    missing = []
#
#    # Check for all executables
#    for element in executableList:
#        isThere = os.path.isfile(element)
#        if (isThere == False):
#            anyWrong = True
#            isPerfect = False
#            missing.append(element)
#    
#    # If any are missing, report. Otherwise, skip
#    if (not anyWrong):
#        i = i + 1
#        # Move back out into kaiju folder
#        os.chdir(home)
#        os.chdir("kaiju")
#        continue
#    
#    else:
#        myText = myText + "*Trying the following module set:*\n"
#        myText = myText + ModuleList[i - 1]
#
#        # Which executables failed?
#        for element in missing:
#            myText = myText + "I couldn't build " + element + "\n"
#            
#        # Move back out into kaiju folder
#        os.chdir(home)
#        os.chdir("kaiju")
#        i = i + 1
#
## If nothing was wrong, change myText
#if (isPerfect == True):
#    myText = ""
#    myText = "Everything built properly!"
#
## If not a test, send message to Slack
#if (not isTest):
#    # Try to send Slack message
#    try:
#        response = client.chat_postMessage(
#            channel="#kaijudev",
#            text=myText,
#        )
#    except SlackApiError as e:
#        # You will get a SlackApiError if "ok" is False
#        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
#
## If it is a test, just print to the command line
#else:
#    print(myText)
