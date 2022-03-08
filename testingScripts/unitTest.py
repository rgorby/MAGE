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

# Get CWD and set main kaiju folder to "home"
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

# Delete everything in the unitTest folder
os.chdir(home)
os.system('rm -r unitTest1')
os.system('rm -r unitTest2')
os.system('mkdir unitTest1')
os.system('mkdir unitTest2')

# Copy pFUnit stuff into Kaiju External
os.chdir(home)
os.system('cp -r /glade/p/hao/msphere/gamshare/pFUnit-4.2.0/ifort-21-MPT/FARGPARSE-1.1 external')
os.system('cp -r /glade/p/hao/msphere/gamshare/pFUnit-4.2.0/ifort-21-MPT/GFTL-1.3 external')
os.system('cp -r /glade/p/hao/msphere/gamshare/pFUnit-4.2.0/ifort-21-MPT/GFTL_SHARED-1.2 external')
os.system('cp -r /glade/p/hao/msphere/gamshare/pFUnit-4.2.0/ifort-21-MPT/PFUNIT-4.2 external')

# Go back to scripts folder
os.chdir(home)
os.chdir("testingScripts")

iteration = 1

# Read in modules.txt and load only the requested modules
file = open('unitModules.txt', 'r')
modules = file.readlines()
#print(modules)

ModuleList = []
myModules = []
tempString = ""

# Create List from separate modules
for line in modules:
    if (line.strip() == "##NEW ENVIRONMENT##"):
        # Set aside what we have already
        ModuleList.append(myModules)
        # Reset
        myModules = []
        iteration += 1
    else:
        myModules.append(line.strip())

# Add the last module set
ModuleList.append(myModules)

for setOfModules in ModuleList:
	for line in setOfModules:
		print(line)

# Create the list of arguments for the first set
arguments = "module purge; module list;"

for line in ModuleList[0]:
	arguments = arguments + "module load " + line + ";"

# BUILD EXECUTABLES AND TESTS
# Move to the correct test folder
os.chdir(home)
os.chdir('unitTest1')
#arguments = arguments + "cd" + home + ";"
#arguments = arguments + "cd kaiju/unitTest1;"
# Invoke cmake
arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON -DENABLE_MPI=ON -DENABLE_MKL=ON;"
# Make gamera, voltron and allTests
arguments = arguments + "make gamera_mpi; make voltron_mpi; make allTests;"
print(arguments)
subprocess.call(arguments, shell=True)

# Create the list of arguments for the second set
# NOT WORKING RIGHT NOW
#arguments = "module purge; module list;"

#for line in ModuleList[1]:
	#arguments = arguments + "module load " + line + ";"

# BUILD EXECUTABLES AND TESTS
# Move to the correct test folder
#os.chdir(home)
#os.chdir('kaiju/unitTest2')
# Invoke cmake
#arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON;"
# Make Gamera, Voltron, and allTests
#arguments = arguments + "make gamera_mpi; make voltron_mpi; make allTests;"
#print(arguments)
#subprocess.call(arguments, shell=True)
#ModuleList.append(tempString)
#arguments = "module purge; module list;"
#print(arguments)
#subprocess.call(arguments, shell=True)

# Submitting the test
# Go to correct directory
os.chdir(home)
os.chdir('tests')
#arguments = "qsub runNonCaseTests.pbs"
#print(arguments)
#submission = subprocess.call(arguments, shell=True, stdout=subprocess.PIPE)
#readString = submission.stdout.read()
#readString = readString.decode('ascii')
#print(submission)

#finalString = readString + "\n"
subprocess.call("cp ../tests/genTestData.pbs ../unitTest1/bin", shell=True)
subprocess.call("cp runNonCaseTests1.pbs ../unitTest1/bin", shell=True)
subprocess.call("cp runNonCaseTests2.pbs ../unitTest1/bin", shell=True)
subprocess.call("cp runCaseTests.pbs ../unitTest1/bin", shell=True)

os.chdir(home)
os.chdir('unitTest1/bin')

# list all modules with spaces between them, to be loaded in the qsub scripts
modset = ""
for line in ModuleList[0]:
    modset = modset + line + " "

# submit job to generate data needed for automated tests
arguments = 'qsub -v MODULE_LIST="' + modset + '" genTestData.pbs'
print(arguments)
submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
readString = submission.stdout.read()
readString = readString.decode('ascii')
print(readString)
dataGenJob = readString.split('.')[0]
print(dataGenJob)

# now submit the three automated testing jobs, all contingent on the data gen job
arguments = 'qsub -v MODULE_LIST="' + modset + '" -W depend=afterok:' + dataGenJob + ' runCaseTests.pbs'
print(arguments)
submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
readString = submission.stdout.read()
readString = readString.decode('ascii')
print(readString)

finalString = readString
firstJob = readString.split('.')[0]
print(firstJob)

arguments = 'qsub -v MODULE_LIST="' + modset + '" -W depend=afterok:' + dataGenJob + ' runNonCaseTests1.pbs'
print(arguments)
submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
readString = submission.stdout.read()
readString = readString.decode('ascii')
print(readString)

finalString = finalString + readString

secondJob = readString.split('.')[0]
print (secondJob)

arguments = 'qsub -v MODULE_LIST="' + modset + '" -W depend=afterok:' + dataGenJob + ' runNonCaseTests2.pbs'
print(arguments)
submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
readString = submission.stdout.read()
readString = readString.decode('ascii')
print(readString)

finalString = finalString + readString

thirdJob = readString.split('.')[0]
print (secondJob)


file = open("jobs.txt", 'w+')
file.write(firstJob + "\n")
file.write(secondJob + "\n")
file.write(thirdJob)

# SUBMIT JOB THAT WILL FOLLOW UP ONCE PREVIOUS JOBS HAVE FINISHED

# HERE IS THE STUFF FOR MOVING TO THE CORRECT FOLDER
# Move to the correct unitTest folder
#        arguments = arguments + "cd unitTest" + str(iteration) + "; "
        # Invoke cmake
#        arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON;"


# Change directory to Kaiju repo
#os.chdir(home)
#os.chdir("kaiju")

# Check build directories for good executables
#myText = ""
#i = 1
#isPerfect = True
#while i <= iteration:
#    isGamera = False
#    isVoltron = False
    
    # Move to next build folder
#    os.chdir("build" + str(i) + "/bin")

#    print(os.getcwd())

#    anyWrong = False
#    missing = []

    # Check for all executables
#    for element in executableList:
#        isThere = os.path.isfile(element)
#        if (isThere == False):
#            anyWrong = True
#            isPerfect = False
#            missing.append(element)
    
    # If any are missing, report. Otherwise, skip
#    if (not anyWrong):
#        i = i + 1
        # Move back out into kaiju folder
#        os.chdir(home)
#        os.chdir("kaiju")
#        continue
    
#    else:
#        myText = myText + "*Trying the following module set:*\n"
#        myText = myText + ModuleList[i - 1]

        # Which executables failed?
#        for element in missing:
#            myText = myText + "I couldn't build " + element + "\n"
            
        # Move back out into kaiju folder
#        os.chdir(home)
#        os.chdir("kaiju")
#        i = i + 1

# If nothing was wrong, change myText
#if (isPerfect == True):
#    myText = ""
#    myText = "Everything built properly!"

# If not a test, send message to Slack
#if (not isTest):
    # Try to send Slack message
#    try:
#        response = client.chat_postMessage(
#            channel="#kaijudev",
#            text=myText,
#        )
#    except SlackApiError as e:
        # You will get a SlackApiError if "ok" is False
#        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
