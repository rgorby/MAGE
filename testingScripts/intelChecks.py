# import os
import sys
# import subprocess
# from os.path import expanduser
# sys.path.insert(1, "./python-slackclient")
# from slack import WebClient
# from slack.errors import SlackApiError
# import logging
# logging.basicConfig(level=logging.DEBUG)
# import time
# import argparse

print(f"Starting {sys.argv[0]}")

# # read arguments
# parser = argparse.ArgumentParser()
# parser.add_argument('-t',action='store_true',default=False, help='Enables testing mode')
# parser.add_argument('-l',action='store_true',default=False, help='Enables loud mode')
# parser.add_argument('-a',action='store_true',default=False, help='Run all tests')
# parser.add_argument('-f',action='store_true',default=False, help='Force the tests to run')
# parser.add_argument('--account',type=str, default='', help='qsub account number')

# args = parser.parse_args()
# isTest = args.t
# beLoud = args.l
# doAll = args.a
# forceRun = args.f
# account = args.account

# # Get Slack API token
# slack_token = os.environ["SLACK_BOT_TOKEN"]
# print(slack_token)
# client = WebClient(token=slack_token)

# # Get CWD and set kaiju to "home"
# calledFrom = os.path.dirname(os.path.abspath(__file__))
# os.chdir(calledFrom)
# orig = os.getcwd()
# os.chdir('..')
# home = os.getcwd()

# # Delete everything in the unitTest folder
# os.chdir(home)
# os.system('rm -r intelChecks')
# os.system('mkdir intelChecks')


# # Go back to scripts folder
# os.chdir(home)
# os.chdir("testingScripts")

# iteration = 1

# # Read in modules.txt and load only the requested modules
# file = open('intelModules.txt', 'r')
# modules = file.readlines()
# #print(modules)

# ModuleList = []
# myModules = []
# tempString = ""

# # Create List from separate modules
# for line in modules:
#     if (line.strip() == "##NEW ENVIRONMENT##"):
#         # Set aside what we have already
#         ModuleList.append(myModules)
#         # Reset
#         myModules = []
#         iteration += 1
#     else:
#         myModules.append(line.strip())

# # Add the last module set
# ModuleList.append(myModules)

# for setOfModules in ModuleList:
# 	for line in setOfModules:
# 		print(line)

# # Create the list of arguments for the first set
# arguments = "module purge; module list;"

# for line in ModuleList[0]:
# 	arguments = arguments + "module load " + line + ";"

# # BUILD EXECUTABLES AND TESTS
# # Move to the correct test folder
# os.chdir(home)
# os.chdir('intelChecks')
# #arguments = arguments + "cd" + home + ";"
# #arguments = arguments + "cd kaiju/unitTest1;"
# # Invoke cmake
# arguments = arguments + "cmake ../ -DALLOW_INVALID_COMPILERS=ON -DDISABLE_DEBUG_BOUNDS_CHECKS=ON -DENABLE_MPI=ON -DENABLE_MKL=ON -DCMAKE_BUILD_TYPE=DEBUG;"
# # Make gamera, voltron and allTests
# arguments = arguments + "make gamera_mpi; make voltron_mpi;"
# print(arguments)
# subprocess.call(arguments, shell=True)

# os.chdir(home)
# os.chdir('testingScripts')
# subprocess.call("cp tinyCase.xml ../intelChecks/bin", shell=True)
# subprocess.call("cp lfmD.h5 ../intelChecks/bin", shell=True)
# subprocess.call("cp bcwind.h5 ../intelChecks/bin", shell=True)
# subprocess.call("cp rcmconfig.h5 ../intelChecks/bin", shell=True)
# subprocess.call("cp intelCheckSubmitMem.pbs ../intelChecks/bin", shell=True)
# subprocess.call("cp intelCheckSubmitThread.pbs ../intelChecks/bin", shell=True)
# subprocess.call("cp memSuppress.sup ../intelChecks/bin", shell=True)
# subprocess.call("cp threadSuppress.sup ../intelChecks/bin", shell=True)

# # SUBMIT INTEL CHECK JOBS
# os.chdir(home)
# os.chdir('intelChecks/bin')

# # list all modules with spaces between them, to be loaded in the qsub scripts
# modset = ""
# for line in ModuleList[0]:
#     modset = modset + line + " "

# # submit memory checker
# arguments = 'qsub -A ' + account + ' -v MODULE_LIST="' + modset + '",KAIJUROOTDIR=' + home + ' intelCheckSubmitMem.pbs'
# print(arguments)
# submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
# readString = submission.stdout.read()
# readString = readString.decode('ascii')
# print(readString)

# firstJobNumber = readString.split('.')[0]
# print(firstJobNumber)

# # submit thread checker
# arguments = 'qsub -A ' + account + ' -v MODULE_LIST="' + modset + '",KAIJUROOTDIR=' + home + ' intelCheckSubmitThread.pbs'
# print(arguments)
# submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
# readString = submission.stdout.read()
# readString = readString.decode('ascii')
# print(readString)

# secondJobNumber = readString.split('.')[0]
# print(secondJobNumber)

# file = open("jobs.txt", 'w+')
# file.write(firstJobNumber + "\n")
# file.write(secondJobNumber)

# # SUBMIT FOLLOW-UP JOB FOR SLACK POSTING
# #os.chdir(home)
# #os.chdir('kaiju/testingScripts')
# #arguments = 'qsub intelCheckReportSubmit.pbs -W depend=after:'
# #arguments = arguments + numberString
# #print(arguments)

# # WAIT ABOUT 1 MINUTE
# #time.sleep(60)

# #report = subprocess.call(arguments, shell=True, stdout=subprocess.PIPE)

# # FINISHED

# # If not a test, send message to Slack
# #if (not isTest):
#     # Try to send Slack message
# #    try:
# #        response = client.chat_postMessage(
# #            channel="#kaijudev",
# #            text=myText,
# #        )
# #    except SlackApiError as e:
#         # You will get a SlackApiError if "ok" is False
# #        assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

print(f"Ending {sys.argv[0]}")
