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

# os.chdir(home)

# # get my current branch
# p = subprocess.Popen("git symbolic-ref --short HEAD", shell=True, stdout=subprocess.PIPE)
# gBranch = p.stdout.read()
# gBranch = gBranch.decode('ascii')
# gBranch = gBranch.rstrip()
# print(gBranch)

# # If the weekly dash base folder doesn't exist, need to generate the restart
# if( not os.path.exists("dashRestarts")):
#     message = "No restart data available for weekly dash on branch " + gBranch + ". Please generate restart data and try again."
#     if(not isTest):
#         try:
#             response = client.chat_postMessage(
#                 channel="#kaijudev",
#                 text=message,
#             )
#         except SlackApiError as e:
#            # You will get a SlackApiError if "ok" is False
#            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
#     else:
#         print(message)
#     exit()

# # Build voltron_mpi.x
# os.system('rm -r weeklyDash')
# os.system('mkdir weeklyDash')
# os.chdir("weeklyDash")

# # Read in modules.txt and load only the requested modules
# file = open('../testingScripts/dashModules.txt', 'r')
# modules = file.readlines()
# #print(modules)

# myModules = []
# tempString = ""

# # Create List from separate modules
# for line in modules:
#     myModules.append(line.strip())

# for line in myModules:
# 	print(line)

# # Create the list of arguments for the first set
# arguments = "module purge; module list;"

# for line in myModules:
# 	arguments = arguments + "module load " + line + ";"

# # BUILD EXECUTABLES
# # Invoke cmake
# arguments = arguments + "cmake ../ -DENABLE_MPI=ON -DENABLE_MKL=ON -DCMAKE_BUILD_TYPE=Release;"
# # Make voltron_mpi
# arguments = arguments + "make voltron_mpi;"
# print(arguments)
# subprocess.call(arguments, shell=True)

# os.chdir("bin")

# # copy additional files from testing folder
# subprocess.call("cp ../../testingScripts/weeklyDashGo.xml .", shell=True)
# subprocess.call("cp ../../testingScripts/weeklyDashGo.pbs .", shell=True)

# # Generate new supporting files
# subprocess.call("genLFM.py -gid Q", shell=True)
# os.system("mv lfmQ.h5 NEWlfmX.h5")
# subprocess.call("cda2wind.py -t0 2016-08-09T02:00:00 -t1 2016-08-09T12:00:00 -o NEWbcwind.h5", shell=True)
# subprocess.call("genRCM.py -o NEWrcmconfig.h5", shell=True)

# # Copy the restart data
# subprocess.call("cp ../../dashRestarts/* .", shell=True)

# # Compare new supporting files to originals
# lfmp = subprocess.Popen("h5diff lfmX.h5 NEWlfmX.h5", shell=True, stdout=subprocess.PIPE)
# lfmp.wait()
# gridDiff = lfmp.stdout.read().decode('ascii').rstrip()
# if(gridDiff != "" or lfmp.returncode != 0):
#     message = "Quad grid for weekly dash has changed on branch " + gBranch + ". Case cannot be run. Please re-generate restart data, and ensure the grid change was intentional."
#     if(not isTest):
#         try:
#             response = client.chat_postMessage(
#                 channel="#kaijudev",
#                 text=message,
#             )
#         except SlackApiError as e:
#            # You will get a SlackApiError if "ok" is False
#            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
#     else:
#         print(message)
#     exit()

# bcp = subprocess.Popen("h5diff bcwind.h5 NEWbcwind.h5", shell=True, stdout=subprocess.PIPE)
# bcp.wait()
# windDiff = bcp.stdout.read().decode('ascii').rstrip()
# if(windDiff != "" or bcp.returncode != 0):
#     message = "solar wind file for weekly dash has changed on branch " + gBranch + ". Case cannot be run. Please re-generate restart data, and ensure the wind data change was intentional."
#     if(not isTest):
#         try:
#             response = client.chat_postMessage(
#                 channel="#kaijudev",
#                 text=message,
#             )
#         except SlackApiError as e:
#            # You will get a SlackApiError if "ok" is False
#            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
#     else:
#         print(message)
#     exit()

# rcmp = subprocess.Popen("h5diff rcmconfig.h5 NEWrcmconfig.h5", shell=True, stdout=subprocess.PIPE)
# rcmp.wait()
# rcmDiff = rcmp.stdout.read().decode('ascii').rstrip()
# if(rcmDiff != "" or rcmp.returncode != 0):
#     message = "rcmconfig for weekly dash has changed on branch " + gBranch + ". Case cannot be run. Please re-generate restart data, and ensure the rcmconfig change was intentional."
#     if(not isTest):
#         try:
#             response = client.chat_postMessage(
#                 channel="#kaijudev",
#                 text=message,
#             )
#         except SlackApiError as e:
#            # You will get a SlackApiError if "ok" is False
#            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
#     else:
#         print(message)
#     exit()

# # Submit the run

# # list all modules with spaces between them, to be loaded in the qsub scripts
# modset = ""
# for line in myModules:
#     modset = modset + line + " "

# # submit weekly dash
# arguments = 'qsub -A ' + account + ' -v MODULE_LIST="' + modset + '",KAIJUROOTDIR=' + home + ' weeklyDashGo.pbs'
# print(arguments)
# submission = subprocess.Popen(arguments, shell=True, stdout=subprocess.PIPE)
# readString = submission.stdout.read()
# readString = readString.decode('ascii')
# print(readString)

# firstJobNumber = readString.split('.')[0]
# print(firstJobNumber)

# file = open("jobs.txt", 'w+')
# file.write(firstJobNumber)

# message = "Run started on branch " + gBranch + " as jobid " + firstJobNumber
# print(message)

print(f"Ending {sys.argv[0]}")
