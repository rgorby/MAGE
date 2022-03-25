import os
import sys
import subprocess
from os.path import expanduser
sys.path.insert(1, "./python-slackclient")
from slack import WebClient
from slack.errors import SlackApiError
import logging
logging.basicConfig(level=logging.DEBUG)
import matplotlib as mpl
mpl.use('Agg')
import h5py
import numpy as np
import datetime
import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv
from astropy.time import Time

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
debugMode = False
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
        elif(str(sys.argv[i]) == '-d'):
            print("Debugging Mode: On")
            debugMode = True
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
    if(debugMode):
        print("Changing branch to development for testing mode")
        gBranch = "development"
    else:
        print("Exitting")
        exit()

# Go to weekly dash folder
os.chdir(home)

if(not os.path.exists("weeklyDash")):
    print("dash results folder does not exist")
    exit()

os.chdir('weeklyDash')
os.chdir('bin')

# Check for jobs.txt
jobsExists = os.path.exists("jobs.txt")

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

if (not os.path.exists(jobFile1)):
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
os.chdir("weeklyDash")

# Get performance data
p = subprocess.Popen('sed --quiet "s/^ \\+UT \\+= \\+\\([0-9-]\\+ [0-9:]\\+\\).*$/\\1/p" ' + home + '/weeklyDash/bin/weeklyDashGo.out', shell=True, stdout=subprocess.PIPE)
utData = p.stdout.read().decode('ascii')
p = subprocess.Popen('sed --quiet "s/^ \\+Running @ *\\([0-9]\\+\\.\\?[0-9]*\\)% of real-time.*$/\\1/p" ' + home + '/weeklyDash/bin/weeklyDashGo.out', shell=True, stdout=subprocess.PIPE)
rtData = p.stdout.read().decode('ascii')

# Split UT and RT data into lists
utData = utData.splitlines()
rtData = rtData.splitlines()

# There is always one extra line of UT data in the front, strip it
utData = utData[1:]

# Make sure the lists are of equal length now (console output not always reliable)
if(len(utData) > len(rtData)):
    utData = utData[:len(rtData)]
elif(len(rtData) > len(utData)):
    rtData = rtData[:len(utData)]
# Convert the rt data to floats
rtData_f = [float(rtData[n]) for n in range(len(rtData))]

# Get simulation data from voltron output file
LW = 0.75 # LineWidth
alpha = 0.25 # Transparency
gCol = "slategrey" # GridColor

fVolt = home + "/weeklyDash/bin/msphere.volt.h5"
print(fVolt)
#Get Dst and CPCP plot info
print("Scraping")
nSteps,sIds = kh5.cntSteps(fVolt)
symh  = kh5.getTs(fVolt,sIds,"SymH")
MJD   = kh5.getTs(fVolt,sIds,"MJD")
BSDst = kh5.getTs(fVolt,sIds,"BSDst")
nCPCP = kh5.getTs(fVolt,sIds,"cpcpN")
sCPCP = kh5.getTs(fVolt,sIds,"cpcpS")
UT = Time(MJD,format='mjd').isot
ut_datetime = [datetime.datetime.strptime(UT[n],'%Y-%m-%dT%H:%M:%S.%f') for n in range(len(UT))] # needed to plot symh

# load old data from an h5 file
masterUT = None
masterRT = None
masterUTsim = None
masterDST = None
masterCPCPn = None
masterCPCPs = None
devpriorUT = None
devpriorRT = None
devpriorUTsim = None
devpriorDST = None
devpriorCPCPn = None
devpriorCPCPs = None
devcurrentUT = None
devcurrentRT = None
devcurrentUTsim = None
devcurrentDST = None
devcurrentCPCPn = None
devcurrentCPCPs = None
if(os.path.exists('previousData.h5')):
    data_object = h5py.File('previousData.h5', 'r')
    if 'masterUT' in data_object:
        masterUT    = [x.decode('utf-8') for x in data_object['masterUT']]
        masterRT    = data_object['masterRT'].value
        masterUTsim = [x.decode('utf-8') for x in data_object['masterUTsim']]
        masterDST   = data_object['masterDST'].value
        masterCPCPn = data_object['masterCPCPn'].value
        masterCPCPs = data_object['masterCPCPs'].value
    if 'devpriorUT' in data_object:
        devpriorUT    = [x.decode('utf-8') for x in data_object['devpriorUT']]
        devpriorRT    = data_object['devpriorRT'].value
        devpriorUTsim = [x.decode('utf-8') for x in data_object['devpriorUTsim']]
        devpriorDST   = data_object['devpriorDST'].value
        devpriorCPCPn = data_object['devpriorCPCPn'].value
        devpriorCPCPs = data_object['devpriorCPCPs'].value
    if 'devcurrentUT' in data_object:
        devcurrentUT    = [x.decode('utf-8') for x in data_object['devcurrentUT']]
        devcurrentRT    = data_object['devcurrentRT'].value
        devcurrentUTsim = [x.decode('utf-8') for x in data_object['devcurrentUTsim']]
        devcurrentDST   = data_object['devcurrentDST'].value
        devcurrentCPCPn = data_object['devcurrentCPCPn'].value
        devcurrentCPCPs = data_object['devcurrentCPCPs'].value
    data_object.close()

# update appropriate data with new data
if(gBranch == 'master'):
    masterUT = utData
    masterRT = rtData_f
    masterUTsim = UT
    masterDST = BSDst
    masterCPCPn = nCPCP
    masterCPCPs = sCPCP
elif(gBranch == 'development'):
    devpriorUT = devcurrentUT
    devpriorRT = devcurrentRT
    devpriorUTsim = devcurrentUTsim
    devpriorDST = devcurrentDST
    devpriorCPCPn = devcurrentCPCPn
    devpriorCPCPs = devcurrentCPCPs
    devcurrentUT = utData
    devcurrentRT = rtData_f
    devcurrentUTsim = UT
    devcurrentDST = BSDst
    devcurrentCPCPn = nCPCP
    devcurrentCPCPs = sCPCP

# Convert date strings into date-time objects
if masterUT is not None:
    masterUTdt = [datetime.datetime.strptime(masterUT[n],'%Y-%m-%d %H:%M:%S') for n in range(len(masterUT))]
    masterUTsimdt = [datetime.datetime.strptime(masterUTsim[n],'%Y-%m-%dT%H:%M:%S.%f') for n in range(len(masterUTsim))]
if devpriorUT is not None:
    devpriorUTdt = [datetime.datetime.strptime(devpriorUT[n],'%Y-%m-%d %H:%M:%S') for n in range(len(devpriorUT))]
    devpriorUTsimdt = [datetime.datetime.strptime(devpriorUTsim[n],'%Y-%m-%dT%H:%M:%S.%f') for n in range(len(devpriorUTsim))]
if devcurrentUT is not None:
    devcurrentUTdt = [datetime.datetime.strptime(devcurrentUT[n],'%Y-%m-%d %H:%M:%S') for n in range(len(devcurrentUT))]
    devcurrentUTsimdt = [datetime.datetime.strptime(devcurrentUTsim[n],'%Y-%m-%dT%H:%M:%S.%f') for n in range(len(devcurrentUTsim))]

# Make Real-Time Performance Plot
fSz = (14,7)
fig = mpl.pyplot.figure(figsize=fSz)
gs = mpl.gridspec.GridSpec(1,1,hspace=0.05,wspace=0.05)
ax=fig.add_subplot(gs[0,0])

if masterRT is not None:
    ax.plot(masterUTdt,masterRT,label="master",linewidth=LW)
if devpriorRT is not None:
    ax.plot(devpriorUTdt,devpriorRT,label="dev prior",linewidth=LW)
if devcurrentRT is not None:
    ax.plot(devcurrentUTdt,devcurrentRT,label="dev current",linewidth=LW)

ax.legend(loc='lower right',fontsize="small")

ax.minorticks_on()
ax.xaxis_date()
xfmt = mpl.dates.DateFormatter('%H:%M \n%Y-%m-%d')
ax.set_ylabel("Percent of Real-Time [%]")
ax.xaxis.set_major_formatter(xfmt)
ax.xaxis.set_minor_locator(mpl.dates.HourLocator())
# mpl.pyplot.grid(True)
# ax.xaxis.grid(True,which='major',linewidth=LW  ,alpha=alpha,color=gCol)
# ax.xaxis.grid(True,which='minor',linewidth=LW/4,alpha=alpha,color=gCol)
# ax.yaxis.grid(True,which='major',linewidth=LW  ,alpha=alpha,color=gCol)
ax.set_title("Real-Time Performance")
fOut = "perfPlots.png"
kv.savePic(fOut)
mpl.pyplot.close('all')


# Make DST Plot
fSz = (14,7)
fig = mpl.pyplot.figure(figsize=fSz)
gs = mpl.gridspec.GridSpec(1,1,hspace=0.05,wspace=0.05)
ax=fig.add_subplot(gs[0,0])

ax.plot(ut_datetime,symh,label="SYM-H",linewidth=2*LW)
if masterDST is not None:
    ax.plot(masterUTsimdt,masterDST,label="master",linewidth=LW)
if devpriorDST is not None:
    ax.plot(devpriorUTsimdt,devpriorDST,label="dev prior",linewidth=LW)
if devcurrentDST is not None:
    ax.plot(devcurrentUTsimdt,devcurrentDST,label="dev current",linewidth=LW)

ax.legend(loc='upper right',fontsize="small")

ax.minorticks_on()
ax.xaxis_date()
xfmt = mpl.dates.DateFormatter('%H:%M \n%Y-%m-%d')
ax.set_ylabel("Dst [nT]")
ax.xaxis.set_major_formatter(xfmt)
ax.xaxis.set_minor_locator(mpl.dates.HourLocator())
# mpl.pyplot.grid(True)
# ax.xaxis.grid(True,which='major',linewidth=LW  ,alpha=alpha,color=gCol)
# ax.xaxis.grid(True,which='minor',linewidth=LW/4,alpha=alpha,color=gCol)
# ax.yaxis.grid(True,which='major',linewidth=LW  ,alpha=alpha,color=gCol)
ax.set_title("BSDst")
fOut = "dstPlots.png"
kv.savePic(fOut)
mpl.pyplot.close('all')
    
# Make CPCP Plot
fSz = (14,7)
fig = mpl.pyplot.figure(figsize=fSz)
gs = mpl.gridspec.GridSpec(1,1,hspace=0.05,wspace=0.05)
ax=fig.add_subplot(gs[0,0])

if masterCPCPn is not None:
    ax.plot(masterUTsimdt,masterCPCPn,color='orange',linestyle='dotted',label="master-North",linewidth=LW)
    ax.plot(masterUTsimdt,masterCPCPs,color='blue',linestyle='dotted',label="master-South",linewidth=LW)
if devpriorCPCPn is not None:
    ax.plot(devpriorUTsimdt,devpriorCPCPn,color='orange',linestyle='dashed',label="dev prior-North",linewidth=LW)
    ax.plot(devpriorUTsimdt,devpriorCPCPs,color='blue',linestyle='dashed',label="dev prior-South",linewidth=LW)
if devcurrentCPCPn is not None:
    ax.plot(devcurrentUTsimdt,devcurrentCPCPn,color='orange',linestyle='solid',label="dev current-North",linewidth=LW)
    ax.plot(devcurrentUTsimdt,devcurrentCPCPs,color='blue',linestyle='solid',label="dev current-South",linewidth=LW)

ax.legend(loc='upper right',fontsize="small")

ax.minorticks_on()
ax.xaxis_date()
xfmt = mpl.dates.DateFormatter('%H:%M \n%Y-%m-%d')
ax.set_ylabel("Dst [nT]")
ax.xaxis.set_major_formatter(xfmt)
ax.xaxis.set_minor_locator(mpl.dates.HourLocator())
# mpl.pyplot.grid(True)
# ax.xaxis.grid(True,which='major',linewidth=LW  ,alpha=alpha,color=gCol)
# ax.xaxis.grid(True,which='minor',linewidth=LW/4,alpha=alpha,color=gCol)
# ax.yaxis.grid(True,which='major',linewidth=LW  ,alpha=alpha,color=gCol)
ax.set_title("CPCP")
fOut = "cpcpPlots.png"
kv.savePic(fOut)
mpl.pyplot.close('all')

# Save the new data as json
with h5py.File('previousData.h5', 'w') as data_object:
    if masterUT is not None:
        data_object.create_dataset('masterUT',    data=[x.encode('utf-8') for x in masterUT])
        data_object.create_dataset('masterRT',    data=masterRT)
        data_object.create_dataset('masterUTsim', data=[x.encode('utf-8') for x in masterUTsim])
        data_object.create_dataset('masterDST',   data=masterDST)
        data_object.create_dataset('masterCPCPn', data=masterCPCPn)
        data_object.create_dataset('masterCPCPs', data=masterCPCPs)
    if devpriorUT is not None:
        data_object.create_dataset('devpriorUT',    data=[x.encode('utf-8') for x in devpriorUT])
        data_object.create_dataset('devpriorRT',    data=devpriorRT)
        data_object.create_dataset('devpriorUTsim', data=[x.encode('utf-8') for x in devpriorUTsim])
        data_object.create_dataset('devpriorDST',   data=devpriorDST)
        data_object.create_dataset('devpriorCPCPn', data=devpriorCPCPn)
        data_object.create_dataset('devpriorCPCPs', data=devpriorCPCPs)
    if devcurrentUT is not None:
        data_object.create_dataset('devcurrentUT',    data=[x.encode('utf-8') for x in devcurrentUT])
        data_object.create_dataset('devcurrentRT',    data=devcurrentRT)
        data_object.create_dataset('devcurrentUTsim', data=[x.encode('utf-8') for x in devcurrentUTsim])
        data_object.create_dataset('devcurrentDST',   data=devcurrentDST)
        data_object.create_dataset('devcurrentCPCPn', data=devcurrentCPCPn)
        data_object.create_dataset('devcurrentCPCPs', data=devcurrentCPCPs)

# If I'm on development, copy latest quick look plots over old ones
if(gBranch == "development"):
    os.system("mv development_qk_msph.png development_qk_msph_old.png")
    os.system("mv development_qk_rcm.png  development_qk_rcm_old.png")
    os.system("mv development_qk_mix.png  development_qk_mix_old.png")

# Move back to simulation data
os.chdir(home)
os.chdir('weeklyDash')
os.chdir('bin')

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
subprocess.call('convert mm.png dmo.png dm.png -append combined_qk_msph.png', shell=True)
subprocess.call('convert mx.png dxo.png dx.png -append combined_qk_mix.png',  shell=True)
subprocess.call('convert mr.png dro.png dr.png -append combined_qk_rcm.png',  shell=True)
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
if(not debugMode):
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
            text="Weekly results complete on branch " + gBranch + ". Latest comparative results attached as replies to this message.\nOr up-to-date results can be viewed on the wiki at https://bitbucket.org/aplkaiju/kaiju/wiki/weeklyDash/dashStatus"
        )
    except SlackApiError as e:
       # You will get a SlackApiError if "ok" is False
       assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'

    if response["ok"]:
        parent_ts = response["ts"]
        try:
            response = client.chat_postMessage(
                channel="#kaijudev",
                text="This was a 4x4x1 (IxJxK) decomposed Quad Resolution Run using 8 nodes for Gamera, 1 for Voltron, and 3 Squish Helper nodes (12 nodes total)",
                thread_ts=parent_ts,
            )
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    
        try:
            response = client.files_upload(
                file='perfPlots.png',
                initial_comment='Real-Time Performance\n\n',
                channels="#kaijudev",
                thread_ts=parent_ts,
                )
            assert response['ok']
            slack_file = response['file']
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    
        try:
            response = client.files_upload(
                file='dstPlots.png',
                initial_comment='DST Plots\n\n',
                channels="#kaijudev",
                thread_ts=parent_ts,
                )
            assert response['ok']
            slack_file = response['file']
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    
        try:
            response = client.files_upload(
                file='cpcpPlots.png',
                initial_comment='CPCP Plots\n\n',
                channels="#kaijudev",
                thread_ts=parent_ts,
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
                thread_ts=parent_ts,
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
                thread_ts=parent_ts,
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
                thread_ts=parent_ts,
                )
            assert response['ok']
            slack_file = response['file']
        except SlackApiError as e:
            # You will get a SlackApiError if "ok" is False
            assert e.response["error"]  # str like 'invalid_auth', 'channel_not_found'
    else:
        print("Failed to post parent message to slack, could not attach replies either.")

else:
    print("weekly run completed successfully on branch " + gBranch)

# Delete jobs.txt
os.chdir(home)
os.chdir('weeklyDash')
os.chdir('bin')
os.remove("jobs.txt")

