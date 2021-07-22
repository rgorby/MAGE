#!/usr/bin/env python

#standard python
import sys
import os
import glob
import argparse
import subprocess
import time
from argparse import RawTextHelpFormatter

#numpy
import numpy as np

#Kaipy and related
from   astropy.time import Time
import h5py
import kaipy.kaiH5 as kaiH5
import kaipy.kaiViz as kv
import kaipy.kaiTools as kaiTools
import spacepy.datamodel as dm



if __name__ == '__main__':
	MainS = """Checks the run for satellites with data avilable and then
	sets up PBS job scripts for running interpolation in parallel."""

	parser = argparse.ArgumentParser(description=MainS,
	formatter_class=RawTextHelpFormatter)
	parser.add_argument('-id',type=str,metavar='runid',default='msphere',
		help='RunID of data (default: %(default)s)')
	parser.add_argument('-path',type=str,metavar='path',default='.',
		help='Path to directory containing files (default: %(default)s)')
	parser.add_argument('-satId',type=str,metavar='Satellite Id',
		default=None,help='Name of Satellite to compare')
	parser.add_argument('-cmd',type=str,metavar='command',default=None,
		help='Full path to sctrack.x command')
	parser.add_argument('-acct',type=str,metavar='acct',default=None,
		help='Account number to use in pbs script')
	parser.add_argument('--keep',action='store_true',
		help='Keep intermediate files')

	args = parser.parse_args()

	fdir = args.path
	ftag = args.id
	cmd = args.cmd
	scRequested = args.satId
	keep = args.keep
	acct = args.acct

	if acct == None:
		print('Must input a valid account to charge, use -acct flag')
		quit()

	if fdir == '.':
		fdir = os.getcwd()

	if None == cmd:
		my_env = os.environ.copy()
		cmd = os.path.join(os.getenv('KAIJUDIR'),'build','bin','sctrack.x')
	if not (os.path.isfile(cmd) and os.access(cmd, os.X_OK)):
		print(cmd,'either not found or not executable')
		sys.exit()

	scIds = kaiTools.getScIds()

	#print(cmddir)
	#print('Extracting GAMERA data along',scId, 'trajectory')
	#cmd = "/Users/wiltbemj/src/kaiju/build/bin/sctrack.x"
	#Pull the timestep information from the magnetosphere files
	(fname,isMPI,Ri,Rj,Rk) = kaiTools.getRunInfo(fdir,ftag)
	nsteps,sIds=kaiH5.cntSteps(fname)
	gamMJD=kaiH5.getTs(fname,sIds,aID='MJD')
	gamT=kaiH5.getTs(fname,sIds,aID='time')
	gamUT = kaiTools.MJD2UT(gamMJD)
	## Deal with startup by using the first non zero time as the inital
	## MJD
	loc = np.argwhere(gamT > 0.0)[0][0]
	t0 = gamUT[loc]
	t1 = gamUT[-1]
	deltaT = np.round(gamT[loc+1]-gamT[loc])
	mjdFileStart = gamMJD[loc]
	secFileStart = gamT[loc]
	numPer = 10
	numSegments=int(np.floor(((t1-t0).total_seconds()/deltaT)/numPer))

	if None == scRequested:
		scToDo = scIds
	else:
		scToDo = []
		scToDo.append(scRequested)

	for scId in scToDo:
		print('Getting spacecraft data for', scId)
		status,data = kaiTools.getSatData(scIds[scId],
			t0.strftime("%Y-%m-%dT%H:%M:%SZ"),
			t1.strftime("%Y-%m-%dT%H:%M:%SZ"),deltaT)

		if status['http']['status_code'] != 200:
			print('No data available for', scId)
		else:
			(scTrackName,xmlFileName) = kaiTools.createInputFiles(data,
				scIds[scId],scId,mjdFileStart,secFileStart,
				fdir,ftag,numSegments)
			lockCmdName = os.path.join(fdir,'makeLock.sh')
			fLock = open(lockCmdName,'w')
			fLock.write("#!/bin/bash\ntouch $1")
			fLock.close()
			os.chmod(lockCmdName,0o775)
			pbsName = kaiTools.genSatCompPbsScript(scId,fdir,cmd,account=acct)
			pbsCmd = ['qsub','-J','1-'+str(numSegments),pbsName]
			results = subprocess.run(pbsCmd,capture_output=True)
			jobId = results.stdout.decode('utf-8').split('.')
			pbsLockName = kaiTools.genSatCompLockScript(scId,fdir,account=acct)
			pbsLockCmd = ['qsub','-W','depend=afterok:'+jobId[0],pbsLockName]
			results = subprocess.run(pbsLockCmd,capture_output=True)
			lockId = results.stdout.decode('utf-8').split('.')
			lockFileName = os.path.join(fdir,scId+'.lock')
			while not os.path.exists(lockFileName):
				time.sleep(10)
			h5name = kaiTools.mergeFiles(scId,fdir,numSegments)
			kaiTools.addGAMERA(data,scIds[scId],h5name)
			kaiTools.matchUnits(data)
			cdfname = os.path.join(fdir, scId + '.comp.cdf')
			if os.path.exists(cdfname):
				print('Deleting %s' % cdfname)
				os.system('rm %s' % cdfname)
			print('Creating CDF file',cdfname,'with',scId,'and GAMERA data')
			dm.toCDF(cdfname,data)
			plotname = os.path.join(fdir,scId+'.png')
			print('Plotting results to',plotname)
			kv.compPlot(plotname,scId,data)
			if not keep:
				h5parts = glob.glob(os.path.join(fdir,scId)+'.*.sc.h5')
				for file in h5parts:
					os.remove(file)
				jobParts = glob.glob(os.path.join(fdir,scId)+
					'.o'+jobId[0].split('[')[0]+'.*')
				for file in jobParts:
					os.remove(file)
				pbsScripts = glob.glob(os.path.join(fdir,scId)+'*pbs')
				for file in pbsScripts:
					os.remove(file)
				lockLog = os.path.join(fdir,scId+'.o'+lockId[0])
				os.remove(lockLog)
				lockFile = os.path.join(fdir,scId+'.lock')
				os.remove(lockFile)
				lockFile = os.path.join(fdir,scId+'.xml')
				os.remove(lockFile)



