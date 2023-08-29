#!/usr/bin/env python

#standard python
import sys
import os
import argparse
from argparse import RawTextHelpFormatter

#numpy
import numpy as np

#Kaipy and related
from   astropy.time import Time
import h5py
import kaipy.kaiH5 as kaiH5
import kaipy.kaiViz as kv
import kaipy.kaiTools as kaiTools
import kaipy.kaijson as kj
import kaipy.satcomp.scutils as scutils
import spacepy.datamodel as dm

if __name__ == '__main__':
	MainS = """Extracts information from satellite trajectory for various
	spacecraft.  Space craft data is pulled from CDAWeb.  Output CDF files
	contain data pulled from CDAWeb along with data extracted from GAMERA.
	Image files of satellite comparisons are also produced.
	"""


	parser = argparse.ArgumentParser(description=MainS,
		formatter_class=RawTextHelpFormatter)
	parser.add_argument('-id',type=str,metavar='runid',default='msphere',
		help='RunID of data (default: %(default)s)')
	parser.add_argument('-path',type=str,metavar='path',default='.',
		help='Path to directory containing REMIX files (default: %(default)s)')
	parser.add_argument('-cmd',type=str,metavar='command',default=None,
		help='Full path to sctrack.x command')
	parser.add_argument('-satId',type=str,metavar='Satellite Id',
		default=None,help='Name of Satellite to compare')
	parser.add_argument('-numSeg',type=int,metavar='Number of segments',
		default=1,help='Number of segments to simulateously process')
	parser.add_argument('--keep',action='store_true',
		help='Keep intermediate files')

	args = parser.parse_args()

	fdir = args.path
	ftag = args.id
	cmd = args.cmd
	scRequested = args.satId
	numSegments = args.numSeg
	keep = args.keep

	if fdir == '.':
		fdir = os.getcwd()

	if None == cmd:
		my_env = os.environ.copy()
		cmd = os.path.join(os.getenv('KAIJUDIR'),'build','bin','sctrack.x')
	if not (os.path.isfile(cmd) and os.access(cmd, os.X_OK)):
		print(cmd,'either not found or not executable')
		sys.exit()

	scIds = scutils.getScIds()

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

	#scToDo =['CLUSTER1']
	if None == scRequested:
		scToDo = scIds
	else:
		scToDo = []
		scToDo.append(scRequested)

	for scId in scToDo:
		print('Getting spacecraft data for', scId)
		status,data = scutils.getSatData(scIds[scId],
			t0.strftime("%Y-%m-%dT%H:%M:%SZ"),
			t1.strftime("%Y-%m-%dT%H:%M:%SZ"),deltaT)

		if status['http']['status_code'] != 200 or data is None:
			print('No data available for', scId)
		else:
			print('Extracting GAMERA data')
			toRe = scutils.extractGAMERA(data,scIds[scId],scId,
				mjdFileStart,secFileStart,fdir,
				ftag,cmd,numSegments,keep)
			scutils.matchUnits(data)
			cdfname = os.path.join(fdir, scId + '.comp.cdf')
			if os.path.exists(cdfname):
				print('Deleting %s' % cdfname)
				os.system('rm %s' % cdfname)
			print('Creating CDF file',cdfname,'with',scId,'and GAMERA data')
			dm.toCDF(cdfname,data)
			plotname = os.path.join(fdir,scId+'.png')
			print('Plotting results to',plotname)
			kv.compPlot(plotname,scId,data)
			print('Computing Errors')
			errname = os.path.join(fdir,scId+'-error.txt')
			scutils.errorReport(errname,scId,data)
			plotname = os.path.join(fdir,scId+'-traj.png')
			print('Plotting trajectory to',plotname)
			kv.trajPlot(plotname,scId,data,toRe)
			
