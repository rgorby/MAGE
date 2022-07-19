import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import dates
from astropy.time import Time
import datetime
import os, sys
import progressbar
import argparse
from argparse import RawTextHelpFormatter
import numpy as np

import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv
import kaipy.kaiTools as kT
import kaipy.satcomp.scutils as scutils
import kaipy.satcomp.scRCM as scRCM

def fmtTKL(AxTKL):
	AxTKL.set_yscale('log')
	AxTKL.tick_params(axis='y', pad=-1)
	AxTKL.yaxis.labelpad = -1
	AxTKL.tick_params(axis='x', pad=-1)
	AxTKL.xaxis.labelpad = -1
	AxTKL.title.set_text(str(ut_tkl[0]))

if __name__=="__main__":
	fdir  = os.getcwd()
	ftag  = "msphere"
	trtag = "RBSP-%s_MAGNETOMETER_1SEC-GSM_EMFISIS-L3.sc.h5"  # Spacecraft trajectory and values along track
	vTag = "H-PAP_RBSPICE"
	tStart = -1
	tEnd = -1
	tStride = 10
	vidOut = "vid_rcm-rbsp-comp"

	tklV_choices = ['press', 'odf']

	jdir = "jstore"

	pIDs = ['cdas','times','track','tkl']

	MainS = """Pulls RBSP data and compares it to synthetic RBSP intensity measurementsfrom the simulation, 
	calculated from extracted RBSP trajectory and PSD files.
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of model data (default: %(default)s)")
	parser.add_argument('-trj',type=str,metavar="scTrk",default=trtag,help="spacecraft trajectory file (default: %(default)s)")
	parser.add_argument('-jdir',type=str,metavar="directory",default=jdir,help="Directory to store and find json files (default: %(default)s)")
	parser.add_argument('-scId', type=str,choices=scRCM.supportedSats[:],default="RBSPB",help="Sat id (default: %(default)s)")
	parser.add_argument('-v', type=str,choices=scRCM.supportedDsets[:],default="Hydrogen_omniflux_RBSPICE",help="Dataset (default: %(default)s)")
	parser.add_argument('-tStart',type=int, default=tStart,help="Starting time step for L vs. E calculation (default: First step in RCM data)")
	parser.add_argument('-tEnd',type=int, default=tEnd,help="Ending time step for L vs. E calculation (default: Last step in RCM data)")
	parser.add_argument('-tStride',type=int, default=tStride,help="Time step stride for L vs. E calculation (default: %(default)s)")
	parser.add_argument('-plotTag',type=str,default="",help="Extra tag for each plot")
	parser.add_argument('-vidOut',type=str,default=vidOut,help="Output directory (relative to -d) for video images (default: %(default)s)")
	parser.add_argument('-tklv', type=str,choices=tklV_choices,default=tklV_choices[0],help="Variable to plot in Lvsk panel (default: %(default)s)")
	parser.add_argument('-forceCalc',type=str,metavar=pIDs,default="",help="Comma-separated process IDs to force recalculation for given process")
	parser.add_argument('-HOPESPICE',action='store_true',help="Combine HOPE and RBSPICE hydrogen data")
	#Finalize parsing
	args  = parser.parse_args()
	fdir  = args.d
	ftag  = args.id
	trtag = args.trj
	jdir  = args.jdir
	scId = args.scId
	vTag = args.v
	tStart = args.tStart
	tEnd = args.tEnd
	tStride = args.tStride
	plotTag = args.plotTag
	vidOut = args.vidOut
	tklv = args.tklv
	fcStr = args.forceCalc
	doHopeSpice = args.HOPESPICE

	#Extract RBSP identifier (A or B)
	scTag = trtag.split('RBSP')[1][:2]
	if '-' in scTag:
		scTag = trtag.split('-')[1][0]
	else:
		scTag = scTag[0]

	if fcStr == "":
		fcIDs = []
	elif "all" in fcStr:
		fcIDs = pIDs
	else:
		fcIDs = fcStr.split(',')


	#======
	#Init data
	#======
	mhdrcm_fname = os.path.join(fdir, ftag+'.mhdrcm.h5')
	rcm_fname    = os.path.join(fdir, ftag+'.rcm.h5'   )
	trackf5      = os.path.join(fdir, trtag            )

	kh5.CheckOrDie(mhdrcm_fname)
	kh5.CheckOrDie(rcm_fname)
	kh5.CheckOrDie(trackf5)
	kh5.CheckDirOrMake(jdir)

	#Sort of start and end times
	rcmNt, rcmSIDs = kh5.cntSteps(rcm_fname)
	rcmSIDs = np.sort(rcmSIDs)
	if tStart == -1:
		tStart = rcmSIDs[0]
	elif tStart < rcmSIDs[0]:
		print("Step '{}' not in RCM times, starting from {}".format(tStart, rcmSIDs[0]))
		tStart = rcmSIDs[0]
	if tEnd == -1:
		tEnd = rcmSIDs[-1]
	elif tEnd > rcmSIDs[-1]:
		print("Step '{}' not in RCM times, ending at {}".format(tEnd, rcmSIDs[-1]))
		tEnd = rcmSIDs[-1]

	#Get start and end times from sctrack file
	isotfmt = '%Y-%m-%dT%H:%M:%S.%f'
	scMJDs = kh5.PullVar(trackf5, 'MJDs')
	#ut = scutils.mjd_to_ut(scMJDs)
	ut = kT.MJD2UT(scMJDs)
	t0r = ut[0].strftime("%Y-%m-%dT%H:%M:%SZ")
	t1r = ut[-1].strftime("%Y-%m-%dT%H:%M:%SZ")

	print("Retrieving RBSPICE dataset")
	if doHopeSpice:
		ephData, scData = scRCM.getSCOmniDiffFlux(scId, "Hydrogen_omniflux_RBSPICE", t0r, t1r, jdir=jdir,forceCalc=('cdas' in fcIDs))
		hope_ephData, hope_scData = scRCM.getSCOmniDiffFlux(scId, "Hydrogen_PAFlux_HOPE", t0r, t1r, jdir=jdir,forceCalc=('cdas' in fcIDs))
	else:
		ephData, scData = scRCM.getSCOmniDiffFlux(scId, vTag, t0r, t1r, jdir=jdir,forceCalc=('cdas' in fcIDs))

	print("\n\nGrabbing RCM time")
	rcmTimes = scRCM.getRCMtimes(rcm_fname,mhdrcm_fname,jdir=jdir,forceCalc=('times' in fcIDs))

	print("\n\nExtracting RCM values along sc trajectory")
	rcmTrack = scRCM.getRCM_scTrack(trackf5, rcm_fname, rcmTimes, jdir=jdir, forceCalc=('track' in fcIDs), scName="RBSP-B")

	print("\n\nConsolidating grids")
	if doHopeSpice:
		scEGrid = scData['energies']
		hope_scEGrid = hope_scData['energies']
		species = scData['species']
		rcmEGrid = rcmTrack[species]['energies']
		eMax = np.max([scEGrid.max(), hope_scEGrid.max(), rcmEGrid.max()])
		eMin = np.min([scEGrid[scEGrid>0].min(), hope_scEGrid[hope_scEGrid>0].min(), rcmEGrid[rcmEGrid>0].min()])
		eGrid = np.logspace(np.log10(eMin), np.log10(eMax), 200, endpoint=True)
		consolData = scRCM.consolidateODFs(scData, rcmTrack, eGrid=eGrid, doPlot=False)
		hope_consolData = scRCM.consolidateODFs(hope_scData, rcmTrack, eGrid=eGrid)
	else:
		eGrid = np.logspace(np.log10(40), np.log10(6E2), 200, endpoint=True)
		consolData = scRCM.consolidateODFs(scData, rcmTrack, eGrid=eGrid)

	print("\n\nCalculating tkl vars(var wedge)")
	#tkldata = scRCM.getIntensitiesVsL('msphere.rcm.h5','msphere.mhdrcm.h5',tStart, tEnd, tStride, jdir=jdir, forceCalc=('tkl' in fcIDs))
	tkldata = scRCM.getVarWedge(rcm_fname, mhdrcm_fname, tStart, tEnd, tStride, 5, jdir=jdir, eGrid=eGrid*1E3, forceCalc=('tkl' in fcIDs))  # This function expects eGrid to be in eV
	
	#Works but very verbose
	#print("\n\nTesting RCM eqlatlon grab")
	rcm_eqlatlon = scRCM.getRCM_eqlatlon(mhdrcm_fname, rcmTimes, tStart, tEnd, tStride)

	print('\n\nPlotting')
	fig = plt.figure(figsize=(20,9))
	gs = gridspec.GridSpec(8,16, wspace=0.8, hspace=0.6)
	cmap_odf = "CMRmap"
	cmap_press = 'viridis'
	cmap_parpress = 'gnuplot2'
	cmap_rcm = "CMRmap"
	
	AxCB_odf = fig.add_subplot(gs[:,0])
	AxSC = fig.add_subplot(gs[0:4,1:8])
	AxRCM = fig.add_subplot(gs[4:8,1:8])

	AxTL = fig.add_subplot(gs[0:2, 8:12])
	AxTKL = fig.add_subplot(gs[2:7,8:12])
	

	AxRCMLatLon = fig.add_subplot(gs[:4, 12:16], projection='polar')
	AxRCMEq = fig.add_subplot(gs[4:7, 12:16])
	
	#If tkl plot will match Diff Flux colorbar, let it use it and make pressure cbar span tkl and eqlatlon plots
	if tklv == 'odf':
		AxCB_press  = fig.add_subplot(gs[7,8:16])
	elif tklv == 'press':
		#TKL will actually show partial pressure, which needs its own cbar
		AxCB_parpress = fig.add_subplot(gs[7,8:12])
		#Original colorbar only spans eqlatlon plots
		AxCB_press  = fig.add_subplot(gs[7,12:16])

	odfnorm = kv.genNorm(1E4, 5E6, doLog=True)
	ut_tkl = kT.MJD2UT(tkldata['MJD'])
	pressnorm = kv.genNorm(1E-2, 75, doLog=False)
	parpressnorm = kv.genNorm(1e-3, 5, doLog=True)
	#Movie time
	outdir = os.path.join(fdir, vidOut)
	kh5.CheckDirOrMake(outdir)
	n_pad = int(np.log10((len(tkldata['MJD'])))) + 1
	ticker = 0

	#Run first iteration manually
	pltmjd = tkldata['MJD'][0]
	if doHopeSpice:
		scRCM.plt_ODF_Comp(AxSC, AxRCM, AxCB_odf, hope_consolData, mjd=pltmjd, norm=odfnorm, cmapName=cmap_odf)
		scRCM.plt_ODF_Comp(AxSC, AxRCM, AxCB_odf, consolData, mjd=pltmjd, norm=odfnorm, cmapName=cmap_odf, forcePop=True)
	else:
		scRCM.plt_ODF_Comp(AxSC, AxRCM, AxCB_odf, consolData, mjd=pltmjd, norm=odfnorm, cmapName=cmap_odf)
	AxSC.title.set_text('RBSP RCM Comparison')
	AxCB_odf.yaxis.set_ticks_position('left')
	AxCB_odf.yaxis.set_label_position('left')
	AxSC.tick_params(axis='y', pad=-1)
	AxSC.yaxis.labelpad = -3
	AxRCM.tick_params(axis='y', pad=-1)
	AxRCM.yaxis.labelpad = -3
	
	scRCM.plt_tl(AxTL, tkldata, AxCB=AxCB_press, mjd=pltmjd, norm=pressnorm, cmapName=cmap_press)
	
	if tklv == 'odf':
		scRCM.plt_tkl(AxTKL, tkldata, vName=tklv, mjd=pltmjd, norm=odfnorm, cmapName=cmap_odf)
	elif tklv == 'press':
		#Actually partial pressure
		scRCM.plt_tkl(AxTKL, tkldata, vName=tklv, AxCB=AxCB_parpress, mjd=pltmjd, norm=parpressnorm, cmapName=cmap_parpress)
	

	AxTL.xaxis.set_ticks_position('top')
	AxTL.xaxis.set_label_position('top')
	AxTL.tick_params(axis='y', pad=-1)
	AxTL.yaxis.labelpad = -1 
	fmtTKL(AxTKL)
	AxCB_press.xaxis.labelpad = -1

	scRCM.plt_rcm_eqlatlon(AxRCMLatLon, AxRCMEq, rcm_eqlatlon, rcmTrack, mjd=pltmjd, norm=pressnorm, cmapName=cmap_press)

	if plotTag != "":
		plt.suptitle(plotTag, fontsize=20)

	ticker += 1
	filename = "{}.{:0>{n}d}.png".format("vid", ticker, n=n_pad)
	ofname = os.path.join(outdir, filename)
	kv.savePic(ofname)
	
	#Now just things that need to update
	bar = progressbar.ProgressBar(max_value=len(tkldata['MJD']))
	for n in range(1,len(tkldata['MJD'])):
		bar.update(n)
		pltmjd = tkldata['MJD'][n]
		
		scRCM.plt_ODF_Comp(AxSC, AxRCM, AxCB_odf, consolData, mjd=pltmjd, norm=odfnorm, cmapName=cmap_odf)
		
		scRCM.plt_tl(AxTL, tkldata, AxCB=AxCB_press, mjd=pltmjd, norm=pressnorm, cmapName=cmap_press)

		if tklv == 'odf':
			scRCM.plt_tkl(AxTKL, tkldata, vName=tklv, mjd=pltmjd, norm=odfnorm, cmapName=cmap_odf, satTrackData=rcmTrack)
		elif tklv == 'press':
			#Actually partial pressure
			scRCM.plt_tkl(AxTKL, tkldata, vName=tklv, AxCB=AxCB_parpress, mjd=pltmjd, norm=parpressnorm, cmapName=cmap_parpress, satTrackData=rcmTrack)
		fmtTKL(AxTKL)
		AxTKL.title.set_text(str(ut_tkl[n]))

		scRCM.plt_rcm_eqlatlon(AxRCMLatLon, AxRCMEq, rcm_eqlatlon, rcmTrack, mjd=pltmjd, norm=pressnorm, cmapName=cmap_press)
		#Doesn't make labels until show is called (i think)
		#tickLabels = [s.get_text() for s in AxTL.get_xticklabels()]
		#AxTL.set_xticklabels(tickLabels)

		ticker += 1
		filename = "{}.{:0>{n}d}.png".format("vid", ticker, n=n_pad)
		ofname = os.path.join(outdir, filename)
		kv.savePic(ofname)
		
	
	plt.show()
