#!/usr/bin/env python
#Make a plot of Dst from Gamera-RCM
import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import h5py
import matplotlib.pyplot as plt
import kaipy.gamera.magsphere as msph
import kaipy.kaiViz as kv
import numpy as np
import kaipy.kaiTools as kt
import datetime
from matplotlib import dates
import matplotlib.gridspec as gridspec
from astropy.time import Time
import os
import kaipy.kaiH5 as kaiH5

if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	ftag = "msphere"
	swfname = "bcwind.h5"
	tpad = 8 #Number of hours beyond MHD to plot
	iMax = -1
	MainS = """Creates simple plot comparing SYM-H from OMNI dataset to Gamera-RCM.
	Need to run or point to directory that has the bcwind and msphere.gam files of interest
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-tpad',type=float,metavar="time padding",default=tpad,help="Time beyond MHD data (in hours) to plot (default: %(default)s)")
	parser.add_argument('-swfile',type=str,metavar='filename',default=swfname,help="Solar wind file name (default: %(default)s)")
	#Finalizing parsing
	args = parser.parse_args()
	fdir = args.d
	tpad = args.tpad
	swfname = args.swfile

	#UT formats for plotting
	isotfmt = '%Y-%m-%dT%H:%M:%S.%f'
	t0fmt   = '%Y-%m-%d %H:%M:%S'
	utfmt   = '%H:%M \n%Y-%m-%d'


	fBC = os.path.join(fdir, swfname)
	kaiH5.CheckOrDie(fBC)
	ut_symh,tD,dstD = kt.GetSymH(fBC)

	fvolt = os.path.join(fdir,ftag+".volt.h5")
	BSDst = kaiH5.getTs(fvolt,sIds=None,aID="BSDst")
	MJD   = kaiH5.getTs(fvolt,sIds=None,aID="MJD")
	I = np.isinf(MJD)
	MJD0 = MJD[~I].min()-1
	MJD[I] = MJD0

	tScl = 1.0/(60.0*60)
	UT = Time(MJD,format='mjd').isot
	
	ut = [datetime.datetime.strptime(UT[n],isotfmt) for n in range(len(UT))]
	if iMax != -1:
		iMax = np.min(len(ut)-1,iMax)
	else:
		iMax = len(ut)-1

	# Remove Restart Step. Tends to cause weird artifacts
	deldt = []
	for it in range(iMax,1,-1):
		dt = ut[it] - ut[it-1]
		dt = dt.total_seconds()
		if dt < 2.:
			deldt.append(it)	

	BSDst = np.delete( BSDst,deldt )
	ut = np.delete( ut,deldt )
	LW = 0.75
	fSz = (14,7)
	fig = plt.figure(figsize=fSz)
	gs = gridspec.GridSpec(1,1,hspace=0.05,wspace=0.05)
	ax=fig.add_subplot(gs[0,0])
	ax.plot(ut_symh,dstD,label="SYM-H",linewidth=2*LW)
	ax.plot(ut,BSDst,label="BSDst",linewidth=LW)
	ax.legend(loc='upper right',fontsize="small",ncol=2)
	ax.axhline(color='magenta',linewidth=0.5*LW)
	ax.xaxis_date()
	xfmt = dates.DateFormatter(utfmt)
	ax.set_ylabel("Dst [nT]")
	ax.xaxis.set_major_formatter(xfmt)
	
	xMinD = np.array(ut_symh).min()
	xMaxD = np.array(ut_symh).max()
	xMinS = np.array(ut).min()
	xMaxS = np.array(ut).max()
	
	if (xMaxD>xMaxS):
		xMax = min(xMaxS+datetime.timedelta(hours=tpad),xMaxD)
	else:
		xMax = xMaxS
	xMin = xMinD

	ax.set_xlim(xMin,xMax)
	kv.savePic("qkdstpic.png")
