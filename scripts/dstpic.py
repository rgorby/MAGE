#!/usr/bin/env python
#Make a plot of Dst from Gamera-RCM
import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import h5py
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
import numpy as np
import kaipy.gamera.dstutils as dstutils
import datetime
from matplotlib import dates
import matplotlib.gridspec as gridspec
from astropy.time import Time
import os

if __name__ == "__main__":
	#Defaults
	Dir = os.getcwd()
	nSk = 4	 #stride of time steps to calculate Dst

	MainS = """Creates simple plot comparing SYM-H from OMNI dataset to Gamera-RCM.
	Need to run or point to directory that has the bcwind and msphere.gam files of interest
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=Dir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-nsk',type=int,metavar="step stride",default=nSk,help="Stride between steps used to calculate Dst (default: %(default)s)")

	#Finalizing parsing
	args = parser.parse_args()
	Dir = args.d
	nSk = args.nsk

	#UT formats for plotting
	isotfmt = '%Y-%m-%dT%H:%M:%S.%f'
	t0fmt = '%Y-%m-%d %H:%M:%S'
	utfmt='%H:%M \n%Y-%m-%d'

	iS = [0,1,3,5,7]#,3]#,5,7]
	X0 = 0.0
	Y0 = 0.0
	Z0 = 0.0

	Dir = os.getcwd()

	NumI = len(iS)

	fBC = "%s/bcwind.h5"%(Dir)
	utD,tD,dstD = dstutils.GetSymH(fBC)
	ut_symh=[]
	[ut_symh.append(datetime.datetime.strptime(utD[n].decode('utf-8'),t0fmt)) for n in range(len(utD))]

	fIn = "%s/msphere.gam.h5"%(Dir)

	s0,sE = dstutils.GetSteps(fIn)

	Ns = len(np.arange(s0,sE,nSk))

	T = np.zeros(Ns)
	MJD = np.zeros(Ns)
	Dst = np.zeros((NumI,Ns))

	Nc = 10

	mp = 0
	for nStp in range(s0,sE,nSk):
		Xc,Yc,Zc,dV = dstutils.GetCC(fIn)
		mjd,t,Bz = dstutils.GetBz(fIn,nStp,Xc,Yc,Zc,dV,X0,Y0,Z0)
		T[mp] = t
		MJD[mp] = mjd
		for i in range(NumI):
			i0 = iS[i]
			Dst[i,mp] = Bz[i0:,:,:].sum()
		if (mp % Nc == 0):
			i0 = np.argmin(np.abs(tD-t))
			dDat = dstD[i0]
			print("Step %d / DST = %f / Data = %f"%(nStp,Dst[0,mp],dDat))
		mp = mp + 1

	tScl = 1.0/(60.0*60)
	UT = Time(MJD+0.04166666666,format='mjd').isot #time shifting by an hour
	ut = [datetime.datetime.strptime(UT[n],isotfmt) for n in range(len(UT))]

	LW = 0.75
	fSz = (14,7)
	fig = plt.figure(figsize=fSz)
	gs = gridspec.GridSpec(1,1,hspace=0.05,wspace=0.05)
	ax=fig.add_subplot(gs[0,0])
	ax.plot(ut_symh,dstD,label="SYM-H",linewidth=2*LW)
	for i in range(NumI):
		iStr = "GameRCM (i>=%d)"%(iS[i]+1)
		ax.plot(ut,Dst[i,:],label=iStr,linewidth=LW)
	ax.legend(loc='lower left',fontsize="small",ncol=2)
	ax.axhline(color='magenta',linewidth=2*LW)
	ax.xaxis_date()
	xfmt = dates.DateFormatter(utfmt)
	ax.set_ylabel("Dst [nT]")
	ax.xaxis.set_major_formatter(xfmt)
	kv.savePic("qkdstpic.png")
