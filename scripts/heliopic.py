#!/usr/bin/env python
#Make a quick figure of a Gamera helio run

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes #do not need for helio if not inserting anything
import matplotlib.gridspec as gridspec
import numpy as np

import kaipy.gamhelio.helioViz as hviz
import kaipy.gamhelio.heliosphere as hsph
import kaipy.gamera.gampp as gampp #serial/MPI gamera output
import os

if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	ftag = "wsa"
	nStp = -1
	fOut = "qkpichelio.png"
	doDen = False
	#noIon = False
	noMPI = False
	#doJy = False
	#doBz = False
	MainS = """Creates simple multi-panel figure for Gamera helio run
	Top Panel - density and radial velocity in equatorial plane
	Bottom Panel - density and radial velocity in meridional plane
	"""

	# For helio we can choose what to show by default and have a couple of options "instead" defaults

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-n' ,type=int,metavar="step" ,default=nStp,help="Time slice to plot (default: %(default)s)")
	#parser.add_argument('-bz'   , action='store_true', default=doBz ,help="Show Bz instead of dBz (default: %(default)s)")
	parser.add_argument('-den'  , action='store_true', default=doDen,help="Show density instead of pressure (default: %(default)s)")
	#parser.add_argument('-jy'   , action='store_true', default=doJy ,help="Show Jy instead of pressure (default: %(default)s)")
	#parser.add_argument('-noion', action='store_true', default=noIon,help="Don't show ReMIX data (default: %(default)s)")
	parser.add_argument('-nompi', action='store_true', default=noMPI,help="Don't show MPI boundaries (default: %(default)s)")
	
	#size of domain - do not need right now. Think how it can be useful.
	#hviz.AddSizeArgs(parser)

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id
	nStp = args.n
	doDen = args.den
	#noIon = args.noion
	noMPI = args.nompi
	doMPI = (not noMPI)
	#doJy = args.jy
	#doBz = args.bz

	#Get domain size
	xyBds = hviz.GetSizeBds()
	print (xyBds)

	#---------------------
	#Do work
	doFast=False
	#doIon = not noIon

	#---------
	#Figure parameters
	figSz = (10,18)
	
	#======
	#Init data
	gsph = hsph.GamsphPipe(fdir,ftag,doFast=doFast)

	if (nStp<0):
		#assume this is working for helo but who knows
		nStp = hsph.sFin
		print("Using Step %d"%(nStp))

	print (nStp)

	#======
	#Setup figure
	fig = plt.figure(figsize=figSz)
	gs = gridspec.GridSpec(6,6,height_ratios=[20,1,20,1,10,1])
	

	AxL0 = fig.add_subplot(gs[0,0:3])
	AxR0 = fig.add_subplot(gs[0,3:])

	AxL1 = fig.add_subplot(gs[2,0:3])
	AxR1 = fig.add_subplot(gs[2,3:])

	AxL2 = fig.add_subplot(gs[4,0:3])
	AxR2 = fig.add_subplot(gs[4,3:])

	#colorbars
	AxC1_0 = fig.add_subplot(gs[1,0:3])
	AxC2_0 = fig.add_subplot(gs[1,3:])

	AxC1_1 = fig.add_subplot(gs[3,0:3])
	AxC2_1 = fig.add_subplot(gs[3,3:])

	AxC1_2 = fig.add_subplot(gs[5,0:3])
	AxC2_2 = fig.add_subplot(gs[5,3:])


	#plotting equatorial magV (left plot in msphere pic)
	hviz.PlotEqMagV(gsph,nStp,xyBds,AxL0,AxC1_0)
	hviz.PlotEqTemp(gsph,nStp,xyBds,AxR0,AxC2_0)

	hviz.PlotMerMagV(gsph,nStp,xyBds,AxL1,AxC1_1)
	hviz.PlotMerDNorm(gsph,nStp,xyBds,AxR1,AxC2_1)
	
	hviz.PlotiSlMagV(gsph,nStp,xyBds,AxL2,AxC1_2)
	hviz.PlotiSlD(gsph,nStp,xyBds,AxR2,AxC2_2)
	
	#Add time (upper left)
	gsph.AddTime(nStp,AxL0,xy=[0.025,0.875],fs="x-large")


	#Add wsa info (lower left) instead of Solar wind params
	#gsph.AddSW(nStp,AxL,xy=[0.625,0.025],fs="small")

	#Add inset RCM plot
	#if (doRCM):
	#	AxRCM = inset_axes(AxL,width="30%",height="30%",loc=3)
	#	rcmpp.RCMInset(AxRCM,rcmdata,nStp,mviz.vP)
	#	#Add some dBz contours
	#	AxRCM.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors=mviz.bz0Col,linewidths=mviz.cLW)
	#	#AxRCM.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(dBz),dbzVals,norm=vDB,cmap=mviz.dbCM,linewidths=0.25)
	#	rcmpp.AddRCMBox(AxL)

	#if (doIon):
	#	gsph.AddCPCP(nStp,AxR,xy=[0.610,0.925])

	#if (doIon):
	#	mviz.AddIonBoxes(gs[0,3:],ion)

	#Add MPI decomp
	#if (doMPI):
	#	hviz.PlotMPI(gsph,AxL0)
	#	hviz.PlotMPI(gsph,AxR0)
	#	hviz.PlotMPI(gsph,AxL1)
	#	hviz.PlotMPI(gsph,AxR1)
	#	hviz.PlotMPI(gsph,AxL2)
	#	hviz.PlotMPI(gsph,AxR2)

	kv.savePic(fOut,bLenX=45)
