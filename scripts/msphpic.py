#!/usr/bin/env python
#Make a quick figure of a Gamera magnetosphere run

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
import numpy as np
import kaipy.gamera.msphViz as mviz
import kaipy.gamera.remixpp as rmpp
import kaipy.gamera.magsphere as msph
import kaipy.gamera.gampp as gampp
import kaipy.gamera.rcmpp as rcmpp
import os

if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	ftag = "msphere"
	nStp = -1
	fOut = "qkpic.png"
	doDen = False
	noIon = False
	noMPI = False

	MainS = """Creates simple multi-panel figure for Gamera magnetosphere run
	Top Panel - Residual vertical magnetic field
	Bottom Panel - Pressure (or density) and hemispherical insets
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-n' ,type=int,metavar="step" ,default=nStp,help="Time slice to plot (default: %(default)s)")
	parser.add_argument('-den', action='store_true', default=doDen,help="Show density instead of pressure (default: %(default)s)")
	parser.add_argument('-noion', action='store_true', default=noIon,help="Don't show ReMIX data (default: %(default)s)")
	parser.add_argument('-nompi', action='store_true', default=noMPI,help="Don't show MPI boundaries (default: %(default)s)")

	mviz.AddSizeArgs(parser)

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id
	nStp = args.n
	doDen = args.den
	noIon = args.noion
	noMPI = args.nompi
	doMPI = (not noMPI)
	
	#Get domain size
	xyBds = mviz.GetSizeBds(args)

	#---------------------
	#Do work
	doFast=False
	doIon = not noIon

	#---------
	#Figure parameters
	figSz = (12,7.5)
	
	#======
	#Init data
	gsph = msph.GamsphPipe(fdir,ftag,doFast=doFast)

	#Check for remix
	rcmChk = fdir + "/%s.mhdrcm.h5"%(ftag)
	doRCM = os.path.exists(rcmChk)
	if (doRCM):
		print("Found RCM data")
		rcmdata = gampp.GameraPipe(fdir,ftag+".mhdrcm")

	if (nStp<0):
		nStp = gsph.sFin
		print("Using Step %d"%(nStp))
	#======
	#Setup figure
	fig = plt.figure(figsize=figSz)
	gs = gridspec.GridSpec(3,4,height_ratios=[20,1.0,1.0],hspace=0.025)

	AxL = fig.add_subplot(gs[0,0:2])
	AxR = fig.add_subplot(gs[0,2:])
	AxC1 = fig.add_subplot(gs[2,0])
	AxC2 = fig.add_subplot(gs[2,3])
	AxC3 = fig.add_subplot(gs[2,1])
	AxC4 = fig.add_subplot(gs[2,2])

	rmpp.cMax = 1.00
	kv.genCB(AxC4,kv.genNorm(rmpp.cMax),"FAC",cM=rmpp.fcMap,Ntk=4)
	rmpp.AddPotCB(AxC3)
	
	mviz.PlotEqB(gsph,nStp,xyBds,AxL,AxC1)
	mviz.PlotMerid(gsph,nStp,xyBds,AxR,doDen,doRCM,AxC2)
	
	gsph.AddTime(nStp,AxL,xy=[0.025,0.89],fs="x-large")
	gsph.AddSW(nStp,AxL,xy=[0.625,0.025],fs="small")

	#Add inset RCM plot
	if (doRCM):
		AxRCM = inset_axes(AxL,width="30%",height="30%",loc=3)
		rcmpp.RCMInset(AxRCM,rcmdata,nStp,mviz.vP)
		rcmpp.AddRCMBox(AxL)

	if (doIon):
		gsph.AddCPCP(nStp,AxR,xy=[0.610,0.925])

	if (doIon):
		dxy = [32.5,32.5]
		gsph.CMIViz(AxR,nStp,dxy=dxy,loc="upper left",doNorth=True)
		gsph.CMIViz(AxR,nStp,dxy=dxy,loc="lower left",doNorth=False)

	#Add MPI decomp
	if (doMPI):
		mviz.PlotMPI(gsph,AxL)
		mviz.PlotMPI(gsph,AxR)

	kv.savePic(fOut,bLenX=45)