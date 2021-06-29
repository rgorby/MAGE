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
import kaipy.remix.remix as remix
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
	noRCM = False
	doJy = False
	doBz = False
	doBigRCM = False
	doSrc = False

	MainS = """Creates simple multi-panel figure for Gamera magnetosphere run
	Top Panel - Residual vertical magnetic field
	Bottom Panel - Pressure (or density) and hemispherical insets
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-n' ,type=int,metavar="step" ,default=nStp,help="Time slice to plot (default: %(default)s)")
	parser.add_argument('-bz'   , action='store_true', default=doBz ,help="Show Bz instead of dBz (default: %(default)s)")
	parser.add_argument('-den'  , action='store_true', default=doDen,help="Show density instead of pressure (default: %(default)s)")
	parser.add_argument('-jy'   , action='store_true', default=doJy ,help="Show Jy instead of pressure (default: %(default)s)")
	parser.add_argument('-noion', action='store_true', default=noIon,help="Don't show ReMIX data (default: %(default)s)")
	parser.add_argument('-nompi', action='store_true', default=noMPI,help="Don't show MPI boundaries (default: %(default)s)")
	parser.add_argument('-norcm', action='store_true', default=noRCM,help="Don't show RCM data (default: %(default)s)")
	parser.add_argument('-bigrcm', action='store_true',default=doBigRCM,help="Show entire RCM domain (default: %(default)s)")
	parser.add_argument('-src'   , action='store_true', default=doSrc ,help="Show source term (default: %(default)s)")

	mviz.AddSizeArgs(parser)

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id
	nStp = args.n
	doDen = args.den
	noIon = args.noion
	#noMPI = args.nompi
	doMPI = (not noMPI)
	doJy = args.jy
	doSrc = args.src
	doBz = args.bz
	noRCM = args.norcm
	doBigRCM = args.bigrcm

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

	if (nStp<0):
		nStp = gsph.sFin
		print("Using Step %d"%(nStp))

	#Check for remix
	rcmChk = fdir + "/%s.mhdrcm.h5"%(ftag)
	rmxChk = fdir + "/%s.mix.h5"%(ftag)
	doRCM = os.path.exists(rcmChk)
	doMIX = os.path.exists(rmxChk)

	if (doRCM):
		print("Found RCM data")
		rcmdata = gampp.GameraPipe(fdir,ftag+".mhdrcm")
		mviz.vP = kv.genNorm(1.0e-2,100.0,doLog=True)
		rcmpp.doEll = not doBigRCM
	if (doMIX):
		print("Found ReMIX data")
		ion = remix.remix(rmxChk,nStp)

	#======
	#Setup figure
	fig = plt.figure(figsize=figSz)
	gs = gridspec.GridSpec(3,6,height_ratios=[20,1,1],hspace=0.025)
	

	AxL = fig.add_subplot(gs[0,0:3])
	AxR = fig.add_subplot(gs[0,3:])

	AxC1 = fig.add_subplot(gs[-1,0:2])
	AxC2 = fig.add_subplot(gs[-1,2:4])
	AxC3 = fig.add_subplot(gs[-1,4:6])


	cbM = kv.genCB(AxC2,kv.genNorm(remix.facMax),"FAC",cM=remix.facCM,Ntk=4)
	AxC2.xaxis.set_ticks_position('top')

	
	Bz = mviz.PlotEqB(gsph,nStp,xyBds,AxL,AxC1,doBz=doBz)

	if (doJy):
		mviz.PlotJyXZ(gsph,nStp,xyBds,AxR,AxC3)
	else:
		mviz.PlotMerid(gsph,nStp,xyBds,AxR,doDen,doRCM,AxC3,doSrc=doSrc)
	

	gsph.AddTime(nStp,AxL,xy=[0.025,0.89],fs="x-large")
	gsph.AddSW(nStp,AxL,xy=[0.625,0.025],fs="small")

	#Add inset RCM plot
	if (not noRCM):
		AxRCM = inset_axes(AxL,width="30%",height="30%",loc=3)
		rcmpp.RCMInset(AxRCM,rcmdata,nStp,mviz.vP)
		#Add some dBz contours
		AxRCM.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors=mviz.bz0Col,linewidths=mviz.cLW)
		#AxRCM.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(dBz),dbzVals,norm=vDB,cmap=mviz.dbCM,linewidths=0.25)
		rcmpp.AddRCMBox(AxL)

	if (doIon):
		gsph.AddCPCP(nStp,AxR,xy=[0.610,0.925])

	if (doIon):
		mviz.AddIonBoxes(gs[0,3:],ion)

	#Add MPI decomp
	if (doMPI):
		mviz.PlotMPI(gsph,AxL)
		mviz.PlotMPI(gsph,AxR)

	kv.savePic(fOut,bLenX=45)
