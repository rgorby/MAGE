#!/usr/bin/env python
#Make a quick figure of a Gamera magnetosphere run

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
import matplotlib.gridspec as gridspec
import numpy as np
import kaipy.gamera.gampp as gampp
import kaipy.gamera.rcmpp as rcmpp
import palettable
import os
import numpy.ma as ma

if __name__ == "__main__":
	#Defaults
	MHDCol = rcmpp.MHDCol
	MHDLW = rcmpp.MHDLW
	fdir = os.getcwd()
	ftag = "msphere"
	nStp = -1
	fOut = "qkrcmpic.png"
	MainS = """Creates simple multi-panel figure for RCM magnetosphere run
	Top Panel - XXX
	Bottom Panel - XXX
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-n' ,type=int,metavar="step" ,default=nStp,help="Time slice to plot (default: %(default)s)")
	parser.add_argument('-big',action='store_true', default=False,help="Plot entire RCM grid (default: %(default)s)")
	parser.add_argument('-wgt', action='store_true',default=False,help="Show wRCM instead of FTE (default: %(default)s)")
	parser.add_argument('-vol', action='store_true',default=False,help="Show FTV instead of FTE (default: %(default)s)")
	parser.add_argument('-kt' , action='store_true',default=False,help="Show temperature instead of FTE (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id + ".mhdrcm"
	nStp = args.n
	doBig = args.big

	doWgt = args.wgt
	doVol = args.vol
	doT   = args.kt

	rcmpp.doEll = not doBig

	#---------
	#Figure parameters
	xTail = -20.0
	xSun = 10.0

	yMax = 15.0
	xyBds = [xTail,xSun,-yMax,yMax]

	figSz = (12,6)
	eCol = "slategrey"
	eLW = 0.15
	cLW = 0.5
	vP = kv.genNorm(1.0e-1,1.0e+2,doLog=True)
	vS = kv.genNorm(0.0,0.25)
	vW = kv.genNorm(0,1)
	vV = kv.genNorm(1.0e-2,1.0,doLog=True)
	vT = kv.genNorm(0,50)

	Nc = 10
	nMin = 1.0
	nMax = 1.0e+3
	vD = kv.genNorm(nMin,nMax,doLog=True)
	cVals = np.logspace(1.0,3.0,Nc)

	pCMap = "viridis"
	sCMap = "terrain"
	dCMap = "cool"
	wCMap = "bwr_r"
	vCMap = "gnuplot2"
	#dCMap = palettable.cmocean.sequential.Algae_20_r.mpl_colormap
	
	#======
	#Init data
	rcmdata = gampp.GameraPipe(fdir,ftag)
	if (nStp<0):
		nStp = rcmdata.sFin
		print("Using Step %d"%(nStp))

	#======
	#Setup figure
	fig = plt.figure(figsize=figSz)
	gs = gridspec.GridSpec(3,3,height_ratios=[20,1.0,1.0],hspace=0.025)

	AxL = fig.add_subplot(gs[0,0])
	AxM = fig.add_subplot(gs[0,1])
	AxR = fig.add_subplot(gs[0,-1])

	AxC1 = fig.add_subplot(gs[-1,0])
	AxC2 = fig.add_subplot(gs[-1,1])
	AxC3 = fig.add_subplot(gs[-1,-1])
	kv.genCB(AxC1,vP,"Pressure [nPa]",cM=pCMap)
	kv.genCB(AxC2,vD,"Density [#/cc]",cM=dCMap)
	if (doWgt):
		kv.genCB(AxC3,vW,r"wRCM",cM=wCMap)
	elif (doVol):
		kv.genCB(AxC3,vV,r"Flux-Tube Volume [Re/nT]",cM=vCMap)
	elif (doT):
		kv.genCB(AxC3,vT,r"Temperature [keV]",cM=vCMap)
	else:	
		kv.genCB(AxC3,vS,r"Flux-Tube Entropy [nPa (R$_{E}$/nT)$^{\gamma}$]",cM=sCMap)

	AxL.clear()
	AxM.clear()
	AxR.clear()

	bmX,bmY = rcmpp.RCMEq(rcmdata,nStp,doMask=True)
	I = rcmpp.GetMask(rcmdata,nStp)
	Ni = (~I).sum()
	
	if (Ni == 0):
		print("No closed field region in RCM, exiting ...")
		exit()

	Prcm  = rcmpp.GetVarMask(rcmdata,nStp,"P"    ,I)
	Nrcm  = rcmpp.GetVarMask(rcmdata,nStp,"N"    ,I)

	Pmhd  = rcmpp.GetVarMask(rcmdata,nStp,"Pmhd" ,I)
	Nmhd  = rcmpp.GetVarMask(rcmdata,nStp,"Nmhd" ,I)
	S     = rcmpp.GetVarMask(rcmdata,nStp,"S"    ,I)
	toMHD = rcmpp.GetVarMask(rcmdata,nStp,"toMHD",I)
	
		
	if (doWgt):
		wRCM  = rcmpp.GetVarMask(rcmdata,nStp,"wIMAG" ,I)
	if (doVol):
		bVol = rcmpp.GetVarMask(rcmdata,nStp,"bVol" ,I)
	if (doBig):
		toRCM = rcmpp.GetVarMask(rcmdata,nStp,"IOpen" ,I)

	AxL.set_title("RCM Pressure")

	AxL.pcolor(bmX,bmY,Prcm,norm=vP,cmap=pCMap)
	AxL.plot(bmX,bmY,color=eCol,linewidth=eLW)
	AxL.plot(bmX.T,bmY.T,color=eCol,linewidth=eLW)
	kv.addEarth2D(ax=AxL)
	kv.SetAx(xyBds,AxL)

	#Handle left
	AxM.set_title("MHD Pressure")
	AxM.pcolor(bmX,bmY,Pmhd,norm=vP,cmap=pCMap)
	
	AxM.contour(bmX,bmY,Nmhd,cVals,norm=vD,cmap=dCMap,linewidths=cLW)
	kv.addEarth2D(ax=AxM)
	kv.SetAx(xyBds,AxM)
	Axs = [AxL,AxM,AxR]

	if (nStp>0):
		for n in range(3):
			Ax = Axs[n]

			CS1 = Ax.contour(bmX,bmY,toMHD,[0.5],colors=MHDCol,linewidths=MHDLW)
			manloc = [(0.0,8.0)]

			fmt = {}
			fmt[0.5] = 'MHD'
			Ax.clabel(CS1,CS1.levels[::2],inline=True,fmt=fmt,fontsize=5,inline_spacing=25,manual=manloc)
			if (doBig):
				CS2 = Ax.contour(bmX,bmY,toRCM,[-0.5],colors=rcmpp.rcmCol,linewidths=MHDLW,linestyles='solid')

	#Handle right
	if (doWgt):
		AxR.set_title("RCM Weight")
		AxR.pcolor(bmX,bmY,wRCM,norm=vW,cmap=wCMap)
	elif (doVol):
		AxR.set_title("Flux-tube Volume")
		AxR.pcolor(bmX,bmY,bVol,norm=vV,cmap=vCMap)
	elif (doT):
		kT = 6.25*Prcm/Nrcm
		AxR.set_title("RCM Temperature")
		AxR.pcolor(bmX,bmY,kT,norm=vT,cmap=vCMap)
	else:	
		AxR.set_title("Flux-Tube Entropy")
		AxR.pcolor(bmX,bmY,S,norm=vS,cmap=sCMap)
	
	AxR.plot(bmX,bmY,color=eCol,linewidth=eLW)
	AxR.plot(bmX.T,bmY.T,color=eCol,linewidth=eLW)
	kv.addEarth2D(ax=AxR)
	kv.SetAx(xyBds,AxR)


	plt.suptitle("Step#%d"%(nStp),fontsize="x-large")

	fOut = "qkrcmpic.png"
	kv.savePic(fOut)
	

