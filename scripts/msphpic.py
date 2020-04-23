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
	doBig = False
	doSmall = False
	doHuge = False
	doBigger = False
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
	parser.add_argument('-big', action='store_true', default=doBig,help="Use larger domain bounds (default: %(default)s)")
	parser.add_argument('-small', action='store_true', default=doSmall,help="Use smaller domain bounds (default: %(default)s)")
	parser.add_argument('-huge', action='store_true', default=doHuge,help="Show full domain (default: %(default)s)")
	parser.add_argument('-bigger', action='store_true', default=doBigger,help="Use larger domain bounds (default: %(default)s)")
	parser.add_argument('-nompi', action='store_true', default=noMPI,help="Don't show MPI boundaries (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id
	nStp = args.n
	doDen = args.den
	noIon = args.noion
	doBig = args.big
	doSmall = args.small
	doHuge = args.huge
	doBigger = args.bigger
	noMPI = args.nompi
	doMPI = (not noMPI)
	

	#---------------------
	#Do work
	doFast=False
	doIon = not noIon

	#---------
	#Figure parameters
	figSz = (12,7.5)
	dbCMap = "RdGy_r"
	pCMap = "viridis"
	dCMap = "viridis"
	
	if (doSmall):
		xTail = -10.0
		xSun = 5.0
	elif (doBig):
		xTail = -100.0
		xSun = 20.0
	elif (doBigger):
		xTail = -200.0
		xSun = 25.0
	elif (doHuge):
		xTail = -350.0
		xSun = 40.0
	else:
		xTail = -40.0
		xSun = 20.0

	yMax = (xSun-xTail)/2.0
	xyBds = [xTail,xSun,-yMax,yMax]

	

	cLW = 0.25
	vcLW = 0.5
	vP = kv.genNorm(1.0e-2,10.0,doLog=True)
	vDB = kv.genNorm(25)
	vDD = kv.genNorm(0,25)
	Nc = 11
	if (doDen):
		cVals = np.linspace(0,25,Nc)
	else:
		cVals = np.logspace(np.log10(1.0),np.log10(10.0),Nc)
	LW = 0.25
	ashd = 0.5

	#======
	#Init data
	gsph = msph.GamsphPipe(fdir,ftag,doFast=doFast)

	#Check for remix
	rcmChk = fdir + "/%s.mhdrcm.h5"%(ftag)
	doRCM = os.path.exists(rcmChk)
	if (doRCM):
		print("Found RCM data")
		rcmdata = gampp.GameraPipe(fdir,ftag+".mhdrcm")
		vDD = kv.genNorm(1.0,1.0e+3,doLog=True)
		if (doDen):
			cVals = np.logspace(0.0,3.0,Nc)

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

	kv.genCB(AxC1,vDB,"Residual Field [nT]",cM=dbCMap,Ntk=7)
	rmpp.cMax = 1.00
	kv.genCB(AxC4,kv.genNorm(rmpp.cMax),"FAC",cM=rmpp.fcMap,Ntk=4)
	rmpp.AddPotCB(AxC3)
	
	if (doDen):
		if (doRCM):
			kv.genCB(AxC2,vDD,"Density [#/cc]",cM=dCMap)
		else:
			kv.genCB(AxC2,vDD,"Density [#/cc]",cM=dCMap,Ntk=7)
		
	else:
		kv.genCB(AxC2,vP,"Pressure [nPa]",cM=pCMap)

	AxL.clear()
	AxR.clear()

	dbz = gsph.DelBz(nStp)
	Bz = gsph.EggSlice("Bz",nStp,doEq=True)

	#Plot left
	AxL.pcolormesh(gsph.xxi,gsph.yyi,dbz,cmap=dbCMap,norm=vDB)
	AxL.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors='magenta',linewidths=cLW)

	if (doDen):
		Dxz = gsph.EggSlice("D",nStp,doEq=False)
		Dxy = gsph.EggSlice("D",nStp,doEq=True)
	else:
		Pxz = gsph.EggSlice("P",nStp,vScl=gsph.pScl,doEq=False)
		Pxy = gsph.EggSlice("P",nStp,vScl=gsph.pScl,doEq=True)
	
	kv.addEarth2D(ax=AxL)
	kv.SetAx(xyBds,AxL)
	gsph.AddTime(nStp,AxL,xy=[0.025,0.89],fs="x-large")
	gsph.AddSW(nStp,AxL,xy=[0.625,0.025],fs="small")
	AxL.set_xlabel('SM-X [Re]')
	AxL.set_ylabel('SM-Y [Re]')

	#Add inset RCM plot
	if (doRCM):
		AxRCM = inset_axes(AxL,width="30%",height="30%",loc=3)
		rcmpp.RCMInset(AxRCM,rcmdata,nStp,vP)
		rcmpp.AddRCMBox(AxL)

	#Add contour to equatorial plot and do right plot
	if (doDen):
		AxR.pcolormesh(gsph.xxi,gsph.yyi,Dxz,cmap=dCMap,norm=vDD)
		AxL.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Dxy),cVals,norm=vDD,cmap=dCMap,linewidths=vcLW)

	else:
		AxR.pcolormesh(gsph.xxi,gsph.yyi,Pxz,cmap=pCMap,norm=vP)
		AxL.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Pxy),cVals,norm=vP,cmap=pCMap,linewidths=vcLW)

	kv.addEarth2D(ax=AxR)
	kv.SetAx(xyBds,AxR)
	AxR.yaxis.tick_right()
	AxR.yaxis.set_label_position('right')
	if (doIon):
		gsph.AddCPCP(nStp,AxR,xy=[0.610,0.925])

	AxR.set_xlabel('SM-X [Re]')
	AxR.set_ylabel('SM-Z [Re]')

	if (doIon):
		dxy = [32.5,32.5]
		gsph.CMIViz(AxR,nStp,dxy=dxy,loc="upper left",doNorth=True)
		gsph.CMIViz(AxR,nStp,dxy=dxy,loc="lower left",doNorth=False)

	#Add MPI decomp
	if (doMPI):
		gCol = "deepskyblue"
		for i in range(gsph.Ri):
			i0 = i*gsph.dNi
			AxL.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],gCol,linewidth=LW,alpha=ashd)
			AxR.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],gCol,linewidth=LW,alpha=ashd)
		if (gsph.Rj>1):
			for j in range(1,gsph.Rj):
				j0 = j*gsph.dNj
				AxL.plot(gsph.xxi[:,j0] ,gsph.yyi[:,j0],gCol,linewidth=LW,alpha=ashd)
				AxL.plot(gsph.xxi[:,j0],-gsph.yyi[:,j0],gCol,linewidth=LW,alpha=ashd)
				AxR.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=LW,alpha=ashd)
				AxR.plot(gsph.xxi[:,j0],-gsph.yyi[:,j0],gCol,linewidth=LW,alpha=ashd)
			#X-axis (+)
			AxL.plot(gsph.xxi[:,0], gsph.yyi[:,0],gCol,linewidth=LW,alpha=ashd)
			AxR.plot(gsph.xxi[:,0], gsph.yyi[:,0],gCol,linewidth=LW,alpha=ashd)
			#X-axis (-)
			j0 = (gsph.Rj)*gsph.dNj
			AxL.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=LW,alpha=ashd)
			AxR.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=LW,alpha=ashd)



	kv.savePic(fOut,bLenX=45)
