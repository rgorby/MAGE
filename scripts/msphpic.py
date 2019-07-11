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
import kaipy.gamera.remixpp as rmpp
import kaipy.gamera.magsphere as msph

import os

if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	ftag = "msphere"
	nStp = 0
	fOut = "qkpic.png"
	doDen = False
	noIon = False
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

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id
	nStp = args.n
	doDen = args.den
	noIon = args.noion

	#---------------------
	#Do work
	doFast=False
	doIon = not noIon

	#---------
	#Figure parameters
	figSz = (12,7.5)
	dbCMap = "RdGy_r"
	bCMap = "inferno"
	pCMap = "viridis"
	dCMap = "viridis"
	vCMap = "bwr"

	doDBZ = True
	doLogP = True
	doMPI = False
	xyBds = [-40,20,-30,30]

	PMin = 0.0 ; PMax = 2.0
	Nc = 11
	vDB = kv.genNorm(25)
	vBB = kv.genNorm(1,250,doLog=True)
	vDD = kv.genNorm(0,25)
	vV  = kv.genNorm(5)
	if (doLogP):
		vP = kv.genNorm(1.0e-3,1.0e+1,doLog=True)
	else:
		vP = kv.genNorm(PMin,PMax)
	

	LW = 0.25
	ashd = 0.25

	#======
	#Init data
	gsph = msph.GamsphPipe(fdir,ftag,doFast=doFast)

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

	if (doDBZ):
		kv.genCB(AxC1,vDB,"Residual Field [nT]",cM=dbCMap,Ntk=7)
	else:
		kv.genCB(AxC1,vBB,"Magnetic Field [nT]",cM=bCMap)
	kv.genCB(AxC4,kv.genNorm(rmpp.cMax),"FAC",cM=rmpp.fcMap,Ntk=5)
	#kv.genCB(AxC3,kv.genNorm(rmpp.pMax),"Potential [kV]",cM=rmpp.pMap,Ntk=7)
	rmpp.AddPotCB(AxC3)

	if (doDen):
		kv.genCB(AxC2,vDD,"Density [#/cc]",cM=dCMap,Ntk=7)
	else:
		if (doLogP):
			kv.genCB(AxC2,vP,"Pressure [nPa]",cM=pCMap)
		else:	
			kv.genCB(AxC2,vP,"Pressure [nPa]",cM=pCMap,Ntk=6)

	AxL.clear()
	AxR.clear()
	if (doDBZ):
		dbz = gsph.DelBz(nStp)
	else:
		MagB = gsph.eqMagB(nStp)
	#Pxy = gsph.EggSlice("P",nStp,vScl=gsph.pScl,doEq=True)

	if (doDen):
		Dxz = gsph.EggSlice("D",nStp,doEq=False)
	else:
		Pxz = gsph.EggSlice("P",nStp,vScl=gsph.pScl,doEq=False)
	if (doDBZ):
		AxL.pcolormesh(gsph.xxi,gsph.yyi,dbz,cmap=dbCMap,norm=vDB)
	else:
		AxL.pcolormesh(gsph.xxi,gsph.yyi,MagB,cmap=bCMap,norm=vBB)
	# if (not doDen):
	# 	AxL.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Pxy),pVals,norm=vP,cmap=pCMap)

	kv.addEarth2D(ax=AxL)
	kv.SetAx(xyBds,AxL)
	gsph.AddTime(nStp,AxL,xy=[0.025,0.935],fs="x-large")
	gsph.AddSW(nStp,AxL,xy=[0.625,0.025],fs="small")
	AxL.set_xlabel('SM-X [Re]')
	AxL.set_ylabel('SM-Y [Re]')

	if (doDen):
		AxR.pcolormesh(gsph.xxi,gsph.yyi,Dxz,cmap=dCMap,norm=vDD)
	else:
		AxR.pcolormesh(gsph.xxi,gsph.yyi,Pxz,cmap=pCMap,norm=vP)

	kv.addEarth2D(ax=AxR)
	kv.SetAx(xyBds,AxR)
	AxR.yaxis.tick_right()
	AxR.yaxis.set_label_position('right')
	if (doIon):
		gsph.AddCPCP(nStp,AxR,xy=[0.65,0.925])

	AxR.set_xlabel('SM-X [Re]')
	AxR.set_ylabel('SM-Z [Re]')


	#Add MPI decomp
	if (doMPI):
		for i in range(gsph.Ri):
			i0 = i*gsph.dNi
			AxL.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],'m',linewidth=LW,alpha=ashd)
			AxR.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],'c',linewidth=LW,alpha=ashd)

	if (doIon):
		fmix = gsph.Gam2Remix(nStp)
		dxy = [32.5,32.5]
		rmpp.CMIViz(AxR,fmix,dxy=dxy)
		rmpp.CMIViz(AxR,fmix,dxy=dxy,loc="lower left",doNorth=False)

	kv.savePic(fOut,bLenX=45)
