#Various helper routines for magnetosphere quick look plots

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
import numpy as np
import kaipy.kaiViz as kv
import kaipy.kaiTools as kt
import kaipy.gamera.magsphere as msph
import kaipy.remix.remix as remix
import os

dbMax = 25.0
dbCM = "RdGy_r"
bzCM = "bwr"
cLW = 0.25
bz0Col = "magenta"
mpiCol = "deepskyblue"

jMax = 10.0 #Max current for contours

eMax = 5.0  #Max current for contours

#Default pressure colorbar
vP = kv.genNorm(vMin=1.0e-2,vMax=10.0,doLog=True)
szStrs = ['small','std','big','dm']
szBds = {}
szBds["std"]      = [-40.0 ,20.0,2.0]
szBds["big"]      = [-100.0,20.0,2.0]
szBds["bigger"]   = [-200.0,25.0,2.0]
szBds["small"]    = [-10.0 , 5.0,2.0]
szBds["dm"]       = [-30.0 ,10.0,40.0/15.0]

#Add different size options to argument
def AddSizeArgs(parser):
	parser.add_argument('-size',type=str,default="std",choices=szStrs,help="Domain bounds options (default: %(default)s)")

#Return domain size from parsed arguments
def GetSizeBds(args):

	szStr = args.size
	szBd = szBds[szStr]

	xTail = szBd[0]
	xSun  = szBd[1]
	yMax = (xSun-xTail)/szBd[2]
	xyBds = [xTail,xSun,-yMax,yMax]

	return xyBds

#Plot equatorial field
def PlotEqB(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True,doBz=False):
	vBZ = kv.genNorm(dbMax)
	vDB = kv.genNorm(dbMax)

	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		if (doBz):
			kv.genCB(AxCB,vBZ,"Vertical Field [nT]",cM=bzCM,Ntk=7)
		else:
			kv.genCB(AxCB,vDB,"Residual Field [nT]",cM=dbCM,Ntk=7)
	#Now do main plotting
	if (doClear):
		Ax.clear()
	Bz = gsph.EggSlice("Bz",nStp,doEq=True)
	if (doBz):
		Ax.pcolormesh(gsph.xxi,gsph.yyi,Bz,cmap=bzCM,norm=vBZ)
	else:
		dbz = gsph.DelBz(nStp)
		Ax.pcolormesh(gsph.xxi,gsph.yyi,dbz,cmap=dbCM,norm=vDB)
	Ax.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors=bz0Col,linewidths=cLW)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM-X [Re]')
		Ax.set_ylabel('SM-Y [Re]')
	return Bz

def PlotMerid(gsph,nStp,xyBds,Ax,doDen=False,doRCM=False,AxCB=None,doClear=True,doDeco=True,doSrc=False):
	CMx = "viridis"
	if (doDen):
		
		if (doRCM):
			vN = kv.genNorm(vMin=1.0,vMax=1.0e+3,doLog=True)
		else:
			vN = kv.genNorm(0,25)
		if (doSrc):
			vID = "SrcD"
			cbStr = "Source Density [#/cc]"
			
		else:
			vID = "D"
			cbStr = "Density [#/cc]"
		Q = gsph.EggSlice(vID,nStp,doEq=False)
	else:
		vN = vP
		if (doSrc):
			vID = "SrcP"
			cbStr = "Source Pressure [nPa]"
		else:
			vID = "P"
			cbStr = "Pressure [nPa]"
		Q = gsph.EggSlice(vID,nStp,doEq=False)
	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vN,cbStr,cM=CMx)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,Q,cmap=CMx,norm=vN)

	kv.SetAx(xyBds,Ax)
	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM-X [Re]')
		Ax.set_ylabel('SM-Z [Re]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')

def PlotJyXZ(gsph,nStp,xyBds,Ax,AxCB=None,jScl=None,doDeco=True):
	if (jScl is None):
		#Just assuming current scaling is nA/m2
		jScl = 1.0
	vJ = kv.genNorm(jMax)
	jCMap = "PRGn"
	Nc = 15
	cVals = np.linspace(-jMax,jMax,Nc)

	if (AxCB is not None):
		AxCB.clear()
		kv.genCB(AxCB,vJ,"Jy [nA/m2]",cM=jCMap)
	Q = jScl*gsph.EggSlice("Jy",nStp,doEq=False)
	#Zero out first shell b/c bad derivative
	print(Q.shape)
	Q[0:2,:] = 0.0
	#Ax.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Q),cVals,norm=vJ,cmap=jCMap,linewidths=cLW)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,Q,norm=vJ,cmap=jCMap)
	kv.SetAx(xyBds,Ax)
	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM-X [Re]')
		Ax.set_ylabel('SM-Z [Re]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')

#Plot equatorial azimuthal electric field
def PlotEqEphi(gsph,nStp,xyBds,Ax,AxCB=None,doClear=True,doDeco=True):
	vE = kv.genNorm(eMax)
	vEMap = "PRGn"
	if (AxCB is not None):
		#Add the colorbar to AxCB
		AxCB.clear()
		kv.genCB(AxCB,vE,r"E$_{phi}$ [mV/m]",cM=vEMap)

	#Now do main plotting
	if (doClear):
		Ax.clear()
	Bx = gsph.EggSlice("Bx",nStp,doEq=True)
	By = gsph.EggSlice("By",nStp,doEq=True)
	Bz = gsph.EggSlice("Bz",nStp,doEq=True)
	Vx = gsph.EggSlice("Vx",nStp,doEq=True)
	Vy = gsph.EggSlice("Vy",nStp,doEq=True)
	Vz = gsph.EggSlice("Vz",nStp,doEq=True)

	# calculating some variables to to plot
	#E=-VxB
	Ex = -(Vy*Bz-Vz*By)*0.001 # [mV/m]
	Ey =  (Vx*Bz-Vz*Bx)*0.001
	Ez = -(Vx*By-Vy*Bx)*0.001

	# coordinate transform
	ppc = np.arctan2(gsph.yyc,gsph.xxc)
	theta = np.pi #eq plane
	Er,Et,Ep = kt.xyz2rtp(ppc,theta,Ex,Ey,Ez)

	Ax.pcolormesh(gsph.xxi,gsph.yyi,Ep,cmap=vEMap,norm=vE)

	kv.SetAx(xyBds,Ax)

	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM-X [Re]')
		Ax.set_ylabel('SM-Y [Re]')
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')

#Add MPI contours
def PlotMPI(gsph,Ax,ashd=0.5):
	gCol = mpiCol
	for i in range(gsph.Ri):
		i0 = i*gsph.dNi
		Ax.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],mpiCol,linewidth=cLW,alpha=ashd)

	if (gsph.Rj>1):
		for j in range(1,gsph.Rj):
			j0 = j*gsph.dNj
			Ax.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)
			Ax.plot(gsph.xxi[:,j0],-gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)
		#X-axis (+)
		Ax.plot(gsph.xxi[:,0], gsph.yyi[:,0],gCol,linewidth=cLW,alpha=ashd)
		#X-axis (-)
		j0 = (gsph.Rj)*gsph.dNj
		Ax.plot(gsph.xxi[:,j0], gsph.yyi[:,j0],gCol,linewidth=cLW,alpha=ashd)

def AddIonBoxes(gs,ion):
	gsRM = gs.subgridspec(20,20)

	wXY = 6
	dX = 1
	dY = 1

	#Northern
	ion.init_vars('NORTH')
	ax = ion.plot('current'  ,gs=gsRM[dY:dY+wXY,dX:dX+wXY],doInset=True)
	#ax.set_title("North",fontsize="small")

	#Southern
	ion.init_vars('SOUTH')
	ax = ion.plot('current'  ,gs=gsRM[-dY-wXY:-dY,dX:dX+wXY],doInset=True)

def plotPlane(gsph,data,xyBds,Ax,AxCB,var='D',vMin=None,vMax=None,doDeco=True,cmap='viridis',doLog=False,midp=None):
	if (AxCB is not None):
		AxCB.clear()
	if (not midp):
		if (vMin is None):
 			vMin = np.min(data)
		if (vMax is None):
			vMax = np.max(data)
	else:
		if ((vMin is None) and (vMax is None)):
			vMax = np.max(np.abs([np.min(data),np.max(data)]))
			vMin = -1.0*vMax

	vNorm = kv.genNorm(vMin,vMax=vMax,doLog=doLog,midP=midp)
	kv.genCB(AxCB,vNorm,cbT=var,cM=cmap,Ntk=7)
	Ax.pcolormesh(gsph.xxi,gsph.yyi,data,cmap=cmap,norm=vNorm)
	kv.SetAx(xyBds,Ax)

	return

def plotXY(gsph,nStp,xyBds,Ax,AxCB,var='D',vMin=None,vMax=None,doDeco=True,cmap='viridis',doLog=False,midp=None):
	data = gsph.EggSlice(var,nStp,doEq=True)
	plotPlane(gsph,data,xyBds,Ax,AxCB,var,vMin=vMin,vMax=vMax,doDeco=doDeco,cmap=cmap,doLog=doLog,midp=midp)
	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM_X [Re]')
		Ax.set_ylabel('SM-Y [Re]')

	return data

def plotXZ(gsph,nStp,xzBds,Ax,AxCB,var='D',vMin=None,vMax=None,doDeco=True,cmap='viridis',doLog=False,midp=None):
	data = gsph.EggSlice(var,nStp,doEq=False)
	plotPlane(gsph,data,xzBds,Ax,AxCB,var,vMin=vMin,vMax=vMax,doDeco=doDeco,cmap=cmap,doLog=doLog,midp=midp)
	if (doDeco):
		kv.addEarth2D(ax=Ax)
		Ax.set_xlabel('SM_X [Re]')
		Ax.set_ylabel('SM-Z [Re]')

	return data
