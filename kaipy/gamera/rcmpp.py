#Various routines for plotting RCM output data
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import matplotlib.pyplot as plt
import h5py
import kaipy.kaiViz as kv
import kaipy.gamera.gampp as gampp
import os
import numpy.ma as ma
import matplotlib.patches as patches
import matplotlib.ticker as plticker

rcBds = [-15,10.0,-12.5,12.5]
pCMap = "viridis"
eCol = "slategrey"
rcmCol = "dodgerblue"
gCol = "cyan"
gLW = 0.15

doXYZ = True #How to calculate "equator"
doCut = False
pCut = 1.0e-8
MHDCol = "red"
eLW = 0.05
MHDLW = 0.5

doEll = True
#Get equatorial coordinates, masked if asked to
def RCMEq(rcmdata,nStp,doMask=False,doXYZ=doXYZ):

	bmX = rcmdata.GetVar("xMin",nStp)
	bmY = rcmdata.GetVar("yMin",nStp)
	if (doXYZ):
		bmZ = rcmdata.GetVar("zMin",nStp)
		bmP = np.arctan2(bmY,bmX)
		bmR = np.sqrt(bmX*bmX + bmY*bmY + bmZ*bmZ)
		bmX = bmR*np.cos(bmP)
		bmY = bmR*np.sin(bmP)
	if (doMask):
		I = GetMask(rcmdata,nStp)
		bmX = ma.masked_array(bmX,mask=I)
		bmY = ma.masked_array(bmY,mask=I)
	return bmX,bmY

def GetVarMask(rcmdata,nStp,Qid="P",I=None):
	if (I is None):
		I = GetMask(rcmdata,nStp)
	Q = rcmdata.GetVar(Qid,nStp)
	Q = ma.masked_array(Q,mask=I)
	return Q

#Calculate mask
#doRCM: Do RCM domain or full closed region
def GetMask(rcmdata,nStp):
	IOpen = rcmdata.GetVar("IOpen",nStp)
	
	if (doEll):
		ioCut = -0.5
	else:
		ioCut = 0.5

	if (doCut):
		Prcm = rcmdata.GetVar("P",nStp)
		I = (IOpen > ioCut) | (Prcm<pCut)
	else:
		I = (IOpen > ioCut)
	return I

#Take axis and rcmdata object and add pressure plot
def RCMInset(AxRCM,rcmdata,nStp,vP):
	if (AxRCM is None):
		AxRCM = plt.gca()

	bmX,bmY = RCMEq(rcmdata,nStp,doMask=True)
	I = GetMask(rcmdata,nStp)
	Ni = (~I).sum()

	if (Ni == 0):
		return

	Prcm  = GetVarMask(rcmdata,nStp,"P"    ,I)
	toMHD = GetVarMask(rcmdata,nStp,"toMHD",I)

	#Start plotting
	AxRCM.pcolor(bmX,bmY,Prcm,norm=vP,cmap=pCMap)
	AxRCM.plot(bmX,bmY,color=eCol,linewidth=eLW)
	AxRCM.plot(bmX.T,bmY.T,color=eCol,linewidth=eLW)

	doCon = (nStp>0) and (toMHD.min()<0.5) and (toMHD.max()>0.5)
	#Add MHD ingestion contour
	if (doCon):
		CS1 = AxRCM.contour(bmX,bmY,toMHD,[0.5],colors=MHDCol,linewidths=MHDLW)
		manloc = [(0.0,8.0)]

		fmt = {}
		fmt[0.5] = 'MHD'
		AxRCM.clabel(CS1,CS1.levels[::2],inline=True,fmt=fmt,fontsize=5,inline_spacing=25,manual=manloc)

	#Add grid
	xS = [-10,-5,0,5]
	yS = [-10,-5,0,5,10]
	for x in xS:
		AxRCM.axvline(x,linewidth=gLW,color=gCol)
	for y in yS:
		AxRCM.axhline(y,linewidth=gLW,color=gCol)

	kv.addEarth2D(ax=AxRCM)
	kv.SetAx(rcBds,AxRCM)
	kv.SetAxLabs(AxRCM,xLab=None,yLab=None)
	AxRCM.spines['bottom'].set_color(rcmCol)
	AxRCM.spines['top'].set_color(rcmCol) 
	AxRCM.spines['right'].set_color(rcmCol)
	AxRCM.spines['left'].set_color(rcmCol)

	AxRCM.set_title("RCM Pressure",fontsize="small",color=rcmCol)

#Add RCM box to other plot
def AddRCMBox(Ax):
	if (Ax is None):
		Ax = plt.gca()
	xy0 = (rcBds[0],rcBds[2])
	H = rcBds[3]-rcBds[2]
	W = rcBds[1]-rcBds[0]
	rcmRec = patches.Rectangle( (xy0),W,H, fill=False,edgecolor=rcmCol)
	Ax.add_patch(rcmRec)
