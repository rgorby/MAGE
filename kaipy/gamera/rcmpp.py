#Various routines for plotting RCM output data
import numpy as np
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import matplotlib.pyplot as plt
import h5py
import kaipy.kaiViz as kv
import kaipy.gamera.gampp as gampp
import os
import numpy.ma as ma
import matplotlib.patches as patches

rcBds = [-15,10.0,-12.5,12.5]
pCMap = "viridis"
eCol = "slategrey"
rcmCol = "darkcyan"
eLW = 0.05

#Take axis and rcmdata object and add pressure plot
def RCMInset(AxRCM,rcmdata,nStp,vP):
	if (AxRCM is None):
		AxRCM = plt.gca()

	bmX = rcmdata.GetVar("xMin",nStp)
	bmY = rcmdata.GetVar("yMin",nStp)
	Prcm = rcmdata.GetVar("P",nStp)
	IOpen = rcmdata.GetVar("IOpen",nStp)
	I = (IOpen > -0.5)
	Ni = (~I).sum()

	if (Ni == 0):
		return
	#If still here we got something to show

	#Do masks
	bmX = ma.masked_array(bmX,mask=I)
	bmY = ma.masked_array(bmY,mask=I)
	Prcm = ma.masked_array(Prcm,mask=I)

	#Start plotting
	AxRCM.pcolor(bmX,bmY,Prcm,norm=vP,cmap=pCMap)
	AxRCM.plot(bmX,bmY,color=eCol,linewidth=eLW)
	AxRCM.plot(bmX.T,bmY.T,color=eCol,linewidth=eLW)
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