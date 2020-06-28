#Various routines for plotting remix output data
#Generally dependant on basemap

import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import matplotlib.pyplot as plt
import h5py
import kaipy.kaiViz as kv

r2deg = 180.0/np.pi
pType = 'npstere'
#pType = 'npaeqd'
#pType = 'nplaea'
cMax = 2.0 #Max current density

dfs = "medium"
fcMap = "bwr" #Colormap for current

#Potential plots
pMax = 60.0 #Max potential [kV]
doAPC = True #Do solid/dash for potential contours
dpC = 5.0 #Spacing [kV] for potential contours
#pMap = "PRGn"
pMax = 35.0 #Max potential [kV]
apMap = "bone" #Absolute value of potential map
pMap = apMap

NumP = 10
pVals = np.linspace(-pMax,pMax,NumP)
cLW = 0.5
cAl = 0.5 #Contour alpha

lon0 = -90
R0 = 6.5/6.38

dLon = 45.0 #Degree spacing for longitude
dLat = 10.0 #Degree spacing for latitude

AxEC = 'silver'
gC = 'dimgrey'

gDash = [1,0]
gLW = 0.125

#Adds inset figure to axis using data from fmix
#dxy = [width%,height%]
def CMIPic(nLat,nLon,llBC,P,C,AxM=None,doNorth=True,loc="upper left",dxy=[20,20]):
	from mpl_toolkits.basemap import Basemap
	if (doNorth):
		tStr = "North"
	else:
		tStr = "South"
		
	vC = kv.genNorm(cMax)
	vP = kv.genNorm(pMax)

	#Now add inset figure to passed axis
	wStr = "%f%%"%(dxy[0])
	hStr = "%f%%"%(dxy[1])

	if (AxM is None):
		AxM = plt.gca()

	Ax = inset_axes(AxM,width=wStr,height=hStr,loc=loc)
	
	for spine in Ax.spines.values():
		spine.set_edgecolor(AxEC)
	Ax.patch.set_alpha(0)
	if (doNorth):
		Ax.set_xlabel(tStr,fontsize=dfs,color=AxEC)
	else:
		Ax.set_title(tStr,fontsize=dfs,color=AxEC)

	
	bmAx = Basemap(ax=Ax,projection=pType,boundinglat=llBC,lon_0=lon0)

	Lat0 = np.linspace(90,llBC,nLat+0)
	Lat1 = np.linspace(90,llBC,nLat+1)
	
	if (doNorth):
		Lon0 = np.linspace(0,360,nLon+0)
		Lon1 = np.linspace(0,360,nLon+1)
	else:
		Lon0 = np.linspace(360,0,nLon+0)
		Lon1 = np.linspace(360,0,nLon+1)
	LonC,LatC = np.meshgrid(Lon0,Lat0)
	LonI,LatI = np.meshgrid(Lon1,Lat1)
	
	#Now do plotting, start w/ gridding
	
	#Set parallels
	gP = np.arange(90-dLat,0.95*llBC,-dLat)
	bmAx.drawparallels(gP,latmax=gP.max(),dashes=gDash,linewidth=gLW,color=gC)
	
	#Set meridians
	gM = np.arange(0,360,dLon)
	for lon in gM:
		bmAx.drawgreatcircle(lon,gP.max(),lon,llBC,linewidth=gLW,color=gC)
	
	#Plot data
	bmAx.pcolormesh(LonI,LatI,C,latlon=True,norm=vC,cmap=fcMap)
	#Do potential contours
	if (doAPC):
		#Do positive contours
		vAPp = kv.genNorm(0,pMax)
		pVals = np.arange(dpC,pMax,dpC)
		if (P.max()>dpC):
			bmAx.contour(LonC,LatC,P,pVals,latlon=True,norm=vAPp,cmap=apMap,alpha=cAl,linewidths=cLW,linestyles='solid')
		#Do negative contours
		vAPm = kv.genNorm(-pMax,0)
		apMapm = apMap+"_r"
		pVals = -pVals[::-1]
		if (P.min()<-dpC):
			bmAx.contour(LonC,LatC,P,pVals,latlon=True,norm=vAPm,cmap=apMapm,alpha=cAl,linewidths=cLW,linestyles='dashed')
	else:
		bmAx.contour(LonC,LatC,P,pVals,latlon=True,norm=vP,cmap=pMap,alpha=cAl,linewidths=cLW)
	
def AddPotCB(Ax,Lab="Potential [kV]",Ntk=7):
	if (doAPC):
		vP = kv.genNorm(0,pMax)
		cm = apMap
	else:
		vP = kv.genNorm(-pMax,pMax)
		cm = pMap
	kv.genCB(Ax,vP,Lab,cM=cm,Ntk=Ntk)

def AddCBs(Ax1,Ax2,Lab1="Potential [kV]",Lab2="FAC",Ntk1=7,Ntk2=5,doFlip=True):
	if (doAPC):
		vP = kv.genNorm(0,pMax)
		cm = apMap
	else:
		vP = kv.genNorm(-pMax,pMax)
		cm = pMap
	kv.genCB(Ax1,vP,Lab1,cM=cm,Ntk=Ntk1)
	kv.genCB(Ax2,kv.genNorm(cMax),Lab2,fcMap,Ntk=Ntk2)
	if (doFlip):
		Ax2.xaxis.tick_top()
		Ax2.xaxis.set_label_position('top')
