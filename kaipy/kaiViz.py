#Various scripts to support visualization of Kaiju data

from kaipy.kdefs import *
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize


#Create 2D equatorial grid (Ni,Nj*2+1) from lfm/egg-style
def genEQGrid(fIn):
	with h5py.File(fIn,'r') as hf:
		Nk,Nj,Ni = hf.get("X")[()].shape
		X3 = hf.get("X")[()].T
		Y3 = hf.get("Y")[()].T
		Z3 = hf.get("Z")[()].T
	print("Ni,Nj,Nk = %d,%d,%d"%(Ni,Nj,Nk))
	#Create grid
	Nr = Ni-1
	Np = 2*(Nj-1)
	xx = np.zeros((Nr+1,Np+1))
	yy = np.zeros((Nr+1,Np+1))

	for i in range(Nr+1):
		for j in range(Nj):
			xx[i,j] = X3[i,j,0]
			yy[i,j] = Y3[i,j,0]
		for j in range(Nj,Np+1):
			jp = Np-j
			xx[i,j] =  X3[i,jp,0]
			yy[i,j] = -Y3[i,jp,0]
	return xx,yy
def getEQVar(fIn,StpN=0,vID="D"):
	gID = "Step#%d"%(StpN)
	with h5py.File(fIn,'r') as hf:
		V = hf.get(gID).get(vID)[()].T
		Ni,Nj,Nk = V.shape #Cell-centered (not node)
	Nr = Ni
	Np = 2*(Nj)
	kp = (Nk)/2
	vv = np.zeros((Nr,Np))
	for i in range(Nr):
		for j in range(Np):
			if (j>=Nj):
				jp = Np-j-1
				vv[i,j] = V[i,jp,kp]
			else:
				jp = j
				vv[i,j] = V[i,jp,0]
	return vv

#Calculate equatorial Bz-D
def getEQBzD(xx,yy,MagM=-EarthM0g):
	Nip,Njp = xx.shape
	Ni = Nip-1
	Nj = Njp-1

	BzD = np.zeros((Ni,Nj))
	xxc = np.zeros((Ni,Nj))
	yyc = np.zeros((Ni,Nj))

	for i in range(Ni):
		for j in range(Nj):
			xc = 0.25*(xx[i,j]+xx[i+1,j]+xx[i,j+1]+xx[i+1,j+1])
			yc = 0.25*(yy[i,j]+yy[i+1,j]+yy[i,j+1]+yy[i+1,j+1])
			xxc[i,j] = xc
			yyc[i,j] = yc
			r = np.sqrt(xc**2.0+yc**2.0)
			rm5 = r**(-5.0)
			BzD[i,j] = -r*r*MagM*rm5
	return xxc,yyc,BzD

#---------------------------------
#Matplotlib helpers

#Set axis bounds w/ aspect ratio
def SetAx(xyBds=[-1,1,-1,1],ax=None,Adj='box'):
	if ax is None:
		ax = plt.gca()
	ax.set_xlim(xyBds[0],xyBds[1])
	ax.set_ylim(xyBds[2],xyBds[3])
	ax.set_aspect('equal',adjustable=Adj)

#Set axis labels and locations
def SetAxLabs(Ax,xLab,yLab,doBot=True,doLeft=True,fs="medium"):
	Ax.set_xlabel(xLab,fontsize=fs)
	Ax.set_ylabel(yLab,fontsize=fs)
	if (not doBot):
		Ax.xaxis.tick_top()
		Ax.xaxis.set_label_position('top')
	if (not doLeft):
		Ax.yaxis.tick_right()
		Ax.yaxis.set_label_position('right')

	#Kill labels if not string is None
	if (xLab is None):
		Ax.xaxis.label.set_visible(False)
		Ax.xaxis.set_visible(False)
		plt.setp(Ax.get_xticklabels(),visible=False)

	if (yLab is None):
		Ax.yaxis.label.set_visible(False)
		Ax.yaxis.set_visible(False)
		plt.setp(Ax.get_yticklabels(),visible=False)

#Set X axis to labels to well formatted date time
def SetAxDate(Ax,fmt='%H:%M \n%Y-%m-%d'):
	Ax.xaxis_date()
	Ax.xaxis.set_major_formatter(mpl.dates.DateFormatter(fmt))

#Adds 2D earth w/ dawn/dusk
def addEarth2D(Re=1, angle=-90, ax=None):
	from matplotlib.patches import Wedge

	if ax is None:
		ax = plt.gca()
	colors=('w','k')
	center = 0
	theta1, theta2 = angle, angle + 180

	w1 = Wedge(center, Re,-90,90,   fc=colors[0],ec='k')
	w2 = Wedge(center, Re, 90, -90, fc=colors[1],ec='k')
	for wedge in [w1, w2]:
		ax.add_artist(wedge)
	return [w1, w2]
#Add inner cutout
def DrawCut(Rin=2.5,ax=None):
	from matplotlib.patches import Wedge
	if ax is None:
		ax = plt.gca()
	w1 = Wedge(0,Rin,0,360,fc=None,ec='k')
	ax.add_artist(w1)
	return w1

#Take cell-centered polar values and add extra phi layer
#Useful for contour plots through +X
def reWrap(V):
	Ni,Nj = V.shape
	Vp = np.zeros((Ni,Nj+1))
	Vp[:,0:Nj] = V
	Vp[:,-1] = V[:,0]
	return Vp

#Image files
#Wrapper to save (and trim) figure
def savePic(fOut,dpiQ=300,doTrim=True,bLenX=20,bLenY=None,doClose=False,doEps=False):
	#Start by saving
	import matplotlib.pyplot as plt
	if (doEps):
		plt.savefig(fOut,dpi=dpiQ,format='eps')
	else:
		plt.savefig(fOut,dpi=dpiQ)
		if (doTrim):
			trimFig(fOut,bLenX,bLenY)
	if (doClose):
		plt.close('all')

#Use imagemagick to trim whitespace off figure
#doEven: Guarantee even number of pixels in X/Y
def trimFig(fName,bLenX=20,bLenY=None,doEven=True):
	import os
	if (bLenY is None):
		bLenY = bLenX

	ComS = 'convert -trim -border %dx%d -bordercolor "#FFFFFF" '%(bLenX,bLenY) + fName + ' ' + fName
	os.system(ComS)

	if (doEven):
		Nx,Ny = picSz(fName)
		#print(Nx,Ny)
		while ( (Nx % 2) != 0 ):
			#print("Shaving X")
			ShaveX(fName)
			Nx,Ny = picSz(fName)
			#print('\t%d,%d'%(Nx,Ny))

		while ( (Ny % 2) != 0 ):
			#print("Shaving Y")
			ShaveY(fName)
			Nx,Ny = picSz(fName)
			#print('\t%d,%d'%(Nx,Ny))

		Nx,Ny = picSz(fName)
		if ( ((Nx % 2) != 0) or ((Ny % 2) != 0) ):
			print("Parity failure on pic sizing")
			print(Nx,Ny)

def picSz(fName):
	from PIL import Image
	with Image.open(fName) as img:
		Nx,Ny = img.size
	return Nx,Ny

def ShaveX(fName):
	import os
	ComS = 'convert -crop -1+0 +repage ' + fName + ' ' + fName
	os.system(ComS)
def ShaveY(fName):
	import os
	ComS = 'convert -crop +0-1 +repage ' + fName + ' ' + fName
	os.system(ComS)

#---------------------------------
#Create colorbar with specified midpoint (grabbed from stack overflow)
class MidpointNormalize(Normalize):
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y))

#Create norm object for MPL
def genNorm(vMin,vMax=None,doLog=False,midP=None):
	from matplotlib.colors import LogNorm
	from matplotlib.colors import Normalize
	if (vMax is None):
		vMin = -np.abs(vMin)
		vMax = np.abs(vMin)
	if (midP is None):
		doMid = False
	else:
		doMid = True

	if (doMid):
		vN = MidpointNormalize(vmin=vMin,vmax=vMax,midpoint=midP)
	elif (doLog):
		vN = LogNorm(vmin=vMin,vmax=vMax)
	else:
		vN = Normalize(vmin=vMin,vmax=vMax)

	return vN

#Create colorbar object into specified axis
def genCB(AxCB,vN,cbT="Title",cM="viridis",doVert=False,cbSz="medium",Ntk=None):
	from matplotlib import ticker

	if (doVert):
		cbOr = "vertical"
	else:
		cbOr = "horizontal"
	cmData = plt.get_cmap(cM)
	cb = mpl.colorbar.ColorbarBase(AxCB,cmap=cmData,norm=vN,orientation=cbOr)
	if (Ntk is not None):
		cb.locator = ticker.MaxNLocator(nbins=Ntk)
		cb.update_ticks()

	cb.set_label(cbT,fontsize=cbSz)
	cb.ax.tick_params(labelsize=cbSz)
	return cb

def labelStr(data, key, vecComp):
	vecLabel=[ 'x', 'y', 'z' ]
	if (vecComp < 3) and (vecComp > -1):
		label=(data['GAMERA_'+key].attrs['AXISLABEL']+
		vecLabel[vecComp]+' ['+
		data['GAMERA_'+key].attrs['UNITS'].decode()+']')
	else:
		label=(data['GAMERA_'+key].attrs['AXISLABEL']+
		' ['+data['GAMERA_'+key].attrs['UNITS'].decode()+']')
	return label


def itemPlot(Ax,data,key,plotNum,numPlots,vecComp=-1):
	#print(key,vecComp)
	if -1 == vecComp:
		maskedData = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
			data[key][:])
		Ax.plot(data['Epoch_bin'],maskedData)
		maskedGamera = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
			data['GAMERA_'+key][:])
		Ax.plot(data['Epoch_bin'],maskedGamera)
	else:
		maskedData = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
			data[key][:,vecComp])
		Ax.plot(data['Epoch_bin'],maskedData)
		maskedGamera = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
			data['GAMERA_'+key][:,vecComp])
		Ax.plot(data['Epoch_bin'],maskedGamera)
	if (plotNum % 2) == 0:
		left = True
	else:
		left = False
	label = labelStr(data, key,vecComp)
	if plotNum == (numPlots-1):
		SetAxLabs(Ax,'UT',label,doLeft=left)
		SetAxDate(Ax)
	else:
		SetAxLabs(Ax,None,label,doLeft=left)
	return

def compPlot(plotname,scId,data):

	numPlots = 0
	keysToPlot = []
	keys = data.keys()
	#print(keys)
	if 'Density' in keys:
		numPlots = numPlots + 1
		keysToPlot.append('Density')
	if 'Pressue' in keys:
		numPlots = numPlots + 1
		keysToPlot.append('Pressue')
	if 'Temperature' in keys:
		numPlots = numPlots + 1
		keysToPlot.append('Temperature')
	if 'MagneticField' in keys:
		numPlots = numPlots + 3
		keysToPlot.append('MagneticField')
	if 'Velocity' in keys:
		numPlots = numPlots + 3
		keysToPlot.append('Velocity')

	figsize = (10,10)
	fig = plt.figure(figsize=figsize)
	gs = fig.add_gridspec(numPlots,1)
	plotNum = 0
	for key in keysToPlot:
		#print('Plotting',key)
		if 'MagneticField' == key or 'Velocity' == key:
			doVecPlot = True
		else:
			doVecPlot = False
		if 0 == plotNum:
			Ax1 = fig.add_subplot(gs[plotNum,0])
			if doVecPlot:
				#print('key',key,'plotNum',plotNum)
				itemPlot(Ax1,data,key,plotNum,numPlots,vecComp=0)
				plotNum = plotNum + 1
				Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=1)
				plotNum = plotNum + 1
				Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=2)
				plotNum = plotNum + 1
			else:
				#print('key',key,'plotNum',plotNum)
				itemPlot(Ax1,data,key,plotNum,numPlots)
				plotNum = plotNum + 1
		else:
			Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
			if doVecPlot:
				#print('key',key,'plotNum',plotNum)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=0)
				plotNum = plotNum + 1
				Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=1)
				plotNum = plotNum + 1
				Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=2)
				plotNum = plotNum + 1
			else:
				#print('key',key,'plotNum',plotNum)
				itemPlot(Ax,data,key,plotNum,numPlots)
				plotNum = plotNum + 1

	Ax1.legend([scId,'GAMERA'],loc='best')
	Ax1.set_title(plotname)
	plt.subplots_adjust(hspace=0)

	savePic(plotname)
