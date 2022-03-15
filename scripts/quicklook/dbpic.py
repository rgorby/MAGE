#!/usr/bin/env python
#Make a quick figure showing ground dB

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
import kaipy.kaiTools as ktools
import matplotlib.gridspec as gridspec
import numpy as np
import kaipy.gamera.gampp as gampp
import kaipy.kaiH5 as kh5
import kaipy.cmaps.kaimaps as kmaps
import sys
import os
import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

if __name__ == "__main__":
	rad2deg = 180.0/np.pi
	bMag = 1000.0
	bLin = 10.0
	jMag = 0.5

	#Defaults
	fdir = os.getcwd()
	ftag = "msphere"
	nStp = -1
	fOut = "qkdbpic.png"
	doJr = False

	k0 = 0 #Vertical slice to use
	MainS = """Creates visualization of ground dB
	NOTE: Assumes ground dB has been calculated using calcdb.x on simulation data.
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-n' ,type=int,metavar="step" ,default=nStp,help="Time slice to plot (default: %(default)s)")
	parser.add_argument('-k0',type=int,metavar="layer" ,default=k0,help="Vertical layer to plot (default: %(default)s)")
	parser.add_argument('-Jr', action='store_true', default=doJr ,help="Show radial component of anomalous current (default: %(default)s)")


	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id + ".deltab"
	nStp = args.n
	k0   = args.k0
	doJr = args.Jr

	#======
	#Init data
	fname = fdir + "/" + ftag + ".h5"
	dbdata = gampp.GameraPipe(fdir,ftag)
	print("---")
	#Get coordinates
	cStr = kh5.PullAtt(fname,"CoordinatesID",s0=None)
	Re   = kh5.PullAtt(fname,"Re",s0=None) #km
	try:
		CoordID = cStr.decode('utf-8')
	except (UnicodeDecodeError, AttributeError):
		CoordID = cStr
	print("Found %s coordinate data ..."%(CoordID))
	if (nStp<0):
		nStp = dbdata.sFin
		print("Using Step %d"%(nStp))

	#Pull data
	dBn = dbdata.GetVar("dBn",nStp)
	Rcc = dbdata.GetVar("Radcc",doVerb=False)[:,:,k0]
	#Get smlat/lon
	smlon = dbdata.GetVar("smlon",nStp,doVerb=False)[:,:,k0]
	smlat = dbdata.GetVar("smlat",nStp,doVerb=False)[:,:,k0]

	NLat,NLon,Nz = dBn.shape
	dBn = dBn[:,:,k0]
	if (k0 >= Nz):
		print("Invalid vertical level, only %d levels found!"%(Nz))
		sys.exit()
	else:
		Z0 = (Rcc.mean()*Re)-Re
		print("Height = %6.2f [km]"%(Z0))
	if (doJr):
		print("Reading Jr ...")
		Jr = dbdata.GetVar("dbJ",nStp,doVerb=False)[:,:,k0]
		Q = Jr
	else:
		Q = dBn

	print(Q.min(),Q.max())
	#Get MJD to UT
	MJD = kh5.tStep(fname,nStp,aID="MJD")
	utS = ktools.MJD2UT([MJD])
	utDT= utS[0]
	
	#=====
	#Do cartopy stuff
	crs = ccrs.PlateCarree()
	
	X0 = dbdata.X[:,:,k0]
	Y0 = dbdata.Y[:,:,k0]
	Z0 = dbdata.Z[:,:,k0]
	R0 = np.sqrt(X0**2.0 + Y0**2.0 + Z0**2.0)
	LonI = np.arctan2(Y0,X0)*rad2deg
	LatI = np.arcsin(Z0/R0) *rad2deg

	LonI = np.linspace(0,360,NLon+1)
	LatI = LatI[:,0]

	LonC = 0.5*(LonI[0:-1] + LonI[1:])
	LatC = 0.5*(LatI[0:-1] + LatI[1:])

	#=====
	#Do figure stuff
	cmap = kmaps.cmDiv
	if (doJr):
		vQ = kv.genNorm(jMag)
		cbStr = "Anomalous current"
	else:
		vQ = kv.genNorm(bMag,doSymLog=True,linP=bLin)
		cbStr = r"$\Delta B_N$ [nT]"
	figSz = (12,6)
	fig = plt.figure(figsize=figSz)
	gs = gridspec.GridSpec(3,1,height_ratios=[20,1.0,1.0],hspace=0.025)

	AxM  = fig.add_subplot(gs [0,0],projection=crs)
	AxCB = fig.add_subplot(gs[-1,0])

	#Make plots
	AxM.pcolormesh(LonI,LatI,Q,norm=vQ,cmap=cmap)
	kv.genCB(AxCB,vQ,cbStr,cM=cmap)

	#Add decoration
	tStr = utDT.strftime("%m/%d/%Y, %H:%M:%S")
	AxM.set_title(tStr,fontsize="x-large")

	AxM.add_feature(Nightshade(utDT,alpha=0.10))
	AxM.coastlines(resolution='110m', color='black', linewidth=0.25)
	
	gl = AxM.gridlines(crs, draw_labels=True,linewidth=1, color='gray', alpha=0.5, linestyle='--')
	gl.top_labels = False
	gl.left_labels = False
	gl.xlines = False
	gl.ylines = False
	gl.xlocator = mticker.FixedLocator([-120, 0, 120])
	gl.ylocator = mticker.FixedLocator([-45, 0, 45])
	AxM.contour(LonC,LatC,smlon,[0.5,90,180,270],linewidths=1.25, linestyles='--',alpha=0.5,colors='grey')

	kv.savePic(fOut)
	