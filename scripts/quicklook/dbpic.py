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
import kaipy.gamera.deltabViz as dbViz
import sys
import os
import cartopy.crs as ccrs

if __name__ == "__main__":
	rad2deg = 180.0/np.pi
	bMag = dbViz.dbMag
	bLin = dbViz.dbLin
	jMag = dbViz.jMag

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
	CoordID,Re = dbViz.GetCoords(fname)
	print("Found %s coordinate data ..."%(CoordID))

	if (nStp<0):
		nStp = dbdata.sFin
		print("Using Step %d"%(nStp))

	#Check vertical level
	Z0 = dbViz.CheckLevel(dbdata,k0,Re)

	#Pull data
	if (doJr):
		print("Reading Jr ...")
		Jr = dbdata.GetVar("dbJ",nStp,doVerb=False)[:,:,k0]
		Q = Jr
	else:
		dBn = dbdata.GetVar("dBn",nStp,doVerb=False)[:,:,k0]
		Q = dBn

	#Get MJD to UT
	MJD = kh5.tStep(fname,nStp,aID="MJD")
	utS = ktools.MJD2UT([MJD])
	utDT= utS[0]
	
	#=====
	#Do cartopy stuff
	crs = ccrs.PlateCarree()
	LatI,LonI,LatC,LonC = dbViz.GenUniformLL(dbdata,k0)

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
	tStr = dbViz.GenTStr(AxM,fname,nStp)
	dbViz.DecorateDBAxis(AxM,crs,utDT)
	kv.savePic(fOut)
	