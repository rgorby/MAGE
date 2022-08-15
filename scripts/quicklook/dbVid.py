#!/usr/bin/env python
#Make video of Gamera magnetosphere run
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

	#Defaults
	fdir = os.getcwd()
	ftag = "msphere"
	oDir = "vid2D"
	k0 = 0 #Vertical slice to use
	nS = 0
	nE = -1
	Nblk = 1 #Number of blocks
	nID = 1 #Block ID of this job

	MainS = """Creates visualization of ground dB
	NOTE: Assumes ground dB has been calculated using calcdb.x on simulation data.
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="directory",default=oDir,help="Subdirectory to write to (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-nS' ,type=int,metavar="Step-Start",default=nS,help="Starting step (default: %(default)s)")
	parser.add_argument('-nE' ,type=int,metavar="Step-End"  ,default=nE,help="Ending step   (default: %(default)s)")
	parser.add_argument('-k0',type=int,metavar="layer" ,default=k0,help="Vertical layer to plot (default: %(default)s)")

	parser.add_argument('-Nblk' ,type=int,metavar="Nblk",default=Nblk,help="Number of job blocks (default: %(default)s)")
	parser.add_argument('-nID' ,type=int,metavar="nID"  ,default=nID,help="Block ID of this job [1-Nblk] (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id + ".deltab"
	nS = args.nS
	nE = args.nE
	k0   = args.k0
	oSub = args.o
	Nblk = args.Nblk
	nID = args.nID

	#======
	#Init data
	fname = fdir + "/" + ftag + ".h5"
	dbdata = gampp.GameraPipe(fdir,ftag)
	print("---")
	#Get coordinates
	CoordID,Re = dbViz.GetCoords(fname)
	print("Found %s coordinate data ..."%(CoordID))
	#Check vertical level
	Z0 = dbViz.CheckLevel(dbdata,k0,Re)

	#Set step bounds
	if (nE<0):
		nE = dbdata.sFin
	if (nS<dbdata.s0):
		nS = dbdata.s0

	#Setup parallel in time stuff
	vO = np.arange(nS,nE+1)
	Nt = len(vO)
	print("Writing %d outputs between minutes %d and %d"%(Nt,nS,nE))

	if (Nblk>1):
		#Figure out work bounds
		dI = (Nt//Nblk)
		i0 = (nID-1)*dI
		i1 = i0+dI
		if (nID == Nblk):
			i1 = Nt #Make sure we get last bit
		print("\tBlock #%d: %d to %d"%(nID,i0,i1))
	else:
		i0 = 0
		i1 = Nt

	#Setup output directory
	oDir = fdir + "/" + oSub
	print("Writing output to %s"%(oDir))

	#Check/create directory if necessary
	if (not os.path.exists(oDir)):
		try:
			print("Creating directory %s"%(oDir))
			os.makedirs(oDir)
		except OSError as exc:
			if exc.errno == errno.EEXIST and os.path.isdir(oDir):
				pass
			else:
				raise

	#=====
	#Do cartopy stuff
	crs = ccrs.PlateCarree()
	LatI,LonI,LatC,LonC = dbViz.GenUniformLL(dbdata,k0)

	#=====
	#Do figure stuff
	cmap = kmaps.cmDiv
	vQ = kv.genNorm(bMag,doSymLog=True,linP=bLin)
	cbStr = r"$\Delta B_N$ [nT]"
	figSz = (12,6)
	fig = plt.figure(figsize=figSz)
	gs = gridspec.GridSpec(3,1,height_ratios=[20,1.0,1.0],hspace=0.025)

	AxM  = fig.add_subplot(gs [0,0],projection=crs)
	AxCB = fig.add_subplot(gs[-1,0])

	kv.genCB(AxCB,vQ,cbStr,cM=cmap)

	#Loop over sub-range
	for i in range(i0,i1):
		nStp = vO[i]
		AxM.clear()
		
		Q = dbdata.GetVar("dBn",nStp,doVerb=False)[:,:,k0]

		#Get MJD to UT
		MJD = kh5.tStep(fname,nStp,aID="MJD")
		utS = ktools.MJD2UT([MJD])
		utDT= utS[0]

		#Do plot
		AxM.pcolormesh(LonI,LatI,Q,norm=vQ,cmap=cmap)

		#Add decoration
		tStr = dbViz.GenTStr(AxM,fname,nStp)
		dbViz.DecorateDBAxis(AxM,crs,utDT)

		#Save
		npl = vO[i]-nS
		fOut = oDir+"/vid.%04d.png"%(npl)
		kv.savePic(fOut,bLenX=45)