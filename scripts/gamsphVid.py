#!/usr/bin/env python
#Make video of Gamera magnetosphere run
import argparse
from argparse import RawTextHelpFormatter
import kaipy.gamera.magsphere as msph
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
import matplotlib.gridspec as gridspec
import numpy as np
import kaipy.gamera.remixpp as rmpp
import os
import errno

if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	ftag = "msphere"
	oDir = "vid2D"

	ts = 0    #[min]
	te = 200  #[min]
	dt = 60.0 #[sec]

	Nblk = 1 #Number of blocks
	nID = 1 #Block ID of this job

	MainS = """Creates simple multi-panel figure for Gamera magnetosphere run
	Left Panel - Residual vertical magnetic field
	Right Panel - Pressure (or density) and hemispherical insets
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="directory",default=oDir,help="Subdirectory to write to (default: %(default)s)")
	parser.add_argument('-ts' ,type=int,metavar="tStart",default=ts,help="Starting time [min] (default: %(default)s)")
	parser.add_argument('-te' ,type=int,metavar="tEnd"  ,default=te,help="Ending time   [min] (default: %(default)s)")
	parser.add_argument('-dt' ,type=int,metavar="dt"    ,default=dt,help="Cadence       [sec] (default: %(default)s)")
	parser.add_argument('-Nblk' ,type=int,metavar="Nblk",default=Nblk,help="Number of job blocks (default: %(default)s)")
	parser.add_argument('-nID' ,type=int,metavar="nID"  ,default=nID,help="Block ID of this job [1-Nblk] (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id
	ts  = args.ts
	te  = args.te
	dt  = args.dt
	oSub = args.o
	Nblk = args.Nblk
	nID = args.nID

	#Setup timing info
	tOut = np.arange(ts*60.0,te*60.0,dt)
	Nt = len(tOut)
	vO = np.arange(0,Nt)

	print("Writing %d outputs between minutes %d and %d"%(Nt,ts,te))
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
	
	#---------
	#Figure parameters
	figSz = (12,7.5)
	dbCMap = "RdGy_r"
	bCMap = "inferno"
	pCMap = "viridis"
	

	doMPI = False
	xyBds = [-40,20,-30,30]

	PMin = 0.0 ; PMax = 2.0
	Nc = 11
	vDB = kv.genNorm(25)
	vP = kv.genNorm(PMin,PMax)
	pVals = np.linspace(PMin,PMax,Nc)

	#======
	#Init data
	gsph = msph.GamsphPipe(fdir,ftag)

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
	kv.genCB(AxC4,kv.genNorm(rmpp.cMax),"FAC",cM=rmpp.fcMap,Ntk=5)
	rmpp.AddPotCB(AxC3)
	kv.genCB(AxC2,vP,"Pressure",cM=pCMap,Ntk=6)

	#Loop over sub-range
	for i in range(i0,i1):
		#Convert time (in seconds) to Step #
		nStp = np.abs(gsph.T-tOut[i]).argmin()+gsph.s0
		print("Minute = %5.2f / Step = %d"%(tOut[i]/60.0,nStp))
		npl = vO[i]

		AxL.clear()
		AxR.clear()
		
		dbz = gsph.DelBz(nStp)
		Pxz = gsph.EggSlice("P",nStp,vScl=gsph.pScl,doEq=False)

		#Start plotting
		AxL.pcolormesh(gsph.xxi,gsph.yyi,dbz,cmap=dbCMap,norm=vDB)

		kv.addEarth2D(ax=AxL)
		kv.SetAx(xyBds,AxL)
		gsph.AddTime(nStp,AxL,xy=[0.025,0.935],fs="x-large")
		gsph.AddSW(nStp,AxL,xy=[0.625,0.025],fs="small")
		AxL.set_xlabel('SM-X [Re]')
		AxL.set_ylabel('SM-Y [Re]')

		AxR.pcolormesh(gsph.xxi,gsph.yyi,Pxz,cmap=pCMap,norm=vP)

		kv.addEarth2D(ax=AxR)
		kv.SetAx(xyBds,AxR)
		AxR.yaxis.tick_right()
		AxR.yaxis.set_label_position('right')
		gsph.AddCPCP(nStp,AxR,xy=[0.65,0.925])

		AxR.set_xlabel('SM-X [Re]')
		AxR.set_ylabel('SM-Z [Re]')


		#Add MPI decomp
		LW = 0.25
		ashd = 0.25
		if (doMPI):
			for i in range(gsph.Ri):
				i0 = i*gsph.dNi
				AxL.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],'m',linewidth=LW,alpha=ashd)
				AxR.plot(gsph.xxi[i0,:],gsph.yyi[i0,:],'c',linewidth=LW,alpha=ashd)


		fmix = gsph.Gam2Remix(nStp)
		dxy = [32.5,32.5]
		rmpp.CMIViz(AxR,fmix,dxy=dxy)
		rmpp.CMIViz(AxR,fmix,dxy=dxy,loc="lower left",doNorth=False)

		fOut = oDir+"/vid.%04d.png"%(npl)
		kv.savePic(fOut,bLenX=45)
