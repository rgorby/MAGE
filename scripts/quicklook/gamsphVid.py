#!/usr/bin/env python
#Make video of Gamera magnetosphere run
import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
import matplotlib.gridspec as gridspec
import numpy as np
import numpy as np
import kaipy.gamera.msphViz as mviz
import kaipy.remix.remix as remix
import kaipy.gamera.magsphere as msph
import kaipy.gamera.gampp as gampp
import kaipy.gamera.rcmpp as rcmpp
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import errno

cLW = 0.25

if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	ftag = "msphere"
	oDir = "vid2D"
	doDen = False
	ts = 0    #[min]
	te = 200  #[min]
	dt = 60.0 #[sec]
	doBig = False #[Use big window]
	noIon = False
	noRCM = False
	doMPI = False #[Add MPI tiling]
	Nblk = 1 #Number of blocks
	nID = 1 #Block ID of this job
	noMPI = False
	doJy = False
	doBz = False
	doBigRCM = False

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
	#parser.add_argument('-nompi', action='store_true', default=noMPI,help="Don't show MPI boundaries (default: %(default)s)")
	parser.add_argument('-bz'   , action='store_true', default=doBz ,help="Show Bz instead of dBz (default: %(default)s)")
	parser.add_argument('-jy'   , action='store_true', default=doJy ,help="Show Jy instead of pressure (default: %(default)s)")
	parser.add_argument('-bigrcm', action='store_true',default=doBigRCM,help="Show entire RCM domain (default: %(default)s)")
	parser.add_argument('-noion', action='store_true', default=noIon,help="Don't show ReMIX data (default: %(default)s)")
	parser.add_argument('-norcm', action='store_true', default=noRCM,help="Don't show RCM data (default: %(default)s)")

	mviz.AddSizeArgs(parser)

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
	doJy = args.jy
	doBz = args.bz
	doBigRCM = args.bigrcm
	
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

	#Get domain size
	xyBds = mviz.GetSizeBds(args)
	
	#---------
	#Figure parameters
	figSz = (12,7.5)


	#======
	#Init data
	gsph = msph.GamsphPipe(fdir,ftag)

	#Check for remix
	rcmChk = fdir + "/%s.mhdrcm.h5"%(ftag)
	rmxChk = fdir + "/%s.mix.h5"%(ftag)
	doRCM = os.path.exists(rcmChk)
	doMIX = os.path.exists(rmxChk)

	if (doRCM and (not args.norcm)):
		print("Found RCM data")
		rcmdata = gampp.GameraPipe(fdir,ftag+".mhdrcm")
		mviz.vP = kv.genNorm(1.0e-2,100.0,doLog=True)
		rcmpp.doEll = not doBigRCM
	if (doMIX and (not args.noion)):
		print("Found ReMIX data")
		
		
	#======
	#Setup figure
	fig = plt.figure(figsize=figSz)
	gs = gridspec.GridSpec(3,6,height_ratios=[20,1,1],hspace=0.025)
	

	AxL = fig.add_subplot(gs[0,0:3])
	AxR = fig.add_subplot(gs[0,3:])

	AxC1 = fig.add_subplot(gs[-1,0:2])
	AxC2 = fig.add_subplot(gs[-1,2:4])
	AxC3 = fig.add_subplot(gs[-1,4:6])

	cbM = kv.genCB(AxC2,kv.genNorm(remix.facMax),"FAC",cM=remix.facCM,Ntk=4)
	AxC2.xaxis.set_ticks_position('top')

	#Loop over sub-range
	for i in range(i0,i1):
		#Convert time (in seconds) to Step #
		nStp = np.abs(gsph.T-tOut[i]).argmin()+gsph.s0
		print("Minute = %5.2f / Step = %d"%(tOut[i]/60.0,nStp))
		npl = vO[i]

		AxL.clear()
		AxR.clear()

		Bz = mviz.PlotEqB(gsph,nStp,xyBds,AxL,AxC1,doBz=doBz)

		if (doJy):
			mviz.PlotJyXZ(gsph,nStp,xyBds,AxR,AxC3)
		else:
			mviz.PlotMerid(gsph,nStp,xyBds,AxR,doDen,doRCM,AxC3)

		gsph.AddTime(nStp,AxL,xy=[0.025,0.89],fs="x-large")
		gsph.AddSW(nStp,AxL,xy=[0.625,0.025],fs="small")

			
		#Add inset RCM plot
		if (doRCM and (not args.norcm)):
			AxRCM = inset_axes(AxL,width="30%",height="30%",loc=3)
			rcmpp.RCMInset(AxRCM,rcmdata,nStp,mviz.vP)
			AxRCM.contour(kv.reWrap(gsph.xxc),kv.reWrap(gsph.yyc),kv.reWrap(Bz),[0.0],colors=mviz.bz0Col,linewidths=mviz.cLW)
			rcmpp.AddRCMBox(AxL)

		if (doMIX and (not args.noion)):
			ion = remix.remix(rmxChk,nStp)
			gsph.AddCPCP(nStp,AxR,xy=[0.610,0.925])
			mviz.AddIonBoxes(gs[0,3:],ion)		

		#Add MPI decomp
		if (doMPI):
			mviz.PlotMPI(gsph,AxL)
			mviz.PlotMPI(gsph,AxR)

		fOut = oDir+"/vid.%04d.png"%(npl)
		kv.savePic(fOut,bLenX=45)
