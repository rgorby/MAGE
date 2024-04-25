#!/usr/bin/env python
#Make video of error between two Gamera cases
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
import kaipy.gamera.magsphere as msph
import kaipy.gamera.rcmpp as rcmpp
from alive_progress import alive_bar
import kaipy.kdefs as kdefs
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import errno
import subprocess
import shutil

cLW = 0.25

def makeMovie(frame_dir,movie_name):
	frame_pattern = frame_dir + "/vid.%04d.png"
	movie_file = os.getcwd() + "/" + movie_name + ".mp4"
	ffmpegExe = "ffmpeg"
	if shutil.which(ffmpegExe) is None:
		ffmpegExe = "ffmpeg4"
		if shutil.which(ffmpegExe) is None:
			print("Could not find any ffmpeg executable. Video will not be generated.")
			return

	cmd = [
	    ffmpegExe, "-nostdin", "-i", frame_pattern,
	    "-vcodec", "libx264", "-crf", "14", "-profile:v", "high", "-pix_fmt", "yuv420p",
	    movie_file,"-y"
	]
	subprocess.run(cmd, check=True)

if __name__ == "__main__":
	#Defaults
	fdir1 = os.getcwd()
	ftag1 = "msphere"
	fdir2 = os.getcwd()
	ftag2 = "msphere"
	oDir = "vid2D"
	ts = 0    #[min]
	te = 200  #[min]
	dt = 60.0 #[sec]
	Nblk = 1 #Number of blocks
	nID = 1 #Block ID of this job
	noMPI = False # Don't add MPI tiling
	noLog = False
	fieldNames = "Bx, By, Bz"
	doVerb = False
	skipMovie = False

	MainS = """Creates simple multi-panel figure for Gamera magnetosphere run
	Left Panel - Residual vertical magnetic field
	Right Panel - Pressure (or density) and hemispherical insets
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d1',type=str,metavar="directory",default=fdir1,help="Directory to read first dataset from (default: %(default)s)")
	parser.add_argument('-id1',type=str,metavar="runid",default=ftag1,help="RunID of first dataset (default: %(default)s)")
	parser.add_argument('-d2',type=str,metavar="directory",default=fdir2,help="Directory to read second dataset from (default: %(default)s)")
	parser.add_argument('-id2',type=str,metavar="runid",default=ftag2,help="RunID of second dataset (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="directory",default=oDir,help="Subdirectory to write to (default: %(default)s)")
	parser.add_argument('-ts' ,type=int,metavar="tStart",default=ts,help="Starting time [min] (default: %(default)s)")
	parser.add_argument('-te' ,type=int,metavar="tEnd"  ,default=te,help="Ending time   [min] (default: %(default)s)")
	parser.add_argument('-dt' ,type=int,metavar="dt"    ,default=dt,help="Cadence       [sec] (default: %(default)s)")
	parser.add_argument('-Nblk' ,type=int,metavar="Nblk",default=Nblk,help="Number of job blocks (default: %(default)s)")
	parser.add_argument('-nID' ,type=int,metavar="nID"  ,default=nID,help="Block ID of this job [1-Nblk] (default: %(default)s)")
	parser.add_argument('-f',type=str,metavar="fieldnames",default=fieldNames,help="Comma-separated fields to plot (default: %(default)s)")
	parser.add_argument('-linear',action='store_true', default=noLog,help="Plot linear line plot instead of logarithmic (default: %(default)s)")
	parser.add_argument('-v',action='store_true', default=doVerb,help="Do verbose output (default: %(default)s)")
	parser.add_argument('-skipMovie',action='store_true', default=skipMovie,help="Skip automatic movie generation afterwards (default: %(default)s)")
	#parser.add_argument('-nompi', action='store_true', default=noMPI,help="Don't show MPI boundaries (default: %(default)s)")

	mviz.AddSizeArgs(parser)

	#Finalize parsing
	args = parser.parse_args()
	fdir1 = args.d1
	ftag1 = args.id1
	fdir2 = args.d2
	ftag2 = args.id2
	ts  = args.ts
	te  = args.te
	dt  = args.dt
	oSub = args.o
	Nblk = args.Nblk
	nID = args.nID
	fieldNames = args.f
	noLog = args.linear
	doVerb = args.v
	#noMPI = args.noMPI
	
	fnList = [item.strip() for item in fieldNames.split(',')]

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
	oDir = os.getcwd() + "/" + oSub
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
	gsph1 = msph.GamsphPipe(fdir1,ftag1)
	gsph2 = msph.GamsphPipe(fdir2,ftag2)
		
	#======
	#Setup figure
	fig = plt.figure(figsize=figSz)
	gs = gridspec.GridSpec(5,2,height_ratios=[20,5,1,5,9],hspace=0.025)
	

	AxTL = fig.add_subplot(gs[0,0])
	AxTR = fig.add_subplot(gs[0,1])
	AxB = fig.add_subplot(gs[-1,0:2])
	AxB2 = AxB.twinx() # second plot on bottom axis

	AxCT = fig.add_subplot(gs[2,0:2])

	print(fig.axes)

	errTimes = []
	errListRel = []
	errListAbs = []
	relColor = "tab:blue"
	absColor = "tab:orange"

	#Loop over sub-range
	titstr = "Comparing '%s' to '%s'"%(fdir1,fdir2)
	with alive_bar(i1-i0,title=titstr.ljust(kdefs.barLab),length=kdefs.barLen,disable=doVerb) as bar:
		for i in range(i0,i1):
			#Convert time (in seconds) to Step #
			nStp = np.abs(gsph1.T-tOut[i]).argmin()+gsph1.s0
			if doVerb:
				print("Minute = %5.2f / Step = %d"%(tOut[i]/60.0,nStp))
			npl = vO[i]

			AxTL.clear()
			AxTR.clear()
			AxB.clear()
			AxB2.clear()

			#plot upper left msph error
			mviz.PlotEqErrRel(gsph1,gsph2,nStp,xyBds,AxTL,fnList,AxCB=AxCT,doVerb=doVerb)
			AxTL.set_title("Equatorial Slice of Relative Error")

			#plot upper right k-axis error
			mviz.PlotLogicalErrRel(gsph1,gsph2,nStp,AxTR,fnList,2,doVerb=doVerb)
			AxTR.set_title("Per-Cell Relative Error along K-Axis")
			if (not noMPI):
				#plot I-MPI decomp on logical plot
				if(gsph2.Ri > 1):
					for im in range(gsph2.Ri):
						i0 = im*gsph2.dNi
						AxTR.plot([i0, i0],[0, gsph2.Nj],"deepskyblue",linewidth=0.25,alpha=0.5)
				#plot J-MPI decomp on logical plot
				if (gsph2.Rj>1):
					for jm in range(1,gsph2.Rj):
						j0 = jm*gsph2.dNj
						AxTR.plot([0, gsph2.Ni],[j0, j0],"deepskyblue",linewidth=0.25,alpha=0.5)

			#plot bottom line plot
			errTimes.append(tOut[i]/60.0)
			errListRel.append(mviz.CalcTotalErrRel(gsph1,gsph2,nStp,fnList,doVerb=doVerb))
			errListAbs.append(mviz.CalcTotalErrAbs(gsph1,gsph2,nStp,fnList,doVerb=doVerb))
			if noLog:
				AxB.plot(errTimes, errListRel,color=relColor)
				AxB2.plot(errTimes, errListAbs,color=absColor)
			else:
				AxB.semilogy(errTimes, errListRel,color=relColor)
				AxB2.semilogy(errTimes, errListAbs,color=absColor)
			AxB.set_xlabel('Time (min)')
			AxB.set_ylabel('Per-Cell Mean Relative Error',color=relColor)
			AxB.tick_params(axis='y',which='both',colors=relColor,left=True,right=True,labelleft=True,labelright=False)
			AxB2.set_ylabel('Per-Cell Mean Absolute Error',color=absColor)
			#AxB2.yaxis.tick_right()
			AxB2.tick_params(axis='y',which='both',colors=absColor,left=True,right=True,labelleft=False,labelright=True)
			AxB.set_title("'" + fieldNames + "' Per-Cell Error Over Time")

			gsph1.AddTime(nStp,AxTL,xy=[0.025,0.84],fs="x-large")

			#Add MPI decomp
			if (not noMPI):
				mviz.PlotMPI(gsph2,AxTL)
    
			fOut = oDir+"/vid.%04d.png"%(npl)
			kv.savePic(fOut,bLenX=45)

			bar()
	makeMovie(oDir,oSub)

