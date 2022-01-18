#!/usr/bin/env python
#Make a quick figure of a Gamera helio run

import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes #do not need for helio if not inserting anything
import matplotlib.gridspec as gridspec
import numpy as np

import kaipy.gamhelio.helioViz as hviz
import kaipy.gamhelio.heliosphere as hsph
import kaipy.gamera.gampp as gampp #serial/MPI gamera output
import os

if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	ftag = "wsa"
	nStp = -1
	fOut = "qkpichelio.png"
	doDen = False
	noMPI = False
	pic = "pic1"
	MainS = """Creates simple multi-panel figure for Gamera helio run
	Top Panel - density and radial velocity in equatorial plane
	Bottom Panel - density and radial velocity in meridional plane
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-n' ,type=int,metavar="step" ,default=nStp,help="Time slice to plot (default: %(default)s)")
	parser.add_argument('-den'  , action='store_true', default=doDen,help="Show density instead of pressure (default: %(default)s)")
	parser.add_argument('-nompi', action='store_true', default=noMPI,help="Don't show MPI boundaries (default: %(default)s)")
	#pic1 is equatorial pic2 is meridional pic3 is 1 AU maps pic4 is Br in rotating frame
	parser.add_argument('-p',type=str,metavar="pictype",default=pic,help="Type of output image (default: %(default)s)")	

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id
	nStp = args.n
	doDen = args.den
	noMPI = args.nompi
	doMPI = (not noMPI)
	pic = args.p

	#Get domain size
	xyBds = hviz.GetSizeBds(pic)
	print (xyBds)

	#---------------------
	#Do work
	doFast=False

	#---------
	#Figure parameters depending on a pic
	if (pic == "pic1" or pic == "pic2"):
		figSz = (10,12.5)
	elif (pic == "pic3"):
		figSz = (10,6.5)
	else:
		figSz = (10,6.)
	#======
	#Init data
	gsph = hsph.GamsphPipe(fdir,ftag,doFast=doFast)

	if (nStp<0):
		nStp = hsph.sFin
		print("Using Step %d"%(nStp))

	#======
	#Setup figure
	fig = plt.figure(figsize=figSz)

	if (pic != "pic4"):	
		gs = gridspec.GridSpec(4,6,height_ratios=[20,1,20,1])
		#plots. Two rows of two plots
		AxL0 = fig.add_subplot(gs[0,0:3])
		AxR0 = fig.add_subplot(gs[0,3:])

		AxL1 = fig.add_subplot(gs[2,0:3])
		AxR1 = fig.add_subplot(gs[2,3:])

		#colorbars. Two rows
		AxC1_0 = fig.add_subplot(gs[1,0:3])
		AxC2_0 = fig.add_subplot(gs[1,3:])

		AxC1_1 = fig.add_subplot(gs[3,0:3])
		AxC2_1 = fig.add_subplot(gs[3,3:])
	else:
		gs = gridspec.GridSpec(2,1,height_ratios=[20,1])
		Ax = fig.add_subplot(gs[0,0])
		AxC = fig.add_subplot(gs[1,0])


	if (pic == "pic1"):
		hviz.PlotEqMagV(gsph,nStp,xyBds,AxL0,AxC1_0)
		hviz.PlotEqD(gsph,nStp,xyBds,AxR0,AxC2_0)

		hviz.PlotEqTemp(gsph,nStp,xyBds,AxL1,AxC1_1)
		hviz.PlotEqBr(gsph,nStp,xyBds,AxR1,AxC2_1)
	elif (pic == "pic2"):
		hviz.PlotMerMagV(gsph,nStp,xyBds,AxL0,AxC1_0)
		hviz.PlotMerDNorm(gsph,nStp,xyBds,AxR0,AxC2_0)

		hviz.PlotMerTemp(gsph,nStp,xyBds,AxL1,AxC1_1)
		hviz.PlotMerBrNorm(gsph,nStp,xyBds,AxR1,AxC2_1)
	elif (pic == "pic3"):
		hviz.PlotiSlMagV(gsph,nStp,xyBds,AxL0,AxC1_0)
		hviz.PlotiSlD(gsph,nStp,xyBds,AxR0,AxC2_0)

		hviz.PlotiSlTemp(gsph,nStp,xyBds,AxL1,AxC1_1)
		hviz.PlotiSlBr(gsph,nStp,xyBds,AxR1,AxC2_1)
	elif (pic == "pic4"):
		hviz.PlotiSlBrRotatingFrame(gsph,nStp,xyBds,Ax,AxC)
	else:
		print ("Pic is empty. Choose pic1 or pic2 or pic3")
	
#	add time (upper left)
#	if (pic == "pic1" or pic == "pic2"):
#		gsph.AddTime(nStp,AxL0,xy=[0.025,0.875],fs="x-large")
#	elif (pic == "pic3"):	
#		gsph.AddTime(nStp,AxL0,xy=[0.015,0.82],fs="small")
#	elif (pic == "pic4"):
#		gsph.AddTime(nStp,Ax,xy=[0.015,0.92],fs="small")
#	else:
#		print ("Pic is empty. Choose pic1 or pic2 or pic3")

	#TO DO: Add wsa info (lower left) instead of Solar wind params
	#gsph.AddSW(nStp,AxL,xy=[0.625,0.025],fs="small")

	#Need to re-write PlotMPI from msphr to helio case
	#Add MPI decomp
	#if (doMPI):
	#	hviz.PlotMPI(gsph,AxL0)
	#	hviz.PlotMPI(gsph,AxR0)
	#	hviz.PlotMPI(gsph,AxL1)
	#	hviz.PlotMPI(gsph,AxR1)
	#	hviz.PlotMPI(gsph,AxL2)
	#	hviz.PlotMPI(gsph,AxR2)

	kv.savePic(fOut,bLenX=45)
