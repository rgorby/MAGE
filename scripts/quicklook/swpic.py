#!/usr/bin/env python
#Creates a time vs distance plot from a 2D slice created by slice.x
import argparse
from argparse import RawTextHelpFormatter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import dates
import kaipy.kaiViz as kv
import kaipy.kaiH5 as kh5
import kaipy.solarWind.swBCplots as swBCplots
import matplotlib.gridspec as gridspec
import numpy as np
import h5py
import os
import errno
import datetime
import matplotlib.dates as mdates


if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	swtag = "bcwind.h5"
	imgtype = 'png'

	MainS = """Creates simple multi-panel figure for the 
	solar wind conditions within the  bcwind file and saves it as swBCplot.pdf
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="swid",default=swtag,help="Solar wind file used (default: %(default)s)")
	parser.add_argument('-type',type=str,metavar="type",default=imgtype,help="Image type (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	swtag = args.id
	imgtype = args.type


	swIn = fdir+'/'+swtag
	kh5.CheckOrDie(swIn)

	# Name the output file the same as the solarwind file with the image extension
	fOut = swtag.split('.')[0]+'.'+imgtype

	# pulling UT variable for plotting
	t0Fmt = "%Y-%m-%d %H:%M:%S"
	utfmt='%H:%M \n%Y-%m-%d'

	UTall  = kh5.PullVar(swIn,"UT")

	utall = []
	for n in range(len(UTall)):
		utall.append(datetime.datetime.strptime(UTall[n].decode('utf-8'),t0Fmt))

	# pulling the solar wind values from the table
	D = kh5.PullVar(swIn,"D")
	Vx = kh5.PullVar(swIn,"Vx")
	Vy = kh5.PullVar(swIn,"Vy")
	Vz = kh5.PullVar(swIn,"Vz")
	Bx = kh5.PullVar(swIn,"Bx")
	By = kh5.PullVar(swIn,"By")
	Bz = kh5.PullVar(swIn,"Bz")
	Temp = kh5.PullVar(swIn,"Temp")
	Tsec = kh5.PullVar(swIn,"T")
	SYMH = kh5.PullVar(swIn,"symh")

	pltInterp = 0*D
	doEps = False
	swBCplots.swQuickPlot(UTall,D,Temp,Vx,Vy,Vz,Bx,By,Bz,SYMH,pltInterp,fOut,doEps=doEps)


