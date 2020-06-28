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

	MainS = """Creates simple multi-panel figure for the 
	solar wind conditions within the  bcwind file and saves it as swBCplot.pdf
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="swid",default=swtag,help="Solar wind file used (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	swtag = args.id


	swIn = fdir+'/'+swtag
	kh5.CheckOrDie(swIn)

	fOut = "swBCplot.pdf"

	# pulling UT variable for plotting
	t0Fmt = "%Y-%m-%d %H:%M:%S"
	utfmt='%H:%M \n%Y-%m-%d'

	UTall  = kh5.PullVar(swIn,"UT")

	utall = []
	for n in range(len(UTall)):
		utall.append(datetime.datetime.strptime(UTall[n].decode('utf-8'),t0Fmt))

	# pulling the solar wind values from the table
	D = kh5.PullVar(swIn,"D")
	Cs = kh5.PullVar(swIn,"Cs")
	Vx = kh5.PullVar(swIn,"Vx")
	Vy = kh5.PullVar(swIn,"Vy")
	Vz = kh5.PullVar(swIn,"Vz")
	Bx = kh5.PullVar(swIn,"Bx")
	By = kh5.PullVar(swIn,"By")
	Bz = kh5.PullVar(swIn,"Bz")
	Va = kh5.PullVar(swIn,"Va")
	Temp = kh5.PullVar(swIn,"Temp")
	Tsec = kh5.PullVar(swIn,"T")
	SYMH = kh5.PullVar(swIn,"symh")
	Mfast = kh5.PullVar(swIn,"Magnetosonic Mach")

	Thr = Tsec/3600.0

	# constants
	gamma = 5/3.0
	mp = 1.67e-27 #Proton mass [kg]

	# calculating the solar wind dynamic pressure 
	Vmag = np.sqrt(Vx**2.+Vy**2.+Vz**2.)
	Pram = mp*D*Vmag**2*1.0e15 # nPa

	#Setup figure
	fSz = (10,14)
	Nr = 6
	Nc = 1
	clrs = ['#7570b3','#1b9e77','#d95f02','black'] # from colorbrewer2.org for colorblind safe
	
	fig = plt.figure(figsize=fSz)
	
	gs = gridspec.GridSpec(Nr,Nc,hspace=0.05,wspace=0.05)
	
	ax11 = fig.add_subplot(gs[0,0])
	ax12 = fig.add_subplot(gs[3,0])
	ax21 = fig.add_subplot(gs[1,0])
	ax22 = fig.add_subplot(gs[4,0])
	ax31 = fig.add_subplot(gs[2,0])
	ax32 = fig.add_subplot(gs[5,0])
    
	smlabel = ['SM-X','SM-Y','SM-Z']
	xvec = np.zeros((len(D),3))+1e9
	
	fig.suptitle("Solar Wind",y=0.92,fontsize=14) 
	Dlim=np.max(D)-np.min(D)
	ax11.plot(utall,D,color=clrs[3])
	for i in range(3):
	    ax11.plot(utall,xvec[:,i],linewidth=4,label=smlabel[i],color=clrs[i])
	kv.SetAxLabs(ax11,"","n [cm^-3]",doBot=True,doLeft=True)
	ax11.set_ylim(np.min(D)-0.05*Dlim,np.max(D)+0.05*Dlim)
	ax11.tick_params(axis="x",direction="in")
	plt.setp(ax11.get_xticklabels(),visible=False)
	ax11.legend(ncol=len(smlabel), bbox_to_anchor=(0.5,1),loc='lower center', fontsize='small')
    
	TScl = 1.0e-6
	ax21.plot(utall,Temp*TScl,color=clrs[3])
	kv.SetAxLabs(ax21,"","T [MK]",doBot=True,doLeft=False)
	ax21.tick_params(axis="x",direction="in")
	plt.setp(ax21.get_xticklabels(),visible=False)
	
	ax31.plot(utall,Pram,color=clrs[3])
	ax31.xaxis_date()
	kv.SetAxLabs(ax31,"","Dynamic P [nPa]",doBot=True,doLeft=True)
	ax31.tick_params(axis="x",direction="in")
	plt.setp(ax31.get_xticklabels(),visible=False)
	
	vScl = 1.0e-3
	secax12 = ax12.twinx()
	ax12.plot(utall,Vx*vScl,color=clrs[0],linewidth=0.95)
	secax12.plot(utall,Vy*vScl,color=clrs[1],linewidth=0.95)
	secax12.plot(utall,Vz*vScl,color=clrs[2],linewidth=0.95)
	secax12.set_ylabel('Vy,z [km/s]')
	kv.SetAxLabs(ax12,"","Vx [km/s]",doBot=True,doLeft=True)
	ax12.tick_params(axis="x",direction="in")
	plt.setp(ax12.get_xticklabels(),visible=False)
	
	ax22.plot(utall,Bx,color=clrs[0],linewidth=0.95)
	ax22.plot(utall,By,color=clrs[1],linewidth=0.95)
	ax22.plot(utall,Bz,color=clrs[2],linewidth=0.95)
	ax22.axhline(y=0.0, color='black', linestyle='--',alpha=0.6,linewidth=0.9)
	kv.SetAxLabs(ax22,"","B [nT]",doBot=True,doLeft=False)
	ax22.tick_params(axis="x",direction="in")
	plt.setp(ax22.get_xticklabels(),visible=False)
    
	ax32.plot(utall,SYMH,color=clrs[3])
	ax32.axhline(y=0.0, color='black', linestyle='--',alpha=0.6,linewidth=0.9)
	ax32.xaxis_date()
	xfmt = dates.DateFormatter(utfmt)
	ax32.xaxis.set_major_formatter(xfmt)
	kv.SetAxLabs(ax32,"UT","SYM/H [nT]",doBot=True,doLeft=True)
	kv.savePic(fOut,doTrim=False)
	plt.close('all')

