#!/usr/bin/env python
"""
this code computes number flux vs energy distribution of RCM diffuse precipitation in units of cm^-2 eV^-1 str^-1 s^-1
"""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import kaipy.kaiViz as kv
import matplotlib.gridspec as gridspec
import numpy as np
import kaipy.kaiH5 as kh5

from matplotlib import rcParams, cycler

import argparse
from argparse import RawTextHelpFormatter

if __name__ == "__main__":
	ftag = "msphere"
	cmap = plt.cm.plasma
	timeStr = ""
	locStr = "-5,0,0"
	numSamples = 6
	
	MainS = """Creates a plot of the differential flux for RCM ions or electrons in units of cm^-2 keV^-1 str^-1
	"""
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-n',type=int,metavar="numSamples",default=numSamples,help="Number of evenly-spaced time samples (default: %(default)s)")
	parser.add_argument('-t',type=str,metavar="times",default=timeStr,help="Comma-separated times (in hours) to plot (example: 1,2,4,6). Ignores numSamples")
	parser.add_argument('-l',type=str,metavar="loc",default=locStr,help="Comma-separated x,y,z values for equatorial location (default: %(default)s)")
	parser.add_argument('-i',action='store_true',default=False,help="Flag to plot ions instead of electrons (default: %(default)s)")
	
	args = parser.parse_args()
	ftag = args.id
	numSamples = args.n
	timeStr = args.t
	locStr = args.l
	doIons = args.i

	fIn = ftag+".rcm.h5"
	x0,y0,z0 = [float(x) for x in locStr.split(',')]

#--Time stuff--
	#Get rcm data timesteps
	nDataSteps,sIds = kh5.cntSteps(fIn)
	dataTimes = kh5.getTs(fIn,sIds,"time")/3600  # [Hrs]
	if timeStr == "":  # If user didn't specify the times they wanted, do evenly spaced samples using numSamples
		tStart = np.amin(np.abs(dataTimes))
		tEnd = np.amax(dataTimes)
		nHrs = np.linspace(tStart, tEnd, num=numSamples, endpoint=True)
	else:  # Otherwise use user-specified times
		nHrs = [float(t) for t in timeStr.split(',')]
	NumStps = len(nHrs)
	rcParams['axes.prop_cycle'] = cycler(color=cmap(np.linspace(0, 1, NumStps)))

	fSz = (14,7)
	fig = plt.figure(figsize=fSz)
	# conversion factor
	massi = 1.67e-27 # mass of ions in kg
	ev = 1.607e-19 # 1ev in J
	nT2T = 1.e-9
	m2cm = 1.e2
	re = 6.380e6 # radius of earth in m

	Legs = []

	eos = False #Flag for reaching end of steps
	for n in range(NumStps):
		
		nStpIdx = np.abs(dataTimes - nHrs[n]).argmin()
		nStp = nStpIdx + sIds.min()
		legStr = "T = +%2.1f Hours"%(dataTimes[nStpIdx])
		if nStpIdx == nDataSteps: # Note that we're about to plot the last timestep and there's no need to plot any more afterwards
			eos=True
		elif nHrs[n] > dataTimes[-1]:
			if not eos:  # If this is the first time we've exceeded the end of the available timesteps, plot the last timestep anyways
				print("%2.1f [hrs] out of time range (%2.1f, %2.1f), using last step time (%2.1f)"%(nHrs[n],dataTimes[0],dataTimes[-1],dataTimes[nStpIdx]))
				eos = True
			else:
				print("%2.1f [hrs] out of time range (%2.1f, %2.1f), skipping."%(nHrs[n],dataTimes[0],dataTimes[-1]))
				continue
		
		#Pull 3D data
		deleeta = kh5.PullVar(fIn,"rcmdeleeta",nStp)
		vm   = kh5.PullVar(fIn,"rcmvm"  ,nStp)
		lamc = kh5.PullVar(fIn,"alamc"  ,nStp)
		bir  = kh5.PullVar(fIn,"bir"    ,nStp)
		sini = kh5.PullVar(fIn,"sini"   ,nStp)
		xeq  = kh5.PullVar(fIn,"rcmxmin"  ,nStp)
		yeq  = kh5.PullVar(fIn,"rcmymin"  ,nStp)
		zeq  = kh5.PullVar(fIn,"rcmzmin"  ,nStp)
		#Pull dt in s
		dtCpl = kh5.tStep(fIn,nStp,aID="dtCpl",aDef=15.0)

		#Do grid calculations on first time
		if (n==0):
			#Set interfaces
			Nlat,Nlon,Nk = deleeta.shape
			kion = (lamc>0).argmax()
			
			if doIons:
				kStart=kion	
				kEnd=Nk
			else:
				kStart=1  # Skip plasmasphere channel
				kEnd=kion

			Nks = kEnd-kStart
			ilamc = lamc[kStart:kEnd]
			ilami = np.zeros(Nks+1)
			for m in range(0,Nks-1):
				nk = m+kStart
				ilami[m+1] = 0.5*(lamc[nk]+lamc[nk+1])
			ilami[Nks] = lamc[-1] + 0.5*(lamc[-1]-lamc[-2])

			#Ensure positive energies in case of electrons
			ilami = np.abs(ilami)
			ilamc = np.abs(ilamc)

		#Find nearest point (everytime)
		dR = np.sqrt( (xeq-x0)**2.0 + (yeq-y0)**2.0 + (zeq-z0)**2.0 )
		i0,j0 = np.unravel_index(dR.argmin(), dR.shape)

		#Get energy bins in eV
		Ki = vm[i0,j0]*ilami 
		Kc = vm[i0,j0]*ilamc
		#Bin interval in eV
		ekscl = vm[i0,j0]*np.diff(ilami)   

		#Get differential flux in 1/(eV cm^2 str s)
		ijEta = deleeta[i0,j0,kStart:kEnd]*bir[i0,j0]*sini[i0,j0]/ekscl/dtCpl*nT2T/(m2cm**2)/(4*np.pi)
		Legs.append(legStr)
		print(legStr)
		print("\tMin/Mean/Max K = %f / %f / %f"%(Kc.min(),Kc.mean(),Kc.max()))
		k0 =ijEta.argmax()
		print("\tMax @ K = %f"%(Kc[k0]))
		print("\tVM = %f"%(vm[i0,j0]))
		plt.loglog(Kc,ijEta)

	sName = "Ions" if doIons else "Electrons"
	titStr = "%s @ (x,y,z) = (%5.2f,%5.2f,%5.2f)"%(sName,x0,y0,z0)
	Ax = plt.gca()
	plt.legend(Legs,prop={'family': 'monospace'},loc='lower left')

	Ax.legend(Legs,loc='lower left')
	Ax.set_xlabel("Energy [eV]")
	Ax.set_ylabel("differential energy flux /cm^2/eV/str/s")
	Ax.set_title(titStr)
	Ax.grid()

	#Ax.set_ylim(1.0e+10,1.0e+19)
	sTag = "_i" if doIons else "_e"
	kv.savePic("qkRcmPrecipSpec%s.png"%sTag)

	#plt.show()
	#
