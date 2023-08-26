#!/usr/bin/env python
#Pulls RBSP data from cdasws and compares to model output
import argparse
from argparse import RawTextHelpFormatter
#import spacepy and cdasws
import spacepy
from spacepy.coordinates import Coords
from spacepy.time import Ticktock
import spacepy.datamodel as dm
import spacepy.plot as splot
from cdasws import CdasWs
#import standard python stuff
import json
import datetime
import os
import errno
#import numpy and matplotlib
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib import dates
from astropy.time import Time
import scipy.interpolate
#Kaipy and related
import kaipy.kaiViz as kv
import kaipy.kaiTools as kaiTools
import kaipy.chimp.kCyl as kc
import kaipy.kaiH5 as kh5
# import kaipy.gamera.gampp as gampp

TINY = 1.0e-8
jScl = 1 # do we need a factor of 4*np.pi somewhere??

def TWin(Ik,Sig):
	import scipy.ndimage.filters as sfil
	import scipy.ndimage as ndimage
	
	IkS = ndimage.gaussian_filter(Ik,sigma=Sig,mode='nearest')

	return IkS

#Given sin^n(alpha) dep. on intensity calculate fraction based on accessible Alpha
def getJScl(Bmag,Beq,en=2.0):
	Na = 360
	A = np.linspace(0,0.5*np.pi,Na)
	da = A[1]-A[0]
	Ia = np.sin(A)**en
	Ic = np.zeros(Ia.shape)
	Nt = len(Bmag)
	I0 = Ia.sum()

	It = np.zeros(Nt)
	for n in range(Nt):
		if (Bmag[n]<TINY):
			It[n] = 0.0
		else:
			Ac = np.arcsin(np.sqrt(Beq[n]/Bmag[n]))
			Ic[:] = Ia[:]
			Icut = (A>Ac)
			Ic[Icut] = 0.0
			It[n] = Ic.sum()/I0
	return It

if __name__ == "__main__":
	#Defaults
	fdir  = os.getcwd()
	ftag  = "eRBpsd.ps.h5"
	trtag = "RBSP-A_MAGNETOMETER_1SEC-GSM_EMFISIS-L3.sc.h5"
	sctag = 'A'
	Ks    = 1000
	R0    = 2.0

	MainS = """Pulls RBSP data and compares it to synthetic RBSP intensity measurementsfrom the simulation, 
	calculated from extracted RBSP trajectory and PSD files.
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of model data (default: %(default)s)")
	parser.add_argument('-k' ,type=float,metavar="energy" ,default=Ks,help="Energy to comparen 1D profile [keV] (default: %(default)s)")
	parser.add_argument('-trj',type=str,metavar="scTrk",default=trtag,help="spacecraft trajectory file (default: %(default)s)")
	parser.add_argument('-sc',type=str,metavar="spacecraft",default=sctag,help="RBSP s/c to plot 'A' or 'B' (default: %(default)s)")
	parser.add_argument('-r0',type=float,metavar="R0",default=R0,help="radius w/in which to mask observations (default: %(default)s)")

	#Finalize parsing
	args  = parser.parse_args()
	fdir  = args.d
	ftag  = args.id
	trtag = args.trj
	sctag = args.sc
	KS    = args.k
	R0    = args.r0

	#======
	#Init data
	fIn = fdir+'/'+ftag
	kh5.CheckOrDie(fIn)
	fTrk = fdir+'/'+trtag
	kh5.CheckOrDie(fTrk)

	isotfmt = '%Y-%m-%dT%H:%M:%S.%f'
	utfmt='%H:%M \n%Y-%m-%d'

	#======
	#Get track data
	#======

	xeq = kh5.PullVar(fTrk,"xeq")
	yeq = kh5.PullVar(fTrk,"yeq")
	scT = kh5.PullVar(fTrk,"T")
	scMJDs = kh5.PullVar(fTrk,"MJDs")

	#Get information for mirror ratio
	Bx = kh5.PullVar(fTrk,"Bx")
	By = kh5.PullVar(fTrk,"By")
	Bz = kh5.PullVar(fTrk,"Bz")
	Bmag = np.sqrt(Bx**2.0 + By**2.0 + Bz**2.0)
	Beq = kh5.PullVar(fTrk,"Beq")

	J0 = getJScl(Bmag,Beq)

	Req = np.sqrt(xeq**2.0 + yeq**2.0)
	Peq = np.arctan2(yeq,xeq)
	Peq[Peq<0] = Peq[Peq<0] + 2*np.pi

	Nsc = len(scT)

	scTi = np.linspace(scT.min(),scT.max(),Nsc+1)
	scMJDi = np.linspace(scMJDs.min(),scMJDs.max(),Nsc+1)

	UTi = Time(scMJDi,format='mjd').isot
	uti = [datetime.datetime.strptime(UTi[n],isotfmt) for n in range(len(UTi))]

	UT = Time(scMJDs,format='mjd').isot
	ut = [datetime.datetime.strptime(UT[n],isotfmt) for n in range(len(UT))]
	
	#======
	#Get RBSP data
	#======

	t0r = uti[0].strftime("%Y-%m-%dT%H:%M:%SZ")
	t1r = uti[-1].strftime("%Y-%m-%dT%H:%M:%SZ")

	if (sctag == 'A' or sctag == 'B'):
		scStr  = "RBSP%s_REL03_ECT-MAGEIS-L2"%(sctag)
		ephStr = "RBSP-%s_MAGNETOMETER_1SEC-GSM_EMFISIS-L3"%(sctag) 
	else:
		print("Unable to find s/c: %s, please set to 'A' or 'B';"%(sctag))
		print("!!Exiting!!")
		quit()

	cdas = CdasWs()

	Qstr = 'FESA'
	status,data = cdas.get_data(scStr,[Qstr],t0r,t1r)
	stat,ephdata = cdas.get_data(ephStr,['coordinates'],t0r,t1r)

	#Get radius (ephemeris data at different cadence)
	Rscr = np.sqrt( ephdata['coordinates'][:,0]**2.0 + ephdata['coordinates'][:,1]**2.0 + ephdata['coordinates'][:,2]**2.0)
	Rscr = Rscr/6380.0
	Tscr = ephdata['Epoch']

	TINY = 1.0e-8

	Estr = 'FEDU_Energy'

	T = data['Epoch']
	E = np.asarray(data[Estr])
	Q = np.transpose(np.asarray(data[Qstr]))

	Q[Q<0] = TINY
	E = E[0,:]
	k0 = (E>0).argmax()

	#Chop out shenanigans
	E = E[k0:]
	Q = Q[k0:,:]
	for n in range(len(T)):
		n0 = np.abs(T[n]-Tscr).argmin() #Closest time in ephemeris
		R = Rscr[n0]
		if (R<=R0):
			Q[:,n] = 0.0

	# get indices for 1D comparison
	k0sc = (np.abs(E-Ks)).argmin()
	Q0 = Q[k0sc,:]
	Q0[Q0<=TINY] = np.nan # dont plot bad data

	#======
	#Get PSD data 
	#======

	xx,yy,Ki,Kc = kc.getGrid(ftag)
	Nk = len(Kc)
	J = np.zeros((Nsc,Nk))

	Nig,Njg = xx.shape
	Ni = Nig-1; Nj = Njg-1
	Nk = len(Kc)

	Nt,sIDs = kh5.cntSteps(ftag)
	psT = kh5.getTs(ftag,sIDs,aID="time")
	psMJDs = kh5.getTs(ftag,sIDs,aID="MJD")

	Ri = xx[:,0] #L interfaces
	Rc = 0.5*(Ri[1:] + Ri[0:-1])
	Pi = np.linspace(0,2*np.pi,Nj+1)
	Pc = 0.5*(Pi[1:] + Pi[0:-1])


	s0 = sIDs.min()
	sE = sIDs.max()
	psTmin = psT.min()
	psTmax = psT.max()
	psMJDMax = psMJDs.max()
	psMJDMin = psMJDs.min()

	psSteps = np.arange(s0,sE)
	Nt = len(psSteps)
	Jrpt = np.zeros((Ni,Nj,Nk,Nt))
	for n in range(Nt):
		nStp = psSteps[n]
		Jrpt[:,:,:,n] = kh5.PullVar(ftag,"jPSD",nStp)

	Lmin = Ri[0]
	Lmax = Ri[-1]
	dphi = 2*np.pi/Nj


	aD = 0.25*0.25
	aF = 0.50*0.25
	aC = 0.50*0.50

	#Create interpolator
	Dims = (Rc,Pc,Kc,psMJDs[0:-1])
	JrptkI = scipy.interpolate.RegularGridInterpolator(Dims,Jrpt,method='linear',bounds_error=False,fill_value=0.0)

	for n in range(Nsc):
		isBad = (Req[n]<=Lmin) or (Req[n]>=Lmax) or (scMJDs[n]<psMJDMin) or (scMJDs[n]>psMJDMax)
		if (isBad):
			J[n,:] = 0.0
		else:
			#TODO: Should be interpolating here
			i0 = (Ri>=Req[n]).argmax() - 1
			j0 = int(np.floor(Peq[n]//dphi))
			t0 = np.abs(psMJDs-scMJDs[n]).argmin()

			for i in range(Nk):
				J[n,i] = JrptkI( [Req[n],Peq[n],Kc[i],scMJDs[n]] )
			J[n,:] = J0[n]*J[n,:]

	k0f = (np.abs(Kc-Ks)).argmin()

	fntSz = "small"

	cmapName = 'gnuplot2'
	vMin=1
	vMax=1.0e6
	norm = kv.genNorm(vMin,vMax,doLog=True)
	fig = plt.figure(figsize=(12,8.5))
	gs = gridspec.GridSpec(2,2,width_ratios=[10,0.25],wspace=0.025,hspace=0.05)
	Ax10 = fig.add_subplot(gs[0,0])
	Ax11 = fig.add_subplot(gs[1,0])
	AxC1 = fig.add_subplot(gs[:,1])

	xBds = [uti[0], uti[-1]]
	jBds = [min(Q0),2*max(Q0)]
	kBds = [min(E),max(E)]

	kv.genCB(AxC1,norm,r'Intensity [$cm^{-2} sr^{-1} s^{-1} keV^{-1}$]',cM=cmapName,doVert=True)

	Ax10.pcolormesh(T,E,Q,norm=norm,cmap=cmapName)
	Ax10.set_facecolor('w')
	Ax10.xaxis.set_ticklabels([])
	Ax10.set_ylabel('RBSP-%s\nEnergy [keV]'%(sctag),fontsize=fntSz,family="monospace")
	Ax10.set_yscale('log')
	Ax10.set_ylim(kBds)
	Ax10.set_xlim(xBds)

	Ax11.pcolormesh(uti,Ki,jScl*J[:,:].T,norm=norm,cmap=cmapName)
	Ax11.set_facecolor('w')
	xfmt = dates.DateFormatter(utfmt)
	Ax11.xaxis.set_major_formatter(xfmt)
	Ax11.set_ylabel('Simulation \nEnergy [keV]',fontsize=fntSz,family="monospace")
	Ax11.set_yscale('log')
	Ax11.set_ylim(kBds)
	Ax11.set_xlim(xBds)
	kv.SetAxDate(Ax11)
	kv.savePic('JktComp.png')

	fig2 = plt.figure(figsize=(15,5))
	gs2 = gridspec.GridSpec(1,1,hspace=0.125)
	Ax20 = fig2.add_subplot(gs2[0,0])

	kStr = "%d [keV]"%(Ks)
	cLW  = 2
	clrs = ['black','#d95f02','#1b9e77','#7570b3'] # from colorbrewer2.org for colorblind safe

	Ax20.plot(T,Q0,color=clrs[0],label='RBSP-%s'%(sctag),linewidth=cLW)
	Ax20.plot(ut,jScl*J[:,k0f],color=clrs[1],label="Simulation",linewidth=cLW)
	Ax20.set_facecolor('w')
	xfmt = dates.DateFormatter(utfmt)
	Ax20.xaxis.set_major_formatter(xfmt)
	Ax20.set_ylabel('Intensity\n'+ r'[$cm^{-2} sr^{-1} s^{-1} keV^{-1}$]',fontsize=fntSz,family="monospace")
	Ax20.set_yscale('log')
	Ax20.set_ylim(jBds)
	Ax20.set_xlim(xBds)
	Ax20.legend(loc='upper right',fontsize=fntSz)
	Ax20.text(0.01,0.95,kStr,color="black",fontsize=fntSz,transform=Ax20.transAxes,family="monospace")
	kv.SetAxDate(Ax20)
	kv.savePic('1dJComp.png')



