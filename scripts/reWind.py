#!/usr/bin/env python

#Converts LFM-style wind file to Gamera

#Reads from ASCII LFM wind file
#Time(min) Density (AMU/cm^-3) Vx(km/s) Vy(km/s) Vz(km/s) Cs(km/s) Bx(nT) By(nT) Bz(nT) B(nT) tilt(rad)

#Writes to HDF5 Gamera wind file
#t,D,V,P,B = [s],[#/cm3],[m/s],[nPa],[nT]

Mp = 1.67e-27 #Proton mass [kg]
gamma = 5/3.0

import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import h5py

if __name__ == "__main__":
	fOut = "bcwind.h5"
	fIn = "SW-SM-DAT"
	nSkip = 3 #Number of header lines to skip

	MainS = """Converts LFM-style solar wind file to Gamera format
	LFM wind parameters (ASCII) / Gamera wind parameters (HDF5)
	LFM   : t,D,V,Cs,B,tilt = min, AMU/cm3, km/s, km/s, nT, rad 
	Gamera: t,D,V,P ,B,tilt = sec, AMU/cm3, m/s , nPa , nT, rad

	"""
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('fIn',type=str,metavar="sw.txt",default=fIn,help="Input LFM wind file (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="wind.h5",default=fOut,help="Output Gamera wind file (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()

	fIn = args.fIn
	fOut = args.o

	print("Reading LFM solar wind from %s"%(fIn))
	lfmD = np.genfromtxt(fIn,skip_header=nSkip)
	Nt,Nv = lfmD.shape
	print("\tFound %d variables and %d lines"%(Nv,Nt))

	#Create holders for Gamera data
	T  = np.zeros(Nt)
	D  = np.zeros(Nt)
	P  = np.zeros(Nt)
	Vx = np.zeros(Nt)
	Vy = np.zeros(Nt)
	Vz = np.zeros(Nt)
	Bx = np.zeros(Nt)
	By = np.zeros(Nt)
	Bz = np.zeros(Nt)
	ThT= np.zeros(Nt) #Tilt

	#Convert LFM time to seconds and reset to start at 0
	T0 = lfmD[:,0].min()
	T = (lfmD[:,0]-T0)*60

	#Density, magnetic field, and tilt don't require scaling
	D   = lfmD[:,1]
	ThT = lfmD[:,10]
	Bx  = lfmD[:,6]
	By  = lfmD[:,7]
	Bz  = lfmD[:,8]

	#Velocity
	vScl = 1.0e+3 #km/s->m/s
	Vx  = vScl*lfmD[:,2]
	Vy  = vScl*lfmD[:,3]
	Vz  = vScl*lfmD[:,4]

	#Now get pressure (nPa) from sound speed (km/s)
	#Convert density: AMU/cm3 -> kg/m3
	Cs = lfmD[:,5]
	Dmks = (D*Mp)*( (1.0e+2)**3.0 ) #kg/m3
	Cmks = Cs*1.0e+3 #m/s
	P = (1.0e+9)*Dmks*Cmks*Cmks/gamma #nPa

print("Writing Gamera solar wind to %s"%(fOut))
with h5py.File(fOut,'w') as hf:
	hf.create_dataset("T" ,data=T)
	hf.create_dataset("D" ,data=D)
	hf.create_dataset("P" ,data=P)
	hf.create_dataset("Vx",data=Vx)
	hf.create_dataset("Vy",data=Vy)
	hf.create_dataset("Vz",data=Vz)
	hf.create_dataset("Bx",data=Bx)
	hf.create_dataset("By",data=By)
	hf.create_dataset("Bz",data=Bz)
	hf.create_dataset("tilt",data=ThT)


