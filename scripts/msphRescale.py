#!/usr/bin/env python
#Up/Down-scales Gamera magnetosphere restart

import argparse
import kaipy.gamera.gamGrids as gg
import kaipy.gamera.magsphereRescale as upscl
from argparse import RawTextHelpFormatter
import numpy as np
import h5py

if __name__ == "__main__":
	rStrs = ['U','D']
	fIn = "msphere.00000.h5"
	fOut = "msphere2x.00000.h5"

	MainS = """Up/Down-scales Gamera magnetosphere restart
	Run types (rx)
	U: Upscale, 2x in each dimension
	D: Downscale, 1/2x in each dimension
	"""	
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-rx',type=str,default="U",choices=rStrs,help="Scaling Specifier (default: %(default)s)")
	parser.add_argument('-i',type=str,metavar="file",default=fIn ,help="Input Restart HDF5 (default: %(default)s)")
	parser.add_argument('-o',type=str,metavar="file",default=fOut,help="Output Restart HDF5 (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	rX = args.rx
	fIn = args.i
	fOut = args.o

	#Now do work
	print("Rescaling %s to %s"%(fIn,fOut))

	#Open output
	oH5 = h5py.File(fOut,'w')

	#Open input
	print("Reading from %s"%(fIn))
	iH5 = h5py.File(fIn,'r')
	Ns,Nv,Nk,Nj,Ni = iH5['Gas'].shape
	G = np.zeros((Ns,Nv,Nk,Nj,Ni))
	M = np.zeros((3,Nk+1,Nj+1,Ni+1))
	G[:,:,:,:,:] = iH5['Gas'][:]
	M[  :,:,:,:] = iH5['magFlux'][:]
	X = iH5['X'][:]
	Y = iH5['Y'][:]
	Z = iH5['Z'][:]
	doGas0 = ('Gas0' in iH5.keys()) #Whether there's a Gas0 array
	if (doGas0):
		G0 = np.zeros((Ns,Nv,Nk,Nj,Ni))
		G0[:,:,:,:,:] = iH5['Gas0'][:]
	#Transfer attributes to output
	for k in iH5.attrs.keys():
		aStr = str(k)
		oH5.attrs.create(k,iH5.attrs[aStr])

	#Close input
	iH5.close()

	if (rX == "U"):
		print("Upscaling ...")
		Xr,Yr,Zr = upscl.upGrid(X,Y,Z)
		Gr = upscl.upGas(X,Y,Z,G,Xr.T,Yr.T,Zr.T)
		FluxR = upscl.upFlux(X,Y,Z,M,Xr,Yr,Zr)
		if (doGas0):
			G0r = upscl.upGas(X,Y,Z,G0,Xr.T,Yr.T,Zr.T)
	else:
		print("Downscaling ...")
		Xr,Yr,Zr = upscl.downGrid(X,Y,Z)
		Gr = upscl.downGas(X,Y,Z,G,Xr.T,Yr.T,Zr.T)
		FluxR = upscl.downFlux(X,Y,Z,M,Xr,Yr,Zr)
		if (doGas0):
			G0r = upscl.downGas(X,Y,Z,G0,Xr.T,Yr.T,Zr.T)
	print("Writing to %s"%(fOut))

	#Write out grid to restart
	oH5.create_dataset("X",data=Xr.T)
	oH5.create_dataset("Y",data=Yr.T)
	oH5.create_dataset("Z",data=Zr.T)

	#Write out gas/flux variables
	oH5.create_dataset("Gas",data=Gr)
	oH5.create_dataset("magFlux",data=FluxR)
	if (doGas0):
		oH5.create_dataset("Gas0",data=G0r)
	#Close output
	oH5.close()


