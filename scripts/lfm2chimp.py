#!/usr/bin/env python
#Turns collection of LFM HDF-4 data to Chimp-style (fake Gamera)


import argparse
import sys
import numpy as np
import os
import gamGrids as gg
import lfm2kaiju as l2k
import h5py
from glob import glob

if __name__ == "__main__":
	#Defaults
	fOut = "ebLFM.h5"

	parser = argparse.ArgumentParser(description="Converts collection of LFM HDF-4 data to Chimp-style (fake Gamera)")
	#Get required info
	parser.add_argument('-o',type=str,metavar="Output file name",default=fOut,help="File to write output grid to (default: %(default)s)")
	#Do only EB
	parser.add_argument('--mhd',dest='doMHD',action='store_true',default=False,help="Do all MHD variables (default: %(default)s)")
	#Do Jupiter
	parser.add_argument('--jupiter',dest='doJupiter',action='store_true',default=False,help="Pull from Jovian grid (default: %(default)s)")

	#Files to interpolate from
	parser.add_argument('hdfs',nargs='+',metavar='lfm.hdf',help="List of files to convert")

	#Finished getting arguments, parse and move on
	args = parser.parse_args()

	hdfs = list(args.hdfs)
	lfmfile = hdfs[0]

	fOut = args.o
	doMHD = args.doMHD
	doJupiter = args.doJupiter


	print("Writing out to %s"%fOut)
	#Get grid from LFM file and write to Gamera-style H5
	if (doJupiter):
		l2k.lfm2gg(lfmfile,fOut,doEarth=False,doJupiter=True)
	else:
		l2k.lfm2gg(lfmfile,fOut,doEarth=True,doJupiter=False)

	#Get time slice information
	Ts = l2k.lfmTimes(hdfs)
	T0 = Ts.min() #Smallest time

	#Loop through sorted timeslices
	lfmSlcs = sorted(zip(Ts-T0,hdfs))
	n = 0

	oH5 = h5py.File(fOut,'r+')

	for lS in lfmSlcs:
		fIn = lS[1]
		#print("Reading %s"%(fIn))

		gID = "Step#%d"%(n)
		#Create group
		oH5.create_group(gID)

		#Copy time
		oH5[gID].attrs.create("time",lS[0])

		#Get LFM V/B fields
		Vx,Vy,Vz,Bx,By,Bz = l2k.lfmFields(fIn)


		oH5[gID].create_dataset("Bx",data=np.single(Bx))
		oH5[gID].create_dataset("By",data=np.single(By))
		oH5[gID].create_dataset("Bz",data=np.single(Bz))

		oH5[gID].create_dataset("Vx",data=np.single(Vx))
		oH5[gID].create_dataset("Vy",data=np.single(Vy))
		oH5[gID].create_dataset("Vz",data=np.single(Vz))

		#Do other MHD variables if requested
		if (doMHD):
			D,P = l2k.lfmFlow(fIn)
			oH5[gID].create_dataset("D",data=np.single(D))
			oH5[gID].create_dataset("P",data=np.single(P))

		n = n+1
	oH5.close()
