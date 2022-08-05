#!/usr/bin/env python
#Simple script to spit out the step information of a Kaiju H5 file
import argparse
from argparse import RawTextHelpFormatter
import os
import kaipy.kaiH5 as kh5
import numpy as np

def MJD2Str(m0):
	from astropy.time import Time
	dtObj = Time(m0,format='mjd').datetime
	tStr = dtObj.strftime("%H:%M:%S") + " " +  dtObj.strftime("%m/%d/%Y") 
	return tStr
if __name__ == "__main__":
	#Defaults

	MainS = """Identifies the domain (in steps and time) of a Kaiju HDF-5 file"""
	
	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('h5F',nargs='+',metavar='Gamera.h5',help="Filename of Gamera HDF5 Output")

	#Finished getting arguments, parse and move on
	args = parser.parse_args()
	#h5F = args.h5F
	for idx, h5F in enumerate(args.h5F):
		print("Reading %s"%(h5F))
		nSteps,sIds = kh5.cntSteps(h5F)
		s0 = sIds.min()
		sE = sIds.max()
		print("\tFound %d steps"%(nSteps))
		print("\tSteps = [%d,%d]"%(s0,sE))
		tMin = kh5.getTs(h5F,np.array([s0]))
		tMax = kh5.getTs(h5F,np.array([sE]))
		print("\tTime = [%f,%f]"%(tMin,tMax))
		MJDs = kh5.getTs(h5F,sIds,"MJD",-np.inf)
		if (MJDs.max()>0):
			MJDMin = MJDs.min()
			MJDMax = MJDs.max()

			print("\tMJD  = [%f,%f]"%(MJDMin,MJDMax))
			tS1 = MJD2Str(MJDMin)
			tS2 = MJD2Str(MJDMax)
			print("\t\tStart: %s"%(tS1))
			print("\t\tStop : %s"%(tS2))

		hStr = kh5.GetHash(h5F)
		bStr = kh5.GetBranch(h5F)

		print("\tGit: Hash = %s / Branch = %s"%(hStr,bStr))
	#---------------------
