#!/usr/bin/env python
#Joins decomposed var gamera output files into single h5 file

import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import h5py
import os
import kaipy.kaiH5 as kh5
import kaipy.gamhelio.heliosphere as hsph
import glob

#Create new file w/ same root vars/attributes as old
def createfile(fIn,fOut):
	print('Creating new output file:',fOut)
	iH5 = h5py.File(fIn,'r')
	oH5 = h5py.File(fOut,'w')
	#Start by scraping all variables from root
	#Copy root attributes
	print("Copying root attributes ...")
	for k in iH5.attrs.keys():
		aStr = str(k)
		oH5.attrs.create(k,iH5.attrs[aStr])
		print("\t%s"%(aStr))
    #Copy root groups
	print("Copying root variables ...")
	for Q in iH5.keys():
		sQ = str(Q)
		#Don't include stuff that starts with "Step"
		if "Step" not in sQ and "timeAttributeCache" not in sQ:
			oH5.create_dataset(sQ,data=iH5[sQ])
			print("\t%s"%(sQ))
	iH5.close()

	return oH5

if __name__ == "__main__":
	dIn = os.getcwd()

	runid = "msphere"
	typeid = "gam"
	varid = "D"
	default_slice = "0:-1:1"
	MainS = """Joins from an MPI gamera run for a single variable, with nslice number of time steps into single file
	
	runid : Run ID

	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-runid',type=str,metavar="runid",default=runid,help="Input run ID (default: %(default)s)")
	parser.add_argument('-typeid',type=str,metavar="typeid",default=typeid,help="Input type ID (default: %(default)s)")
	parser.add_argument('-varid',type=str,metavar="var",default=varid,help="Input var (default: %(default)s)")
	parser.add_argument(
        "-d", type=str, metavar="directory", default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
	parser.add_argument(
        "--nslice", type=lambda n: [int(item) for item in n.split(':')], metavar="step slice", default=default_slice,
        help="Slice for range of time slice(s) to plot (default: %(default)s)"
    )    
	parser.add_argument(
        "-p", "--parallel", action="store_true", default=False,
        help="Read from HDF5 in parallel (default: %(default)s)."
    )
	parser.add_argument(
        "-nw", "--nworkers", type=int, metavar="nworkers", default=4,
        help="Number of parallel workers (default: %(default)s)"
    )
	#Finalize parsing
	args = parser.parse_args()
	runid = args.runid
	typeid = args.typeid
	varid = args.varid
	slices = args.nslice
	fdir = args.d
	doFast = False
	doParallel = args.parallel
	nWorkers = args.nworkers

	fOut = "%s.%s.%s.h5"%(runid,varid,typeid)

	fIn = glob.glob('%s*.%s.h5'%(runid,typeid))
	fIn.sort()

	N = len(fIn)

	print("Found %d files, writing output to %s"%(N,fOut))
	if (N == 0):
		print("No files found, exiting")
		exit()
	#Create file w/ attributes and root variables as first file
	oH5 = createfile(fIn[0],fOut)

	gsph = hsph.GamsphPipe(fdir, runid, doFast=doFast, doParallel=doParallel, nWorkers=nWorkers)
	s0 = 0
	nS = slices[0]
	nE = slices[1]
	nStride = slices[2]
	dN = nE-nS+1
	print("Reading steps %d to %d every %d steps"%(nS,nE,nStride))
	print("\tWriting steps from %d to %d"%(s0,s0+dN-1))
	iH5 = h5py.File(fIn[0],'r')

	#Loop over steps in the input file
	for s in range(nS,nE+1,nStride):
		#Input
		igStr = "Step#%d"%(s)
		ogStr = "Step#%d"%(s0)
		
		oH5.create_group(ogStr)
		print("Copying %s to %s"%(igStr,ogStr))

		#Group atts
		for k in iH5[igStr].attrs.keys():
			aStr = str(k)
			oH5[ogStr].attrs.create(k,iH5[igStr].attrs[aStr])
			#print(aStr)
		#Group var
		print("\tCopying %s"%(varid))
		Q = gsph.GetVar(varid,s)
		oH5[ogStr].create_dataset(varid,data=Q)
		#Update s0
		s0 = s0 + 1

	iH5.close()
	#Done
	oH5.close()
