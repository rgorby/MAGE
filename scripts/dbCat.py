#!/usr/bin/env python
#Joins decomposed DB files into single

import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import h5py
import os
import kaipy.kaiH5 as kh5
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
		if "Step" not in sQ:
			oH5.create_dataset(sQ,data=iH5[sQ])
			print("\t%s"%(sQ))
	iH5.close()

	return oH5

if __name__ == "__main__":
	dIn = os.getcwd()

	runid = "msphere"
	typid = "deltab"
	MainS = """Joins blocks created by calcdb.x (or similar CHIMP routines) into single file
	
	runid : Run ID

	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-runid',type=str,metavar="runid",default=runid,help="Input run ID (default: %(default)s)")
	parser.add_argument('-typeid',type=str,metavar="typeid",default=typeid,help="Input type ID (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	runid = args.runid
	typeid = args.typeid
	
	dbIns = glob.glob('%s.????.%s.h5'%(runid,typeid))
	dbIns.sort()
	
	fOut = "%s.%s.h5"%(runid,typeid)

	N = len(dbIns)
	print("Found %d files, writing output to %s"%(N,fOut))
	if (N == 0):
		print("No files found, exiting")
		exit()
	#Create file w/ attributes and root variables as first file
	oH5 = createfile(dbIns[0],fOut)

	s0 = 0 #Current step

	#Now loop over files
	for n in range(N):
		fIn = dbIns[n]
		
		Ns,sIDs = kh5.cntSteps(fIn)
		nS = sIDs.min()
		nE = sIDs.max()
		dN = nE-nS+1
		print("Reading steps %d to %d from %s"%(nS,nE,fIn))
		print("\tWriting to %d to %d"%(s0,s0+dN-1))
		iH5 = h5py.File(fIn,'r')

		#Loop over steps in the input file
		for s in range(nS,nE+1):
			#Input
			igStr = "Step#%d"%(s)
			ogStr = "Step#%d"%(s-nS+s0)

			oH5.create_group(ogStr)
			#print("Copying %s to %s"%(igStr,ogStr))

			#Group atts
			for k in iH5[igStr].attrs.keys():

				aStr = str(k)
				oH5[ogStr].attrs.create(k,iH5[igStr].attrs[aStr])
				#print(aStr)
			#Group vars
			for Q in iH5[igStr].keys():
				sQ = str(Q)
				#print("\tCopying %s"%(sQ))
				oH5[ogStr].create_dataset(sQ,data=iH5[igStr][sQ])

		iH5.close()
		#Update s0
		s0 = s0 + dN
	#Done
	oH5.close()
