#!/usr/bin/env python
# Joins decomposed eb files generated using "Parallel In Time" into a single file
# Example: bit.ly/3OQg71F

import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import h5py
import os
import kaipy.kaiH5 as kh5
import glob
import kaipy.kdefs as kd

tEps = 1.0e-3 #Small time

#Create new file w/ same root vars/attributes as old
def createfile(fIn,fOut,doLink=False):
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
		#Skip cache, we add it later
		if kd.grpTimeCache in sQ:
			continue
		#Don't include stuff that starts with "Step"
		if "Step" not in sQ:
			if doLink:
				oH5[sQ] = h5py.ExternalLink(fIn, sQ)
			else:
				oH5.create_dataset(sQ,data=iH5[sQ])
			print("\t%s"%(sQ))
	iH5.close()

	return oH5

if __name__ == "__main__":
	dIn = os.getcwd()

	runid = "msphere"
	typeid = "deltab"
	MainS = """Joins blocks created by calcdb.x (or similar CHIMP routines) into single file
	
	runid : Run ID

	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-runid',type=str,metavar="runid",default=runid,help="Input run ID (default: %(default)s)")
	parser.add_argument('-typeid',type=str,metavar="typeid",default=typeid,help="Input type ID (default: %(default)s)")
	parser.add_argument('--link',action='store_true',help="Create links to existing files rather than copy data (default: %(default)s)")

	#Finalize parsing
	args = parser.parse_args()
	runid = args.runid
	typeid = args.typeid
	doLink = args.link

	globStr = '%s.????.%s.h5'%(runid,typeid)
	dbIns = glob.glob(globStr)
	dbIns.sort()

	if doLink:
		fOut = "%s.%s.link.h5"%(runid,typeid)
	else:
		fOut = "%s.%s.h5"%(runid,typeid)

	N = len(dbIns)
	print("Found %d files, writing output to %s"%(N,fOut))
	if (N == 0):
		print("No files found, exiting")
		exit()
	#Create file w/ attributes and root variables as first file
	oH5 = createfile(dbIns[0],fOut, doLink)

	# Store and concat timeAttributeCache to add at the very end
	timeCacheVars = {}
	with h5py.File(dbIns[0], 'r') as tacF:
		for k in tacF[kd.grpTimeCache].keys():
			timeCacheVars[k] = np.array([], dtype=tacF[kd.grpTimeCache][k].dtype)

	s0 = 0 #Current step
	nowTime = 0.0
	oldTime = -np.inf

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

		# Grow timeAttributeCache
		for k in iH5[kd.grpTimeCache].keys():
			data = iH5[kd.grpTimeCache][k][:]
			if k == 'step':
				data += s0  # Cache for merged h5 file needs to remap original steps to their position in merged file
			timeCacheVars[k] = np.append(timeCacheVars[k], data, axis=0)

		#Loop over steps in the input file
		for s in range(nS,nE+1):
			#Input
			igStr = "Step#%d"%(s)
			ogStr = "Step#%d"%(s0)

			#Check if this is too close to last value
			nowTime = kh5.tStep(fIn,s)
			#print(nowTime,oldTime)
			if ( np.abs(nowTime-oldTime)<=tEps):
				print("\tSkipping step %d"%(s))
				continue
			else:
				#Good value, update old time
				oldTime = nowTime

			if doLink:
				oH5[ogStr] = h5py.ExternalLink(fIn, igStr)
			else:
				oH5.create_group(ogStr)
				print("Copying %s to %s"%(igStr,ogStr))

				#Group atts
				for k in iH5[igStr].attrs.keys():

					aStr = str(k)
					oH5[ogStr].attrs.create(k,iH5[igStr].attrs[aStr])
					#print(aStr)
				#Group vars
				for Q in iH5[igStr].keys():
					sQ = str(Q)
					oH5[ogStr].create_dataset(sQ,data=iH5[igStr][sQ])
			#Update s0
			s0 = s0 + 1

		iH5.close()

	# Write timeAttributeCache to output file
	print("Writing " + kd.grpTimeCache)
	tag = oH5.create_group(kd.grpTimeCache)
	for k in timeCacheVars:
		tag.create_dataset(k, data=timeCacheVars[k], dtype=timeCacheVars[k].dtype)
	
	#Done
	oH5.close()
