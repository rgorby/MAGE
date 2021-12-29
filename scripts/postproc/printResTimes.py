#!/usr/bin/env python3
"""
Overengineered script to print the time of each restart file by parsing h5dump output
Uses only packages available on Cheyenne
"""
import subprocess
import glob
import argparse
from argparse import RawTextHelpFormatter

def sortFn(elem): #Used to sort final list in order of nRes
    	return int(elem['nRes'])
		
def getAttrKeyValue(lineList):
	for line in lineList:
		if 'ATTRIBUTE' in line:
			key = line.split('"')[1]
		if '(0)' in line:
			value = line.split('(0):')[1]
	return key, value

#Return dictionary of attrs for a single restart file
def getKVsFromFile(fName):
	spOutput = subprocess.check_output(['h5dump','-a','nRes','-a','t','-a','DATETIME',fName])
	output = spOutput.decode('utf-8').split('\n')

	#Parse attributes (this could be better)
	attrLocs = []
	for i in range(len(output)):
		if 'ATTRIBUTE' in output[i]:
			attrLocs.append(i)
	
	attrs = {} 
	for i in range(len(attrLocs)):
		if i == len(attrLocs)-1:
			k,v = getAttrKeyValue(output[attrLocs[i]:])
		else:
			k,v = getAttrKeyValue(output[attrLocs[i]:attrLocs[i+1]])
		attrs[k] = v
	return attrs

if __name__=='__main__':
	idStr_noMPI = ".gam.Res.*.h5"
	idStrMPI = "_0*_0*_0*_0000_0000_0000.gam.Res.*.h5"

	ftag = "msphere"
	timeFmt = "m"
	MainS = """Overengineered script to print the time of each restart file by parsing h5dump output
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-f',type=str,metavar="timeFmt",default=timeFmt,help="Time format [s,m,h] (default: %(default)s)")
	parser.add_argument('-dt',action='store_true',default=False,help="Print the datetime instead of simulation time (default: %(default)s)")
	args = parser.parse_args()
	ftag = args.id
	timeFormat = args.f
	doDatetime = args.dt

	if timeFormat not in ['s','m','h']:
		print('Unrecognized value "%s" for time format. Using "%s"'%(timeFormat, timeFmt))
		timeFormat = timeFmt
	if timeFormat == 's':
		timeFormat = 'sec'
		timeMult = 60
	elif timeFormat == 'h':
		timeFormat = 'hr'
		timeMult = 1./60
	else:
		timeFormat = 'min'
		timeMult = 1
	
	#Get list of files
	globStr = ftag + idStr_noMPI
	print("Looking for nonMPI restarts...",end='')
	fileList = glob.glob(globStr)
	if len(fileList) == 0:
		globStr = ftag + idStrMPI
		print("Not found\nLooking for MPI restarts...",end='')
		fileList = glob.glob(globStr)
	if len(fileList) == 0:
		print("Not found\nCheck id (globStr = %s)"%(globStr))
		quit()
	print("Found")

	attrList = []
	#Build list of attr dicts from list of files
	for fStr in fileList:
		if 'XXXXX' in fStr:
			continue
		fAttrs = getKVsFromFile(fStr)
		fAttrs['fname'] = fStr
		attrList.append(fAttrs)

	#Print list (time from Gam restarts only, for now)
	attrList.sort(key=sortFn)
	for entry in attrList:
		if doDatetime:
			formatString = " {}: {}".format(entry['fname'], entry['DATETIME'])
		else:
			formatString = " {}: {:4.2f} [{}]".format(entry['fname'], float(entry['t'])*timeMult, timeFormat)
		print(formatString)


