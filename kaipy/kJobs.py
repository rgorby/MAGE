#Various routines to automate generation of run scripts
#XML input decks, PBS scripts
import xml.etree.ElementTree as et
import xml.dom.minidom

XV0 = 23 #Magic number for XML header length

#Various machine defaults, for now just assuming Cheyenne

Nth = 72 #Thread number
Nc  = 72 #CPU number
stkSz = 128
qStr = "regular"
wStr = "12:00:00"
aStr = "UJHB0010"


#----------------
#XML routines

#Create XML block
#Root name (ie Gamera)
#List of block names (sim,prob,xxx)
#List of dictionaries for each block
def genXBlk(rStr,bStrs,bDicts,doHeader=False):
	rXML = et.Element(rStr)
	Nd = len(bStrs)
	for i in range(Nd):
		et.SubElement(rXML,bStrs[i],bDicts[i])

	oStr = xml2str(rXML,doHeader=doHeader)
	oStr = oStr[:-1] #Remove trailing newline
	
	return oStr

#Adds label/dictionary to list
def addXBlk(bStr,bDic,bStrs,bDics):
	bStrs.append(bStr)
	bDics.append(bDic)
#----------------
#PBS routines


#Create OMP job script
#ComX = command string, ie "pusher.x"
#inStr = input deck name (no trailing .xml)
#mStr = module to load, ie module restore XXX


#Will run,
#ComX inStr.xml > inStr.out

#For job array runs
#Submit w/ qsub -J 1-2 XXX.pbs
def OMPJob(fOut,ComX,inStr,Nthz=None,mod="kaiju",doArray=True):

	pbsT = """
#!/bin/bash
#PBS -A %s
#PBS -N %s
#PBS -j oe
#PBS -q %s
#PBS -l walltime=%s
#PBS -l select=%s

source ~/.bashrc
"""
	#Generate core/thread string
	if (Nthz is None):
		Nthz = Nth
	cStr = "1:ncpus=%d:ompthreads=%d"%(Nc,Nthz)
	pbsHD = pbsT % (aStr,inStr,qStr,wStr,cStr)

	if (mod is not None):
		mStr = "module restore %s"%(mod)
		pbsHD = pbsHD + mStr + "\n"
	pbsB = """
module list
hostname
date
export OMP_NUM_THREADS=%d
export OMP_STACKSIZE=%dM
export JNUM=${PBS_ARRAY_INDEX:-0}
%s
date
"""
	if (doArray):
		xStr = "omplace -nt $OMP_NUM_THREADS %s %s.xml ${JNUM} > %s.${JNUM}.out"%(ComX,inStr,inStr)
	else:		
		xStr = "omplace -nt $OMP_NUM_THREADS %s %s.xml > %s.out"%(ComX,inStr,inStr)
	pbsX = pbsB % (Nthz,stkSz,xStr)

	oStr = pbsHD+pbsX
	with open(fOut,"w") as fID:
		fID.write(oStr)

	#print(oStr)


#----------------
#Utilities

#Turn XML data to string, optional strip XML version header
def xml2str(xIn,doHeader=False):
	xmlStr = xml.dom.minidom.parseString(et.tostring(xIn)).toprettyxml(indent="    ")
	xOut = xmlStr
	if (not doHeader):
		xOut = xmlStr[23:]

	return xOut

