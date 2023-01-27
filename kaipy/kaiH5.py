import h5py
import numpy as np
import os, sys, subprocess
import kaipy.kdefs as kdefs

#Generate MPI-style name
def genName(bStr,i,j,k,Ri,Rj,Rk,nRes=None):
	n = k + j*Rk + i*Rj*Rk
	if (nRes is None):
		fID = bStr + "_%04d_%04d_%04d_%04d_%04d_%04d.gam.h5"%(Ri,Rj,Rk,i,j,k)
	else:
		fID = bStr + "_%04d_%04d_%04d_%04d_%04d_%04d"%(Ri,Rj,Rk,i,j,k)+".gam.Res.%05d.h5"%(nRes)
	return fID
#Generate old-style MPI name
def genNameOld(bStr,i,j,k,Ri,Rj,Rk,nRes=None):
	n = k + j*Rk + i*Rj*Rk
	if (nRes is None):
		fID = bStr + "_%04d_%04d_%04d_%04d_%04d_%04d_%012d.h5"%(Ri,Rj,Rk,i,j,k,n)
	else:
		fID = bStr + "_%04d_%04d_%04d_%04d_%04d_%04d_%012d"%(Ri,Rj,Rk,i,j,k,n)+".Res.%05d.h5"%(nRes)
	return fID

#Quick check and exit routine
def CheckOrDie(fname):
	import os.path
	isExist = os.path.exists(fname)
	if (not isExist):
		if (len(fname) == 0):
			fname = "<EMPTY>"
		sys.exit("Unable to find file: %s"%(fname))

#Check directory exists and make it if not
def CheckDirOrMake(fdir):
	import os
	isDir = os.path.isdir(fdir)
	if (not isDir):
		print("Creating %s"%(fdir))
		os.makedirs(fdir)

# Stamp file with git hash (using this script's location)
def StampHash(fname):
	with h5py.File(fname, 'a') as f5:
		cwd = os.path.dirname(os.path.realpath(__file__))
		try:
			gh = subprocess.check_output(['git', '-C', cwd, 'rev-parse', 'HEAD']).decode('ascii').strip()
		except:
			print("ERROR: Couldn't grab git hash")
			gh = "XXXXX"
		f5.attrs['GITHASH'] = gh

# Stamp file with git branch (using this script's location)
def StampBranch(fname):
	with h5py.File(fname, 'a') as f5:
		cwd = os.path.dirname(os.path.realpath(__file__))
		try:
			gb = subprocess.check_output(['git', '-C', cwd, 'rev-parse', '--abbrev-ref', 'HEAD']).decode('ascii').strip()
		except:
			print("ERROR: Couldn't grab git branch")
			gb = "XXXXX"
		f5.attrs['GITBRANCH'] = gb

#Get git hash from file if it exists
def GetHash(fname):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		rStr = hf.attrs.get("GITHASH","NONE")
	try:
		hStr = rStr.decode('utf-8')
	except (UnicodeDecodeError, AttributeError):
		hStr = rStr
	return hStr

#Get git branch from file if it exists
def GetBranch(fname):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		rStr = hf.attrs.get("GITBRANCH","NONE")
	try:
		hStr = rStr.decode('utf-8')
	except (UnicodeDecodeError, AttributeError):
		hStr = rStr
	return hStr

#Return time at step n
def tStep(fname,nStp=0,aID="time",aDef=0.0):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		gID = "Step#%d"%(nStp)
		#t = hf[gID].attrs.get(aID,aDef)
		if aID in hf[gID].attrs:
			t = hf[gID].attrs[aID]
		else:
			t = aDef
	return t
	
def cntSteps(fname,doTryRecover=True,s0=0):
	from alive_progress import alive_it

	try:
		CheckOrDie(fname)
		with h5py.File(fname,'r') as hf:
			Steps = [grp for grp in alive_it(hf.keys(),title="#-Steps".ljust(kdefs.barLab),length=kdefs.barLen) if "Step#" in grp]
		sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=np.int)
		sIds.sort()
		nSteps = len(Steps)
		
		return nSteps,sIds
	except (ValueError, IndexError) as e:
		print("!!Warning: h5 file contains unreadable steps")

	if not doTryRecover:
		print("  Can try again with 'cntSteps(fname,doTryRecover=True)'")
		raise ValueError("Unreadable steps")
	else:
		print("Trying to read again while skipping bad values")
		with h5py.File(fname,'r') as hf:
			badCounts = 0
			s = s0
			sIds = []
			while badCounts < 2000:
				try:
					sName = 'Step#' + str(s)
					tryReadGrp = hf[sName]
					#If still here, step is readable
					sIds.append(s)
					print("Read step " + str(s))
				except ValueError:
					badCounts += 1
					print("Bad count on step " + str(s))
				s+=1
			nSteps = len(sIds)
			sIds = np.array(sIds)
			sIds.sort()
		return nSteps,sIds

#More general version of cntSteps, useful for Step#X/Line#Y
def cntX(fname,gID=None,StrX="/Step#"):
	with h5py.File(fname,'r') as hf:
		if (gID is not None):
			grps = hf[gID].values()
		else:
			grps = hf.values()
		grpNames = [str(grp.name) for grp in grps]
		#Steps = [stp if "/Step#" in stp for stp in grpNames]
		Steps = [stp for stp in grpNames if StrX in stp]
		nSteps = len(Steps)

		sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=np.int)
		sIds.sort()
		return nSteps,sIds

def getTs(fname,sIds=None,aID="time",aDef=0.0):
	from alive_progress import alive_bar
	
	if (sIds is None):
		nSteps,sIds = cntSteps(fname)
	Nt = len(sIds)
	T = np.zeros(Nt)
	CheckOrDie(fname)
	titStr = "Time series: %s"%(aID)

	with h5py.File(fname,'r') as hf, alive_bar(Nt,title=titStr.ljust(kdefs.barLab),length=kdefs.barLen) as bar:
		for idx, n in enumerate(sIds):
			gId = "/Step#%d"%(n)
			T[idx] = hf[gId].attrs.get(aID,aDef)
			bar()
	return T

#Get shape/dimension of grid
def getDims(fname,vID="X",doFlip=True):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		Dims = hf["/"][vID].shape
	Ds = np.array(Dims,dtype=np.int)
	if (doFlip):
		Ds = np.flip(Ds,axis=0)
	return Ds

#Get root variables
def getRootVars(fname):
	CheckOrDie(fname)

	with h5py.File(fname,'r') as hf:
		vIds = [str(k) for k in hf.keys() if "Step" not in str(k)]
	#Remove coordinates from list of root variables
	xyzS = ["X","Y","Z"]
	vIds = [v for v in vIds if v not in xyzS]

	return vIds

#Get variables in initial Step
def getVars(fname,smin):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		gId = "/Step#%d"%(smin)
		stp0 = hf[gId]
		vIds = [str(k) for k in stp0.keys()]
	return vIds

#Get variable data
def PullVar(fname,vID,s0=None):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		if (s0 is None):
			V = hf[vID][()].T
		else:
			gId = "/Step#%d"%(s0)
			V = hf[gId][vID][()].T
	return V

#Get attribute data from Step#s0 or root (s0=None)
def PullAtt(fname,vID,s0=None):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		if (s0 is None):
			Q = hf.attrs[vID]
			
		else:
			gID = "/Step#%d"%(s0)
			Q = hf[gID].attrs[vID]
	return Q
