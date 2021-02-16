import h5py
import numpy as np

#Generate MPI-style name
def genName(bStr,i,j,k,Ri,Rj,Rk,nRes=None):
	n = k + j*Rk + i*Rj*Rk
	if (nRes is None):
		fID = bStr + "_%04d_%04d_%04d_%04d_%04d_%04d.gam.h5"%(Ri,Rj,Rk,i,j,k)
	else:
		fID = bStr + "_%04d_%04d_%04d_%04d_%04d_%04d"%(Ri,Rj,Rk,i,j,k)+".Res.%05d.h5"%(nRes)
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
		print("Unable to find file: %s"%(fname))
		print("!!Exiting!!")
		quit()

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
#Return time at step n
def tStep(fname,nStp=0,aID="time",aDef=0.0):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		gID = "Step#%d"%(nStp)
		t = hf[gID].attrs.get(aID,aDef)
	return t
	
def cntSteps(fname):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		Steps = [grp for grp in hf.keys() if "Step#" in grp]
	sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=np.int)
	nSteps = len(Steps)
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
		return nSteps,sIds

def getTs(fname,sIds=None,aID="time",aDef=0.0):
	if (sIds is None):
		nSteps,sIds = cntSteps(fname)
	Nt = len(sIds)
	T = np.zeros(Nt)
	i0 = sIds.min()
	i1 = sIds.max()
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		for n in range(i0,i1+1):
			gId = "/Step#%d"%(n)
			T[n-i0] = hf[gId].attrs.get(aID,aDef)
	return T

#Get shape/dimension of grid
def getDims(fname,doFlip=True):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		Dims = hf["/"]["X"].shape
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
def getVars(fname,s0):
	CheckOrDie(fname)
	with h5py.File(fname,'r') as hf:
		gId = "/Step#%d"%(s0)
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
